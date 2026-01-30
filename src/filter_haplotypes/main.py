"""
Main entry point for the FilterHaplotypes command-line tool.
This module orchestrates the entire de-duplication pipeline, from parsing input files
to generating the final interactive report and filtered assembly.
"""

import argparse
import logging
import multiprocessing
import os
import sys
from pathlib import Path
from typing import List, Dict, Any, Tuple

from src.filter_haplotypes.parsers.fasta_parser import parse_fasta
from src.filter_haplotypes.parsers.paf_parser import parse_paf, get_primary_targets
from src.filter_haplotypes.parsers.mash_parser import parse_mash, build_mash_lookup, get_mash_distance
from src.filter_haplotypes.parsers.busco_parser import parse_busco, get_busco_counts
from src.filter_haplotypes.core.models import ContigSummary, Status
from src.filter_haplotypes.core.filtering import (
    calculate_initial_redundancy,
    tile_and_score_contig,
    get_overlapping_pairs,
    estimate_distance_threshold,
    run_tournament_on_target,
    screen_unaligned_contig
)
from src.filter_haplotypes.utils.logging import setup_logging, worker_configurer
from src.filter_haplotypes.utils.stats import calculate_assembly_stats
from src.filter_haplotypes.visualization.report_generator import generate_report, write_filtered_fasta

def main():
    parser = argparse.ArgumentParser(
        description="FilterHaplotypes: Reference-based de-duplication of genome assemblies.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Mandatory
    parser.add_argument("-p", "--paf", required=True, help="Alignment file (PAF) generated via minimap2 -c")
    parser.add_argument("-m", "--mash", required=True, help="Mash distances TSV (Query vs Query)")
    parser.add_argument("-f", "--fasta", required=True, help="Original assembly FASTA file")
    
    # Optional
    parser.add_argument("-b", "--busco", help="Optional BUSCO full_table.tsv")
    parser.add_argument("-o", "--output", default="./output", help="Output directory for results")
    
    # Configurable
    parser.add_argument("--min-mq", type=int, default=10, help="Mapping Quality threshold")
    parser.add_argument("--min-overlap", type=int, default=1, help="Minimum overlap bases to trigger overlap condition")
    parser.add_argument("--min-size-safeguard", type=float, default=0.50, help="Ratio for size-based retention")
    parser.add_argument("--distance-threshold", type=float, help="Mash distance threshold (Overrides estimator if supplied)")
    parser.add_argument("--threads", type=int, default=max(1, multiprocessing.cpu_count() - 1), help="Number of CPU cores for parallel processing")
    parser.add_argument("--max-tournament-iterations", type=int, default=100000, help="Maximum iterations for tournament loop")

    args = parser.parse_args()
    
    output_dir = Path(args.output)
    log_queue, log_listener = setup_logging(output_dir)
    
    logger = logging.getLogger(__name__)
    try:
        logger.info("Starting FilterHaplotypes pipeline...")

        # Phase 1: Pre-processing
        logger.info("Phase 1: Parsing and filtering PAF...")
        paf_df = parse_paf(args.paf, args.min_mq)
        primary_paf_df = get_primary_targets(paf_df)
        
        # Phase 2: Initialization
        logger.info("Phase 2: Initializing ContigSummary objects...")
        # GC calculation is parallelized inside parse_fasta
        fasta_data = parse_fasta(args.fasta, args.threads)
        busco_map = parse_busco(args.busco) if args.busco else {}
        
        summary_list = []
        query_to_summary = {}
        for q_id, (gc, length) in fasta_data.items():
            cs = ContigSummary(
                query_id=q_id,
                query_length=length,
                gc_content=gc,
                busco_genes=busco_map.get(q_id, set())
            )
            summary_list.append(cs)
            query_to_summary[q_id] = cs
            
        # Mark aligned contigs
        aligned_query_ids = primary_paf_df['query_id'].unique()
        for q_id in aligned_query_ids:
            if q_id in query_to_summary:
                summary = query_to_summary[q_id]
                summary.status = Status.ALIGNED_RETAINED
                # primary_target_id
                target_id = primary_paf_df[primary_paf_df['query_id'] == q_id]['target_id'].iloc[0]
                summary.target_id = target_id

        # Phase 3: Alignment Filtering and Scoring
        logger.info("Phase 3: Tiling and scoring alignments...")
        calculate_initial_redundancy(summary_list, primary_paf_df)
        
        # Parallelize Tiling
        query_groups = [group for _, group in primary_paf_df.groupby('query_id')]
        with multiprocessing.Pool(args.threads, initializer=worker_configurer, initargs=(log_queue,)) as pool:
            tiling_results = pool.starmap(tile_and_score_contig, [(group, args.min_overlap) for group in query_groups])
        
        for q_id, intervals, score, max_as, tiled_out in tiling_results:
            if q_id in query_to_summary:
                cs = query_to_summary[q_id]
                cs.intervals = intervals
                cs.sum_normalized_score = score
                cs.max_alignment_score = max_as
                cs.tiled_out_count = tiled_out

        # Phase 4: Mash Distance Threshold Calculation
        logger.info("Phase 4: Calculating Mash distance threshold...")
        mash_df = parse_mash(args.mash)
        mash_lookup = build_mash_lookup(mash_df)
        
        overlapping_pairs = get_overlapping_pairs(summary_list, args.min_overlap)
        overlap_distances = []
        for id1, id2 in overlapping_pairs:
            dist = get_mash_distance(mash_lookup, id1, id2)
            if dist is not None:
                overlap_distances.append(dist)
        
        logger.info(f"Overlapping pairs: {len(overlap_distances)}")
        
        if args.distance_threshold is not None:
            dist_threshold = args.distance_threshold
            method = "User-supplied"
        else:
            dist_threshold, method = estimate_distance_threshold(overlap_distances)
            
        logger.info(f"Final distance threshold: {dist_threshold:.4f} (Method: {method})")

        # Phase 5: Redundancy Resolution (Iterative Tournament)
        logger.info("Phase 5: Running iterative tournament...")
        # Group by target_id
        target_groups = {}
        for c in summary_list:
            if c.status == Status.ALIGNED_RETAINED:
                if c.target_id not in target_groups:
                    target_groups[c.target_id] = []
                target_groups[c.target_id].append(c)
        
        # Parallelize by target group
        target_items = list(target_groups.values())
        with multiprocessing.Pool(args.threads, initializer=worker_configurer, initargs=(log_queue,)) as pool:
            tournament_results = pool.starmap(run_tournament_on_target, [
                (group, mash_lookup, dist_threshold, args.min_overlap, args.min_size_safeguard, args.max_tournament_iterations)
                for group in target_items
            ])
        
        for group in tournament_results:
            for updated_cs in group:
                original = query_to_summary[updated_cs.query_id]
                original.status = updated_cs.status
                original.disqualifier = updated_cs.disqualifier
                original.discarded_reason = updated_cs.discarded_reason
                original.retained_reason = updated_cs.retained_reason

        # Phase 6: Unaligned Contig Handling
        logger.info("Phase 6: Screening unaligned contigs...")
        unaligned = [c for c in summary_list if c.status == Status.UNALIGNED_RETAINED]
        # Sort by length descending as requested
        unaligned.sort(key=lambda x: x.query_length, reverse=True)
        
        retained_so_far = [c for c in summary_list if c.status == Status.ALIGNED_RETAINED]
        
        current_retained = list(retained_so_far)
        
        with multiprocessing.Pool(args.threads, initializer=worker_configurer, initargs=(log_queue,)) as pool:
            screened_unaligned = pool.starmap(screen_unaligned_contig, [
                (u, current_retained, mash_lookup, dist_threshold) for u in unaligned
            ])
        
        # Now update summary_list and handle U-U redundancy for those still RETAINED
        unaligned_retained_after_aligned_check = [u for u in screened_unaligned if u.status == Status.UNALIGNED_RETAINED]
        
        final_unaligned_retained = []
        for u in unaligned_retained_after_aligned_check:
            is_redundant = False
            for r in final_unaligned_retained:
                dist = get_mash_distance(mash_lookup, u.query_id, r.query_id)
                if dist is not None and dist < dist_threshold:
                    u.status = Status.UNALIGNED_DISCARDED
                    u.disqualifier = r.query_id
                    u.discarded_reason["Mash_Redundancy"] = True
                    is_redundant = True
                    break
            if not is_redundant:
                final_unaligned_retained.append(u)
                
        # Update the summary list with results from Phase 6
        for u in screened_unaligned:
            original = query_to_summary[u.query_id]
            original.status = u.status
            original.disqualifier = u.disqualifier
            original.discarded_reason = u.discarded_reason

        num_unaligned_total = len(unaligned)
        num_unaligned_retained = len(final_unaligned_retained)
        if num_unaligned_total > 0:
            perc_unaligned_retained = (num_unaligned_retained / num_unaligned_total) * 100
            logger.info(f"Unaligned contigs retained: {num_unaligned_retained} ({perc_unaligned_retained:.2f}%)")
            logger.info(f"Unaligned contigs discarded: {num_unaligned_total - num_unaligned_retained} ({100 - perc_unaligned_retained:.2f}%)")

        # Phase 7: Reporting
        logger.info("Phase 7: Generating reports...")
        stats_initial = calculate_assembly_stats([c.query_length for c in summary_list])
        busco_initial = get_busco_counts(args.busco, set(query_to_summary.keys())) if args.busco else {}
        
        retained_ids = {c.query_id for c in summary_list if c.status in [Status.ALIGNED_RETAINED, Status.UNALIGNED_RETAINED]}
        busco_filtered = get_busco_counts(args.busco, retained_ids) if args.busco else {}
        
        generate_report(
            summary_list,
            overlap_distances,
            dist_threshold,
            stats_initial,
            busco_initial,
            busco_filtered,
            output_dir,
            args.threads
        )
        
        write_filtered_fasta(Path(args.fasta), output_dir / 'filtered_assembly.fasta', retained_ids)
        
        logger.info(f"Pipeline complete. Results saved in {output_dir}")
    except Exception as e:
        logger.error(f"Critical failure: {e}")
        sys.exit(1)
    finally:
        log_listener.stop()

if __name__ == "__main__":
    main()
