"""
Core filtering logic for FilterHaplotypes.
Includes alignment tiling, score finalization, Mash distance threshold estimation,
and the iterative tournament algorithm for redundancy resolution.
"""

import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
from typing import List, Dict, Tuple, Any, Optional
from src.filter_haplotypes.core.models import ContigSummary, Status
import logging

logger = logging.getLogger(__name__)

def calculate_initial_redundancy(summary_list: List[ContigSummary], paf_df: pd.DataFrame):
    """
    Phase 3 Step 3: Calculate initial redundancy for aligned contigs.
    
    :param summary_list: List of ContigSummary objects.
    :param paf_df: Filtered PAF DataFrame (primary targets only).
    """
    total_contigs_with_overlaps = 0
    total_overlapping_bases = 0
    total_bases_in_assembly = sum(c.query_length for c in summary_list)

    aligned_summaries = {c.query_id: c for c in summary_list if c.status == Status.ALIGNED_RETAINED}
    
    for query_id, group in paf_df.groupby('query_id'):
        if query_id not in aligned_summaries:
            continue
            
        summary = aligned_summaries[query_id]
        intervals = sorted(zip(group['target_start'], group['target_end']))
        
        if not intervals:
            continue
            
        # Calculate overlapping bases
        overlapping_bases = 0
        current_max = -1
        
        # Simple interval overlap calculation for total covered bases
        # and then subtract from total sum of interval lengths to get redundant bases?
        # "total number of bases covered by more than one alignment"
        
        # To find bases covered by >1 alignment:
        # We can use a coordinate count approach or sorted intervals
        
        events = []
        for start, end in intervals:
            events.append((start, 1))
            events.append((end, -1))
        
        events.sort()
        
        coverage = 0
        last_pos = events[0][0]
        redundant_bases = 0
        
        for pos, delta in events:
            if coverage > 1:
                redundant_bases += (pos - last_pos)
            coverage += delta
            last_pos = pos
            
        summary.initial_overlapping_bases = redundant_bases
        if redundant_bases > 0:
            total_contigs_with_overlaps += 1
            total_overlapping_bases += redundant_bases

    # Logging stats about redundant alignments at [DEBUG] level
    num_aligned = len(aligned_summaries)
    if num_aligned > 0:
        perc_contigs = (total_contigs_with_overlaps / len(summary_list)) * 100 if summary_list else 0
        perc_bases = (total_overlapping_bases / total_bases_in_assembly) * 100 if total_bases_in_assembly > 0 else 0
        
        logger.debug(f"Query contigs containing alignments with overlaps: {total_contigs_with_overlaps}")
        logger.debug(f"Percentage of contigs in the assembly they represent: {perc_contigs:.2f}%")
        logger.debug(f"Total percentage of bases covered by multiple alignments out of all bases in the assembly: {perc_bases:.2f}%")

def tile_and_score_contig(query_group: pd.DataFrame, overlap_tolerance: int = 10) -> Tuple[str, List[Tuple[int, int]], float, int, int]:
    """
    Phase 3 Step 4 & 5: Tiling algorithm and score finalization for a single contig.
    
    :param query_group: DataFrame containing PAF records for a single query_id (on its primary target).
    :param overlap_tolerance: Maximum allowed overlap on target (bp).
    :return: Tuple (query_id, intervals, sum_normalized_score, max_alignment_score, tiled_out_count).
    """
    query_id = query_group['query_id'].iloc[0]
    query_length = query_group['query_len'].iloc[0]
    
    # Sorting: AS descending, then alignment length descending
    sorted_group = query_group.sort_values(by=['AS', 'aln_len'], ascending=[False, False])
    
    accepted_intervals: List[Tuple[int, int]] = []
    total_as = 0
    max_as = 0
    tiled_out_count = 0
    
    for row in sorted_group.itertuples(index=False):
        start, end, score = row.target_start, row.target_end, row.AS
        
        # Check overlap with already accepted intervals
        is_overlapping = False
        for a_start, a_end in accepted_intervals:
            # Overlap length = min(end1, end2) - max(start1, start2)
            overlap = min(end, a_end) - max(start, a_start)
            if overlap > overlap_tolerance:
                is_overlapping = True
                break
        
        if not is_overlapping:
            accepted_intervals.append((start, end))
            total_as += score
            if score > max_as:
                max_as = score
        else:
            tiled_out_count += 1
            
    sum_normalized_score = total_as / query_length if query_length > 0 else 0.0
    
    return query_id, accepted_intervals, sum_normalized_score, max_as, tiled_out_count

def estimate_distance_threshold(overlap_distances: List[float]) -> Tuple[float, str]:
    """
    Phase 4 Step 6: Estimate Mash distance threshold based on distribution of overlapping pairs.
    
    :param overlap_distances: List of Mash distances for overlapping contig pairs.
    :return: A tuple (distance_threshold, method_used).
    """
    if not overlap_distances or len(overlap_distances) < 1000:
        return 0.05, "Default (Insufficient pairs < 1000)"
    
    data = np.array(overlap_distances)
    mean_dist = np.mean(data)
    
    # Kernel Density Estimation
    try:
        kde = gaussian_kde(data, bw_method='scott')
        x_eval = np.linspace(0, 0.2, 500)
        kde_vals = kde.evaluate(x_eval)
        
        # Find local minima (valleys)
        minima_indices = (np.diff(np.sign(np.diff(kde_vals))) > 0).nonzero()[0] + 1
        
        if len(minima_indices) > 0:
            # Check if bimodal by checking for at least two peaks separated by a valley
            # A simple heuristic: if there's a local minimum in [0, 0.2] it might be bimodal
            valley_idx = minima_indices[0]
            valley_x = x_eval[valley_idx]
            return float(valley_x), "KDE Valley (Bimodal)"
            
    except Exception as e:
        logger.warning(f"KDE estimation failed: {e}")

    if mean_dist <= 0.1:
        threshold = np.percentile(data, 95)
        return float(threshold), "95th Percentile (Unimodal Mean <= 0.1)"
    else:
        return 0.05, "Default (Unimodal Mean > 0.1)"

def get_overlapping_pairs(summary_list: List[ContigSummary], min_overlap: int = 1) -> List[Tuple[str, str]]:
    """
    Identify pairs of contigs with overlapping intervals on the same target.
    
    :param summary_list: List of ContigSummary objects.
    :param min_overlap: Minimum overlap bases.
    :return: List of tuples (query_id1, query_id2).
    """
    aligned = [c for c in summary_list if c.status == Status.ALIGNED_RETAINED]
    target_groups = {}
    for c in aligned:
        if c.target_id not in target_groups:
            target_groups[c.target_id] = []
        target_groups[c.target_id].append(c)
        
    overlapping_pairs = []
    for target_id, contigs in target_groups.items():
        for i in range(len(contigs)):
            c1 = contigs[i]
            for j in range(i + 1, len(contigs)):
                c2 = contigs[j]
                
                # Check for any overlapping intervals
                has_overlap = False
                for s1, e1 in c1.intervals:
                    for s2, e2 in c2.intervals:
                        overlap = min(e1, e2) - max(s1, s2)
                        if overlap >= min_overlap:
                            has_overlap = True
                            break
                    if has_overlap:
                        break
                
                if has_overlap:
                    overlapping_pairs.append((c1.query_id, c2.query_id))
                    
    return overlapping_pairs

def run_tournament_on_target(
    contigs: List[ContigSummary],
    mash_lookup: Dict[str, Dict[str, float]],
    distance_threshold: float,
    min_overlap: int = 1,
    min_size_safeguard: float = 0.5,
    max_tournament_iterations: int = 100000
) -> List[ContigSummary]:
    """
    Phase 5 Step 7 & 8: Perform tournament for a group of contigs aligned to the same target.
    
    :param contigs: List of ContigSummary objects for a specific target_id.
    :param mash_lookup: Nested dictionary for Mash distance lookups.
    :param distance_threshold: Mash distance threshold.
    :param min_overlap: Minimum overlap bases.
    :param min_size_safeguard: Ratio for size-based retention.
    :param max_tournament_iterations: Maximum iterations for the tournament loop.
    :return: Updated list of ContigSummary objects.
    """
    if not contigs:
        return []

    # Helper to check competition rules
    def check_competition(C: ContigSummary, O: ContigSummary) -> bool:
        # 1. Shared Target (assumed by grouping)
        # 2. Overlap Condition
        has_overlap = False
        for sC, eC in C.intervals:
            for sO, eO in O.intervals:
                if min(eC, eO) - max(sC, sO) >= min_overlap:
                    has_overlap = True
                    break
            if has_overlap: break
        if not has_overlap: return False

        # 3. Active Status (O must be RETAINED)
        if O.status != Status.ALIGNED_RETAINED: return False

        # 4. Superior Score
        if O.sum_normalized_score > C.sum_normalized_score:
            pass
        elif O.sum_normalized_score == C.sum_normalized_score:
            # Tie-break: earlier in sorted list wins. 
            # We will sort 'contigs' at the start of tournament.
            # So if O appears before C in 'contigs', O wins.
            idx_C = contigs.index(C)
            idx_O = contigs.index(O)
            if idx_O < idx_C:
                logger.warning(f"Score tie between {C.query_id} and {O.query_id}. {O.query_id} wins by position.")
            else:
                return False
        else:
            return False

        # 5. Similarity
        dist = mash_lookup.get(C.query_id, {}).get(O.query_id)
        if dist is None or dist >= distance_threshold:
            return False

        # 6. Size Safeguard
        if O.query_length < min_size_safeguard * C.query_length:
            return False

        return True

    # Sorting for Pass 1: target_id (already grouped), then min start coordinate
    contigs.sort(key=lambda x: min(i[0] for i in x.intervals) if x.intervals else 0)

    # Pass 1: Initial Sweep
    # This pass identifies contigs that are immediately redundant against better-scoring contigs.
    for C in contigs:
        if C.status != Status.ALIGNED_RETAINED: continue
        
        disqualifier_found = False
        for O in contigs:
            if O.query_id == C.query_id: continue
            # If O is superior to C, mark C for discard
            if check_competition(C, O):
                C.status = Status.ALIGNED_DISCARDED
                C.disqualifier = O.query_id
                C.discarded_reason["Round1"] = True
                disqualifier_found = True
                break
        
        if not disqualifier_found:
            # If not discarded, determine the specific reasons why it was kept
            any_overlap = False
            for O in contigs:
                if O.query_id == C.query_id: continue
                # Identify if there's any overlap with other contigs on this target
                ovl = False
                for sC, eC in C.intervals:
                    for sO, eO in O.intervals:
                        if min(eC, eO) - max(sC, sO) >= min_overlap:
                            ovl = True; break
                    if ovl: break
                if not ovl: continue
                
                any_overlap = True
                # If it overlaps with a retained contig but isn't discarded, record why
                if O.status == Status.ALIGNED_RETAINED:
                    if C.sum_normalized_score > O.sum_normalized_score:
                        C.retained_reason["Score"] = True
                    
                    dist = mash_lookup.get(C.query_id, {}).get(O.query_id)
                    # Divergent enough to be kept despite lower score?
                    if dist is not None and dist > distance_threshold and C.sum_normalized_score < O.sum_normalized_score:
                        C.retained_reason["Mash"] = True
                    
                    # Too large to be discarded by a much smaller contig?
                    if O.query_length < min_size_safeguard * C.query_length and C.sum_normalized_score < O.sum_normalized_score:
                        C.retained_reason["Size"] = True
            
            # If no overlaps at all, it's unique to this locus
            if not any_overlap:
                C.retained_reason["Unique"] = True

    # Pass 2+: Iterative Loop
    # This loop handles 'orphans' and ensures stability across status changes.
    status_changed = True
    iterations = 0
    while status_changed and iterations < max_tournament_iterations:
        status_changed = False
        iterations += 1
        
        # Identify Orphans: Discarded contigs whose disqualifier was itself discarded later
        orphans = []
        id_to_contig = {c.query_id: c for c in contigs}
        for c in contigs:
            if c.status == Status.ALIGNED_DISCARDED:
                dq = id_to_contig.get(c.disqualifier)
                if dq and dq.status == Status.ALIGNED_DISCARDED:
                    orphans.append(c)
        
        # Sort orphans by genomic position to maintain deterministic behavior
        orphans.sort(key=lambda x: min(i[0] for i in x.intervals) if x.intervals else 0)
        
        for C in orphans:
            # The Challenge: See if the orphan can be promoted or if it's still disqualified
            
            # 1. Challenge existing winners
            for R in contigs:
                if R.status != Status.ALIGNED_RETAINED: continue
                # Does orphan C discard current winner R?
                old_status = C.status
                C.status = Status.ALIGNED_RETAINED # Temporarily treat as active for competition check
                if check_competition(R, C):
                    R.status = Status.ALIGNED_DISCARDED
                    R.disqualifier = C.query_id
                    R.discarded_reason["OrphanOverride"] = True
                    # Reset retained reasons as it's now discarded
                    for k in ["OrphanRecovery", "Score", "Mash", "Size"]:
                        R.retained_reason[k] = False
                    status_changed = True
                C.status = old_status

            # 2. Check if C is still disqualified by any remaining winner R
            disqualified_by_R = None
            for R in contigs:
                if R.status != Status.ALIGNED_RETAINED: continue
                if check_competition(C, R):
                    disqualified_by_R = R
                    break
            
            if disqualified_by_R:
                # Still disqualified, but potentially by a new winner
                C.disqualifier = disqualified_by_R.query_id
                C.discarded_reason["OrphanOverride"] = True
                for k in ["OrphanRecovery", "Score", "Mash", "Size"]:
                    C.retained_reason[k] = False
            else:
                # Promote C to RETAINED as no active winner disqualifies it
                C.status = Status.ALIGNED_RETAINED
                C.disqualifier = None
                C.retained_reason["OrphanRecovery"] = True
                
                # Update competitive retained reasons
                any_overlap = False
                for R in contigs:
                    if R.query_id == C.query_id: continue
                    ovl = False
                    for sC, eC in C.intervals:
                        for sR, eR in R.intervals:
                            if min(eC, eR) - max(sC, sR) >= min_overlap:
                                ovl = True; break
                        if ovl: break
                    if not ovl: continue
                    
                    any_overlap = True
                    if R.status == Status.ALIGNED_RETAINED:
                        if C.sum_normalized_score > R.sum_normalized_score:
                            C.retained_reason["Score"] = True
                        dist = mash_lookup.get(C.query_id, {}).get(R.query_id)
                        if dist is not None and dist > distance_threshold and C.sum_normalized_score < R.sum_normalized_score:
                            C.retained_reason["Mash"] = True
                        if R.query_length < min_size_safeguard * C.query_length and C.sum_normalized_score < R.sum_normalized_score:
                            C.retained_reason["Size"] = True
                
                if not any_overlap:
                    C.retained_reason["Unique"] = True
                    
                status_changed = True

    target_id = contigs[0].target_id if contigs else "unknown"
    if iterations >= max_tournament_iterations:
        logger.warning(f"Tournament limit reached ({max_tournament_iterations}) for target {target_id}.")
    else:
        logger.debug(f"Tournament for {target_id} stabilized after {iterations} iterations.")

    return contigs

def screen_unaligned_contig(
    unaligned_contig: ContigSummary,
    retained_contigs: List[ContigSummary],
    mash_lookup: Dict[str, Dict[str, float]],
    distance_threshold: float
) -> ContigSummary:
    """
    Phase 6 Step 9: Identify if an unaligned contig is redundant against retained contigs.
    
    :param unaligned_contig: The unaligned ContigSummary to check.
    :param retained_contigs: List of currently RETAINED contigs (ALIGNED or UNALIGNED).
    :param mash_lookup: Nested dictionary for Mash distance lookups.
    :param distance_threshold: Mash distance threshold.
    :return: Updated unaligned_contig.
    """
    u_id = unaligned_contig.query_id
    
    for c in retained_contigs:
        if c.query_id == u_id: continue
        
        dist = mash_lookup.get(u_id, {}).get(c.query_id)
        if dist is not None and dist < distance_threshold:
            unaligned_contig.status = Status.UNALIGNED_DISCARDED
            unaligned_contig.disqualifier = c.query_id
            unaligned_contig.discarded_reason["Mash_Redundancy"] = True
            break
            
    return unaligned_contig
