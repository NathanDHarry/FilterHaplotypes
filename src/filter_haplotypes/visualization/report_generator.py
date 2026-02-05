"""
Report generation module for FilterHaplotypes.
Generates the final TSV summary, filtered FASTA, and interactive HTML dashboard.
"""

from Bio import SeqIO
import plotly.graph_objects as go
import multiprocessing
from jinja2 import Environment, FileSystemLoader
from pathlib import Path
from typing import List, Dict, Any, Set
from src.filter_haplotypes.core.models import ContigSummary, Status
from src.filter_haplotypes.utils.stats import calculate_assembly_stats, calculate_l_curve
import pandas as pd

def process_contig_metrics(c: ContigSummary) -> Dict[str, Any]:
    """
    Calculate individual metrics for a contig in parallel.
    
    :param c: ContigSummary object.
    :return: Dictionary of metrics.
    """
    matching_length = sum(end - start for start, end in c.intervals)
    return {
        'query_id': c.query_id,
        'query_length': c.query_length,
        'gc_content': c.gc_content,
        'matching_length': matching_length,
        'status': c.status.value,
        'target_id': c.target_id,
        'sum_normalized_score': c.sum_normalized_score,
        'disqualifier': c.disqualifier,
        'Round1_discard': c.discarded_reason['Round1'],
        'OrphanOverride_discard': c.discarded_reason['OrphanOverride'],
        'Mash_Redundancy_discard': c.discarded_reason['Mash_Redundancy'],
        'Score_retain': c.retained_reason['Score'],
        'Mash_retain': c.retained_reason['Mash'],
        'Size_retain': c.retained_reason['Size'],
        'OrphanRecovery_retain': c.retained_reason['OrphanRecovery'],
        'Unique_retain': c.retained_reason['Unique']
    }

def generate_report(
    summary_list: List[ContigSummary],
    overlap_distances: List[float],
    distance_threshold: float,
    stats_initial: Dict[str, Any],
    busco_initial: Dict[str, Any],
    busco_filtered: Dict[str, Any],
    output_dir: Path,
    threads: int = 1,
    run_parameters: Dict[str, Any] = None
):
    """
    Phase 7 Step 10: Generate all output files and the interactive HTML report.
    
    :param summary_list: Final list of ContigSummary objects.
    :param overlap_distances: List of Mash distances used for thresholding.
    :param distance_threshold: Final Mash distance threshold.
    :param stats_initial: Initial assembly stats.
    :param busco_initial: Initial BUSCO stats.
    :param busco_filtered: Filtered BUSCO stats.
    :param output_dir: Directory to save outputs.
    :param threads: Number of threads for parallel metric calculation.
    :param run_parameters: Dictionary of configurable parameters used for the run.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. Calculate individual contig metrics in parallel (Phase 7 Requirement)
    with multiprocessing.Pool(processes=threads) as pool:
        summary_data = pool.map(process_contig_metrics, summary_list)
    
    df_summary = pd.DataFrame(summary_data)
    df_summary.to_csv(output_dir / 'summary_report.tsv', sep='\t', index=False, encoding='utf-8')
    
    # 2. Assembly Stats for Filtered
    retained_lengths = [c.query_length for c in summary_list if c.status in [Status.ALIGNED_RETAINED, Status.UNALIGNED_RETAINED]]
    stats_filtered = calculate_assembly_stats(retained_lengths)
    
    # 3. Filtering Decision Logic Sums
    retained_counts = {
        'Score': sum(1 for c in summary_list if c.retained_reason['Score']),
        'Mash': sum(1 for c in summary_list if c.retained_reason['Mash']),
        'Size': sum(1 for c in summary_list if c.retained_reason['Size']),
        'OrphanRecovery': sum(1 for c in summary_list if c.retained_reason['OrphanRecovery']),
        'Unique': sum(1 for c in summary_list if c.retained_reason['Unique'])
    }
    discarded_counts = {
        'Round1': sum(1 for c in summary_list if c.discarded_reason['Round1']),
        'OrphanOverride': sum(1 for c in summary_list if c.discarded_reason['OrphanOverride']),
        'Mash_Redundancy': sum(1 for c in summary_list if c.discarded_reason['Mash_Redundancy'])
    }
    
    # 4. Plots
    # Mash Plot
    fig_mash = go.Figure()
    if overlap_distances:
        fig_mash.add_trace(go.Histogram(x=list(overlap_distances), name='Overlap Distances', nbinsx=50))
        fig_mash.add_vline(x=distance_threshold, line_width=3, line_dash="dash", line_color="red", annotation_text=f"Threshold: {distance_threshold:.4f}")
    fig_mash.update_layout(title="Mash Distances Distribution", xaxis_title="Distance", yaxis_title="Count")
    mash_plot_json = fig_mash.to_json()
    
    # L-Curve Plot
    initial_lengths = [c.query_length for c in summary_list]
    xi, yi = calculate_l_curve(initial_lengths)
    xf, yf = calculate_l_curve(retained_lengths)
    
    fig_l = go.Figure()
    fig_l.add_trace(go.Scatter(x=list(xi), y=list(yi), mode='lines', name='Initial Assembly'))
    fig_l.add_trace(go.Scatter(x=list(xf), y=list(yf), mode='lines', name='Filtered Assembly'))
    fig_l.update_layout(title="L-Curve (Cumulative Assembly Size)", xaxis_title="Contig Rank", yaxis_title="Cumulative Bases")
    l_curve_json = fig_l.to_json()
    
    # GC BlobPlot
    fig_gc = go.Figure()
    # Define colors for statuses for consistency
    status_colors = {
        Status.ALIGNED_RETAINED: 'blue',
        Status.UNALIGNED_RETAINED: 'green',
        Status.ALIGNED_DISCARDED: 'red',
        Status.UNALIGNED_DISCARDED: 'orange'
    }
    
    for status in Status:
        status_contigs = [c for c in summary_list if c.status == status]
        if not status_contigs:
            continue
        fig_gc.add_trace(go.Scatter(
            x=[c.query_length for c in status_contigs],
            y=[c.gc_content for c in status_contigs],
            mode='markers',
            name=status.value,
            marker=dict(color=status_colors.get(status, 'gray')),
            text=[c.query_id for c in status_contigs],
            hoverinfo='text+x+y'
        ))
    
    fig_gc.update_layout(
        title="GC-content vs. Contig Length",
        xaxis_type="log",
        xaxis_title="Contig Length (bp)",
        yaxis_title="GC-content (%)"
    )
    gc_plot_json = fig_gc.to_json()
    
    # 5. Render HTML
    template_dir = Path(__file__).parent / 'templates'
    env = Environment(loader=FileSystemLoader(str(template_dir)))
    template = env.get_template('report.html')
    
    html_content = template.render(
        stats_initial=stats_initial,
        stats_filtered=stats_filtered,
        busco_initial=busco_initial,
        busco_filtered=busco_filtered,
        retained_counts=retained_counts,
        discarded_counts=discarded_counts,
        mash_plot_json=mash_plot_json,
        l_curve_json=l_curve_json,
        gc_plot_json=gc_plot_json,
        distance_threshold=distance_threshold,
        run_parameters=run_parameters if run_parameters else {}
    )
    
    with open(output_dir / 'report.html', 'w', encoding='utf-8') as f:
        f.write(html_content)

def write_filtered_fasta(
    input_fasta: Path,
    output_fasta: Path,
    retained_ids: Set[str]
):
    """
    Write sequences of retained contigs to a new FASTA file.
    
    :param input_fasta: Path to the original FASTA.
    :param output_fasta: Path to the output FASTA.
    :param retained_ids: Set of query_ids to retain.
    """
    retained_records = []
    for record in SeqIO.parse(input_fasta, "fasta"):
        if record.id in retained_ids:
            retained_records.append(record)
            
    with open(output_fasta, "w", encoding='utf-8') as f:
        SeqIO.write(retained_records, f, "fasta")
