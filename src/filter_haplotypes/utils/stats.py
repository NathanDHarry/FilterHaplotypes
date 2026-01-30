"""
Assembly statistics calculation utilities.
Includes N50 and L-curve data generation.
"""

import numpy as np
from typing import List, Dict, Any, Tuple

def calculate_assembly_stats(lengths: List[int]) -> Dict[str, int]:
    """
    Calculate assembly statistics (N50, N60, N70, N80, N90, N100) and total bases.
    
    :param lengths: List of contig lengths.
    :return: Dictionary with stats.
    """
    if not lengths:
        return {f"N{i}": 0 for i in range(50, 110, 10)} | {"Total Bases": 0, "Num Contigs": 0}
        
    lengths_sorted = sorted(lengths, reverse=True)
    total_bases = sum(lengths_sorted)
    num_contigs = len(lengths_sorted)
    
    stats = {
        "Total Bases": total_bases,
        "Num Contigs": num_contigs
    }
    
    cumulative_sum = 0
    nx_targets = {i: total_bases * (i / 100.0) for i in range(50, 110, 10)}
    nx_values = {}
    nx_counts = {}
    
    current_nx = 50
    for i, length in enumerate(lengths_sorted):
        cumulative_sum += length
        while current_nx <= 100 and cumulative_sum >= nx_targets[current_nx]:
            nx_values[f"N{current_nx}"] = length
            nx_counts[f"N{current_nx}_count"] = i + 1
            current_nx += 10
            
    for i in range(50, 110, 10):
        stats[f"N{i}"] = nx_values.get(f"N{i}", 0)
        stats[f"N{i}_count"] = nx_counts.get(f"N{i}_count", 0)
        
    return stats

def calculate_l_curve(lengths: List[int]) -> Tuple[List[int], List[int]]:
    """
    Calculate data for L-curve (Cumulative assembly size).
    
    :param lengths: List of contig lengths.
    :return: Tuple of (x_counts, y_cumulative_bases).
    """
    lengths_sorted = sorted(lengths, reverse=True)
    y = np.cumsum(lengths_sorted).tolist()
    x = list(range(1, len(lengths_sorted) + 1))
    return x, y
