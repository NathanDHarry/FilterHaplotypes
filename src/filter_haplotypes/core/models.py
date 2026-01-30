"""
Data models for FilterHaplotypes.
Defines the ContigSummary class and filtering Status enum.
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import List, Tuple, Set, Optional, Dict

class Status(Enum):
    """
    Enum representing the current filtering status of a contig.
    """
    ALIGNED_RETAINED = "ALIGNED_RETAINED"
    ALIGNED_DISCARDED = "ALIGNED_DISCARDED"
    UNALIGNED_RETAINED = "UNALIGNED_RETAINED"
    UNALIGNED_DISCARDED = "UNALIGNED_DISCARDED"

@dataclass
class ContigSummary:
    """
    Data class representing a summary of a contig, its alignments, and its filtering status.
    """
    # General Attributes
    query_id: str
    query_length: int
    gc_content: float = 0.0
    busco_genes: Set[str] = field(default_factory=set)
    status: Status = Status.UNALIGNED_RETAINED
    discarded_reason: Dict[str, bool] = field(default_factory=lambda: {
        "Round1": False,
        "OrphanOverride": False,
        "Mash_Redundancy": False
    })
    retained_reason: Dict[str, bool] = field(default_factory=lambda: {
        "Score": False,
        "Mash": False,
        "Size": False,
        "OrphanRecovery": False,
        "Unique": False
    })
    disqualifier: Optional[str] = None

    # Alignment-Specific Attributes
    target_id: Optional[str] = None
    intervals: List[Tuple[int, int]] = field(default_factory=list)
    sum_normalized_score: float = 0.0
    max_alignment_score: int = 0
    initial_overlapping_bases: int = 0
    tiled_out_count: int = 0  # To store the count of discarded (tiled-out) alignment events
