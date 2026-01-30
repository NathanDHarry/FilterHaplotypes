"""
FASTA file parser for FilterHaplotypes.
Handles sequence reading and parallel GC content calculation.
"""

from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import multiprocessing
from typing import Dict, Tuple

def calculate_gc(sequence_record) -> Tuple[str, float, int]:
    """
    Calculate GC content and length for a single sequence record.

    :param sequence_record: A Biopython SeqRecord object.
    :return: A tuple containing query_id, gc_content (percentage), and length.
    """
    gc = gc_fraction(sequence_record.seq) * 100
    return str(sequence_record.id), float(gc), len(sequence_record.seq)

def parse_fasta(fasta_path: str, threads: int = 1) -> Dict[str, Tuple[float, int]]:
    """
    Parse a FASTA file and calculate GC content and length for each contig in parallel.

    :param fasta_path: Path to the FASTA file.
    :param threads: Number of threads to use for parallel processing.
    :return: A dictionary mapping query_id to a tuple of (gc_content, length).
    """
    records = list(SeqIO.parse(fasta_path, "fasta"))
    
    with multiprocessing.Pool(processes=threads) as pool:
        results = pool.map(calculate_gc, records)
    
    return {query_id: (gc, length) for query_id, gc, length in results}
