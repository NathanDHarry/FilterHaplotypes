"""
BUSCO table parser for FilterHaplotypes.
Handles reading full_table.tsv and assessing assembly completeness.
"""

import pandas as pd
from typing import Dict, Set
import logging

logger = logging.getLogger(__name__)

def parse_busco(busco_path: str) -> Dict[str, Set[str]]:
    """
    Parse a BUSCO full_table.tsv file and extract complete BUSCO IDs for each contig.

    :param busco_path: Path to the BUSCO table.
    :return: A dictionary mapping sequence/contig ID to a set of complete BUSCO IDs.
    """
    try:
        # BUSCO table often has many comment lines starting with #
        df = pd.read_csv(busco_path, sep='\t', comment='#', header=None, encoding='utf-8')
    except Exception as e:
        logger.error(f"Failed to read BUSCO file {busco_path}: {e}")
        return {}

    if df.empty:
        logger.warning(f"BUSCO file {busco_path} is empty or only contains comments.")
        return {}

    # Map standard BUSCO columns
    # 0: Busco id, 1: Status, 2: Sequence
    columns = {0: 'busco_id', 1: 'status', 2: 'sequence'}
    df = df.rename(columns=columns)
    
    # Filter for 'Complete' or 'Duplicated'
    complete_statuses = ['Complete', 'Duplicated']
    df_complete = df[df['status'].isin(complete_statuses)]
    
    # Group by sequence and collect busco_ids
    busco_map = df_complete.groupby('sequence')['busco_id'].apply(set).to_dict()
    
    return busco_map

def get_busco_counts(busco_path: str, retained_contigs: Set[str]) -> Dict[str, int]:
    """
    Calculate BUSCO completeness for a set of retained contigs.

    :param busco_path: Path to the BUSCO table.
    :param retained_contigs: Set of retained contig IDs.
    :return: Dictionary with 'complete_single' and 'duplicated' counts.
    """
    try:
        df = pd.read_csv(busco_path, sep='\t', comment='#', header=None, encoding='utf-8')
    except Exception as e:
        logger.error(f"Failed to read BUSCO file {busco_path}: {e}")
        return {'complete_single': 0, 'duplicated': 0}

    # 0: Busco id, 1: Status, 2: Sequence
    df = df.rename(columns={0: 'busco_id', 1: 'status', 2: 'sequence'})
    
    # Filter for complete/duplicated genes found in retained contigs
    df_retained = df[df['sequence'].isin(retained_contigs)]
    
    # BUSCO logic: 
    # A gene is 'Complete' if it's found at least once.
    # It's 'Duplicated' if it's found more than once across all retained contigs.
    
    # Count occurrences of each BUSCO ID in retained contigs
    counts = df_retained[df_retained['status'].isin(['Complete', 'Duplicated'])]['busco_id'].value_counts()
    
    complete_single = (counts == 1).sum()
    duplicated = (counts > 1).sum()
    
    return {
        'complete_single': int(complete_single),
        'duplicated': int(duplicated)
    }
