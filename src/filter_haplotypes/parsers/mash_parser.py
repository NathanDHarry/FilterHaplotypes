"""
Mash distance file parser for FilterHaplotypes.
Handles reading Mash TSV output and building high-performance lookup structures.
"""

import pandas as pd
from typing import Dict, Tuple, Optional
import logging

logger = logging.getLogger(__name__)

def parse_mash(mash_path: str) -> pd.DataFrame:
    """
    Parse a Mash distance file (TSV) and filter by P-value.

    :param mash_path: Path to the Mash distance file.
    :return: A pandas DataFrame containing Mash records with P-value < 0.05.
    """
    columns = ['id1', 'id2', 'distance', 'p_value', 'hashes']
    try:
        df = pd.read_csv(mash_path, sep='\t', header=None, names=columns, encoding='utf-8')
    except Exception as e:
        logger.error(f"Failed to read Mash file {mash_path}: {e}")
        raise

    # Filter by P-value < 0.05
    initial_count = len(df)
    df_filtered = df[df['p_value'] < 0.05].copy()
    filtered_count = len(df_filtered)
    
    if initial_count > filtered_count:
        logger.info(f"Filtered out {initial_count - filtered_count} Mash records with P-value >= 0.05")
    
    return df_filtered

def build_mash_lookup(df: pd.DataFrame) -> Dict[str, Dict[str, float]]:
    """
    Build a high-performance nested dictionary for Mash distance lookups.
    Stores distance for both (id1, id2) and (id2, id1).

    :param df: Filtered Mash DataFrame.
    :return: A nested dictionary {id1: {id2: distance}}.
    """
    lookup = {}
    for row in df.itertuples(index=False):
        id1, id2, dist = row.id1, row.id2, row.distance
        
        if id1 not in lookup:
            lookup[id1] = {}
        lookup[id1][id2] = dist
        
        if id2 not in lookup:
            lookup[id2] = {}
        lookup[id2][id1] = dist
        
    return lookup

def get_mash_distance(lookup: Dict[str, Dict[str, float]], id1: str, id2: str) -> Optional[float]:
    """
    Retrieve Mash distance between two IDs from the lookup dictionary.

    :param lookup: The nested dictionary {id1: {id2: distance}}.
    :param id1: Query ID 1.
    :param id2: Query ID 2.
    :return: Mash distance or None if not found.
    """
    if id1 == id2:
        return 0.0
    return lookup.get(id1, {}).get(id2)
