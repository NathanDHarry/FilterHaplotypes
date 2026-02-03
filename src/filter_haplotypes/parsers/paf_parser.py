"""
PAF alignment file parser for FilterHaplotypes.
Handles reading PAF records, extracting alignment scores (AS),
and identifying primary target loci for each query contig.
"""

import pandas as pd
import numpy as np
from typing import List, Dict, Any, Optional
import logging

logger = logging.getLogger(__name__)

def parse_paf(paf_path: str, min_mq: int = 10) -> pd.DataFrame:
    """
    Parse a PAF file and filter by mapping quality and extract AS:i tag.

    :param paf_path: Path to the PAF file.
    :param min_mq: Minimum mapping quality threshold.
    :return: A pandas DataFrame containing the filtered PAF records with an 'AS' column.
    """
    # PAF mandatory columns
    columns = [
        'query_id', 'query_len', 'query_start', 'query_end', 'strand',
        'target_id', 'target_len', 'target_start', 'target_end',
        'n_match', 'aln_len', 'mq'
    ]
    
    try:
        # PAF can have varying number of columns due to tags
        df = pd.read_csv(paf_path, sep='\t', header=None, low_memory=False, encoding='utf-8')
    except Exception as e:
        logger.error(f"Failed to read PAF file {paf_path}: {e}")
        raise

    if df.empty:
        logger.warning(f"PAF file {paf_path} is empty.")
        return pd.DataFrame(columns=columns + ['AS'])

    # Map mandatory columns
    num_mandatory = len(columns)
    df_mandatory = df.iloc[:, :num_mandatory].copy()
    df_mandatory.columns = columns
    
    # Extract AS:i tag from remaining columns
    def extract_as_tag(row):
        for val in row[num_mandatory:]:
            if isinstance(val, str) and val.startswith('AS:i:'):
                try:
                    return int(val.split(':')[-1])
                except ValueError:
                    continue
        return None

    df_mandatory['AS'] = df.apply(extract_as_tag, axis=1)
    
    # Validation: Check for AS tag
    if df_mandatory['AS'].isnull().any():
        num_missing = df_mandatory['AS'].isnull().sum()
        logger.warning(f"{num_missing} records missing 'AS:i:' tag in {paf_path}")
    
    # Filter by MQ
    df_filtered = df_mandatory[df_mandatory['mq'] >= min_mq].copy()
    
    # Drop records without AS tag if necessary, or default to 0
    df_filtered['AS'] = df_filtered['AS'].fillna(0).astype(int)
    
    return df_filtered

def get_primary_targets(df: pd.DataFrame) -> pd.DataFrame:
    """
    Identify the primary target locus for each query contig.
    Selection: Highest 90th percentile of AS, then longest alignment length, then alphabetical target_id.

    :param df: Filtered PAF DataFrame.
    :return: DataFrame with only primary target records for each query_id.
    """
    if df.empty:
        return df

    # Calculate 90th percentile of AS per target for each query
    # and also get the max aln_len per target (for tie-breaking)
    target_stats = df.groupby(['query_id', 'target_id']).agg(
        p90_as=('AS', lambda x: np.percentile(x, 90)),
        max_aln_len=('aln_len', 'max')
    ).reset_index()

    # Sort to find the best target per query
    # Primary: 90th percentile AS, Secondary: max aln_len, Tertiary: alphabetical target_id
    target_stats_sorted = target_stats.sort_values(
        by=['query_id', 'p90_as', 'max_aln_len', 'target_id'],
        ascending=[True, False, False, True]
    )
    
    # For each query_id, define primary target
    primary_targets = target_stats_sorted.groupby('query_id').first().reset_index()[['query_id', 'target_id']]
    
    # Filter original df to keep only records aligning to primary targets
    merged = df.merge(primary_targets, on=['query_id', 'target_id'], how='inner')
    
    return merged
