import pytest
import pandas as pd
from src.filter_haplotypes.core.models import ContigSummary, Status
from src.filter_haplotypes.core.filtering import (
    tile_and_score_contig,
    estimate_distance_threshold,
    run_tournament_on_target,
    screen_unaligned_contig,
    count_unique_single_copy_orthologs,
    get_adjusted_score
)

def test_tile_and_score_contig():
    # Create a mock query group
    data = {
        'query_id': ['Q1', 'Q1', 'Q1'],
        'query_len': [1000, 1000, 1000],
        'target_start': [100, 500, 150],
        'target_end': [300, 700, 350],
        'AS': [200, 300, 250],
        'aln_len': [200, 200, 200]
    }
    df = pd.DataFrame(data)
    
    # Sorted by AS: [300 (500-700), 250 (150-350), 200 (100-300)]
    # 300 accepted.
    # 250 overlaps with 300? No. (150-350 and 500-700)
    # 250 accepted.
    # 200 overlaps with 250 (100-300 vs 150-350). Overlap = 150. > 10.
    # 200 discarded.
    
    q_id, intervals, score, max_as, tiled_out = tile_and_score_contig(df, overlap_tolerance=10)
    
    assert q_id == 'Q1'
    assert len(intervals) == 2
    assert (500, 700) in intervals
    assert (150, 350) in intervals
    assert max_as == 300
    assert score == (300 + 250) / 1000
    assert tiled_out == 1

def test_estimate_distance_threshold():
    # Case < 1000 pairs
    dist, method = estimate_distance_threshold([0.01] * 500)
    assert dist == 0.05
    assert "Insufficient" in method
    
    # Case bimodal (mocking valley)
    # This is harder to mock without many points, but let's try a simple set
    dists = [0.01]*600 + [0.04]*100 + [0.15]*600
    dist, method = estimate_distance_threshold(dists)
    assert dist > 0
    # KDE might not find valley if too few points or distribution is weird
    # but it should return something in [0, 0.2] or fallback
    assert 0 <= dist <= 0.2

def test_run_tournament_on_target():
    c1 = ContigSummary(query_id='C1', query_length=1000, target_id='T1', 
                       intervals=[(100, 500)], sum_normalized_score=0.8, status=Status.ALIGNED_RETAINED)
    c2 = ContigSummary(query_id='C2', query_length=1000, target_id='T1', 
                       intervals=[(200, 600)], sum_normalized_score=0.9, status=Status.ALIGNED_RETAINED)
    
    mash_lookup = {'C1': {'C2': 0.01}, 'C2': {'C1': 0.01}}
    dist_threshold = 0.05
    
    # C2 should discard C1
    results = run_tournament_on_target([c1, c2], mash_lookup, dist_threshold, all_contigs=[c1, c2])
    
    c1_res = next(c for c in results if c.query_id == 'C1')
    c2_res = next(c for c in results if c.query_id == 'C2')
    
    assert c2_res.status == Status.ALIGNED_RETAINED
    assert c1_res.status == Status.ALIGNED_DISCARDED
    assert c1_res.disqualifier == 'C2'
    assert c1_res.discarded_reason['Round1'] is True

def test_screen_unaligned_contig():
    u = ContigSummary(query_id='U1', query_length=1000, status=Status.UNALIGNED_RETAINED)
    r = ContigSummary(query_id='R1', query_length=1000, status=Status.ALIGNED_RETAINED)
    
    mash_lookup = {'U1': {'R1': 0.01}, 'R1': {'U1': 0.01}}
    dist_threshold = 0.05
    
    res = screen_unaligned_contig(u, [r], mash_lookup, dist_threshold)
    assert res.status == Status.UNALIGNED_DISCARDED
    assert res.disqualifier == 'R1'
    assert res.discarded_reason['Mash_Redundancy'] is True

def test_count_unique_single_copy_orthologs():
    # Test basic counting of unique orthologs
    c1 = ContigSummary(query_id='C1', query_length=1000, status=Status.ALIGNED_RETAINED,
                       busco_genes={'BUSCO1', 'BUSCO2', 'BUSCO3'})
    c2 = ContigSummary(query_id='C2', query_length=1000, status=Status.ALIGNED_RETAINED,
                       busco_genes={'BUSCO2', 'BUSCO4'})
    c3 = ContigSummary(query_id='C3', query_length=1000, status=Status.ALIGNED_RETAINED,
                       busco_genes={'BUSCO5'})
    
    all_contigs = [c1, c2, c3]
    
    # Count unique orthologs for C1, excluding C1 and C2
    # C1 has {BUSCO1, BUSCO2, BUSCO3}
    # Other contigs (only C3) have {BUSCO5}
    # So unique orthologs in C1 not in C3: BUSCO1, BUSCO2, BUSCO3 = 3
    exclude_ids = {'C1', 'C2'}
    count = count_unique_single_copy_orthologs(c1, all_contigs, exclude_ids)
    assert count == 3
    
    # Count unique orthologs for C2, excluding C1 and C2
    # C2 has {BUSCO2, BUSCO4}
    # Other contigs (only C3) have {BUSCO5}
    # So unique orthologs in C2 not in C3: BUSCO2, BUSCO4 = 2
    count = count_unique_single_copy_orthologs(c2, all_contigs, exclude_ids)
    assert count == 2
    
    # Test with no BUSCO genes
    c_empty = ContigSummary(query_id='C_empty', query_length=1000, status=Status.ALIGNED_RETAINED)
    count = count_unique_single_copy_orthologs(c_empty, all_contigs, {'C_empty'})
    assert count == 0

def test_busco_ortholog_rule_prevents_disqualification():
    # C has more unique orthologs than O, so O cannot disqualify C
    # even though O has better score
    c1 = ContigSummary(query_id='C1', query_length=1000, target_id='T1',
                       intervals=[(100, 500)], sum_normalized_score=0.7, 
                       status=Status.ALIGNED_RETAINED,
                       busco_genes={'BUSCO1', 'BUSCO2', 'BUSCO3'})
    c2 = ContigSummary(query_id='C2', query_length=1000, target_id='T1',
                       intervals=[(200, 600)], sum_normalized_score=0.9,
                       status=Status.ALIGNED_RETAINED,
                       busco_genes={'BUSCO4'})
    c3 = ContigSummary(query_id='C3', query_length=1000, target_id='T2',
                       status=Status.ALIGNED_RETAINED,
                       busco_genes={'BUSCO4'})
    
    mash_lookup = {'C1': {'C2': 0.01}, 'C2': {'C1': 0.01}}
    dist_threshold = 0.05
    
    # C1 has 3 unique orthologs (BUSCO1, BUSCO2, BUSCO3 not in C3)
    # C2 has 0 unique orthologs (BUSCO4 is in C3)
    # So C2 cannot disqualify C1 despite better score
    results = run_tournament_on_target([c1, c2], mash_lookup, dist_threshold, 
                                       all_contigs=[c1, c2, c3])
    
    c1_res = next(c for c in results if c.query_id == 'C1')
    c2_res = next(c for c in results if c.query_id == 'C2')
    
    assert c1_res.status == Status.ALIGNED_RETAINED
    assert c2_res.status == Status.ALIGNED_RETAINED

def test_busco_ortholog_rule_allows_disqualification():
    # C has fewer or equal unique orthologs than O, so O can disqualify C
    c1 = ContigSummary(query_id='C1', query_length=1000, target_id='T1',
                       intervals=[(100, 500)], sum_normalized_score=0.7,
                       status=Status.ALIGNED_RETAINED,
                       busco_genes={'BUSCO1'})
    c2 = ContigSummary(query_id='C2', query_length=1000, target_id='T1',
                       intervals=[(200, 600)], sum_normalized_score=0.9,
                       status=Status.ALIGNED_RETAINED,
                       busco_genes={'BUSCO2', 'BUSCO3'})
    c3 = ContigSummary(query_id='C3', query_length=1000, target_id='T2',
                       status=Status.ALIGNED_RETAINED,
                       busco_genes={'BUSCO4'})
    
    mash_lookup = {'C1': {'C2': 0.01}, 'C2': {'C1': 0.01}}
    dist_threshold = 0.05
    
    # C1 has 1 unique ortholog (BUSCO1)
    # C2 has 2 unique orthologs (BUSCO2, BUSCO3)
    # So C2 can disqualify C1
    results = run_tournament_on_target([c1, c2], mash_lookup, dist_threshold,
                                       all_contigs=[c1, c2, c3])
    
    c1_res = next(c for c in results if c.query_id == 'C1')
    c2_res = next(c for c in results if c.query_id == 'C2')
    
    assert c2_res.status == Status.ALIGNED_RETAINED
    assert c1_res.status == Status.ALIGNED_DISCARDED
    assert c1_res.disqualifier == 'C2'


def test_get_adjusted_score_without_all_contigs():
    # Test that when all_contigs is None, all BUSCO genes are counted
    c1 = ContigSummary(query_id='C1', query_length=1_000_000, 
                       sum_normalized_score=1.0,
                       status=Status.ALIGNED_RETAINED,
                       busco_genes={'BUSCO1', 'BUSCO2', 'BUSCO3'})
    
    # Without all_contigs, should use all 3 BUSCO genes
    # busco_density = 3 / 1.0 = 3.0
    # adjusted_score = 1.0 * (1 + 3.0 * 0.1) = 1.0 * 1.3 = 1.3
    score = get_adjusted_score(c1, busco_bonus_factor=0.1, all_contigs=None)
    assert score == 1.3
    
    # With busco_bonus_factor=0.0, should return base score
    score = get_adjusted_score(c1, busco_bonus_factor=0.0, all_contigs=None)
    assert score == 1.0

def test_get_adjusted_score_with_unique_buscos():
    # Test that when all_contigs is provided, only unique BUSCO genes are counted
    c1 = ContigSummary(query_id='C1', query_length=1_000_000,
                       sum_normalized_score=1.0,
                       status=Status.ALIGNED_RETAINED,
                       busco_genes={'BUSCO1', 'BUSCO2', 'BUSCO3'})
    c2 = ContigSummary(query_id='C2', query_length=1_000_000,
                       sum_normalized_score=0.9,
                       status=Status.ALIGNED_RETAINED,
                       busco_genes={'BUSCO2', 'BUSCO4'})
    c3 = ContigSummary(query_id='C3', query_length=1_000_000,
                       status=Status.ALIGNED_RETAINED,
                       busco_genes={'BUSCO5'})
    
    all_contigs = [c1, c2, c3]
    
    # For C1: BUSCO1 and BUSCO3 are unique (BUSCO2 is in C2)
    # unique_busco_count = 2
    # busco_density = 2 / 1.0 = 2.0
    # adjusted_score = 1.0 * (1 + 2.0 * 0.1) = 1.0 * 1.2 = 1.2
    score = get_adjusted_score(c1, busco_bonus_factor=0.1, all_contigs=all_contigs)
    assert score == 1.2
    
    # For C2: BUSCO4 is unique (BUSCO2 is in C1)
    # unique_busco_count = 1
    # busco_density = 1 / 1.0 = 1.0
    # adjusted_score = 0.9 * (1 + 1.0 * 0.1) = 0.9 * 1.1 = 0.99
    score = get_adjusted_score(c2, busco_bonus_factor=0.1, all_contigs=all_contigs)
    assert score == 0.99
    
    # For C3: BUSCO5 is unique
    # unique_busco_count = 1
    # busco_density = 1 / 1.0 = 1.0
    # adjusted_score = 0.0 * (1 + 1.0 * 0.1) = 0.0
    score = get_adjusted_score(c3, busco_bonus_factor=0.1, all_contigs=all_contigs)
    assert score == 0.0

def test_get_adjusted_score_all_shared_buscos():
    # Test when all BUSCO genes are shared with other contigs
    c1 = ContigSummary(query_id='C1', query_length=1_000_000,
                       sum_normalized_score=1.0,
                       status=Status.ALIGNED_RETAINED,
                       busco_genes={'BUSCO1', 'BUSCO2'})
    c2 = ContigSummary(query_id='C2', query_length=1_000_000,
                       sum_normalized_score=0.9,
                       status=Status.ALIGNED_RETAINED,
                       busco_genes={'BUSCO1', 'BUSCO2', 'BUSCO3'})
    
    all_contigs = [c1, c2]
    
    # For C1: BUSCO1 and BUSCO2 are both in C2, so unique_busco_count = 0
    # busco_density = 0 / 1.0 = 0.0
    # adjusted_score = 1.0 * (1 + 0.0 * 0.1) = 1.0
    score = get_adjusted_score(c1, busco_bonus_factor=0.1, all_contigs=all_contigs)
    assert score == 1.0

def test_get_adjusted_score_no_busco_genes():
    # Test contig with no BUSCO genes
    c1 = ContigSummary(query_id='C1', query_length=1_000_000,
                       sum_normalized_score=1.0,
                       status=Status.ALIGNED_RETAINED)
    c2 = ContigSummary(query_id='C2', query_length=1_000_000,
                       sum_normalized_score=0.9,
                       status=Status.ALIGNED_RETAINED,
                       busco_genes={'BUSCO1'})
    
    all_contigs = [c1, c2]
    
    # C1 has no BUSCO genes, should return base score
    score = get_adjusted_score(c1, busco_bonus_factor=0.1, all_contigs=all_contigs)
    assert score == 1.0
