import pytest
import pandas as pd
from src.filter_haplotypes.core.models import ContigSummary, Status
from src.filter_haplotypes.core.filtering import (
    tile_and_score_contig,
    estimate_distance_threshold,
    run_tournament_on_target,
    screen_unaligned_contig
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
    results = run_tournament_on_target([c1, c2], mash_lookup, dist_threshold)
    
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
