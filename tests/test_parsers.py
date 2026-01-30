import pytest
from pathlib import Path
from src.filter_haplotypes.parsers.fasta_parser import parse_fasta
from src.filter_haplotypes.parsers.paf_parser import parse_paf, get_primary_targets
from src.filter_haplotypes.parsers.mash_parser import parse_mash, build_mash_lookup
from src.filter_haplotypes.parsers.busco_parser import parse_busco

DATA_DIR = Path("data")

def test_parse_fasta():
    fasta_path = DATA_DIR / "example_contigs.fasta"
    if not fasta_path.exists():
        pytest.skip("Example FASTA not found")
    
    results = parse_fasta(str(fasta_path), threads=1)
    assert len(results) > 0
    # Check one contig
    first_id = list(results.keys())[0]
    gc, length = results[first_id]
    assert 0 <= gc <= 100
    assert length > 0

def test_parse_paf():
    paf_path = DATA_DIR / "example_alignments.paf"
    if not paf_path.exists():
        pytest.skip("Example PAF not found")
        
    df = parse_paf(str(paf_path), min_mq=10)
    assert not df.empty
    assert 'AS' in df.columns
    assert (df['mq'] >= 10).all()

def test_get_primary_targets():
    paf_path = DATA_DIR / "example_alignments.paf"
    if not paf_path.exists():
        pytest.skip("Example PAF not found")
        
    df = parse_paf(str(paf_path), min_mq=10)
    primary_df = get_primary_targets(df)
    
    # Each query_id should only have one target_id
    assert primary_df.groupby('query_id')['target_id'].nunique().max() == 1

def test_parse_mash():
    mash_path = DATA_DIR / "example_mash.dist"
    if not mash_path.exists():
        pytest.skip("Example Mash not found")
        
    df = parse_mash(str(mash_path))
    assert not df.empty
    assert (df['p_value'] < 0.05).all()
    
    lookup = build_mash_lookup(df)
    assert len(lookup) > 0

def test_parse_busco():
    busco_path = DATA_DIR / "example_BUSCO_full_table.tsv"
    if not busco_path.exists():
        pytest.skip("Example BUSCO not found")
        
    busco_map = parse_busco(str(busco_path))
    assert len(busco_map) > 0
    for q_id, genes in busco_map.items():
        assert isinstance(genes, set)
        assert len(genes) > 0
