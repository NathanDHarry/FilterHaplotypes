import pytest
import subprocess
import sys
from pathlib import Path
import shutil

DATA_DIR = Path("data")
OUTPUT_DIR = Path("test_output")

def test_full_pipeline():
    if OUTPUT_DIR.exists():
        shutil.rmtree(OUTPUT_DIR)
    
    paf = DATA_DIR / "example_alignments.paf"
    mash = DATA_DIR / "example_mash.dist"
    fasta = DATA_DIR / "example_contigs.fasta"
    busco = DATA_DIR / "example_BUSCO_full_table.tsv"
    
    if not all(p.exists() for p in [paf, mash, fasta, busco]):
        pytest.skip("Integration test data missing")
        
    cmd = [
        sys.executable, "-m", "src.filter_haplotypes.main",
        "-p", str(paf),
        "-m", str(mash),
        "-f", str(fasta),
        "-b", str(busco),
        "-o", str(OUTPUT_DIR),
        "--threads", "2"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    print(result.stdout)
    print(result.stderr)
    
    assert result.returncode == 0
    assert (OUTPUT_DIR / "filtered_assembly.fasta").exists()
    assert (OUTPUT_DIR / "summary_report.tsv").exists()
    assert (OUTPUT_DIR / "report.html").exists()
    assert (OUTPUT_DIR / "log.txt").exists()
