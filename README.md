# FilterHaplotypes: Reference-Based De-duplication

FilterHaplotypes is a command-line tool designed to remove redundant contigs (haplotigs) from highly duplicated genome assemblies (e.g., polyploid or pooled samples). It leverages a reference genome to identify the most representative contig for each genomic region by combining alignment scores (Minimap2) with sequence similarity clustering (Mash 2.0).

## Algorithmic Phases

### Phase 1: Pre-processing
Retains quality alignments (`MappingQuality` >= 20) and identifies the primary `target` locus for each query contig based on the highest Alignment Score (`AS:i`).
- Potential problem: If a contig aligns to multiple targets, it may be difficult to determine the primary target locus. The highest alignment score is probably not the best choice. Look into using the highest 90th percentile `AS:i` score (among all targets' alignments) instead.

### Phase 2: Contig Summary Initialization
Creates a central data structure (`ContigSummary`) for every contig in the de-novo assembly. Aligned contigs are linked to their primary targets. At the end of this phase, a GC content filter is applied: the median GC% is calculated from contigs in the 90th percentile for length, and any contig whose GC% deviates by more than 5 percentage points from this median is discarded as a potential contaminant.

### Phase 3: Alignment Filtering and Scoring
Uses a greedy tiling algorithm to resolve overlapping alignments for each contig on its primary target. A final normalized alignment score is calculated using only these non-redundant alignments.

### Phase 4: Mash Distance Threshold Calculation
Analyzes the distribution of Mash distances among overlapping contig pairs. If not provided by the user, it automatically estimates a biologically meaningful threshold using Kernel Density Estimation (KDE) to distinguish between redundant and unique sequences.

### Phase 5: Redundancy Resolution (Iterative Tournament)
Runs a multi-pass tournament within each genomic locus. Contigs compete based on score, sequence similarity (Mash), and size. An iterative rescue loop ensures that "orphaned" contigs (those whose disqualifiers were themselves discarded) get a fair challenge to maintain optimal coverage.
- May want to include some consideration of BUSCO genes on competing contigs.

### Phase 6: Unaligned Contig Handling
Identifies contigs that did not align to the reference and screens them for redundancy against all retained contigs using Mash distances. Can be bypassed with the `--aligned-only` flag, which discards all unaligned contigs.

### Phase 7: Reporting
Generates the de-duplicated FASTA assembly, a detailed TSV summary of all filtering decisions, and an interactive HTML report with assembly statistics, BUSCO completeness, and diagnostic plots (Mash distribution, BlobPlot, L-curve).
-Note: Add the user's command prompt to the HTML report output.

## Installation

```bash
# Clone the repository
git clone https://github.com/NathanDHarry/FilterHaplotypes.git
cd FilterHaplotypes

# Install dependencies
pip install -e .
```

### Dependencies
- Python >= 3.8
- Biopython
- Pandas
- Plotly
- Jinja2
- SciPy
- NumPy

## Usage

```bash
filter_haplotypes -p alignments.paf -m mash_distances.dist -f contigs.fasta -o output_dir
```

### Mandatory Arguments
- `-p, --paf`: Alignment file (PAF) generated via `minimap2 -c`.
- `-m, --mash`: Mash distance TSV (Query vs Query).
- `-f, --fasta`: Original assembly FASTA file.

### Optional Arguments
- `-b, --busco`: Optional BUSCO `full_table.tsv` for completeness reporting.
- `-o, --output`: Directory for results (Default: `./output`).
- `--min-mq`: Mapping Quality threshold (Default: 20).
- `--aligned-only`: Discard all unaligned contigs (bypass Phase 6).
- `--min-overlap`: Minimum overlap bases to trigger competition (Default: 1).
- `--min-size-safeguard`: Ratio for size-based retention (Default: 0.50).
- `--distance-threshold`: Mash distance threshold (Overrides estimator if supplied).
- `--threads`: Number of CPU cores for parallel processing (Default: Max - 1).
- `--max-tournament-iterations`: Maximum iterations for tournament loop (Default: 100000).

## Example Command

```bash
filter_haplotypes --paf data/example_alignments.paf --mash data/example_mash.dist --fasta data/example_contigs.fasta --busco data/example_BUSCO_full_table.tsv --output my_results
```
