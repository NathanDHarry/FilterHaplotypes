
# FilterHaplotypes Phase 1: Basic Contig Selection

## Purpose

To filter redundant contigs from highly duplicated genome assemblies by selecting the highest-quality, non-overlapping contig per genomic region using reference-based alignment analysis.

## Problem Statement

Genome assemblies from polyploid or pooled samples often contain multiple haplotypes or redundant contigs for single-copy regions. While existing tools like `purge_dups` and `Redundans` have limitations, a reference-based approach can effectively identify and retain only the best-quality contig per genomic region.

**Example use case**: A HiFi assembly from a pooled sample of five full-sibling individuals sequenced on PacBio Revio produces highly duplicated contigs for single-copy orthologs. Using a previous reference genome from the same species, Phase 1 filters out lower-quality haplotypes while retaining non-redundant, high-quality contigs suitable for further scaffolding.

## Solution Overview

**Phase 1** provides a streamlined, focused approach: given a pre-computed minimap2 alignment file (PAF format) of query contigs aligned to a reference genome, the tool:

1. Divides the reference into fixed-size windows
2. Selects the best-scoring contig per window using composite quality-coverage metrics
3. Resolves overlaps on the query sequence using greedy selection
4. Outputs only the final non-redundant contig set

This lightweight implementation serves as a foundation for future enhancements while remaining immediately useful for basic contig filtering tasks.

---

## Workflow Overview

### Step-by-Step Process

1. **Input Validation**
   - Verify PAF alignment file format and integrity
   - Verify query contig FASTA file format and integrity
   - Confirm both files contain expected data

2. **Window Definition**
   - Divide the reference genome into fixed-size windows (default: 50,000 bp)
   - Window boundaries are determined by position on the reference sequence

3. **Per-Window Candidate Selection (Phase A)**
   - For each reference window, identify all query contigs with alignments in that window
   - For each candidate contig, compute:
     - **Alignment Score**: Raw minimap2 alignment score (from PAF column 12)
     - **Window Coverage**: Percentage of the reference window covered by the contig's alignment
     - **Orientation**: Direction of alignment (forward/reverse) for consistency checking
   - Select the single best-scoring candidate per window
   - Result: Candidate contig set (may contain overlaps on query sequence)

4. **Composite Priority Scoring (Phase B Setup)**
   - For all candidates, calculate composite priority scores that balance quality with coverage:
     ```
     Priority Score = Alignment Score + (Coverage Weight × Normalized Query Length)
     ```
   - Normalize query alignment lengths relative to the longest candidate in the set
   - Apply configured coverage weight (default: 0.1)

5. **Greedy Overlap Resolution (Phase B)**
   - Sort candidates by composite priority (highest to lowest)
   - Iteratively select contigs:
     - Initialize empty kept-contig set
     - For each candidate (in priority order):
       - If it overlaps any contig already in the kept set on the query sequence: **discard**
       - If it doesn't overlap any kept contig: **add to kept set**
   - Result: Final non-redundant contig set with maximized coverage

6. **Output Generation**
   - Write selected contigs to output FASTA file
   - Generate summary log with per-window statistics
   - Report filtering decisions and metrics

---

## Contig Selection and Overlap Resolution

### Phase A: Per-Window Selection

Each reference window independently identifies candidate contigs based on alignment quality:

- **Selection Criterion**: Highest alignment score in the window
- **Result**: One candidate per window (candidates may overlap on query sequence)
- **Metrics Computed**:
  - Alignment score (primary metric)
  - Window coverage percentage
  - Alignment orientation (forward/reverse)

### Phase B: Greedy Overlap Resolution with Composite Scoring

The tool employs a **greedy coverage-maximization strategy** combined with **composite priority scoring** to preserve non-redundant contigs while maximizing total sequence coverage.

#### Composite Priority Metric

The tool prioritizes contigs using a composite metric that balances alignment quality with alignment length:
```text
Priority Score = Alignment Score + (Coverage Weight × Normalized Query Length)
```

Where:
- **Alignment Score**: minimap2's raw alignment score (column 12 in PAF; primary quality indicator)
- **Normalized Query Length**: The query contig alignment length (query end – query start) normalized relative to the longest candidate alignment in the full candidate set
- **Coverage Weight**: A configurable weight (default: 0.1) that determines how much alignment length influences the final priority relative to score quality

This composite metric ensures that:
- **High-quality alignments are prioritized** — Alignment score remains the dominant factor
- **Comparable-quality contigs don't lose coverage** — When alignment scores are similar, longer contigs are favored
- **Balanced trade-off** — A dramatically longer but slightly lower-quality contig won't unconditionally displace a short high-quality one

#### Algorithm Steps

1. **Collect all per-window candidates** into a single set
2. **Calculate composite priority scores** for all candidates
3. **Sort candidates by composite priority** (highest to lowest)
4. **Greedily select contigs**:
   - Initialize an empty set of kept contigs
   - For each candidate contig (in descending priority order):
     - Check if it overlaps any previously-kept contig on the query sequence
     - If **no overlap**: add to kept set
     - If **overlap exists**: discard candidate
5. **Return**: Final non-redundant set with maximized coverage

#### Overlap Definition

Two contigs are considered **overlapping** if their query sequence intervals intersect on the same strand. Formally, for two contigs A and B:
- Overlapping if: `max(queryStart_A, queryStart_B) < min(queryEnd_A, queryEnd_B)`
- Non-overlapping otherwise

#### Key Features

- **Quality-First with Coverage Awareness**: Prioritizes alignment quality while recognizing the value of additional genomic sequence
- **Prevents Pathological Loss**: Avoids discarding significantly longer contigs when their quality metrics are only marginally lower
- **Configurable Balance**: The coverage weight can be tuned to adjust the relative importance of length vs. quality
- **Transparent Scoring**: All contigs are evaluated using the same composite metric, ensuring consistent, reproducible results
- **Non-Redundancy Guarantee**: No two kept contigs overlap on the query sequence
- **Computational Efficiency**: Single pass through sorted contigs; O(n log n) complexity dominated by sorting

#### Detailed Example

**Initial candidates after per-window selection:**
- Contig_A: query positions 1,000–30,000 (30,000 bp), alignment score 95
- Contig_B: query positions 20,000–50,000 (30,000 bp), alignment score 92
- Contig_C: query positions 35,000–75,000 (40,000 bp), alignment score 91
- Contig_D: query positions 80,000–95,000 (15,000 bp), alignment score 90

**Composite score calculation** (Coverage Weight = 0.1, max query length = 40,000 bp):
- Contig_A: `95 + (0.1 × 30000/40000) = 95 + 0.75 = 95.75`
- Contig_B: `92 + (0.1 × 30000/40000) = 92 + 0.75 = 92.75`
- Contig_C: `91 + (0.1 × 40000/40000) = 91 + 1.0 = 92.0`
- Contig_D: `90 + (0.1 × 15000/40000) = 90 + 0.375 = 90.375`

**Overlap relationships on query sequence:**
- Contig_A and Contig_B: overlap at positions 20,000–30,000
- Contig_B and Contig_C: overlap at positions 35,000–50,000
- Contig_A and Contig_C: **no overlap** (gap at 30,000–35,000)
- Contig_D: **no overlap** with any other contig (gap at 75,000–80,000)

**Greedy selection process** (sorted by composite priority, descending):
1. Contig_A (priority 95.75) → No kept contigs yet → **Keep**
2. Contig_B (priority 92.75) → Overlaps Contig_A (20,000–30,000) → **Discard**
3. Contig_C (priority 92.0) → Does not overlap Contig_A (gap at 30,000–35,000) → **Keep**
4. Contig_D (priority 90.375) → Does not overlap Contig_A or Contig_C → **Keep**

**Final output:**
- Contig_A: query positions 1,000–30,000 (30,000 bp)
- Contig_C: query positions 35,000–75,000 (40,000 bp)
- Contig_D: query positions 80,000–95,000 (15,000 bp)

**Total coverage**: 85,000 bp of non-redundant, high-quality contigs

**Note**: Contig_C was retained despite slightly lower alignment score than Contig_B because:
1. It provides substantially more unique coverage (40,000 bp vs. 30,000 bp for B's unique portion)
2. Its composite priority score (92.0) is comparable to Contig_B's (92.75)
3. It does not overlap Contig_A, allowing both to coexist

#### Why This Approach

- **vs. Score-Only Greedy**: Balances quality with coverage; prevents pathological cases where a dramatically longer contig with marginally lower quality is discarded
- **vs. Single Best Per Overlap Cluster**: Maximizes coverage by keeping independent contigs across different regions, rather than discarding all but one contig from each overlapping group
- **vs. Length-First Selection**: Ensures quality remains the primary criterion, avoiding selection of low-identity or fragmented contigs simply because they are long
- **vs. Exhaustive Optimization**: Greedy selection is computationally tractable and deterministic, avoiding exponential complexity while maintaining transparent, reproducible decisions

#### Configurable Coverage Weight

- **Parameter**: `--coverage-weight` (default: 0.1)
- **Effect**: Controls how much alignment length influences the final priority relative to alignment score
  - **Higher values** (e.g., 0.2–0.5): Emphasize coverage; longer contigs become more competitive with slightly higher-scoring contigs
  - **Lower values** (e.g., 0.01–0.05): Emphasize quality; alignment score dominates almost entirely, length matters only for near-identical scores
  - **Recommended tuning**: Start with 0.1; increase if you observe loss of valuable longer contigs, decrease if shorter contigs are prioritized too aggressively

#### Limitations and Assumptions

- Assumes alignment score adequately captures quality differences between contigs
- Does not consider alignment fragmentation (number of alignment blocks per window)
- Does not evaluate internal contig structure (duplications, repeats, inversions)
- Does not assess read coverage depth for selected contigs
- May discard lower-scoring contigs with markedly different fragmentation patterns if the score difference is significant

---

## Usage

### Installation

Clone the repository and set up the virtual environment:
```bash
bash git clone <repository_url> 
cd FilterHaplotypes 
python -m venv .venv 
source .venv/bin/activate # On Windows: .venv\Scripts\activate pip install -r requirements.txt
```

### Command-Line Interface

#### Basic Usage
```bash
bash python src/main.py --alignment <path_to_paf_file> --contigs <path_to_fasta_file> --output <path_to_output_fasta>
```

#### Complete Example

```bash
python src/main.py \
  --alignment alignments.paf \
  --contigs query_contigs.fasta \
  --output filtered_contigs.fasta \
  --window-size 50000 \
  --coverage-weight 0.1 \
  --summary summary.log \
  --log-level INFO
```
### Required Arguments

| Argument | Short | Type | Description |
| --- | --- | --- | --- |
| `--alignment` | `-a` | FILE | Path to PAF alignment file (minimap2 output). File must be readable and in valid PAF format. |
| `--contigs` | `-c` | FILE | Path to query contig FASTA file. File must be readable and contain valid FASTA sequences. |
| `--output` | `-o` | FILE | Path for output FASTA file. Will be created or overwritten. Parent directory must be writable. |
### Optional Arguments

| Argument | Short | Type | Default | Description |
| --- | --- | --- | --- | --- |
| `--window-size` | `-w` | INT | 50000 | Fixed reference window size in base pairs. Determines the resolution at which per-window contig selection occurs. Must be > 0. |
| `--coverage-weight` | `-cw` | FLOAT | 0.1 | Weight controlling how alignment length influences priority scoring. Higher values favor longer contigs; lower values favor higher scores. Must be ≥ 0.0. |
| `--summary` | `-s` | FILE | (none) | Path for output summary log file. If not specified, no summary log is generated. TSV format. |
| `--log-level` | `-l` | CHOICE | INFO | Logging verbosity level. Valid values: DEBUG, INFO, WARNING, ERROR. DEBUG produces detailed workflow trace; INFO shows major steps; WARNING shows only warnings and errors. |
| `--memory-limit` | `-m` | INT | (system) | User-specified memory limit in MB. If system memory usage would exceed this limit, the tool will attempt to process files in chunks (future enhancement). If not specified, defaults to 80% of available system memory. |
| `--help` | `-h` | — | — | Display help message and exit. |
| `--version` | `-v` | — | — | Display version information and exit. |
### Input File Specifications
#### PAF Alignment File
- **Format**: PAF (Pairwise Alignment Format)
- **Required Columns**:
    - Column 1: Query sequence name (contig ID)
    - Column 2: Query sequence length
    - Column 3: Query alignment start (0-based)
    - Column 4: Query alignment end (0-based)
    - Column 5: Relative strand (+/-)
    - Column 6: Target sequence name (reference contig ID)
    - Column 7: Target sequence length
    - Column 8: Target alignment start (0-based)
    - Column 9: Target alignment end (0-based)
    - Column 10: Number of residue matches
    - Column 11: Alignment block length
    - Column 12: Mapping quality (0–60)
    - **Column 12 (AS tag)**: Alignment score (SAM-style, used for scoring)

- **Generation**: `minimap2 -cx asm10 reference.fasta query.fasta > alignments.paf`
- **Note**: Must contain at least alignment score; coordinates must be valid (start < end)

#### FASTA Contig File
- **Format**: FASTA
- **Requirements**:
    - Standard FASTA format (header line with `>` prefix, sequences on following lines)
    - Sequence IDs must match query names in PAF file
    - No empty sequences
    - ASCII text encoding

- **Notes**: Sequences can span multiple lines; will be read and concatenated

### Output File Specifications
#### Filtered Contigs FASTA (`--output`)
- **Format**: FASTA
- **Content**: Only contigs selected by the greedy overlap resolution
- **Order**: Sorted by query start position (ascending) across all windows
- **Format**: Standard FASTA (header line with `>` prefix, 80-character wrapped sequences)
- **Example**:
```text
>contig_A
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
>contig_B
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
 ```
#### Summary Log (`--summary`)
- **Format**: Tab-separated values (TSV)
- **Header**: `RefWindow_Start RefWindow_End SelectedContig QueryStart QueryEnd QueryLength AlignmentScore CompositeScore WindowCoverage Orientation DiscardedContigs DiscardReason`
- **Rows**: One row per reference window, containing:
    - **RefWindow_Start**: 0-based start position of reference window
    - **RefWindow_End**: 0-based end position of reference window
    - **SelectedContig**: Name of selected contig (or "NONE" if no alignment in window)
    - **QueryStart**: 0-based start position of contig alignment on query
    - **QueryEnd**: 0-based end position of contig alignment on query
    - **QueryLength**: Length of contig alignment on query (end – start)
    - **AlignmentScore**: Alignment score from PAF
    - **CompositeScore**: Calculated composite priority score
    - **WindowCoverage**: Percentage of reference window covered by alignment (0–100)
    - **Orientation**: Forward (+) or reverse (–) alignment
    - **DiscardedContigs**: Comma-separated list of contigs discarded due to overlap
    - **DiscardReason**: Reason for discard (e.g., "lower_composite_priority", "overlaps_selected_contig")

- **Example**:
```text
  RefWindow_Start	RefWindow_End	SelectedContig	QueryStart	QueryEnd	QueryLength	AlignmentScore	CompositeScore	WindowCoverage	Orientation	DiscardedContigs	DiscardReason
  0	50000	contig_A	1000	30000	29000	95	95.75	100	+	contig_X	lower_composite_priority
  50000	100000	contig_C	35000	75000	40000	91	92.0	100	+	contig_Y	overlaps_selected_contig
```
### Validation and Error Handling
The tool performs comprehensive input validation:
- **File Existence**: Verifies input files exist and are readable
- **Format Validation**:
    - PAF: Checks for minimum required columns, valid coordinates (start < end), numeric fields
    - FASTA: Checks for valid format (header lines, non-empty sequences)

- **Consistency Checks**: Verifies contig IDs in PAF match those in FASTA
- **Data Sanity**: Checks for reasonable window sizes, coverage weights, coordinate ranges
- **Disk Space**: Verifies output directory is writable and has sufficient space
- **Error Messages**: Provides detailed, user-friendly error messages explaining validation failures

On error, the tool:
1. Logs detailed error message to console and log file
2. Exits with non-zero status code
3. Does not generate partial output files

## Requirements
### System Requirements
- **Operating System**: Linux, macOS, or Windows (with bash/WSL for command-line usage)
- **Python**: 3.13 or higher
- **RAM**: Minimum 2 GB (depends on input file size)
- **Disk Space**: At least 2× the combined size of input files for processing

### Python Dependencies
- Standard library only (no external dependencies for Phase 1)
    - : CLI argument parsing `argparse`
    - `logging`: Event logging
    - : File path handling `pathlib`
    - : Temporary file management `tempfile`
    - , : System operations `sys``os`

## Project Structure
```text
FilterHaplotypes/
├── src/
│   ├── file_handler.py       # Input parsing and validation
│   ├── metrics.py            # Per-window selection and composite scoring
│   ├── evaluation.py         # Overlap detection and greedy resolution
│   ├── output_handler.py     # Output file generation and logging
│   ├── utils.py              # Shared utilities (CLI, logging)
│   └── main.py               # Workflow orchestration and entry point
├── tests/
│   ├── test_file_handler.py
│   ├── test_metrics.py
│   ├── test_evaluation.py
│   ├── test_output_handler.py
│   └── test_integration.py
├── data/
│   ├── example_alignments.paf
│   ├── example_contigs.fasta
│   └── expected_output.fasta
├── environment.yaml          # Conda environment specification
├── requirements.txt          # pip dependencies (minimal for Phase 1)
├── README_phase1.md          # This file
└── Outline_phase1.md         # Technical design outline
```
## Implementation Details
### Architecture
The implementation uses a modular, object-oriented design with clear separation of concerns:
- **`file_handler.py`**: Input file parsing and validation
    - Reads and validates PAF alignment data
    - Reads and validates FASTA contig data
    - Returns structured data objects for downstream processing

- **`metrics.py`**: Per-window selection and composite scoring
    - Divides reference into fixed-size windows
    - Groups alignments by window
    - Selects best contig per window
    - Calculates composite priority scores for all candidates

- **`evaluation.py`**: Overlap detection and greedy resolution
    - Detects overlaps between contigs on query sequence
    - Implements greedy selection algorithm
    - Returns final non-redundant contig set

- **`output_handler.py`**: Output file generation
    - Writes selected contigs to FASTA file
    - Generates summary log with statistics
    - Formats output data according to specifications

- **`utils.py`**: Shared utilities
    - CLI argument parsing and help message generation
    - Logger initialization and configuration
    - File and directory utilities

- **`main.py`**: Workflow orchestration
    - Coordinates execution of above modules
    - Handles error propagation
    - Manages temporary files and cleanup

### Data Flow
1. **Input**: User specifies PAF file, FASTA file, output path, and optional parameters
2. **Parsing**: FileHandler reads and validates input files, returns data structures
3. **Windowing**: MetricsCalculator divides reference into windows
4. : For each window, select best-scoring contig **Per-Window Selection**
5. **Overlap Resolution**: Evaluate all candidates, calculate composite scores, greedily select non-overlapping set
6. **Output**: Write filtered contigs to FASTA, generate summary log

### Error Handling Strategy
- Input validation before processing begins
- Informative error messages at each step
- Graceful failure with cleanup of partial outputs
- Detailed logging for debugging

### Performance Considerations
- **Time Complexity**: O(n log n) for sorting, O(n²) worst-case for overlap detection (O(n) for greedy selection with efficient overlap checking)
- **Space Complexity**: O(n) for storing alignment data in memory
- **Scalability**: Current implementation assumes data fits in memory; future versions may support chunked processing for very large alignments

## Examples
Example input and output files are provided in the `data/` directory:
- **`example_alignments.paf`**: Sample minimap2 alignment output
- **`example_contigs.fasta`**: Sample query contigs
- **`expected_output.fasta`**: Expected filtered contig output

To test the tool:
```bash
python src/main.py \
  --alignment data/example_alignments.paf \
  --contigs data/example_contigs.fasta \
  --output /tmp/test_output.fasta \
  --summary /tmp/test_summary.log \
  --log-level DEBUG
```
Compare the output to `expected_output.fasta` for validation.

## Limitations and Future Work
### Current Scope (Phase 1)
**Implemented:**
- ✅ Selects best contig per reference window based on alignment score
- ✅ Computes composite priority scores balancing quality and coverage
- ✅ Resolves overlaps using greedy selection with non-redundancy guarantee
- ✅ Outputs non-redundant contig set in FASTA format
- ✅ Generates detailed summary log
- ✅ Comprehensive input validation and error handling
- ✅ Configurable window size and coverage weight
- ✅ Detailed logging and debugging support

**Not Implemented (Out of Scope):**
- ❌ Scaffolding of selected contigs
- ❌ Evaluation of contig internal structure (duplications, repeats, inversions)
- ❌ Read coverage depth analysis
- ❌ Dynamic window sizing based on alignment density
- ❌ Support for SAM/BAM alignment formats
- ❌ Fragmentation scoring (alignment block count and spacing)
- ❌ Parallelization of window processing

### Future Enhancements (Phase 2+)
**Planned improvements:**
- **Advanced Metrics**: Fragmentation score (number and spacing of alignment blocks)
- **Structural Analysis**: Detection of internal duplications, palindromic repeats, inversions
- **Coverage Analysis**: Alignment of HiFi reads back to assembly for depth-of-coverage metrics
- **Dynamic Windows**: Adaptive window sizing based on alignment density in reference regions
- **Format Support**: SAM/BAM input in addition to PAF
- **Performance**: Parallelization for multi-reference processing
- **Scaffolding**: Integration with scaffolding tools to assemble selected contigs
- **Visualization**: Summary plots and interactive reports

## Implementation Notes for Development
When implementing the code based on this specification, keep in mind:
1. **PAF Parsing**: Carefully parse all required columns; many minimap2 outputs include optional fields. Alignment score may be in column 12 (SAM AS tag format) or as an optional tag.
2. **Coordinate Systems**: Ensure consistent use of 0-based, half-open intervals (start inclusive, end exclusive) throughout.
3. **Overlap Detection**: Implement efficient overlap checking; consider interval tree or sorted interval data structure for large candidate sets.
4. **Floating-Point Precision**: When computing composite scores and normalized lengths, use appropriate precision to avoid comparison artifacts.
5. **Memory Efficiency**: For large PAF files (100k+ alignments), consider streaming processing rather than loading entire file into memory.
6. **Logging**: Provide detailed debug output during per-window selection and overlap resolution for troubleshooting; users should be able to trace why specific contigs were selected or discarded.
7. **FASTA Output**: Ensure standard FASTA format compliance; wrap sequences to 80 characters per line; preserve sequence order.
8. **Summary Log**: TSV format is preferred for compatibility with downstream analysis tools (R, Python pandas, etc.); ensure proper escaping of special characters.

## Troubleshooting
### Common Issues
**Issue**: "PAF file contains invalid coordinates (start >= end)"
- **Solution**: Verify PAF file was generated with minimap2 using correct syntax. Check for file corruption or improper editing.

**Issue**: "Contig IDs in FASTA do not match PAF query names"
- **Solution**: Ensure FASTA and PAF files are consistent. Check for whitespace or special character differences in sequence headers. Regenerate alignment if files were modified.

**Issue**: No contigs selected (empty output file)
- **Solution**: Check log file for details. Possible causes: window size too small, coverage weight set too high, or alignment file is empty. Try increasing window size or decreasing coverage weight.

**Issue**: All contigs from input appear in output (no filtering occurred)
- **Solution**: Check that overlap detection is working correctly. Verify that candidate contigs actually overlap on query sequence. Try with smaller window size to increase candidate set diversity.

## Citation and Attribution
If you use FilterHaplotypes Phase 1 in your research, please cite:
- [Full citation to be added]

## License
[License to be specified]
## Contact
[Contact information to be added]
