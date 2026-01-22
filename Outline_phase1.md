
# FilterHaplotypes Phase 1: Technical Design and Implementation Outline

## Overview

This document provides a detailed technical design and architectural outline for FilterHaplotypes Phase 1, which implements basic contig filtering based on reference-guided alignment analysis. The tool processes minimap2 PAF alignment files and selects non-redundant, high-quality contigs using per-window selection combined with greedy overlap resolution and composite quality-coverage scoring.

---

## Workflow Architecture

### High-Level Data Flow

**Input Phase:**
- PAF alignment file → Parse PAF records
- FASTA contig file → Load contig sequences

**Processing Phase:**
- Reference windows → Group alignments by window
- Per-window selection → Find best contig per window
- Composite scoring → Calculate priority scores
- Greedy resolution → Select non-overlapping set

**Output Phase:**
- Write filtered FASTA file
- Generate summary log (TSV)

---

## Module Specifications

### 1. **`file_handler.py`** – Input Validation and Parsing

#### Responsibilities
- Validate and parse PAF alignment files
- Validate and parse FASTA contig files
- Perform format checks and consistency validation
- Return structured data objects for downstream processing

#### Data Structures

**`PAFRecord` (NamedTuple or dataclass)**
Fields:
- query_name: str # Column 1: Query sequence identifier
- query_length: int # Column 2: Query sequence total length
- query_start: int # Column 3: 0-based query alignment start
- query_end: int # Column 4: 0-based query alignment end
- strand: str # Column 5: Relative strand (+/-)
- target_name: str # Column 6: Target (reference) sequence name
- target_length: int # Column 7: Target sequence total length
- target_start: int # Column 8: 0-based target alignment start
- target_end: int # Column 9: 0-based target alignment end
- num_matches: int # Column 10: Number of residue matches
- alignment_length: int # Column 11: Alignment block length
- mapping_quality: int # Column 12: Mapping quality (0–60)
- alignment_score: int # Alignment score (from SAM AS tag or column 12)

**`ContigRecord` (NamedTuple or dataclass)**
Fields:
- name: str # Contig sequence identifier (from FASTA header)
- sequence: str # Nucleotide sequence
- length: int # Sequence length

**`WindowRecord` (NamedTuple or dataclass)**
Fields:
- window_id: int # Sequential window number
- ref_name: str # Reference contig name
- window_start: int # 0-based window start position
- window_end: int # 0-based window end position
- alignments: List[PAFRecord] # All alignments in this window


#### Key Methods

**`validate_file(file_path: str) -> bool`**
- Verify file exists and is readable
- Return: True if valid, False otherwise
- Raise: FileNotFoundError, PermissionError

**`read_paf(file_path: str) -> List[PAFRecord]`**
- Parse PAF format file line-by-line
- Validate each record: 
  - Minimum 12 columns present
  - Coordinates valid (start < end)
  - Numeric fields are integers
- Return: List of PAFRecord objects sorted by target position
- Raise: ValueError on format error, IOError on read failure

**`read_fasta(file_path: str) -> Dict[str, ContigRecord]`**
- Parse FASTA format file
- Handle multi-line sequences (concatenate across lines)
- Validate:
  - Headers start with `>`
  - No empty sequences
  - Sequence IDs are unique
- Return: Dictionary mapping contig name to ContigRecord
- Raise: ValueError on format error, IOError on read failure

**`validate_consistency(paf_records: List[PAFRecord], contig_dict: Dict[str, ContigRecord]) -> bool`**
- Verify all PAF query names exist in FASTA file
- Check for contig length consistency between PAF column 2 and FASTA
- Return: True if consistent, False otherwise
- Raise: ValueError with detailed mismatch information

---

### 2. **`metrics.py`** – Per-Window Selection and Composite Scoring

#### Responsibilities
- Divide reference into fixed-size windows
- Group alignments by window
- Select best-scoring contig per window
- Calculate composite priority scores for all candidates
- Provide detailed metrics for each candidate

#### Data Structures

**`CandidateRecord` (NamedTuple or dataclass)**
Fields:
- window_id: int # Reference window ID
- ref_name: str # Reference contig name
- ref_start: int # Window start position
- ref_end: int # Window end position
- query_name: str # Selected contig name
- query_start: int # Query alignment start
- query_end: int # Query alignment end
- query_length: int # Query alignment length
- alignment_score: int # Alignment score
- window_coverage: float # Coverage % of reference window (0–100)
- orientation: str # Forward (+) or reverse (–)
- composite_score: float # Calculated composite priority score

**`MetricsData` (NamedTuple or dataclass)**
Fields:
- ref_name: str # Reference contig name
- ref_length: int # Reference contig total length
- window_size: int # Fixed window size in bp
- total_windows: int # Number of windows created
- candidates: List[CandidateRecord] # All per-window candidates
- max_query_length: int # Maximum query alignment length (for normalization)

#### Key Methods

**`parse_windows(ref_name: str, ref_length: int, window_size: int) -> List[WindowRecord]`**
- Divide reference into contiguous, non-overlapping windows
- Windows: [0, window_size), [window_size, 2×window_size), ..., [last, ref_length)
- Return: List of WindowRecord objects with empty alignment lists
- Validate: window_size > 0

**`group_alignments_by_window(paf_records: List[PAFRecord], windows: List[WindowRecord]) -> List[WindowRecord]`**
- For each PAF record, assign to all windows that overlap its target alignment
- An alignment overlaps window if: `max(align_start, window_start) < min(align_end, window_end)`
- Return: WindowRecord objects with populated alignment lists
- Side effect: Modifies window objects

**`select_best_contig_per_window(windows: List[WindowRecord]) -> List[CandidateRecord]`**
- For each window with alignments:
  - Find alignment with highest alignment_score
  - If tie: use secondary criterion (e.g., highest coverage, then earliest target start)
- For each window without alignments:
  - Create candidate with query_name = "NONE"
- Return: List of CandidateRecord (one per window)
- Calculate metrics:
  - `window_coverage = (min(align_end, window_end) - max(align_start, window_start)) / (window_end - window_start) × 100`
  - `orientation = "+" if strand == "+" else "−"`

**`calculate_composite_scores(candidates: List[CandidateRecord], coverage_weight: float) -> List[CandidateRecord]`**
- Find max_query_length among candidates
- For each candidate:
  - If query_name == "NONE": composite_score = 0
  - Else: `composite_score = alignment_score + (coverage_weight × (query_length / max_query_length))`
- Return: Updated candidate list with composite_score populated
- Validate: coverage_weight ≥ 0.0

**`get_candidate_metrics(candidate: CandidateRecord) -> Dict[str, Any]`**
- Return: Dictionary with candidate metadata for logging/debugging
- Fields: window_id, query_name, alignment_score, composite_score, window_coverage, query_length

---

### 3. **`evaluation.py`** – Overlap Detection and Greedy Resolution

#### Responsibilities
- Detect overlaps between contigs on query sequence
- Implement greedy overlap resolution algorithm
- Track discard decisions and reasons
- Return final non-redundant contig set

#### Data Structures

**`OverlapInfo` (NamedTuple or dataclass)**
Fields:
- contig_a: str # First contig name
- contig_b: str # Second contig name
- overlap_start: int # Overlap region start
- overlap_end: int # Overlap region end
- overlap_length: int # Length of overlap

**`DiscardDecision` (NamedTuple or dataclass)**
Fields:
- discarded_contig: str # Contig that was discarded
- kept_contig: str # Contig that was kept
- reason: str # Reason for discard
- overlap_info: Optional[OverlapInfo] # Overlap details if applicable

#### Key Methods

**`detect_overlap(contig_a: CandidateRecord, contig_b: CandidateRecord) -> Optional[OverlapInfo]`**
- Check if two contigs overlap on query sequence
- Overlap condition: `max(a.query_start, b.query_start) < min(a.query_end, b.query_end)`
- If overlapping:
  - overlap_start = `max(a.query_start, b.query_start)`
  - overlap_end = `min(a.query_end, b.query_end)`
  - overlap_length = overlap_end − overlap_start
  - Return: OverlapInfo
- Else: Return None

**`check_overlap_with_set(candidate: CandidateRecord, kept_set: List[CandidateRecord]) -> Optional[CandidateRecord]`**
- Check if candidate overlaps any contig in kept_set
- Return: First overlapping contig from kept_set, or None if no overlap
- Optimization: Sort kept_set by query position for early termination

**`greedy_select_contigs(candidates: List[CandidateRecord]) -> Tuple[List[CandidateRecord], List[DiscardDecision]]`**
- Input: Candidates sorted by composite_score (descending)
- Algorithm:
  ```
  kept_set = []
  discard_decisions = []
  
  for each candidate:
    if candidate.query_name == "NONE":
      continue  # Skip no-alignment windows
    
    overlapping_contig = check_overlap_with_set(candidate, kept_set)
    
    if overlapping_contig is None:
      kept_set.append(candidate)
    else:
      discard_decisions.append(DiscardDecision(
        discarded_contig = candidate.query_name,
        kept_contig = overlapping_contig.query_name,
        reason = "overlaps_selected_contig",
        overlap_info = detect_overlap(candidate, overlapping_contig)
      ))
  
  return kept_set, discard_decisions
  ```
- Return: Tuple of (kept candidates, discard decisions list)

**`resolve_overlaps(all_candidates: List[CandidateRecord]) -> Tuple[List[CandidateRecord], List[DiscardDecision]]`**
- Sort candidates by composite_score (descending)
- Call: greedy_select_contigs(sorted_candidates)
- Return: Final non-redundant set and discard decisions

---

### 4. **`output_handler.py`** – Output File Generation and Logging

#### Responsibilities
- Write selected contigs to FASTA file
- Generate summary log (TSV format)
- Format output according to specifications
- Handle file I/O and error conditions

#### Key Methods

**`write_fasta(selected_contigs: Dict[str, str], output_path: str, contig_order: List[str]) -> None`**
- Write selected contigs to FASTA file
- Order: Sorted by query start position (use contig_order list)
- Format: Standard FASTA
  - Header: `>contig_name`
  - Sequence: Wrapped to 80 characters per line
- Validation:
  - Output directory is writable
  - Sufficient disk space available
  - No existing file, or overwrite confirmed
- Raise: IOError, OSError on write failure

**`write_summary_log(candidates: List[CandidateRecord], discard_decisions: List[DiscardDecision], summary_path: str) -> None`**
- Write TSV summary log file
- Header: `RefWindow_Start\tRefWindow_End\tSelectedContig\tQueryStart\tQueryEnd\tQueryLength\tAlignmentScore\tCompositeScore\tWindowCoverage\tOrientation\tDiscardedContigs\tDiscardReason`
- Rows: One per window
  - For windows with selected contig: Full candidate information
  - For windows with discarded contigs: Additional discard info
  - For windows with no alignment: query_name = "NONE", other fields empty/0
- Data formatting:
  - Numeric fields: Integers for positions/lengths, floats (2 decimals) for scores/coverage
  - Discard list: Comma-separated contig names
  - Discard reason: One of: "lower_composite_priority", "overlaps_selected_contig", "no_alignment"
- Validation:
  - Output directory is writable
  - Sufficient disk space available
- Raise: IOError, OSError on write failure

**`generate_summary_stats(candidates: List[CandidateRecord], discard_decisions: List[DiscardDecision]) -> Dict[str, Any]`**
- Calculate summary statistics
- Return: Dictionary with:
  - total_windows: int
  - windows_with_alignments: int
  - windows_without_alignments: int
  - selected_contigs_count: int
  - total_selected_bp: int
  - total_discarded_contigs: int
  - average_composite_score: float
  - median_query_length: float

---

### 5. **`utils.py`** – Shared Utilities

#### Responsibilities
- CLI argument parsing and validation
- Logger initialization and configuration
- File system utilities
- Common constants and helper functions

#### Key Classes and Functions

**`CLIParser` (class)**
- Responsibilities:
  - Parse command-line arguments
  - Validate argument values
  - Generate help message
  
- Key Methods:
  - `__init__()`: Initialize parser with all arguments
  - `parse(args: List[str]) -> argparse.Namespace`: Parse arguments and return namespace
  - Methods should include:
    - Validation for file paths (existence, readability, writability)
    - Validation for numeric parameters (window-size > 0, coverage-weight ≥ 0)
    - Validation for log-level (one of: DEBUG, INFO, WARNING, ERROR)

**`LoggerConfig` (class)**
- Responsibilities:
  - Initialize logger with specified level
  - Configure console and file logging
  
- Key Methods:
  - `__init__(log_level: str, log_file: Optional[str] = None)`: Initialize logger
  - `get_logger() -> logging.Logger`: Return configured logger

**Utility Functions**

**`validate_window_size(window_size: int) -> bool`**
- Verify window_size > 0
- Return: True if valid

**`validate_coverage_weight(coverage_weight: float) -> bool`**
- Verify coverage_weight ≥ 0.0
- Return: True if valid

**`get_output_directory(output_path: str) -> str`**
- Extract directory from output file path
- Create directory if not exists
- Return: Directory path
- Raise: OSError if creation fails or path not writable

**`estimate_file_size(file_path: str) -> int`**
- Return: File size in bytes

**`check_disk_space(directory: str, required_bytes: int) -> bool`**
- Check if directory has sufficient free space
- Return: True if sufficient, False otherwise

---

### 6. **`main.py`** – Workflow Orchestration

#### Responsibilities
- Parse command-line arguments
- Orchestrate module execution
- Handle errors and logging
- Manage temporary files and cleanup

#### Main Workflow Function

**`main(args: List[str]) -> int`**
- Entry point for CLI execution
- Execution sequence:
  
  ```
  1. Parse and validate command-line arguments
  2. Initialize logger
  3. Log start of execution
  4. Try:
       a. Read and validate PAF file (file_handler.read_paf)
       b. Read and validate FASTA file (file_handler.read_fasta)
       c. Validate consistency (file_handler.validate_consistency)
       d. Infer reference contigs and lengths from PAF records
       e. Create windows for each reference contig (metrics.parse_windows)
       f. Group alignments by window (metrics.group_alignments_by_window)
       g. Select best contig per window (metrics.select_best_contig_per_window)
       h. Calculate composite scores (metrics.calculate_composite_scores)
       i. Sort candidates by composite score (descending)
       j. Resolve overlaps (evaluation.resolve_overlaps)
       k. Extract selected contig sequences
       l. Write output FASTA (output_handler.write_fasta)
       m. Write summary log if requested (output_handler.write_summary_log)
       n. Log completion statistics
       o. Return: 0 (success)
  5. Except ValidationError:
       a. Log detailed error message
       b. Print user-friendly message to stderr
       c. Return: 1 (failure)
  6. Except FileError:
       a. Log detailed error message
       b. Print user-friendly message to stderr
       c. Return: 2 (file I/O failure)
  7. Except Exception:
       a. Log full stack trace
       b. Print generic error message to stderr
       c. Return: 255 (unexpected error)
  8. Finally:
       a. Clean up temporary files
       b. Close all file handles
  ```

**`run_workflow(paf_path: str, fasta_path: str, output_path: str, window_size: int, coverage_weight: float, summary_path: Optional[str], logger: logging.Logger) -> Tuple[int, List[str]]`**
- Execute filtering workflow
- Return: Tuple of (exit_code, output_file_list)

**Error Handling Strategy**
- Use custom exception classes:
  - `ValidationError`: Invalid input format or values
  - `FileError`: File I/O problems
  - `ProcessingError`: Logic errors during processing
- Catch exceptions at module boundaries
- Provide informative error messages with context
- Use logger.exception() for full stack traces in DEBUG mode

---

## Data Flow Details

### Input Parsing Phase
1. CLI arguments → CLIParser → validated parameters
2. PAF file → FileHandler.read_paf() → List[PAFRecord]
3. FASTA file → FileHandler.read_fasta() → Dict[str, ContigRecord]
4. Consistency check → FileHandler.validate_consistency() → bool

### Windowing and Selection Phase
1. PAF records → Infer reference contigs and lengths
2. For each reference contig:
   - Create windows → List[WindowRecord]
   - Group alignments into windows → Updated WindowRecord list
   - Select best per window → List[CandidateRecord]
3. All candidates → Calculate composite scores → Updated CandidateRecord list

### Resolution and Output Phase
1. All candidates (sorted by composite score) → Greedy overlap resolution
2. Kept contigs → Extract sequences from FASTA
3. Kept contigs + decisions → Write FASTA and summary log
4. Statistics → Log and return

---

## Implementation Considerations

### Performance Optimization
- **Alignment Indexing**: Use dictionary (target_name → alignments) for fast window lookup
- **Overlap Detection**: Sort candidates by query position; use two-pointer technique for O(n) worst-case
- **Memory Management**: Stream large PAF files if needed (Phase 1: assume in-memory)
- **File I/O**: Use buffered reading/writing for large files

### Robustness
- **Validation**: Check all inputs at module boundaries
- **Error Messages**: Provide context (file name, line number, reason)
- **Logging**: DEBUG level for trace details, INFO for major steps, WARNING/ERROR for issues
- **Testing**: Unit tests for each module; integration tests for full workflow

### User Experience
- **Help Message**: Clear description of all arguments, defaults, examples
- **Error Recovery**: Graceful failure with meaningful messages
- **Diagnostics**: Optional DEBUG logging for troubleshooting
- **Output Clarity**: Summary statistics and discard reasons logged

### Extensibility
- **Modular Design**: Each module has single responsibility
- **Data Structures**: Use standard types (lists, dicts) for easy manipulation
- **Configuration**: Parametrize window size and coverage weight
- **Plugin Points**: Future: add custom scoring functions, filter strategies

---

## Testing Strategy

### Unit Tests (in `tests/` directory)

**`test_file_handler.py`**
- Test PAF parsing: valid files, invalid formats, edge cases
- Test FASTA parsing: valid files, invalid formats, multi-line sequences
- Test consistency validation: matching/mismatching IDs

**`test_metrics.py`**
- Test window creation: boundary conditions, edge cases
- Test alignment grouping: correct window assignment
- Test per-window selection: correct best-contig identification
- Test composite scoring: correct score calculation, normalization

**`test_evaluation.py`**
- Test overlap detection: overlapping/non-overlapping pairs
- Test greedy selection: correct priority order, overlap handling
- Test discard tracking: correct reasons and contig names

**`test_output_handler.py`**
- Test FASTA writing: format compliance, sequence wrapping
- Test TSV log writing: format compliance, field values
- Test summary statistics: correct aggregation

**`test_utils.py`**
- Test CLI parsing: valid/invalid arguments
- Test logger configuration: correct setup
- Test file utilities: path handling, space checking

### Integration Tests (in `tests/test_integration.py`)
- Test full workflow with example data
- Compare output to expected files
- Test with various parameter combinations
- Test error handling and recovery

### Example Data (in `data/` directory)
- `example_alignments.paf`: Small, valid PAF file with multiple contigs and windows
- `example_contigs.fasta`: Corresponding FASTA file
- `expected_output.fasta`: Expected filtered output
- `expected_summary.log`: Expected summary log

---

## Deployment and Distribution

### Environment Setup
- **Python Version**: 3.13+
- **Dependencies**: Standard library only (no external packages)
- **Virtual Environment**: venv (virtualenv)

### Distribution Files
- `environment.yaml`: Conda specification (if using conda)
- `requirements.txt`: pip dependencies (minimal for Phase 1)
- Installation script or README instructions

### CI/CD Considerations
- Unit tests must pass
- Integration tests with example data
- Code style checks (if applicable)
- Documentation generation

---

## Future Extensibility Points

### Phase 2 Enhancements
1. **Fragmentation Scoring**: Add alignment block count/spacing analysis
2. **Structural Assessment**: Internal contig duplicate detection
3. **Coverage Analysis**: Read alignment depth metrics
4. **Dynamic Windows**: Adaptive sizing based on alignment density

### Code Organization for Growth
- Clear module boundaries allow adding new scoring modules
- Central metrics calculation point (metrics.py) for new metrics
- Evaluation module (evaluation.py) can incorporate new resolution strategies
- Output module easily extended with new report formats

### Performance Improvements
- Chunked processing for very large PAF files
- Parallel window processing with multiprocessing
- Caching of computed metrics
- Indexed structures for faster lookups

---

## Conclusion

This outline provides a complete technical specification for Phase 1 implementation. Each module has clearly defined responsibilities, data structures, and methods. The modular architecture allows for testing at multiple levels and provides clear extension points for future enhancements while maintaining a clean, maintainable codebase.

The workflow is designed to be efficient (O(n log n) dominated by sorting), scalable to reasonably large inputs (100k+ alignments), and robust in handling various input formats and error conditions.
























