# Object-Oriented Outline for FilterHaplotypes

## Overview

To create a maintainable, object-oriented (OO) codebase for `FilterHaplotypes`, the outline should be designed around **classes** and **objects** that encapsulate specific functionality. Each module and workflow component will correspond to a group of cohesive classes and methods within its dedicated namespace. This approach emphasizes separation of concerns, reusability, and scalability.

---

## Object-Oriented Design Workflow

### Suggested Workflow Updates (OO Approach)

1. **Input Validation and Preprocessing**
   - Encapsulate input parsing, validation, and preprocessing tasks in a `FileHandler` class.
   - Create specialized classes to handle the reference genome and HiFi assembly contigs.

2. **Alignment**
   - Use an `Aligner` class to coordinate external alignment tools like `minimap2` and manage alignment results.
   - Store alignment data in dedicated `Alignment` objects for processing by other classes.

3. **Metrics Computation**
   - Design classes like `MetricsCalculator` for alignment statistics, `StructuralAnalyzer` for HiFi contig features (duplications, repeats, etc.), and `CoverageAnalyzer` for read depth calculations.

4. **Evaluation and Filtering**
   - Build a `ContigEvaluator` class to manage scoring, ranking, and rejection of contigs.
   - Resolve contig overlaps using an `OverlapResolver` class to ensure high-quality scaffold input.

5. **Scaffolding**
   - Design a `Scaffolder` class to handle the assembly of selected contigs into scaffolds using the reference genome.

6. **Utilities and CLI**
   - Create a `CommandLineInterface` class for argument parsing, help message generation, and memory management.
   - Include general utilities like `Logger` for detailed logging and `TempFileManager` for intermediate file handling.

---

## Proposed Classes and Responsibilities

Here is a breakdown of the suggested OO design, organized into modules with their respective classes and descriptions.

### 1. **`file_handler.py`** – File Input/Output and Validation

#### **Classes**
- **`FileHandler`**
  - **Responsibilities:**
    - Validate input files (e.g., FASTA/FASTQ) for format/consistency.
    - Support parsing of reference and HiFi assembly inputs.
    - Handle output file generation and summary logs.

  - **Key Methods:**
    - `validate_file(file_path)`: Ensure file is in the correct format and readable.
    - `split_reference_genome(output_dir)`: Break the reference genome scaffolds into contigs.
    - `write_results(output_path)`: Export scaffolds and filtered contigs.

- **`GenomeParser`**
  - **Responsibilities:**
    - Encapsulate operations for parsing genome sequences (e.g., detecting stretches of `Ns`).
  - **Key Methods:**
    - `parse_genome(fasta_file)`: Parse genome assembly into contigs for further use.

---

### 2. **`alignment.py`** – Alignment Handling

#### **Classes**
- **`Aligner`**
  - **Responsibilities:**
    - Interface with alignment tools (e.g., `minimap2`).
    - Manage and store alignment outputs for processing.
  - **Key Methods:**
    - `run_alignment(hifi_file, ref_file, output_file)`: Align HiFi contigs to reference genome.
    - `parse_alignment(output_file)`: Extract statistics and details from alignment results.

- **`AlignmentResult`**
  - **Responsibilities:**
    - Store alignment data for HiFi contigs against reference windows in an organized structure.
  - **Attributes:**
    - `alignment_blocks`: Number and spacing of alignments.
    - `orientation`: Orientation consistency of contigs.

---

### 3. **`metrics.py`** – Metrics Computation

#### **Classes**
- **`MetricsCalculator`**
  - **Responsibilities:**
    - Compute alignment-based scores (e.g., fragmentation score).
    - Aggregate other metrics (e.g., coverage depth, structural quality).
  - **Methods:**
    - `calculate_fragmentation_score(alignment_result)`: Compute fragmentation.
    - `aggregate_metrics()`: Combine multiple metrics into a total score.

- **`CoverageAnalyzer`**
  - **Responsibilities:**
    - Calculate read coverage depth for HiFi assembly.
  - **Methods:**
    - `compute_read_coverage(hifi_reads, contig_assembly)`: Align reads and compute coverage per contig.

- **`StructuralAnalyzer`**
  - **Responsibilities:**
    - Identify internal features of contigs such as duplications, repeats, and inversions.
  - **Methods:**
    - `detect_internal_features(contig_sequence)`: Identify structural anomalies.
    - `apply_penalty()`: Apply score penalties for problematic contigs.

---

### 4. **`evaluation.py`** – Evaluation and Filtering

#### **Classes**
- **`ContigEvaluator`**
  - **Responsibilities:**
    - Manage the ranking and filtering of contigs based on computed metrics.
  - **Methods:**
    - `evaluate_contigs(scores)`: Generate final contig rankings.
    - `apply_filters()`: Remove problematic contigs.

- **`OverlapResolver`**
  - **Responsibilities:**
    - Resolve overlapping contigs by comparing scores and retaining the highest scoring one.
  - **Methods:**
    - `resolve_overlaps(contig_list)`: Iteratively resolve overlapping candidates.

---

### 5. **`scaffolding.py`** – Scaffolding

#### **Classes**
- **`Scaffolder`**
  - **Responsibilities:**
    - Assemble high-quality contigs into scaffolds using reference genome.
  - **Methods:**
    - `build_scaffold(filtered_contigs, reference_genome)`: Align and link contigs into scaffolds.
    - `output_scaffolds(output_dir)`: Write scaffolded results to output.

---

### 6. **`utils.py`** – Utilities and Helper Functions

#### **Classes**
- **`CommandLineInterface`**
  - **Responsibilities:**
    - Handle argument parsing and CLI interaction.
  - **Methods:**
    - `parse_arguments()`: Parse command-line arguments.
    - `generate_help_message()`: Provide detailed usage information.

- **`Logger`**
  - **Responsibilities:**
    - Log steps and events during workflow execution.
  - **Methods:**
    - `log(message)`: Log standard message or warnings.
    - `log_error(message)`: Log critical errors.

- **`TempFileManager`**
  - **Responsibilities:**
    - Manage temporary files during workflow execution.
  - **Methods:**
    - `create_temp_dir()`: Create a temporary working directory.
    - `cleanup()`: Remove all temporary files post-run.

---

## Revised Project Structure

### Suggested Directory Layout

- **`FilterHaplotypes/`**
  - **`src/`**
    - `file_handler.py` – Input/output operations.
    - `alignment.py` – Alignment handling.
    - `metrics.py` – Metrics computation.
    - `evaluation.py` – Contig evaluation and filtering.
    - `scaffolding.py` – Contig scaffolding logic.
    - `utils.py` – Shared helper utilities.
    - `main.py` – Workflow entry point.
  - **`tests/`**
    - Unit tests for all modules.
    - Example input/output for regression testing.
  - **`data/`**
    - Example datasets (e.g., FASTA input files).
  - **`config/`**
    - Configuration files for customizable parameters.
  - **`README.md`** – Project overview and documentation.
  - **`environment.yml`** – Dependency management for miniconda.

---

## Implementation Notes

1. **Encapsulation:** Each module handles a specific responsibility. Classes encapsulate related functionality where applicable.
2. **Extensibility:** New alignment tools, analysis methods, or metrics can be integrated by extending existing classes or adding new ones without disrupting the workflow.
3. **Inter-Class Communication:** Classes pass objects (e.g., `AlignmentResult`, `MetricsCalculator`) to share data, avoiding redundant recalculations.
4. **Error Handling:** Input validators and error logging ensure robustness and user-friendly diagnostics.

This OO structure promotes readability, flexibility, and maintainability, making the project easier to debug, test, and enhance.