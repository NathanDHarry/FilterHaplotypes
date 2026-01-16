FilterHaplotypes

Purpose: To filter out assembled contigs that are likely alternate haplotypes in 
genome assemblies which are made from pooled samples that are polyploid.

Problem: I have sequenced a sample of five full-sib individuals' genomes on 
the PacBio Revio genome sequencing platform. This has yielded over 20 million individual
HiFi long reads. I have assembled these reads using the Hifiasm assembler with maximum 
haplotig purging parameters, and still the assembly is highly duplicated for single-copy
orthologs. Tools like purge_dups and Redundans are not able to correct this problem.

Solution: I can use a previous reference genome from the same species to filter out contigs
in the highly duplicated hifi assembly that are likely alternate haplotigs, and then 
scaffold the remaining contigs into scaffolds using that same reference genome.

Approach: The workflow involves aligning each HiFi contig to the reference genome using minimap2, after breaking the reference genome scaffolds into contigs at stretches of N's indicating scaffolding. Each reference contig will be divided into dynamically sized windows, starting with a default length based on the median lengths of the reference and HiFi assembly contigs. These windows will then be adjusted stepwise based on the number of HiFi contig alignments per window.

To evaluate contigs, several key metrics will be computed. First, alignment metrics will be gathered, including the fragmentation score, which assesses the number of alignment blocks per reference-window pair. This score helps differentiate clean haplotypes from chimeras, SV-rich contigs, and misassemblies by examining:

- The number of alignment blocks per window.
- The spacing between alignment blocks.
- Orientation consistency across alignment blocks.

Additionally, the original HiFi reads will be aligned back to the assembly to assess read coverage depth. Internal self-alignments of HiFi contigs will be performed to detect features such as internal duplications, palindromic repeats, large internal repeats, and inversions. Contigs with such features will receive a negative score modifier.

Once all metrics are computed, each reference window will be evaluated to select the best-matching HiFi contig based on alignment metrics, read coverage depth, and internal contig structure. Overlapping contigs in the selected set will be resolved by choosing the contig with the better score based on these factors.

By integrating these analyses, particularly the fragmentation score, the workflow effectively identifies and filters out problematic contigs, enabling the selection of high-quality contigs for scaffolding.

Implementation: The project will be coded in Python, and is intended to be distributed using miniconda3 and should be easily usable from a standard Unix command line using argument flags and user-supplied input data streams. It should have comprehensive help messaging, robust input file format validations, alternative flagging detection (-h / --help), user-supplied memory limits, and post-run temp file cleanup.
