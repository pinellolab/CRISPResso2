# CRISPResso2
CRISPResso2 is a software pipeline for the analysis of genome editing experiments. It is designed to enable rapid and intuitive interpretation of results produced by amplicon sequencing. 

Briefly, CRISPResso2:
- aligns sequencing reads to a reference sequence
- quantifies insertions, mutations and deletions to determine whether a read is modified or unmodified by genome editing
- summarizes editing results in intuitive plots and datasets

## What can I do with CRISPResso2?
CRISPResso2 can be used to analyze genome editing outcomes using cleaving nucleases (e.g. Cas9 or Cpf1) or noncleaving nucleases (e.g. base editors). The following operations can be automatically performed:
- filtering of low-quality reads
- adapter trimming
- alignment of reads to one or multiple reference sequences (in the case of multiple alleles)
- quantification of HDR and NHEJ outcomes (if the HDR sequence is provided)
- quantification frameshift/inframe mutations and identification affected splice sites (if an exon sequence is provided)
- visualization of the indel distribution and position (for cleaving nucleases)
- visualization of distribution and position of substitutions (for base editors)
- visualization of alleles and their frequencies

## CRISPResso processing
#### Quality filtering
Input reads are first filtered based on the quality score (phred33) in order to remove potentially false positive indels. The filtering based on the phred33 quality score can be modulated by adjusting the optimal parameters (see additional notes below).
#### Adapter trimming
Next, adapters are trimmed from the reads. If no adapter are present, select 'No Trimming' under the 'Trimming adapter' heading in the optional parameters. If reads contain adapter sequences that need to be trimmed, select the adapters used for trimming under the ‘Trimming adapter’ heading in the optional parameters. Possible adapters include Nextera PE, TruSeq3 PE, TruSeq3 SE, TruSeq2 PE, and TruSeq2 SE. The adapters are trimmed from the reads using Trimmomatic.
#### Read merging
If paired-end reads are provided, reads are merged using FLASh . This produces a single read for alignment to the amplicon sequence, and reduces sequencing errors that may be present at the end of sequencing reads.
#### Alignment
The preprocessed reads are then aligned to the reference sequence with a global sequence alignment algorithm that takes into account our biological knowledge of nuclease function. If multiple alleles are present at the editing site, each allele can be passed to CRISPResso2 and sequenced reads will be assigned to the reference sequence or origin.
#### Visualization and analysis
Finally, after analyzing the aligned reads,	a set of informative graphs are generated, allowing for the quantification and visualization of the position and type of outcomes within the amplicon sequence.

## How is CRISPResso2 different from CRISPResso?
CRISPResso2 introduces four key innovations for the analysis of genome editing data:
1) Comprehensive analysis of sequencing data from base editors. We have added additional analysis and visualization capabilities especially for experiments using base editors.
2) Allele specific quantification of heterozygous references. If the targeted editing region has more than one allele, reads arising from each allele can be deconvoluted. 
3) A novel biologically-informed alignment algorithm. This algorithm incorporates knowledge about the mutations produced by gene editing tools to create more biologically-likely alignments.
4) Ultra-fast processing time.

## Installation
CRISPResso2 can be used via the Docker virtualization system. To run CRISPResso2, first download and install docker: https://docs.docker.com/engine/installation/
	

Then type the command:

	docker pull lucapinello/crispresso

See an example on how to run CRISPResso from a Docker image in the section **TESTING CRISPResso** below.


## Usage
CRISPResso2 is designed be run on a single amplicon. For experiments involving multiple amplicons in the same fastq, see the instructions for CRISPRessoPooled or CRISPRessoWGS below. 

CRISPresso2 requires only two parameters: input sequences in the form of fastq files, and the amplicon sequence to align to. For example: 
```CRISPResso -r1 reads.fastq.gz -a AATGTCCCCCAATGGGAAGTTCATCTGGCACTGCCCACAGGTGAGGAGGTCATGATCCCCTTCTGGAGCTCCCAACGGGCCGTGGTCTGGTTCATCATCTGTAAGAATGGCTTCAAGAGGCTCGGCTGTGGTT
```

Below is a complete list of parameters:

### Parameter list
-h or --help: show a help message and exit.

-a or --amplicon_seq: The amplicon sequence used for the experiment.

-r1 or --fastq_r1: The first fastq file.

-r2 or --fastq_r2 FASTQ_R2: The second fastq file for paired end reads.

-an or --amplicon_name: A name for the reference amplicon can be given. If multiple amplicons are given, multiple names can be specified here. (default: Reference)

-amas or --amplicon_min_alignment_score: Amplicon Minimum Alignment Score; score between 0 and 100; sequences must have at least this homology score with the amplicon to be aligned (can be comma-separated list of multiple scores, corresponding to amplicon sequences given in --amplicon_seq) After reads are aligned to a reference sequence, the homology is calculated as the number of bp they have in common. If the aligned read has a homology less than this parameter, it is discarded. This is useful for filtering erroneous reads that do not align to the target amplicon, for example arising from alternate primer locations. (default: 60)

--default_min_aln_score or --min_identity_score: Default minimum homology score for a read to align to a reference amplicon (default: 60)

--expand_ambiguous_alignments: If more than one reference amplicon is given, reads that align to multiple reference amplicons will count equally toward each amplicon. Default behavior is to exclude ambiguous alignments. (default: False)

-g or --guide_seq: sgRNA sequence, if more than one, please separate by commas. Note that the sgRNA needs to be input as the guide RNA sequence (usually 20 nt) immediately adjacent to but not including the PAM sequence (5' of NGG for SpCas9). If the PAM is found on the opposite strand with respect to the Amplicon Sequence, ensure the sgRNA sequence is also found on the opposite strand. The CRISPResso convention is to depict the expected cleavage position using the value of the parameter '--quantification_window_center' nucleotides from the 3' end of the guide. In addition, the use of alternate nucleases besides SpCas9 is supported. For example, if using the Cpf1 system, enter the sequence (usually 20 nt) immediately 3' of the PAM sequence and explicitly set the '--cleavage_offset' parameter to 1, since the default setting of -3 is suitable only for SpCas9. (default: )

-e or --expected_hdr_amplicon_seq: Amplicon sequence expected after HDR. The expected HDR amplicon sequence can be provided to quantify the number of reads showing a successful HDR repair. Note that the entire amplicon sequence must be provided, not just the donor template. CRISPResso2 will quantify identified instances of NHEJ, HDR, or mixed editing events. (default: )

-c or --coding_seq: Subsequence/s of the amplicon sequence covering one or more coding sequences for frameshift analysis. Sequences of exons within the amplicon sequence can be provided to enable frameshift analysis and splice site analysis by CRISPResso2. If more than one (for example, split by intron/s), please separate by commas. Users should provide the subsequences of the reference amplicon sequence that correspond to coding sequences (not the whole exon sequence(s)!). (default: )

-q or --min_average_read_quality: Minimum average quality score (phred33) to keep a read (default: 0)

-s or --min_single_bp_quality: Minimum single bp score (phred33) to keep a read (default: 0)

--min_bp_quality_or_N: Bases with a quality score (phred33) less than this value will be set to "N" (default: 0)

--file_prefix: File prefix for output plots and tables (default: )

-n or --name: Output name of the report (default: the names is obtained from the filename of the fastq file/s used in input) (default: )

-o or --output_folder: Output folder to use for the analysis (default: current folder)

--split_paired_end: Splits a single fastq file containing paired end reads in two files before running CRISPResso (default: False)

--trim_sequences: Enable the trimming of Illumina adapters with Trimmomatic (default: False)

--trimmomatic_options_string: Override options for Trimmomatic (default: ILLUMINACLIP:/Users/luca/anaconda/lib/python2.7/site-packages/CRISPResso-0.8.0-py2.7.egg/CRISPResso/data/NexteraPE-PE.fa:0:90:10:0:true). This parameter is useful to specify different adaptor sequences used in the experiment if you need to trim them.

--min_paired_end_reads_overlap: Parameter for the FLASH read merging step. Minimum required overlap length between two reads to provide a confident overlap. (default: None)

--max_paired_end_reads_overlap: Parameter for the FLASH merging step. Maximum overlap length expected in approximately 90% of read pairs.  Please see the FLASH manual for more information.  (default: None)

-w or --quantification_window_size or --window_around_sgrna: Defines the size of the quantification window(s) centered around the position specified by the "-- cleavage_offset" or "--quantification_window_center" parameter in relation to the provided guide RNA sequence (--sgRNA). Indels overlapping this quantification window are included in classifying reads as modified or unmodified. A value of 0 disables this window and indels in the entire amplicon are considered. (default: 1)

-wc or --quantification_window_center or --cleavage_offset: Center of quantification window to use within respect to the 3' end of the provided sgRNA sequence. Remember that the sgRNA sequence must be entered without the PAM. For cleaving nucleases, this is the predicted cleavage position. The default is -3 and is suitable for the Cas9 system. For alternate nucleases, other cleavage offsets may be appropriate, for example, if using Cpf1 this parameter would be set to 1. For base editors, this could be set to -17. (default: -3)

--exclude_bp_from_left: Exclude bp from the left side of the amplicon sequence for the quantification of the indels (default: 15)

--exclude_bp_from_right: Exclude bp from the right side of the amplicon sequence for the quantification of the indels (default: 15)

--ignore_substitutions: Ignore substitutions events for the quantification and visualization (default: False)

--ignore_insertions: Ignore insertions events for the quantification and visualization (default: False)

--ignore_deletions: Ignore deletions events for the quantification and visualization (default: False)

--discard_indel_reads: Discard reads with indels in the quantification window from analysis (default: False)

--needleman_wunsch_gap_open: Gap open option for Needleman-Wunsch alignment (default: -20)

--needleman_wunsch_gap_extend: Gap extend option for Needleman-Wunsch alignment (default: -2)

--needleman_wunsch_gap_incentive: Gap incentive value for inserting indels at cut sites (default: 1)

--needleman_wunsch_aln_matrix_loc: Location of the matrix specifying substitution scores in the NCBI format (see ftp://ftp.ncbi.nih.gov/blast/matrices/) (default: EDNAFULL)

--keep_intermediate: Keep all the intermediate files (default: False)

--dump: Dump numpy arrays and pandas dataframes to file for debugging purposes (default: False)

--plot_window_size or --offset_around_cut_to_plot: Window around quantification window center to plot.  Plots alleles centered at each guide. (default: 40)

--min_frequency_alleles_around_cut_to_plot: Minimum % reads required to report an allele in the alleles table plot. (default: 0.2)

--max_rows_alleles_around_cut_to_plot: Maximum number of rows to report in the alleles table plot. (default: 50)

--conversion_nuc_from: For base editor plots, this is the nucleotide targeted by the base editor (default: C)

--conversion_nuc_to: For base editor plots, this is the nucleotide produced by the base editor (default: T)

--base_editor_output: Outputs plots and tables to aid in analysis of base editor studies. If base editor output is selected, plots showing the frequency of substitutions in the quantification window are generated. The target and result bases can also be set to measure the rate of on-target conversion at bases in the quantification window. (default: False)

-qwc or --quantification_window_coordinates: Bp positions in the amplicon sequence specifying the quantification window. This parameter overrides values of the "--quantification_window_center", "-- cleavage_offset", "--window_around_sgrna" or "-- window_around_sgrna" values. Any indels outside this window are excluded. Ranges are separted by the dash sign like "start-stop", and multiple ranges can be separated by the underscore (\_). A value of 0 disables this filter. (can be comma-separated list of values, corresponding to amplicon sequences given in --amplicon_seq e.g. 5-10,5-10_20-30 would specify the 5th-10th bp in the first reference and the 5th-10th and 20th-30th bp in the second reference) (default: None)

--crispresso1_mode: Parameter usage as in CRISPResso 1. In particular, if this flag is set the plot window length will be twice the value 'plot_window_size' (whereas in CRISPResso2 the window length will be plot_window_size) and the old output files 'Mapping_statistics.txt', and 'Quantification_of_editing_frequency.txt' are created, and the new files 'nucleotide_frequency_table.txt' and 'substitution_frequency_table.txt' and figure 2a and 2b are suppressed, and the files 'selected_nucleotide_percentage_table.txt' are not produced hen the flag `--base_editor_output` is set (default: False)

--auto: Infer amplicon sequence from most common reads (default: False)

--debug: Show debug messages (default: False)

--no_rerun: Don't rerun CRISPResso2 if a run using the same parameters has already been finished. (default: False)

--suppress_report: Suppress output report (default: False)

## Troubleshooting
Please check that your input file(s) are in FASTQ format (compressed fastq.gz also accepted).

If you get an empty report, please double check that your amplicon sequence is correct and in the right orientation. It can be helpful to inspect the first few lines of your FASTQ file - the start of the amplicon sequence should match the start of your sequences. If not, check to see if the files are trimmed (see point below).

It is important to check if your reads are trimmed or not. CRISPResso2 assumes that the reads ARE ALREADY TRIMMED! If reads are not already trimmed, select the adapters used for trimming under the ‘Trimming Adapter’ heading under the ‘Optional Parameters’. This is FUNDAMENTAL to CRISPResso analysis. Failure to trim adaptors may result in false positives. This will result in a report where you will observe an unrealistic 100% NHEJ in the pie chart in figure 2 and a sharp peak at the edges of the reference amplicon in figure 4.

It is possible to use CRISPResso with single end reads. In this case, select ‘Single end reads’ under the ‘Experimental Design’ heading.

The quality filter assumes that your reads uses the Phred33 scale, and it should be adjusted for each user’s specific application. A reasonable value for this parameter is 30.

If your amplicon sequence is longer than your sequenced read length, the R1 and R2 reads should overlap by at least 10bp. For example, if you sequence using 150bp reads, the maximum amplicon length should be 290 bp.

Especially in repetetive regions, multiple alignments may have the best score. If you want to investigate alternate best-scoring alignments, you can view all alignments using this tool: http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Gotoh. As input, sequences from the 'Alleles_frequency_table.txt' can be used. Specifically, for a given row, the value in the 'Aligned_Sequence' should be entered into the 'Sequence a' box after removing any dashes, and the value in the 'Reference_Sequence' should be entered into the 'Sequence b' box after removing any dashes. The alternate alignments can be selected in the 'Results' panel in the Output section.

## CRISPResso2 output
The output of CRISPResso2 consists of a set of informative graphs that allow for the quantification and visualization of the position and type of outcomes within an amplicon sequence.

### Data file descriptions
*CRISPResso2_report.html* is a summary report that can be viewed in a web browser containing all of the output plots and summary statistics. 

*Alleles_frequency_table.txt* is a tab-separated text file that shows all reads and alignments to references. The first column shows the aligned sequence of the sequenced read. The second column shows the aligned sequence of the reference sequence. Gaps in each of these columns represent insertions and deletions. The next column 'Reference_Name' shows the name of the reference that the read aligned to. The fourth column, 'Read_Status' shows whether the read was modified or unmodified. The fifth through seventh columns ('n_deleted', 'n_inserted', 'n_substituted') show the number of bases deleted, inserted, and substituted as compared to the reference sequence. The eighth column shows the number of reads having that sequence, and the ninth column shows the percentage of all reads having that sequence.

*CRISPResso_mapping_statistics.txt* is a tab-delimited text file showing the number of reads in the input ('READS IN INPUTS') the number of reads after filtering, trimming and merging (READS AFTER PREPROCESSING), the number of reads aligned (READS ALIGNED) and the number of reads for which the alignment had to be computed vs read from cache. 

*CRISPResso_quantification_of_editing_frequency.txt* is a tab-delimited text file showing the number of reads aligning to each reference amplicon, as well as the status (modified/unmodified, number of insertions, deletions, and/or substitions) of those reads. 

*CRISPResso_RUNNING_LOG.txt* shows a log of the CRISPResso run.

The remainder of the files are produced for each amplicon, and each file is prefixed by the name of the amplicon.

*Alleles_frequency_table_around_cut_site_for_NNNNN.txt* is a tab-separated text file that shows alleles and alignments to the specified reference for a subsequence around the predicted cut site (here, shown by 'NNNNN'). This data report is produced for each amplicon when a guide is found in the amplicon sequence.

*quantification_window_substitution_frequency_table.txt* is a tab-separated text file that shows the frequency of substitutions in the amplicon sequence in the quantification window. The first row shows the reference sequence. The following rows show the number of substitutions to each base. For example, the first numeric value in the second row (marked ‘A’) shows the number of bases that have a substitution resulting in an A at the first basepair of the amplicon sequence. The number of unmodified bases at each position is now shown in this table (because they aren’t substitutions). Thus, if the first basepair of the amplicon sequence is an A, the first value in the first row will show 0.

*substitution_frequency_table.txt* is a tab-separated text file that shows the frequency of substitutions in the amplicon sequence across the entire amplicon. The first row shows the reference sequence. The following rows show the number of substitutions to each base. For example, the first numeric value in the second row (marked ‘A’) shows the number of bases that have a substitution resulting in an A at the first basepair of the amplicon sequence. The number of unmodified bases at each position is now shown in this table (because they aren’t substitutions). Thus, if the first basepair of the AMPLICON sequence is an A, the first value in the first row will show 0.

*insertion_histogram.txt* is a tab-separated text file that shows a histogram of the insertion sizes in the amplicon sequence in the quantification window. Insertions outside of the quantification window are not included. The ins_size column shows the insertion length, and the fq column shows the number of reads having that insertion size.

*deletion_histogram.txt* is a tab-separated text file that shows a histogram of the deletion sizes in the amplicon sequence in the quantification window. Deletions outside of the quantification window are not included. The del_size column shows length of the deletion, and the fq column shows the number of reads having that number of substitutions.

*substitution_histogram.txt* is a tab-separated text file that shows a histogram of the number of substitutions in the amplicon sequence in the quantification window. Substitutions outside of the quantification window are not included. The sub_count column shows the number of substitutions, and the fq column shows the number of reads having that number of substitutions.

*effect_vector_insertion.txt* is a tab-separated text file with a one-row header that shows the percentage of reads with an insertion at each base in the reference sequence. The first column shows the 1-based position of the amplicon, and the second column shows the percentage of reads with a insertion at that location.

*effect_vector_deletion.txt* is a tab-separated text file with a one-row header that shows the percentage of reads with a deletion at each base in the reference sequence. The first column shows the 1-based position of the amplicon, and the second column shows the percentage of reads with a deletion at that location.

*effect_vector_substitution.txt* is a tab-separated text file with a one-row header that shows the percentage of reads with a substitution at each base in the reference sequence. The first column shows the 1-based position of the amplicon, and the second column shows the percentage of reads with a substitution at that location.

*effect_vector_combined.txt* is a tab-separated text file with a one-row header that shows the percentage of reads with any modification (insertion, deletion, or substitution) at each base in the reference sequence. The first column shows the 1-based position of the amplicon, and the second column shows the percentage of reads with a modification at that location.

*modification_count_vectors.txt* is a tab-separated file showing the number of modifications for each position in the amplicon. The first row shows the amplicon sequence, and successive rows show the number of reads with insertions (row 2), insertions_left (row 3), deletions (row 4), substitutions (row 5) and the sum of all modifications (row 6). Additionally, the last row shows the number of reads aligned.

If an insertion occurs between bases 5 and 6, the insertions vector will be incremented at bases 5 and 6. However, the insertions_left vector will only be incremented at base 5 so the sum of the insertions_left row represents an accurate count of the number of insertions, whereas the sum of the insertions row will yield twice the number of insertions.

*quantification_window_modification_count_vectors.txt* is a tab-separated file showing the number of modifications for positions in the quantification window of the amplicon. The first row shows the amplicon sequence in the quantification window, and successive rows show the number of reads with insertions (row 2), insertions_left (row 3), deletions (row 4), substitutions (row 5) and the sum of all modifications (row 6). Additionally, the last row shows the number of reads aligned.

*nucleotide_frequency_table.txt* is a tab-separated file showing the number of each residue at each position in the amplicon. The first row shows the amplicon sequence, and successive rows show the number of reads with an A (row 2), C (row 3), G (row 4), T (row 5), N (row 6), or a deletion (-) (row 7) at each position.

*quantification_window_nucleotide_frequency_table.txt* is a tab-separated file showing the number of each residue at positions in the quantification window of the amplicon. The first row shows the amplicon sequence in the quantification window, and successive rows show the number of reads with an A (row 2), C (row 3), G (row 4), T (row 5), N (row 6), or a deletion (-) (row 7) at each position.

*nucleotide_percentage_table.txt* is a tab-separated file showing the percentage of each residue at each position in the amplicon. The first row shows the amplicon sequence, and successive rows show the percentage of reads with an A (row 2), C (row 3), G (row 4), T (row 5), N (row 6), or a deletion (-) (row 7) at each position.

*quantification_window_nucleotide_percentage_table.txt* is a tab-separated file showing the percentage of each residue at positions in the quantification window of the amplicon. The first row shows the amplicon sequence in the quantification window, and successive rows show the percentage of reads with an A (row 2), C (row 3), G (row 4), T (row 5), N (row 6), or a deletion (-) (row 7) at each position.

The following report files are produced when the base editor mode is enabled:

*quantification_window_selected_nucleotide_percentage_table.txt* is a tab-separated text file that shows the percentage of each base at selected nucleotides in the quantification window. If the base editing experiment targets cytosines (as set by the --base_editor_from parameter), each C in the quantification window will be numbered (e.g. C5 represents the cytosine at the 5th position in the reference). The percentage of each base at these selected target cytosines is reported, with the first row showing the numbered cytosines, and the remainder of the rows showing the percentage of each nucleotide present at these locations.

*quantification_window_selected_nucleotide_frequency_table.txt* is a tab-separated text file that shows the frequency of each base at selected nucleotides in the quantification window. If the base editing experiment targets cytosines (as set by the --base_editor_from parameter), each C in the quantification window will be numbered (e.g. C5 represents the cytosine at the 5th position in the reference). The number of reads containing each base at these selected target cytosines is reported, with the first row showing the numbered cytosines, and the remainder of the rows showing the frequency of each nucleotide present at these locations.

The following report files are produced when the amplicon contains a coding sequence:

*frameshift_analysis.txt* is a text file describing the number of noncoding, in-frame, and frameshift mutations. This report file is produced when the amplicon contains a coding sequence.

*splice_sites_analysis.txt* is a text file describing the number of splicing sites that are unmodified and modified. This file report is produced when the amplicon contains a coding sequence.

*effect_vector_insertion_noncoding.txt* is a tab-separated text file with a one-row header that shows the percentage of reads with a noncoding insertion at each base in the reference sequence. The first column shows the 1-based position of the amplicon, and the second column shows the percentage of reads with a noncoding insertion at that location. This report file is produced when amplicon contains a coding sequence.

*effect_vector_deletion_noncoding.txt* is a tab-separated text file with a one-row header that shows the percentage of reads with a noncoding deletion at each base in the reference sequence. The first column shows the 1-based position of the amplicon, and the second column shows the percentage of reads with a noncoding deletion at that location. This report file is produced when amplicon contains a coding sequence.

*effect_vector_substitution_noncoding.txt* is a tab-separated text file with a one-row header that shows the percentage of reads with a noncoding substitution at each base in the reference sequence. The first column shows the 1-based position of the amplicon, and the second column shows the percentage of reads with a noncoding substitution at that location. This report file is produced when amplicon contains a coding sequence.

