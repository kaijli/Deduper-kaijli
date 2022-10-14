Deduper Part 1: Pseudo Code

Write up a strategy for writing a Reference Based PCR Duplicate Removal tool. That is, given a sorted sam file of uniquely mapped reads, remove all PCR duplicates (retain only a single copy of each read). Develop a strategy that avoids loading everything into memory. You should not write any code for this portion of the assignment. Be sure to:

    Define the problem
    Write examples:
        Include a properly formated sorted input sam file
        Include a properly formated expected output sam file
    Develop your algorithm using pseudocode
    Determine high level functions
        Description
        Function headers
        Test examples (for individual functions)
        Return statement

For this portion of the assignment, you should design your algorithm for single-end data, with 96 UMIs. UMI information will be in the QNAME, like so: NS500451:154:HWKTMBGXX:1:11101:15364:1139:GAACAGGT. Discard any UMIs with errors (or error correct, if you're feeling ambitious).

# Define the problem
PCR duplicates can create bias in data, from factors such as GC content and insert length, that would lead to erroneous conclusions from the downstream data analysis. 
The goal of Deduper is to remove PCR duplicates from SAM files after alignment of sequencing data to the reference genome or transcriptome. Using the SAM columns for chromosome, position, and strand (RNAME, POS, FLAG, respectively) we would be able to observe the alignment position of the reads to start identifying the duplicates. The next step would be to use the CIGAR string to identify soft clipping, an indication of whether 