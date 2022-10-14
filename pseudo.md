# Deduper Part 1: Pseudo Code

## Define the problem
PCR duplicates can create bias in data, from factors such as GC content and insert length, that would lead to erroneous conclusions from the downstream data analysis. 

The goal of Deduper is to remove PCR duplicates from SAM files after alignment of sequencing data to the reference genome or transcriptome. Using the SAM columns for chromosome, position, and strand (RNAME, POS, FLAG, respectively) we would be able to observe the alignment position of the reads to start identifying the duplicates. The next step would be to use the CIGAR string to identify soft clipping, an indication of whether the ends of the read are slightly changed or not.   

UMIs group all PCR duplicates together. Each starting strand in PCR has its own UMI, or at least should. The UMI should match the first portion of the sequencing read, and as sequencing goes on, the UMIs will be the first to degrade. Assuming no errors in UMI attachment, duplication, or sequencing, reads that are tagged by different UMIs but have the same sequence are thsus not PCR duplicates, but instead biological duplicates. If there are UMI errors, those could be corrected by observing UMI read abundance and trying to match them based on similarity. 

## Example SAM files
to sort by chromosome and read position  
    `grep -v '@' test.sam | sort -k3,3n -k4,4n`

or, to use samtoolos, can do  
    `samtools sort test.sam`  
but not super sure.   


## Pseudocode Algorithm

Samtools sort
```
sort by chromosome/general location (RNAME, col 3)
sort by position (POS, col 4)
sort by strand (FLAG, col 2)
    bit 4 tells us if the read is mapped or not
        remove the unmapped reads (if the bit is flipped to 1)
    bit 16 tells us if the read is from plus or minus strand        # assuming stranded kit 
        0 is plus strand            
        1 is minus strand
```

Adjust for soft clipping
```
get cigar string (CIGAR, col 6)
check for S 
    if S is found, then it takes the number before S (probs regex for number before S that stops at any other characters)
        if S is at beginning vs end of CIGAR string, take note
            adds or subtracts from mapping position depending on S position
                S strand at beginning means subtracting from 1-based leftmost position for both plus and minus strands 
```

Single-end reads
```
check for same location, position 
will be running through plus strand first
    compare sequences, note the most common
    once it hits minus strand, compare with rev comp
probably will wittle down using UMI, location, position, seq match
make new file?
```

Known UMIs
```
regex for tail of QNAME as UMI (clipped from seq read)
is UMI in known set of 96 UMIs?
does the UMI match the beginning of the read sequence?
    if not, we can check if it's soft clipped by taking tail end of UMI or however much matches in the read
    or we can check if it's just wrong
    without a matching UMI, probably toss it out
start grouping UMIs using dictionary or new file
    the most commonly seen seq read will probably be the actual seq
```

## High Level Functions
- Will maybe be using a lot of SAMTOOLs, looking for mapped reads using 'samtools flags'  
- not sure if it's worth splitting main SAM file into multiple while parsing and grouping. 
- will be using rev comp from bioinfo.py
- will probably use phred score averager from bioinfo.py if needed
- thinking of using Snakemake to integrate samtools and python scripts
- regex will be used to look for the UMIs
- bash will be used to sort SAM files. 

Don't think there's any real new high level functions that will be needed so far. 
