#!/usr/bin/env python
'''
python script to parse over blastp results and 
retain the best hits as a new file. 
'''
import re
import argparse

def get_args():
    parser = argparse.ArgumentParser(description="A program that selects the first instance of PCR duplicates from a sorted, \
        uniquely aligned SAM file and writes a new SAM file without the PCR duplicates.")
    parser.add_argument("-f", "--file", help="Absolute file path to sorted SAM file for deduping")
    parser.add_argument("-o", "--outputfile", help="Absolute file path to deduped SAM file")
    parser.add_argument("-u", "--umi", help="Name of file containing list of UMIs")
    return parser.parse_args()

args = get_args()
input_name = args.file
output_name = args.outputfile
umi_name = args.umi

# umi_name = "STL96.txt"

print("Files imported.")

def cigar_parse(flag: int, cigar: str, pos: int) -> int:
    '''
    Takes a flag integer from sam file and determines whether 16th bit is flipped.
    If 16th bit is flipped (1), then it's the reverse strand.
    If 16th bit is not flipped (0), then it's the forward strand.

    Takes CIGAR string, which strand the read is on, and the left-most starting position.
    Looks for front and end softclipping in CIGAR string using 'S'.
    Looks for deletions from reference, skipped regions from reference, and alignment matches using D,N,M.
    Determines which clipping value to subtract from starting position based on strand.
    Using regex, 'findall' collects list of string numbers for each CIGAR variable. 
    Turn string numbers into integers using map().
    Uses sum of lists to modify 5' starting position of reads. 
    '''
    if (flag & 16) != 16:                                           # forward strand
        front = sum(list(map(int,re.findall(r'^(\d+)S', cigar))))  
        pos = pos - front     
    else:                                                           # reverse strand
        end = sum(list(map(int,re.findall(r'(\d+)S$', cigar))))  
        deletion = sum(list(map(int,re.findall(r'(\d+)D', cigar)))) 
        skip = sum(list(map(int,re.findall(r'(\d+)N', cigar)))) 
        match = sum(list(map(int,re.findall(r'(\d+)M', cigar)))) 
        pos = pos + end + deletion + skip + match
    return pos

def count_uniq_chrom(chrom: str, uniq_chroms: dict) -> dict:
    if chrom in uniq_chroms.keys():
        uniq_chroms[chrom] += 1
    else:
        uniq_chroms[chrom] = 1
    return uniq_chroms

umi_set = set()
wrong_umi,removed = 0,0
uniq_chroms = dict()
with open(umi_name, "r") as fh:
    for line in fh:
        line = line.strip()
        umi_set.add(line)
print("UMIs extracted")

with open(output_name, "w") as fw:          #open file to write out deduped sam file
    with open(input_name, "r") as fr:       #open file to read in sorted, aligned sam file
        entries = dict()
        for line in fr:
            line = line.strip()
            if line[0] == "@":
                fw.write(f"{line}\n")
            else:
                record = line.split()
                qname = record[0]
                umi = re.findall(r'([A-Z]+)$', qname)[0]
                if umi in umi_set:
                    flag = int(record[1])
                    chrom = record[2]
                    pos = int(record[3])
                    cigar = record[5]
                    start = cigar_parse(flag, cigar, pos)
                    identifier = (umi,chrom,start)
                    if identifier in entries.keys():
                        removed += 1
                    else:
                        uniq_chroms = count_uniq_chrom(chrom, uniq_chroms)
                        entries[identifier] = 1
                        fw.write(f"{line}\n")
                else:
                    wrong_umi += 1
print(f"{wrong_umi=}")  
print(f"{removed=}")  
print(f"{uniq_chroms=}")                  
print("finished")
