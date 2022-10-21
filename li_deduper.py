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
    parser.add_argument("-f", "--file", help="Absolute file path to deduped SAM file")
    parser.add_argument("-o", "--outputfile", help="Absolute file path to sorted SAM file")
    parser.add_argument("-u", "--umi", help="Name of file containing list of UMIs")
    return parser.parse_args()

args = get_args()
input_name = args.file
output_name = args.outputfile
umi_name = args.umi

# test files
input_name = "/Users/kaitlynli/Documents/UFO/bioinformatics/Bi624/Deduper-kaijli/test.sam"
umi_name = "/Users/kaitlynli/Documents/UFO/bioinformatics/Bi624/Deduper-kaijli/STL96.txt"
output_name = "/Users/kaitlynli/Documents/UFO/bioinformatics/Bi624/Deduper-kaijli/test_out.sam"

def check_strand(flag: int) -> bool:
    '''
    Takes a flag integer and determines whether 16th bit is flipped.
    If 16th bit is flipped (1), then it's the reverse strand.
    If 16th bit is not flipped (0), then it's the forward strand. 
    '''
    if (flag & 16) != 16:   # check if 16th bit is flipped
        fwd = True          # 16th bit is not flipped, 0
    else:
        fwd = False         # 16th bit is flipped, 1
    return fwd



def cigar_parse(cigar: str, fwd: bool, pos: int) -> int:
    '''
    Takes CIGAR string, which strand the read is on, and the left-most starting position.
    Looks for front and end softclipping in CIGAR string using 'S'.
    Looks for deletions from reference, skipped regions fromreference, and alignment matches using D,N,M.
    Determines which clipping value to subtract from starting position based on strand.
    Using regex, 'findall' collects list of string numbers for each CIGAR variable. 
    Turn string numbers into integers using map().
    Uses sum of lists to modify 5' starting position of reads. 
    '''
    if fwd:
        front = sum(list(map(int,re.findall(r'^(\d+)S', cigar))))  
        pos = pos - front     
    else:                                                            # reverse strand
        end = sum(list(map(int,re.findall(r'(\d+)S$', cigar))))  
        deletion = sum(list(map(int,re.findall(r'(\d+)D', cigar)))) 
        skip = sum(list(map(int,re.findall(r'(\d+)N', cigar)))) 
        match = sum(list(map(int,re.findall(r'(\d+)M', cigar)))) 
        pos = pos + end + deletion + skip + match
    # print(f"{pos=}")
    return pos



umi_set = set()
with open(umi_name, "r") as fh:
    for line in fh:
        line = line.strip()
        umi_set.add(line)

with open(output_name, "w") as fw:
    with open(input_name, "r") as fr:
        entries = dict()
        for line in fr:
            line = line.strip()
            if line[0] == "@":
                fw.write(f"{line}\n")
            else:
                record = line.split()
                qname = record[0]
                flag = int(record[1])
                chrom = record[2]
                pos = int(record[3])
                cigar = record[5]
                u = re.findall(r'([A-Z]+)$', qname)
                umi = u[0]
                if umi in umi_set:
                    fwd = check_strand(flag)
                    start = cigar_parse(cigar, fwd, pos)
                    identifier = (umi,fwd)
                    if identifier in entries.keys():
                        entries[identifier] +=1
                    else:
                        entries[identifier] = 1
                        fw.write(f"{line}\n")
print("finished")
