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
file_name = "/Users/kaitlynli/Documents/UFO/bioinformatics/Bi624/Deduper-kaijli/sorted_test.sam"
umi_name = "/Users/kaitlynli/Documents/UFO/bioinformatics/Bi624/Deduper-kaijli/STL96.txt"

def read_umi_file(filename: str) -> set:
    '''
    read text file containing known UMIs in the reads
    ''' 
    umi_set = set()
    with open(filename, "r") as fh:
        for line in fh:
            line = line.strip()
            umi_set.add(line)
    return umi_set

with open(input_name, "r") as fh:
    header = list()
    collection = dict()
    # qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual = "","","","","","","","","","",""
    for line in fh:
        line = line.strip()
        if line[0] == "@":
            header.append(line)
        else:
            (qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual)=line.split()
            
        
    pass
    