# -*- coding: utf-8 -*-
"""
Created on Sat Jun  7 14:27:49 2014

@author: RDT
"""

from Bio import SeqIO
import csv
import sys

def get_seqs(seqfile, matches, outfile):
    '''Open fasta infile and return iterator of SeqRecords with protein sequences.'''
    records = SeqIO.parse(seqfile, 'fasta')
    seqlist=[]
    matchlist=[]
    with open(matches,'r') as f:
        for row in f:
            matchlist.append(row.strip())
        print matchlist
    for rec in records:
        for item in matchlist:
            if item in rec.id:
                seqlist.append(rec)
            else:
                continue
    
    with open(outfile, 'w') as f:
        SeqIO.write(seqlist, f, 'fasta')


def main():
	assert len(sys.argv) == 3, "usage: python pullseqs_fromtrinity.py <Sequences.fasta> <blast_hits.txt>"
	infile = sys.argv[1]
	blast_hits = sys.argv[2]
	outfile = infile[:infile.index('.')]+blast_hits[:blast_hits.index('.')]+'.fasta'
	get_seqs(infile,blast_hits,outfile)
	
main()
