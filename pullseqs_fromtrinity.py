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
    with open(matches,'rU') as f:
    	reader = csv.reader(f,delimiter='\t')
        for row in reader:
        	print row
	        matchlist.append([row[0],row[-2]])
        print matchlist
    for rec in records:
        for item in matchlist:
            if item[0] == rec.id:
            	rec.description+=' '
            	rec.description+=item[1]
                seqlist.append(rec)
            else:
                continue
    
    with open(outfile, 'w') as f:
        SeqIO.write(seqlist, f, 'fasta')


def main():
	assert len(sys.argv) == 3, "usage: python pullseqs_fromtrinity.py <sequences.fasta> <blast_hits.txt>"
	infile = sys.argv[1]
	print "infile is", infile
	blast_hits = sys.argv[2]
	print "blast file is", blast_hits
	outfile = blast_hits[:blast_hits.index('.')]+'.fasta'
	get_seqs(infile,blast_hits,outfile)
	print "outfile is", outfile
	
main()

