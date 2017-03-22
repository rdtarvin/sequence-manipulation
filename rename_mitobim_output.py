# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 2017

@author: RDT
"""

from Bio import SeqIO
import csv
import sys
import datetime

def get_seqs(seqfile, name, iter, outfile):
    '''Open fasta infile and return iterator of SeqRecords with protein sequences.'''
    records = SeqIO.parse(seqfile, 'fasta')
    seqlist=[]
    matchlist=[]
    for rec in records:
	rec.description = rec.id[:rec.id.index('_bb')]+'_'+iter+'_'+str(datetime.date.today())
        rec.id = name+'_'+rec.id[:rec.id.index('_bb')]
	seqlist.append(rec)
    
    with open(outfile, 'w') as f:
        SeqIO.write(seqlist, f, 'fasta')


def main():
	assert len(sys.argv) == 4, "usage: python rename_mitobim_output.py <sequences.fasta> <iter> <seqname>"
	infile = sys.argv[1]
	iteration = sys.argv[2]
	name= sys.argv[3]
	print "infile is", infile
	outfile = infile[:infile.index('.unpadded')]+'.unpadded.renamed.fasta'
	get_seqs(infile,name,iteration,outfile)
	print "outfile is", outfile
	
main()

