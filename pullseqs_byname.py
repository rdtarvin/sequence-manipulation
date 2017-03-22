# -*- coding: utf-8 -*-
"""
Created on Sat Jun  7 14:27:49 2014

@author: RDT
"""

from Bio import SeqIO
import csv
import sys

def get_seqs(seqfile, seqs_to_pull, outfile):
    '''Open fasta infile and return iterator of SeqRecords with protein sequences.'''
    records = SeqIO.parse(seqfile, 'fasta')
    seqlist=[]
    matchlist=[]
    organism=seqfile.split('_')[0]
#    with open(seqs_to_pull,'rU') as f:
#	for line in f:
#		matchlist.append(line)
#        print matchlist
    matchlist=seqs_to_pull
    for rec in records:
#        for item in matchlist:
#	    print item
        if matchlist == rec.id:
	    if len(rec) > 160:
                rec.description+=' '
                rec.description+=seqfile
	        rec.id = organism+'_'+rec.id
                seqlist.append(rec)
	        print "-- %s appended" %matchlist
	    else:
		print "%s not appended because length (%d) was too short" %(matchlist, len(rec))
        else:
            continue
    
    with open(outfile, 'a') as f:
        SeqIO.write(seqlist, f, 'fasta')


def main():
	assert len(sys.argv) == 4, "usage: python pullseqs_byname.py <sequences.fasta> <sequence_names.txt> <genename>"
	infile = sys.argv[1]
	print "infile is", infile
	wanted_sequences = sys.argv[2]
	print "sequences that will be pulled are", wanted_sequences
	genename = sys.argv[3]
	outfile = genename+'.fasta'
	get_seqs(infile,wanted_sequences,outfile)
	print "appending to ", outfile
	
main()

