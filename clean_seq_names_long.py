# -*- coding: utf-8 -*-
"""
Created on Mon Apr 14 20:27:07 2014

@author: RDT
"""

'''This script reorganizes the rec.id and rec.description of a fasta file.'''

from Bio import SeqIO
import sys

def clean_seqs(infile,genename,database):
    '''Open fasta infile and return iterator of SeqRecords with protein sequences.'''
    records = SeqIO.parse(infile, 'fasta')
    blasthit=str(infile[:-6])
    newrecords=[]
    print infile
    for rec in records:
        items=(rec.description).split(' ') # turns description into a list
        print items
	print len(items)
        rec.id = genename+'_'+items[1] # adds gene name to sequence ID
        newitem=''
#        print length
	if len(items) > 2:
        	for i in range(3,len(items)):
            		newitem='%s %s ' %(newitem,items[i]) # concatenates paths
        	if items[1][0] == '_':
            		items[1] = items[1][1:]
        	rec.description="%s %s %s %s" %(database, items[1], items[2], newitem) # rewrites description
#        print rec.id
#        print rec.description
        	newrecords.append(rec)
	else:
		rec.description="%s %s" %(database, items[1]) # rewrites description
    		newrecords.append(rec)
    outfile=genename+'_clean.fasta'
    with open(outfile,'w') as f:
        SeqIO.write(newrecords, f, 'fasta')    

if __name__ == '__main__':
    infile = sys.argv[1]
    gene = sys.argv[2]
    database= sys.argv[3]
    clean_seqs(infile,gene,database)
