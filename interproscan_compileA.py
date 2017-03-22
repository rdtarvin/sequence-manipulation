"""This file should be executed as such:
python interproscan_compile.py input_file.txt sequence_file.fasta blast_file.outfmt6
"""


#!/usr/bin/env python

# import pylauncher
import sys
import os
import csv
import time
from Bio import SeqIO



## global variables ##
input_file=sys.argv[1]
start_date=(time.strftime("%d/%m/%Y"))
date_formatted=time.strftime("%Y-%m-%d")
start_time=(time.strftime("%H:%M:%S"))
name=input_file.strip('.txt')+'_logfileA.txt'
logfile=open(name,'w')

def shorten_filename(input_file,new_ending):
	"""This function takes two positional arguments:
		1. the input file, with a full path
		2. the new ending for the file
	It returns a tuple including the new file without the path that starts with today's 
	date and with the new ending, as well as the shortened version of the filename."""
	
	fileshort=(input_file.split('/')[-1]).split('.')[0]
	return (date_formatted+"_"+fileshort+new_ending,fileshort)
	
	
def parsefile(input_file,sep='\t'):
    """This function takes one positional argument, the input file, and one optional argument,
    the file separator. The default is tab-delimited but can be specified for other types. The
    function parses the input file into a list of lists and returns this as a list."""
    
    start = time.clock()
    with open(input_file,'rU') as f:
        reader=csv.reader(f,delimiter=sep)
        d=list(reader)
	print >>logfile,"%s was parsed in %fs." %(input_file,(time.clock() - start))
# 	print >>logfile, "### Parsed file is:\n",d,"\n"
    return d
    
#parsefile('/Users/RDT/Documents/Research/NextGenPilot_Harvard/annotations/EtriChaBM_95_interproscan_merged.txt')

def makelist(parsed_input_file,seq_column): 
    """This function takes two positional arguments:
    	1. the parsed input file, output from the parsefile() function
    	2. the column that you want to extract, here the one with the sequence name
    It returns a list of unique items found in that column and the number of sequences in that list."""
    
    start = time.clock()
    seqlist=[]
    for line in parsed_input_file:
        if line[seq_column] not in seqlist: 
            seqlist.append(line[seq_column]) #adds unique sequences to a list
    length=len(seqlist)
    assert length == len(set(seqlist)), "sequences selected incorrectly"
#     print >>logfile, "### Sequence list is:\n",seqlist,"\n"
    
    print >>logfile,"Sequence list created in %fs." %(time.clock() - start)
    
    output_file,short_name=shorten_filename(input_file,'_seqlist.txt')
    with open(output_file,'w') as f:
    	print >>f, seqlist
#     return seqlist, length

# input_file='testfile.txt'
# parsed_input_file=parsefile(input_file)
# makelist(parsed_input_file,input_file,0)
#parsed_input_file=parsefile('2014-03-14_tcdb_blast.csv')
#(genelist,length)=makelist(parsed_input_file,12)



def main():
	parsed_input_file=parsefile(input_file)
	makelist(parsed_input_file,0)


	end_date=(time.strftime("%d/%m/%Y"))
	end_time=(time.strftime("%H:%M:%S"))

	print >>logfile,"Start date and time was %s %s" %(start_date, start_time)
	print >>logfile,"End date and time was %s %s" %(end_date, end_time)
	logfile.close()

main()


#### addendum to makelist(), to print list to file:
#     output_file=shorten_filename(input_file,'_seqlist.txt')
#     with open(output_file, 'w') as f:
#         for item in seqlist:
#         	if '\t' in item:
#         		item.strip('\t')
#         	f.writelines('%s\n' %item)


#### alternative to GOlist printing
# 	with open(output_file, 'w') as f:
# 		for k in dict_perframe:
# 			f.write("%s\t%s\t%s\t%s\t%s\n" %(short_name,k[0],k[1],k[2],'\t'.join(dict_perframe[k])))
#     f.close()