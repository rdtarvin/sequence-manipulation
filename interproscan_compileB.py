
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
name=input_file.strip('.txt')+'_logfileB.txt'
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
	print >>logfile,"%s was parsed in %.3fs." %(input_file,(time.clock() - start))
# 	print >>logfile, "### Parsed file is:\n",d,"\n"
    return d

def rowfinder(seqlist,parsed_input_file,ID_column): 
    '''This function takes three positional arguments:
    	1. a sequence list, the output created by makelist()
    	2. the parsed database created by parsefile()
    	3. and the column from the input file with the GO identifiers, here 4
    It returns a two dictionaries:
    	1. dict_perframe: protein identifiers (values) for each frame/length sequence combination (key)
    	2. dict_perseq: protein identifiers (values) for each sequence (key), with frames combined for each seq
    It then prints a file ('GOidentifier_perframe.csv') with identifiers listed for each sequence frame.'''
    
    # creates a new dictionary with keys as a 3-item tuple (sequence ID, frame, sequence length)
    # and values as GO and protein identifiers
    start = time.clock()
    with open(seqlist,'r') as f:
    	seqlist=f.read()
    dict_perframe={}
    for seq in seqlist: #for each IDd sequence in the input file
		print seq
		IDlist=[]
		for line in parsed_input_file: 
			if seq in line[0]: #find all lines in the parsed database matching the sequence
				IDlist.append(line[ID_column]) #pulls the identifiers column, adds to a list
				length=line[2] #saves sequences length value
			dict_perframe[(seq[:seq.rfind('_')],seq[seq.rfind('_')+1:],length)]=IDlist #creates new dictionary entry

	# prints to new file each sequence, frame, length, and GO / protein numbers
    
    # reduces above dictionary to have one entry per sequence
    dict_perseq={}
    for key in dict_perframe:
    	if key[0] not in dict_perseq: # adds a new dictionary entry if not already existent
    		assert key[0] not in dict_perseq.keys(), "%s is already in the dictionary!" %key[0]
    		dict_perseq[key[0]]=list(set(dict_perframe[key])) # ensures that there are no repeat identifiers
    	else: # otherwise adds new unique values to the key
    		for value in dict_perframe[key]:
    			if value not in dict_perseq[key[0]]:
    				dict_perseq[key[0]].append(value)
    print >>logfile,"Dictionaries created and identifier file printed in %.3fs." %(time.clock() - start)
#     print >>logfile, "### Dictionary is:\n",dict_perframe,"\n"
#     print >>logfile, "### New dictionary is:\n",dict_perseq,"\n"
	
    return dict_perframe,dict_perseq   
    
def check_dicts(perframe,perseq):
	"""This function just checks that all of the sequences are in both dictionaries.
	"""
	start = time.clock()
	framelist=[]
	seqlist=[]
	for item in perframe:
		framelist.append(item[0])
	framelist=list(set(framelist))
	for item in perseq:
		seqlist.append(item)
	seqlist=list(set(seqlist))
	assert len(framelist) == len(seqlist), "not the same length"
	for item in seqlist:
		if item in framelist:
			continue
		else:
			print >>logfile, "%s not in framelist" %item
			break
	for item in framelist:
		if item in seqlist:
			continue
		else:
			print >>logfile, "%s not in seqlist" %item
			break
			
	print >>logfile,"Dictionaries checked in %.3fs." %(time.clock() - start)
	
def print_dict_perframe(dict_perframe):
	start = time.clock()
	
	output_file,short_name=shorten_filename(input_file,'_GOidentifier_perframe.csv')
	with open(output_file, 'w') as f:
		GOwriter= csv.writer(f, delimiter=',')
		for k,v in dict_perframe.items():
			k_list=list(k)
			print k_list
			for item in v:
				k_list.append(item)
			GOwriter.writerow(k_list)
	print >>logfile,"GO dictionary printed in %fs." %(time.clock() - start)

def save_dicts(perframe,perseq):
	start = time.clock()
	output_fileF,short_name=shorten_filename(input_file,'_perframe.txt')
	output_fileS,short_name=shorten_filename(input_file,'_perseq.txt')
	with open(output_fileF,'w') as fF:
		print >>fF, perframe
	with open(output_fileS,'w') as fS:
		print >>fS, perseq
    	
	print >>logfile,"Dictionaries saved in %.3fs." %(time.clock() - start)

def main():
	parsed_input_file=parsefile(input_file)
	seqlist,short_name=shorten_filename(input_file,'_seqlist.txt')
	perframe,perseq=rowfinder(seqlist,parsed_input_file,4)
	check_dicts(perframe,perseq)
	print_dict_perframe(perframe)
	save_dicts(perframe,perseq)

	end_date=(time.strftime("%d/%m/%Y"))
	end_time=(time.strftime("%H:%M:%S"))

	print >>logfile,"Start date and time was %s %s" %(start_date, start_time)
	print >>logfile,"End date and time was %s %s" %(end_date, end_time)
	logfile.close()

main()
