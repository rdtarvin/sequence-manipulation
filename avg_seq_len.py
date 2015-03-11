#!//anaconda/bin/python
""" 
	This script takes a trinity file, pulls the len=[] variable and averages the contig
	length for every component (ie all subcomponents (c) and reads (seq) for each component
	(comp). It should be run 'python avg_seq_len.py Trinity_tissue.fasta'. There is an extra 
	function to check that all reads are included in the output. This will be run only if 
	you execute 'python avg_seq_len.py Trinity_tissue.fasta test'; if 'test' is not 
	included, this will not run. It is slow so I don't recommend running this unless you
	suspect an issue.
	
	Written 11 March 2015 by Rebecca D Tarvin rdtarvin@gmail.com
"""

from Bio import SeqIO
import re
import sys
import time

def parse(trinity_file):
	'''
		Accepts a trinity assembly, returns a dictionary of component:length.
	'''
	
	print "Parsing file..."
	lengths={}
	for record in SeqIO.parse(open(trinity_file, 'rU'),'fasta'): 
		lengths[record.name]=(len(record.seq))
	return lengths


def average_contigs(length_dict, output):
	'''
		Accepts the lengths dictionary and the output filename, returns a dictionary with
		component:[lengths of all subcomponents and reads] (for the test function), and
		prints a csv file with [component name, # reads, avg length].
	'''
	
	# make a dictionary of all lengths for each unique component (comp00000:[lengths])
	print "Compiling dictionary..."
	comp_dict={}
	for item in length_dict:
		component=item.split('_')[0]
		component+='_'
		if component not in comp_dict:
			comp_dict[component]=[length_dict[item]]
		else:
			comp_dict[component].append(length_dict[item])
	
	# write average length and read number for each component to a file
	print "Printing to file..."
	with open(output,'w') as f:
		for item in comp_dict:
			average=(sum(comp_dict[item])/float(len(comp_dict[item])))
			f.write('%s,%s,%s\n' %(item[:-1],str(len(comp_dict[item])),str((round(average,2)))))
	print "Average contig lengths printed to %s" %output
	return comp_dict

# test run only when sys.argv[2] == 'test'
# THIS TEST IS TIME CONSUMING 
def output_test(trinity_file,comp_dict):
	start=time.clock()
	seqIDs=[]
	for record in SeqIO.parse(open(trinity_file,'rU'),'fasta'):
		seqIDs.append(record.id)
	for item in comp_dict:
		counts=len(re.findall(item,str(seqIDs)))
		assert counts == len(comp_dict[item]), "not equal %s: %s, %s" %(item[:-1],counts,len(comp_dict[item]))
		print "all reads accounted for %s: %s in trinity file, %s in dictionary" %(item[:-1],counts,len(comp_dict[item]))
	print "Test passed. Completed in %ss" %(time.clock()-start)



def main():
	trinity_file=sys.argv[1]
	output=trinity_file[trinity_file.index('_')+1:trinity_file.index('.')]+'_count.csv'
	length_dict=parse(trinity_file)
	comp_dict=average_contigs(length_dict,output)
	try:
		sys.argv[2] == 'test'
		output_test(trinity_file,comp_dict)
	except:
		print "No tests were run."

main()