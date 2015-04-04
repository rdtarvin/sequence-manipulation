## this is a python script to rename sequence files from MACROGEN
## execute
## written by RDT on Tuesday March 31 2015 rdtarvin@gmail.com

import csv
import sys
import re
import os

def read_file(csv_input):
	with open(csv_input,'rU') as f:
		reader = csv.reader(f)
		d = list(reader)
	values=[]
	keys=[]
	for item in d:
		keys.append(item[0])
		values.append(item[1])
	dictionary = dict(zip(keys,values))
	return dictionary
    
def move_files(dictionary, directory):
	files = os.listdir(directory)
	ab1 = []
	cmd_list = []
	for item in files:
		if item.endswith('.ab1'):
			ab1.append(item)
	for item in ab1:
		well = item.split('_')[1]
		well = well[:well.index('.')]
		newfilename = dictionary[well] + '.ab1'
		command = 'mv ' + directory + '/' + item + ' ' + directory + '/' + newfilename
		cmd_list.append(command)
	for cmd in cmd_list:
		os.system(cmd)
		
def main():
	assert len(sys.argv) == 3, "\n correct input is: python rename_macrogen.py <csv_file> <directory> \n where <csv_file> is a csv file with two columns, well number then sample name \n and <directory> is the location of the .ab1 files from macrogen"
	dictionary = read_file(sys.argv[1])
	move_files(dictionary,sys.argv[2])

main()

	