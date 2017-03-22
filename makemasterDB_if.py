import pylauncher
import csv
import os
import sys

def findfiles(directory,pattern):
    '''This makes a list of all files with a given ending in a given directory'''
    allfiles=os.listdir(directory) #creates a variable that is a list of files in a directory
    files=[]
    for FILE in allfiles:    
        if pattern in FILE: #if file ends in specified pattern, adds to a list
            files.append(FILE)
    return files #returns a list of files

#files=findfiles('/Users/RDT/Documents/Research/NextGenPilot_Harvard/BLASToutput/text','.txt')
#files

import time
def organize(files):
    '''Takes a list of files and outputs one clean database.'''
    data=[]
    for FILE in files:
        FILE.strip("'") #removes ' ' so python can recognize file name
        with open(FILE,'rU') as f: #opens file
            fname=FILE
            tplace=fname.index('_') #finds location of _ in file name
            eplace=fname.index('.') #finds location of . in file name
            tissue=fname[:tplace] #finds tissue name in file name
            database=fname[(tplace+1):eplace] #finds database name in file name
            print "\nProcessing blast hits of %s from %s database\n" %(tissue,database)
            for line in f:
                line2=[]
                line=line.split('\t')  # separates row into columns      
                line[13]=str(line[13]) #makes the gene name row a string
                line2=line[0:12] #appends columns 1 to 12 to line variable
                line2.append(line[13]) #appends gene name to line variable
                line2.append(tissue) #appends tissue name
                line2.append(database) #appends database name
                data.append(line2) #adds line to master database list
    date=time.strftime("%Y-%m-%d") #finds today's date
    filename=date+'_blast.csv' #creates a file name
    with open(filename,'w') as f: #writes master database to file
        writer=csv.writer(f)
        writer.writerows(data)
    print "DONE"

if __name__ == '__main__':
    print 'It seems to be working'
    directory=sys.argv[1]
    pattern=sys.argv[2]
    files=findfiles(directory,pattern)
    organize(files)
    print (time.strftime("%d/%m/%Y"))
    print (time.strftime("%H:%M:%S"))
