#!/usr/bin/env python
 
'''These functions can be used together to make lists of sequences for all genes
in a master file of blast outputs. The master file should be created using the
script makemasterDB.py, which concatenates a bunch of blast datasets into one
csv file. The file MUST be csv. The gene_parser.py script
must be executed from within a folder that contains the blastDatabase.txt file 
and all trinity files that correspond to that master blast file. As is, this 
program has to be executed from the python interpreter. If you rerun this file
it will replace any fasta files with the same output names.'''

import pylauncher
import sys
import csv
import time
print (time.strftime("%d/%m/%Y"))
print (time.strftime("%H:%M:%S"))

def parsefile(filename):
    '''This function parses a tab-delimited file from R into a list of lists.'''
    with open(filename,'rU') as f:
        reader=csv.reader(f)
        d=list(reader)
    return d

#parsefile('testfile.txt') #to test

def makelist(filename,column): 
    '''This function takes a parsed database and a column number and returns a 
    list of unique items found in that column. It is intended to be used to 
    make a list of unique genes in the database.'''
    genelist=[]
    for line in filename:
        #print line
        if line[column] not in genelist: 
            genelist.append(line[column])
    length=len(genelist)
    #print genelist
    date=time.strftime("%Y-%m-%d")
    output_file=date+'_genelist.txt'
    with open(output_file, 'w') as f:
        for item in genelist:
            if '\t' in item:
                item.strip('\t')
            f.writelines('%s\n' %item)
    return genelist, length
    
#filelist=parsefile('testfile.txt')
#makelist(filelist,12)
#filelist=parsefile('2014-03-14_tcdb_blast.csv')
#(genelist,length)=makelist(filelist,12)

def rowfinder(genelist,database,column): 
    '''This function takes a genelist (from makelist()), the parsed database that 
    genelist was made from, and the column where the genename was and 
    returns a dictionary containing subsetted databases (value) for each 
    gene (key)'''
    allDB={}
    for gene in genelist:
        #print gene
        newDB=[]
        for line in database:
            if gene == line[column]:
                newDB.append(line)
        #print newDB
        allDB[gene]=newDB
    print 'The genes in this file are:\n'
    for item in allDB:
        print item
    return allDB   

#filelist=parsefile('testfile.txt')
#(genelist,length)=makelist(filelist,12)
#rowfinder(genelist,filelist,12)


def makedict(allDB):
    '''This function takes a dictionary of genes (keys) and subsetted databases
    corresponding to each gene (values). It returns a dictionary that pairs 
    each gene with a list of sequences that need to be pulled from trinity files.'''
    trinities={}
    for gene in allDB: #for each gene in the file
        filePairs=[]
        for line in allDB[gene]: #for line in each subsetted database, find query and filename
            filename='Trinity_'+line[13].strip('"')+'.fasta' #writes filename of trinity file
            sequence=line[0].strip('"')
            #print filename
            #print sequence,gene
            filePairs.append((filename,sequence)) #make list of trinity and sequence pairs for each gene
        #print 'The sequences found for the gene %s are %s' %(gene,filePairs)
        #print
        trinities[gene]=filePairs #add list of pairs to dictionary, with gene name as key
        #print 'The trinities for the gene %s are %s' %(gene,trinities)
    print '\n***I have determined which sequences should be pulled for each gene.***'
    print 'The sequences found for each of the following genes are: \n'
    for item in trinities: #prints dictionary
        print item+':'
        for seqs in trinities[item]:
            print seqs
        print '%d sequences total' %(len(trinities[item]))
        print
    return trinities
        
#filelist=parsefile('testfile.txt')
#(genelist,length)=makelist(filelist,12)
#allDB=rowfinder(genelist,filelist,12)
#newDict=makedict(allDB)

from Bio import SeqIO
def pullSequences(dictionary):
    '''This function takes the dictionary which lists the filename and sequences
    to be pulled for each gene. It creates a new fasta file for each gene with 
    lists of sequences pulled from all blast files that aligned to that gene.'''
    print '***Now pulling sequences.***'    
    for gene in dictionary:
        fastaList=[]
        print '--------------------------------------------------'
        print 'Pulling sequences for %s: \n' %(gene)
        for item in dictionary[gene]: #for each sequence listed for a gene
            filename=item[0]
            sequence=item[1]
            #print '    Pulling %s from %s' %(sequence,filename)
            #print
            for record in SeqIO.parse(open(filename, 'rU'),'fasta'): #all files must be in folder
                if record.id == sequence: #find the sequence in the right file
                    record.description = record.description+' '+filename
                    fastaList.append(record) #copy to a list
        print 'The %d sequences pulled for the gene %s are:\n' %(len(fastaList),gene)
        for item in fastaList:
            print item
            print
        while '/' in gene: #this removes '/' from file names 
            loc=gene.index('/')
            gene1=gene[:loc]
            gene2=gene[(loc+1):]
            gene=gene1+'--'+gene2
            if len(gene) > 220:   #cuts file name if its too long
                gene=gene[:220]
        output_file=gene+'.fasta'
        with open(output_file, 'w') as f:
            SeqIO.write(fastaList, f, "fasta") #export to a file
            f.close()


#if __name__ == '__main__':

def main():
    print 'It seems to be working'
    filename=sys.argv[1]
    column=int(sys.argv[2])
    filelist=parsefile(filename)
    (genelist,length)=makelist(filelist,column)
    allDB=rowfinder(genelist,filelist,column)
    newDict=makedict(allDB)
    pullSequences(newDict)
    print (time.strftime("%d/%m/%Y"))
    print (time.strftime("%H:%M:%S"))

main()
