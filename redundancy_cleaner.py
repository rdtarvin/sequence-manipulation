
import sys
from Bio import SeqIO
def sequence_cleaner(fasta_file,min_length=100):
    #create our hash table to add the sequences
    sequences={}
	
	#Using the biopython fasta parse we can read our fasta input
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        #Take the current sequence
        sequence=str(seq_record.seq).upper()
        #Check if the current sequence is according to the user parameters
        if len(sequence)>=min_length:
       # If the sequence passed in the test "is It clean?" and It isnt in the hash table , the sequence and Its id are going to be in the hash
            if sequence not in sequences:
                sequences[sequence]=seq_record.description
       #If It is already in the hash table, We're just gonna concatenate the ID of the current sequence to another one that is already in the hash table
            else:
                sequences[sequence]+="\t"+seq_record.description
 
 
    #Write the clean sequences
 
    #Create a file in the same directory where you ran this script
    output_file=open("clear_"+fasta_file,"w+")
    #Just Read the Hash Table and write on the file as a fasta format
    for sequence in sequences:
            output_file.write(">"+sequences[sequence]+"\n"+sequence+"\n")
    output_file.close()
 
    print "CLEAN!!!\nPlease check clear_"+fasta_file
 
#sequence_cleaner(infile)
 
if __name__ == '__main__':
    infile = sys.argv[1]
    sequence_cleaner(infile)
