#This Python script requires Biopython 1.51 or later
## source: https://news.open-bio.org/2009/12/14/interleaving-paired-fastq-files-with-biopython/
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import itertools
import sys

def interleave(prefix):
	#Setup variables
	file_f = prefix+"_1.fastq"
	file_r = prefix+"_2.fastq"
	file_out = prefix+"_interleaved.fastq"

	handle = open(file_out, "w")
	count = 0

	f_iter = FastqGeneralIterator(open(file_f,"rU"))
	r_iter = FastqGeneralIterator(open(file_r,"rU"))
	for (f_id, f_seq, f_q), (r_id, r_seq, r_q) in itertools.izip(f_iter,r_iter):
		assert f_id.split(' ')[0] == r_id.split(' ')[0]
		count += 2
		#Write out both reads with "/1" and "/2" suffix on ID
		handle.write("@%s/1\n%s\n+\n%s\n@%s/2\n%s\n+\n%s\n" % (f_id.split(' ')[0], f_seq, f_q, r_id.split(' ')[0], r_seq, r_q))
	handle.close()
	print "%i records written to %s" % (count, file_out)

def main():
	prefix = sys.argv[1]
	interleave(prefix)
	
main()
	
