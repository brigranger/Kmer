#!/usr/bin/python
#author:Brian Granger
#date: 4/12/11
#use: split the sequences in the Euks folder into sequences of a certain size (10k or 1k)
#and select a certain amount of sequences from each of the species, for feature... ization
#and use in the GSOM to test for classification

import glob
from Bio import SeqIO
length = 1000 # size of sequences to grab
flist = glob.glob('./Euks/*.fna') #go in the Euks folder and grab every fna (will ignore zipped folders, but those eventually will hopefully be fnas as well...?
for fname in flist: #for each file
	print fname
	fh = open(fname) #open it.
	out = open(fname[0:len(fname)-4] + '_selected.fasta', 'w')
	for record in SeqIO.parse(fh, 'fasta'): #grab each sequence
		seq = record.seq
		for i in range(0, len(seq), length):
			out.write('>' + record.id + '_'+fname+'#' + str(i) + '\n')
			out.write(str(seq)[i:i+length] + '\n')

	fh.close()
	out.close()

#re-look at the k-mer papers: they took THE "460,000" nonoverlapping sequences of 1000 bp, so that's what I attempted to do above.

