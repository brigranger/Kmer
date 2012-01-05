#! /usr/bin/python
from Bio import SeqIO
from degen import krep
from degen import regen

#author: Brian Granger
#date: 4/12/11 (for adding comments at least) Also hoping to add in partial counts for degenerate nucleotides

nucs = ['A', 'C', 'G', 'T'] #list of possible nucleotides
tetras = {} #dictionary to store count of tetramer occurances (?)
tetlist = [] #list of possible 4-mers
for n in nucs: #iterate through 4 instances of nucs. 1 for each position in the 4-mer
	for u in nucs:
		for c in nucs:
			for s in nucs:
				tetlist.append(n+u+c+s) #add the 4-mer to the list of possibles
				tetras[n+u+c+s] = 0 #put the 4-mer in the dictionary of tetramers


fh = open('contigs.fan') #opens our contig fasta (nucleotide) file
out = open('k-mer_count.txt', 'w') #open ouput file
dat = open('k-mer_dat.txt', 'w')
out.write('contig_id') #construct the header to the output file
for tet in tetlist:
	out.write('\t' + tet)
out.write('\n') # header finished
count = 0 #initialize a count
for record in SeqIO.parse(fh, 'fasta'): #have SeqIO parse the contigs file, one record at a time
	for i in range(0, len(record.seq)-4): #slide the window along the sequence
		tet = str(record.seq[i:i+4]) #grab a window of size 4
		if tet in tetras: #it it's in our dictionary (no degeneracy)
			tetras[tet] = tetras[tet]+1 #increase the count
		else: #we have degenerate nucleotides
			regenerates = krep(tet) #use regen.py to regenerate the possibly nucleotides
			for reg in regenerates:
				tetras[reg] = tetras[reg] + (1.0/len(regenerates))
			
	count = count + 1 #a record count to keep track of when to output a "progress update"
	if (count%5000 == 0): #every 1000 sequences...
		print record.id #print the name of the sequence
		print tetras #and the tetramer feature vector obtained
#save the feature vector to a file
	out.write(record.id) #sequence name first
	for t in tetlist: #then output the counts for each of the tetramers
		out.write('\t' + str(tetras[t]))
		dat.write(str(tetras[t]) + '\t')
	out.write('\n')
	dat.write('\n')
	for t in tetlist: #reset the feature vector to zeros so we can run it again on the next sequence
		tetras[t] = 0
fh.close()
