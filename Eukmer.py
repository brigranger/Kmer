#! /usr/bin/python
from Bio import SeqIO
from degen import krep
from degen import regen
import glob

#author: Brian Granger
#date: 4/12/11 (for adding comments at least) Also hoping to add in partial counts for degenerate nucleotides
#degenerate nucs should be handled by degen.py

speclass = {'AT':1 , 'CE':2 , 'DM':3, 'hs':4, 'PF':5, 'SC':6}
nucs = ['A', 'C', 'G', 'T'] #list of possible nucleotides
tetras = {} #dictionary to store count of tetramer occurances (?)
tetlist = [] #list of possible 4-mers
for n in nucs: #iterate through 4 instances of nucs. 1 for each position in the 4-mer
	for u in nucs:
		for c in nucs:
			for s in nucs:
				tetlist.append(n+u+c+s) #add the 4-mer to the list of possibles
				tetras[n+u+c+s] = 0 #put the 4-mer in the dictionary of tetramers

out = open('Euk-mer_count.txt', 'w')
dat = open('Euk-mer_dat.txt','w')
cla = open('Euk-mer_class.txt', 'w')
flist = glob.glob('./Euks/*_selected.fasta') #go in the Euks folder and grab every fna (will ignore zipped folders, but those eventually will hopefully be fnas as well$
out.write('contig_id') #construct the header to the output file
for tet in tetlist:
	out.write('\t' + tet)
out.write('\n') # header finished
count = 0 #initialize a count
for fname in flist: #for each file
        print fname
	spec = fname.split('/')
	spec = spec[len(spec)-1]
	spec = spec.split('_')
	spec = spec[0]
	spec = speclass[spec]
        fh = open(fname) #open it.
	for record in SeqIO.parse(fh, 'fasta'): #have SeqIO parse the contigs file, one record at a time
		if len(record.seq) > 999:
			for i in range(0, len(record.seq)-4): #slide the window along the sequence
				tet = str(record.seq[i:i+4]) #grab a window of size 4
				if tet in tetras: #it it's in our dictionary (no degeneracy)
					tetras[tet] = tetras[tet]+1 #increase the count
				else: #we have degenerate nucleotides
					regenerates = krep(tet) #use regen.py to regenerate the possibly nucleotides
					for reg in regenerates:
						tetras[reg] = tetras[reg] + (1.0/len(regenerates))
			
			count = count + 1 #a record count to keep track of when to output a "progress update"
			if (count%10000 == 0): #every 10000 sequences...
				print record.id #print the name of the sequence
				print tetras #and the tetramer feature vector obtained
#save the feature vector to a file
			out.write(record.id) #sequence name first
			for t in tetlist: #then output the counts for each of the tetramers
				out.write('\t' + str(tetras[t]))
				dat.write('\t' + str(tetras[t]))
			out.write('\n')
			dat.write('\n')
			cla.write(str(spec) + '\t')
			for t in tetlist: #reset the feature vector to zeros so we can run it again on the next sequence
				tetras[t] = 0
fh.close()
