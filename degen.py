#!/usr/bin/python
#author: Brian Granger
#date: 4/12/11
#this code is to provide a function to translate a nucleotide sequence with degenerate nucleotides into all it's possible forms. i.e. the sequence ACBT should result in a list of sequences: ACCT, ACGT, and ACTT.

def regen(n):
	if n in nucs:
		return [n]
	else:
		if n == 'N':
			return ['A','C','G','T'] # N is any nucleotide (placed first since most likely)
		if n == 'B':
			return ['C','G','T'] #B is not A
		if n == 'D':
			return ['A','G','T'] #D is not C
		if n == 'H':
			return ['A','C','T'] #H is not G
		if n == 'K':
			return ['G','T'] #K is... G or T
		if n == 'M':
			return ['A','C'] #M is A or C
		if n == 'R':
			return ['A',"G"] #R is A or G
		if n == "S":
			return ["C","G"] #S is C or G
		if n == "V":
			return ["A","C","G"] #V is not T, couldn't use U since that's uracil
		if n == "W":
			return ["A","T"] #W is A or T
		if n == "Y":
			return ["C","T"] #Y is C or T
		else:
			return -1

def krep(seq): #
	nseqs = ['']
	for n in seq:
		if n not in nucs:
			x = regen(n)
			if x == -1:
				return x
			nuseqs = []
		else:
			nuseqs = []
			x = [n]

		for s in nseqs:
			for u in x:
				nuseqs.append(s+u)
		nseqs = nuseqs
	return nseqs
nucs = ['A','C','G','T']

