'''This script contains fucntions for reading fasta files and writing csv 
analysis results'''

import math
import os
from collections import OrderedDict

def version():
	print("\n\tNumbat v0.1 (12/27/17) is a package for identifying viral \
genomes and genes by quantifying\n\tk-mers and analyzing codon frequencies.")
	print("\n\tCopyright 2017 by Shawn Rupp, Varsani Lab, Biodesign Institute, \
Arizona State University.")
	print("\tThis program comes with ABSOLUTELY NO WARRANTY.\n\tThis is free \
software, and you are welcome to redistribute it under certain conditions.\n")
	quit()

#-----------------------------------------------------------------------------

def GetCodons(seq):
	# Returns list of codons
	cdef list codons = []
	cdef int i
	cdef list nuc = ["A", "T", "C", "G"]
	seq = seq.upper()
	for i in range(0, len(seq), 3):
		c = seq[i:i+3]
		if len(c) == 3:
			if c[0] in nuc and c[1] in nuc and c[2] in nuc:
				codons.append(c)
	return codons

def ReadFasta(infile, done=[], stdout=True):
	# Returns fasta file in list of lists [header, sequence]
	cdef int l = 0
	cdef int total = 0
	cdef str seq = ""
	cdef str line
	cdef str head = ""
	contigs = []
	if stdout == True:
		print("\tReading fasta input...")
	# Read fasta entries into bins
	with open(infile, "r") as fasta:
		for line in fasta:
			line = line.strip()
			if line:
				if line[0] == ">":
					l += 1
					if seq:
						contigs.append([head, seq])
						seq = ""
					head = line[1:]
					if head in done:
						# Skip completed contigs
						head = ""
				else:
					total += len(line)
					if head:
						seq += line.upper()
	if head:
		contigs.append([head, seq])
	return contigs, l, total

def CheckOutput(outfile, header):
	# Makes output file if needed and reads in any completed queries
	cdef int first = 1
	cdef str line
	done = set()
	if os.path.isfile(outfile):
		print("\tReading previous output...")
		with open(outfile, "r") as output:
			for line in output:
				if first == 1:
					# Skip header
					first = 0
				else:
					# Save query names
					done.add(line.split(",")[0])
	else:
		print("\tGenerating new output file...")
		with open(outfile, "w") as output:
			# Initialize file and write header
			output.write(header)
	return list(done)

#-----------------------------------------------------------------------------

def AppendCSV(outfile, name, results):
	# Appends from list of lists of contig/gene results to csv
	cdef str line
	with open(outfile, "a") as output:
		for i in results:
			# Add results line by line
			line = name
			for j in i:
				line += "," + str(j)
			output.write(line + "\n")

#-----------------------------------------------------------------------------

# Initialize codon dictionary as OrderedDict
Codons = OrderedDict([("GCA",0.0), ("GCC",0.0), ("GCG",0.0), ("GCT",0.0), 
("AAC",0.0), ("AAT",0.0), ("GAC",0.0), ("GAT",0.0), ("TGC",0.0), ("TGT",0.0), ("GAA",0.0), 
("GAG",0.0), ("TTC",0.0), ("TTT",0.0), ("GGA",0.0), ("GGC",0.0), ("GGG",0.0), ("GGT",0.0), 
("CAC",0.0), ("CAT",0.0), ("ATA",0.0), ("ATC",0.0), ("ATT",0.0), ("AAA",0.0), ("AAG",0.0), 
("CTA",0.0), ("CTC",0.0), ("CTG",0.0), ("CTT",0.0), ("TTA",0.0), ("TTG",0.0), ("ATG",0.0), 
("CCA",0.0), ("CCC",0.0), ("CCG",0.0), ("CCT",0.0), ("CAA",0.0), ("CAG",0.0), ("AGA",0.0), 
("AGG",0.0), ("CGA",0.0), ("CGC",0.0), ("CGG",0.0), ("CGT",0.0), ("AGC",0.0), ("AGT",0.0), 
("TCA",0.0), ("TCC",0.0), ("TCG",0.0), ("TCT",0.0), ("ACA",0.0), ("ACC",0.0), ("ACG",0.0), 
("ACT",0.0), ("GTA",0.0), ("GTC",0.0), ("GTG",0.0), ("GTT",0.0), ("TGG",0.0), ("NNN",0.0), 
("TAC",0.0), ("TAT",0.0), ("TAA",0.0), ("TAG",0.0), ("TGA",0.0), ("Length",0.0)])
