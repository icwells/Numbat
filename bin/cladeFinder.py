'''This script will quantify the frequencies of codons in gene sequences
and compare between refernece and target sequences.'''

import argparse
import os
import math
from sys import stdout
from datetime import datetime
from collections import OrderedDict
from fastaIO import *

def compareFrequencies(ref, freq):
	# Returns accession of genome with closest codon frequencies
	low = -1.0
	for i in ref.keys():
		diff = 0.0
		for j in ref[i].keys():
			if j in freq.keys():
				if freq[j] == 0 or ref[i][j] == 0:
					# Skip clades if they have exclusive codons
					if freq[j] == 0 and ref[i][j] != 0:
						diff = 100.0
						break
					elif freq[j] != 0 and ref[i][j] == 0:
						diff = 100.0
						break
				# Record absolute value of difference in frequencies
				diff += abs(freq[j]-ref[i][j])
		if diff < low:
			low = diff
			hit = i
		elif low == -1.0:
			low = diff
			hit = i
	return hit, low

def sortSpecies(fasta):
	# Calculates codon frequencies by virus species
	species = {}
	taxa = {}
	for i in fasta:
		acc = i[0].split("-")[0]
		species[i[0]] = acc
	for i in fasta:
		if i[0] in species.keys():
			# accession: [seqs]
			if species[i[0]] in taxa.keys():
				taxa[species[i[0]]].append(i[1])
			else:
				taxa[species[i[0]]] = [i[1]]
	return taxa

def identifyClades(outfile, mx, ref, fasta):
	# Quantifies codon frequencies and compares to reference
	taxa = sortSpecies(fasta)
	quant = getCodonFrequencies(taxa)
	keys = list(quant.keys())
	length = len(keys)
	for idx,i in enumerate(keys):
		# Convert to proportions
		row = quant[i]
		k = list(row.keys())
		freq = dict.fromkeys(k)
		l = row["Length"]
		if l != 0:
			for j in row.keys():
				if j != "Length":
					freq[j] = row[j]/l
			besthit, score = compareFrequencies(ref, freq)
			if 0.0 <= score <= mx:
				AppendCSV(outfile, i, [[besthit, score]])
			stdout.write(("\r\t{:.0%} of query species have finished.").format(idx/length))
	# Clear carriage return
	print()

#-----------------------------------------------------------------------------

def getRefFrequencies(infile):
	# Reads codon frequencies from file
	head = True
	ref = {}
	key = list(Codons.keys())
	key.remove("Length")
	with open(infile, "r") as freq:
		for line in freq:
			if head == False:
				splt = line.strip().split(",")
				# Start new entry (dereference by redeclaring type)
				ref[splt[0]] = OrderedDict(Codons)
				del ref[splt[0]]["Length"]
				for idx,i in enumerate(splt[1:]):
					# Add frequency to dict
					ref[splt[0]][key[idx]] = float(i)
			elif head == True:
				head = False
	return ref

def writeFrequencies(outfile, quant):
	# Calculates frequencies and writes them to file
	c = list(Codons.keys())
	print("\tWriting frequencies to file...")
	# Build header (skipping length key)
	header = "Accession"
	for i in c:
		if i != "Length":
			header += "," + i
	with open(outfile, "w") as output:
		output.write(header + "\n")
		for i in quant.keys():
			row = quant[i]
			string = i
			l = row["Length"]
			if l != 0:
				for j in row.keys():
					if j != "Length":
						# Append frequency as a proportion
						string += "," + str(row[j]/l)
				output.write(string + "\n")


def getCodonFrequencies(taxa):
	# Quantifies codon frequencies by genome accession
	quant = {}
	print("\tQuantifying codons for each clade...")
	for i in taxa.keys():
		for j in taxa[i]:
			codons = GetCodons(j)
			l = len(codons)
			if i not in quant.keys():
				# Start new entry (dereference by redeclaring type)
				quant[i] = OrderedDict(Codons)
			# Record number of codons and occurances
			quant[i]["Length"] += l
			for c in codons:
				if c in Codons.keys():
					quant[i][c] += 1
	return quant

def sortTaxa(clades, fasta):
	# Sorts input sequences by host taxonomy
	taxa = {}
	for i in fasta:
		h = i[0].split("__")[1]
		h = h.replace("-", " ")
		if h != "NA":
			if h in clades.keys():
				# Clade name: [seqs]
				if clades[h] in taxa.keys():
					taxa[clades[h]].append(i[1])
				else:
					taxa[clades[h]] = [i[1]]
	return taxa

def getTaxonomy(ref, genus):
	# Get reference assignments
	clades = {}
	# Get column number
	c = 5
	if genus == True:
		c = 6
	with open(ref, "r") as f:
		for line in f:
			splt = line.split(",")
			clades[splt[0]] = splt[c]
	return clades

#-----------------------------------------------------------------------------

def writeHosts(outfile, hosts):
	# Writes host list to file
	with open(outfile, "w") as output:
		for i in hosts:
			output.write(i)

def getHosts(infile):
	# Extracts unique host names from fasta
	hosts = set()
	with open(infile, "r") as f:
		for line in f:
			if line[0] == ">":
				h = line.split("__")[1]
				if h != "NA":
					hosts.add(h.replace("-", " "))
	return list(hosts)

def main():
	starttime = datetime.now()
	parser = argparse.ArgumentParser(description = "This script will quantify \
the frequencies of codons in gene sequences and compare between reference \
and target sequences.")
	parser.add_argument("-v", action = "store_true", 
help = "Prints version info and exits.")
	parser.add_argument("--extract", action = "store_true",
help = "Extracts unique host names from reference fasta created with extractORFs.py.")
	parser.add_argument("--quantify", action = "store_true",
help = "Quantifies codons in reference genes.")
	parser.add_argument("--genus", action = "store_true", default = False,
help = "Quantify codon frequencies at the genus level (quantifies at family level by default).")
	parser.add_argument("-r", 
help = "Path to reference taxonomy file/quantified codon frequencies file.")
	parser.add_argument("-i", help = "Path to input fasta file.")
	parser.add_argument("-o", help = "Path to output file.")
	parser.add_argument("-m", type = float, default = 0.1,
help = "Maximum allowable difference between query and hit frequencies (default = 0.1).")
	args = parser.parse_args()
	if args.v:
		version()
	elif args.extract:
		print("\n\tExtracting host names from reference...")
		hosts = getHosts(args.i)
		writeHosts(args.o, hosts)
	elif args.quantify == True:
		print("\n\tCalculating codon frequencies in reference sequences...")
		name = args.i[args.i.rfind("/")+1:]
		idx = min(name.find("_"), name.find("."))
		fasta,_,_ = ReadFasta(args.i)
		taxonomy = getTaxonomy(args.r, args.genus)
		taxa = sortTaxa(taxonomy, fasta)
		quant = getCodonFrequencies(taxa)
		writeFrequencies(args.o, quant)
	else:
		print("\n\tComparing codon frequencies between target and reference sequences...")
		header = "Query,BestHit,AverageDifference\n"
		with open(args.o, "w") as output:
			# Initialize file and write header
			output.write(header)
		fasta,_,_ = ReadFasta(args.i)
		ref = getRefFrequencies(args.r)
		identifyClades(args.o, args.m, ref, fasta)
	print(("\tTotal runtime: {}\n").format(datetime.now()-starttime))

if __name__ == "__main__":
	main()
