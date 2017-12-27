'''This script will extract open reading frames from NCBI flat files and 
return a fasta file with virus name and accession, and host family, 
genus, and species.'''

import argparse
import os
from datetime import datetime
from sys import stdout
from glob import glob
from shutil import rmtree
from fastaIO import version

def getSeq(dna, line):
	# Extracts coordinates from line and subsets dna sequence
	f = True
	exons = []
	seqs = ""
	coor = []
	# Isolate coordinates and strand
	line = line.split()[-1].strip().replace("(", "")
	line = line.replace(")", "")
	line = line.replace("<", "")
	line = line.replace(">", "")
	if "complement" in line:
		f = False
		line = line.replace("complement", "")
	if "join" in line:
		line = line.replace("join", "")
		exons = line.split(",")
	if not exons:
		# Get single exon
		exons = [line]
	for i in exons:
		s = i.split("..")
		if len(s) == 2:
			start = int(s[0].replace(".", ""))
			stop = int(s[1].replace(".", ""))
			if stop < start:
				# Swap values (coordinates, etc. will be printed properly)
				s = stop
				stop = start
				start = s
			seqs += dna[start-1:stop]
			if f == False:
				coor.append("{}..{}".format(stop, start))
			else:	
				coor.append("{}..{}".format(start, stop))
	if f == False:
		# Reverse sequences
		seqs = seqs[::-1]
	if len(coor) > 1:
		c = coor[0]
		for i in coor[1:]:
			coor += "," + i
	else:
		c = coor[0]
	return seqs, c

def extract(entry):
	# Extracts coding sequences from genome
	dna = ""
	cds = []
	orfs = []
	host = "NA"
	hosts = set()
	for i,line in enumerate(entry):
		if "ORIGIN" in line:
			origin = i
			break
	# Isolate genome
	for line in entry[origin+1:]:
		line = line.strip().upper()
		line = line[line.find(" ")+1:]
		dna += line.replace(" ", "")
	# Get starting index for each gene
	for i,line in enumerate(entry[:origin]):
		if "CDS " in line:
			cds.append(i)
		elif "/host=" in line:
			# Get host name
			host = line[line.find('"')+1:line.rfind('"')]
			hosts.add(host)
			if "," in host:
				host = host[:host.find(",")]
			if host.count(" ") >= 1:
				# Keep only genus and species names
				splt = host.split()
				host = splt[0] + "-" + splt[1]
	for idx,i in enumerate(cds):
		if ".." in entry[i]:
			prod = "NA"
			pid = "NA"
			end = idx + 1
			if idx == len(cds)-1:
				end = origin
			# Subset sequences
			coor = entry[i].replace("CDS", "")
			coor = coor.strip()
			next = 1
			while coor[-1] == ",":
				# Add second line
				coor += entry[i+next].strip()
				next += 1
			seqs, coordinates = getSeq(dna, coor)
			# Only check lines before next entry
			for line in entry[i:end]:
				if "/product=" in line:
					prod = line[line.find('"')+1:line.rfind('"')].replace(" ", "-")
				elif "/protein_id=" in line:
					pid = line[line.find('"')+1:line.rfind('"')]
			# Insert tilda as placeholder for accession
			string = (">~-{}-{}-[{}]__{}\n{}\n").format(pid, prod, coordinates, host, seqs)
			orfs.append(string)
	return orfs

def getORFs(tmpdir, outfile):
	# Subsets open reading frames and write to output file
	orfs = []
	infiles = glob(tmpdir + "*")
	l = len(infiles)
	with open(outfile, "w") as output:
		for idx,ff in enumerate(infiles):
			# Read in file and extract reading frames
			acc = ff[ff.rfind("/")+1:ff.find(".")]
			with open(ff, "r") as flatfile:
				entry = flatfile.readlines()
			orfs = extract(entry)
			for i in orfs:
				# Write orfs with accession
				output.write(i.replace("~", acc))
			stdout.write(("\r\tExtracted ORFs from {} of {} entries.").format(idx+1, l))

def splitFlatFile(infile):
	# Splits flat file into one file per genome
	entry = []
	name = ""
	print("\n\tSplitting flat file...")
	# Get output directory
	outdir = infile[:infile.rfind("/")+1] + "SplitFlatFile/"
	if not os.path.isdir(outdir):
		os.mkdir(outdir)
	with open(infile, "r") as flatfile:
		for line in flatfile:
			if line.strip():
				if line.strip() != "//":
					entry.append(line)
					if line.split()[0] == "ACCESSION":
						name = line.replace("ACCESSION", "").strip()
						name = name.replace("-", "_")
						if "." in name:
							name = name[name.rfind(".")+1:]
				elif line.strip() == "//":
					# Write entry and reset
					with open(outdir + name + ".ff", "w") as output:
						for i in entry:
							output.write(i)
					entry = []
					name = ""
	return outdir

def main():
	starttime = datetime.now()
	parser = argparse.ArgumentParser(description = "This script will extract \
open reading frames from NCBI flat files in fasta format.")
	parser.add_argument("-v", action = "store_true", 
help = "Prints version info and exits.")
	parser.add_argument("-i", help = "Path to input flat file.")
	args = parser.parse_args()
	if args.v:
		version()
	name = args.i[args.i.rfind("/")+1:]
	outfile = args.i[:args.i.rfind("/")+1] + name[:name.find(".")+1] + "ORFs.fa"
	tmpdir = splitFlatFile(args.i)
	getORFs(tmpdir, outfile)
	rmtree(tmpdir)
	print(("\n\tFinished. Runtime: {}\n").format(datetime.now()-starttime)) 

if __name__ == "__main__":
	main()
