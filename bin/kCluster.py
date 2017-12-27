'''This script will perform a K-mer analysis between given fasta contigs 
and NCBI RefSeqs.'''

import argparse
import os
import sys
from datetime import datetime
from functools import partial
from multiprocessing import Pool, cpu_count
from fastaIO import *
from kclust import *

def main():
	starttime = datetime.now()
	parser = argparse.ArgumentParser(description = "kCluster will perform k-means \
clustering on a fasta file of reference sequences and perform a k-mer based local \
alignment between a fasta of input reads and clustered reference sequences.")
	parser.add_argument("-v", action = "store_true", 
help = "Prints version info and exits.")
	parser.add_argument("--cluster", action = "store_true",
help = "Quantifies k-mers in reference genes.")
	parser.add_argument("-k", type = int, default = 15,
help = "Length of k-mers (default = 15; must be a multiple of three).")
	parser.add_argument("-n", type = int, default = 1000,
help = "Approximate number of genes per cluster (default = 1000).")
	parser.add_argument("-m", help = "Path to trained cluster directory.")
	parser.add_argument("-i", help = "Path to input fasta file.")
	parser.add_argument("-o", help = "Path to output file/directory.")
	parser.add_argument("-a", type = float, default = 0.05,
help = "Alpha value for e-value (default = 0.05).")
	parser.add_argument("-t", type = int, default = 1,
help = "Number of processors to use for analysis (default = 1).")
	args = parser.parse_args()
	if args.v:
		version()
	if args.cluster:
		if args.k % 3 != 0:
			print("\n\tError: K is not a multiple of three. Exiting.")
			quit()
		# Make output directory
		if args.o[-1] != "/":
			args.o += "/"
		if not os.path.isdir(args.o):
			os.mkdir(args.o)
		print(("\n\tCounting k-mers of length {} in reference file.").format(args.k))
		fasta,_,total = ReadFasta(args.i)
		model = KMeans(args.k, args.n, total)
		model.Cluster(fasta)
		model.WriteKmers(args.o)
	else:
		print(("\n\tIdentifying genes with clustered model.").format(args.k))
		cpu = args.t
		if cpu > cpu_count():
			cpu = cpu_count()
		header = "Query,k-merLength,queryLength,Hit,HitLength,PercentIdentical,E-Value\n"
		done = CheckOutput(args.o, header)
		model = KMeans()
		model.GetModel(args.m)
		# Compile k-mers in parallel
		fasta, l, _ = ReadFasta(args.i, done)
		l = float(l)
		pool = Pool(processes = cpu)
		func = partial(IdentifyGenes, args.o, args.a, model)
		print(("\n\tIdentifying genes with {} threads....").format(cpu))
		for i,_ in enumerate(pool.imap_unordered(func, fasta), 1):
				sys.stderr.write("\r\t{0:.0%} of contigs have finished".format(i/l))
		pool.close()
		pool.join()
		print()
	print(("\tTotal runtime: {}\n").format(datetime.now()-starttime))

if __name__ == "__main__":
	main()
