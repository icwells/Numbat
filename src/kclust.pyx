'''This script contains functions for counting and comparing k-mers.'''

import math
from sys import stdout
from math import floor
from random import randint
from copy import deepcopy
from glob import glob
from fastaIO import ReadFasta,AppendCSV

class KMER():
	def __init__(self, k, name, seq):
		self.K = k
		self.Name = name.replace(",", "")
		self.Seq = seq
		self.L = len(seq)
		self.Total = -1
		self.Unique = -1
		self.Kmers = []
		self.__countKmers__()
		
	def __countKmers__(self):
		# Splits sequence into kmers of codons and creates list of occurances
		cdef int i
		cdef str sub
		cdef list codons = []
		cdef int c = self.K/3
		if len(self.Seq)%self.K != 0:
			# Subset to a length divisible by k (drop remainder of dividing length by k)
			self.Seq = self.Seq[:floor(self.L/self.K)*self.K]
		# Convert to codons
		for i in range(0, len(self.Seq), 3):
			codons.append(self.Seq[i:i+3])
			i += 3
		for i in range(len(codons[:-self.K+1])):
			# Subset next k nucleotides
			sub = "".join(codons[i:i+c])
			self.Kmers.append(sub)
		self.Total = len(self.Kmers)
		self.Unique = len(set(self.Kmers))

#--------------KMeans---------------------------------------------------------

class KMeans():
	def __init__(self, k=0, n=1000, t=0, tolerance=0.05, max_iter=300):
		self.Total = t
		self.K = k
		self.N = 0
		self.__n__ = n
		self.Dir = ""
		self.Tolerance = tolerance
		self.Max = max_iter
		self.Centroids = {}
		self.Previous = {}
		self.Clusters = {}

	def __getDist__(self, x1, y1, x2, y2):
		# Caclulates Euclidean distance
		return math.sqrt((x1-y1)**2+(x2-y2)**2)

	def __avgDist__(self):
		# Returns average change in centroids
		cdef float d = 0.0
		cdef int i
		for i in self.Centroids.keys():
			d += self.__getDist__(self.Centroids[i][0], self.Previous[i][0], 
									self.Centroids[i][1], self.Previous[i][1])
		return d/self.N

	def __getCentroids__(self):
		# Returns centroid of new cluster
		cdef int l = 0
		cdef int u = 0
		cdef int idx
		cdef list avg
		# Get average of centroids
		for idx in self.Centroids.keys():
			l += self.Centroids[idx][0]
			u += self.Centroids[idx][1]
		avg = [l/self.N, u/self.N]
		for idx in self.Clusters.keys():
			length = len(self.Clusters[idx])
			if length > 0:
				l = 0
				u = 0
				for i in self.Clusters[idx]:
					l += i.L
					u += i.Unique
				self.Centroids[idx] = [l/length, u/length]
			elif length == 0:
				# Assign pseudo-random centroid by getting average of random centroid and means
				grab = self.Centroids[randint(0, self.N-1)]
				self.Centroids[idx] = [(grab[0]+avg[0])/2, (grab[1]+avg[1])/2]

	def GetCluster(self, gene, exclude=[]):
		# Returns index of closest centroid
		cdef float dist = -1.0
		cdef int idx = 0
		cdef float d
		for i in self.Centroids.keys():
			if i not in exclude:
				d = self.__getDist__(gene.L, self.Centroids[i][0], gene.Unique, self.Centroids[i][1])
				if d != 0.0:
					if d < dist or dist == -1.0:
						# Keep lower distance or assign starting distance
						dist = d
						idx = i
		return idx

	def __kMeansCluster__(self, kmers):
		# Sorts kmers into n clusters using k-means algorithm
		cdef float pd = -1.0
		cdef float dist
		cdef int idx
		cdef int i
		cdef int c
		cdef int ind
		# Assign initial centroids at random
		for idx in range(0,self.N):
			# Store centroids as list
			self.Centroids[idx] = [kmers[idx].L, kmers[idx].Unique]
		for i in range(self.Max):
			# Iterate up to maximum number of iterations
			stdout.write(("\r\tClustering k-mers. Iteration: {}").format(i+1))
			self.Clusters.clear()
			for idx in range(self.N):
				self.Clusters[idx] = []
			for gene in kmers:
				c = self.GetCluster(gene)
				self.Clusters[c].append(gene)
			# Use deep copy to avoid pointing to same dict
			self.Previous = deepcopy(self.Centroids)
			self.__getCentroids__()
			# Compare average change in distance of centroids
			dist = self.__avgDist__()
			if pd > -1.0:
				if abs(dist-pd)/pd <= self.Tolerance:
					break
				else:
					pd = dist
			else:
				pd = dist
		print("\n\tFinished clustering k-mers.")

	def Cluster(self, fasta):
		# Clusters refseqs and returns clustered kmers
		cdef list target = []
		cdef float count = 1.0
		cdef float l
		cdef list i
		l = float(len(fasta))
		# Get number of bins
		self.N = math.ceil(l/self.__n__)
		print("\tIdentifying k-mers in reference sequences...")
		for i in fasta:
			# Count kmers and add to list
			gene = KMER(self.K, i[0], i[1])
			# Erase kmers list to reduce memory load
			gene.Kmers.clear()
			target.append(gene)
			stdout.write(("\r\tCounted k-mers in {:.0%} of sequences.").format((count/l)))
			count += 1.0
		print("\n\tFinished counting k-mers.")
		self.__kMeansCluster__(target)

#-----------------I/O---------------------------------------------------------

	def WriteKmers(self, outdir):
		# Saves clustered reference squences to file
		print("\tWriting clustered sequences to file...")
		cdef str outfile
		cdef int i
		outfile = outdir + str(self.K) + "mer-Cluster{}.fa"
		with open(outdir + str(self.K) + "mer-Centroids.csv", "w") as output:
			# Record kmer length and centroids
			output.write(("TotalSequenceLength,{}\n").format(self.Total))
			output.write(("kmerLength,{}\n").format(self.K))
			for i in self.Centroids.keys():
				output.write(("{},{},{}\n").format(i, self.Centroids[i][0], self.Centroids[i][1]))
		for i in self.Clusters.keys():
			# Save each cluster to own file
			with open(outfile.format(i), "w") as output:
				for j in self.Clusters[i]:
					# Write clusters to invididual files
					output.write((">{}\n{}\n").format(j.Name, j.Seq))

	def __checkInfiles__(self, indir):
		# Esnures only one centroid file is present
		infiles = glob(indir + "*mer-Centroids.csv")
		if len(infiles) < 1:
			print("\t[ERROR] Centroid file not found. Exiting.\n")
			quit()
		elif len(infiles) > 1:
			print("\t[ERROR] Multiple Centroid files found. Exiting.\n")
		return infiles[0]

	def GetModel(self, indir):
		# Reads in clustered reference sequences
		cdef int first = 0
		cdef str line
		cdef list splt
		cdef int key
		cdef list val
		if indir[-1] != "/":
			indir += "/"
		infile = self.__checkInfiles__(indir)
		self.Dir = indir
		print("\tReading clustered centroids...")
		with open(infile, "r") as f:
			for line in f:
				line = line.strip()
				splt = line.split(",")
				if first == 2:
					key = int(splt[0])
					val = [float(splt[1]), float(splt[2])]
					self.Centroids[key] = val
				elif first == 1:
					self.K = int(splt[1])
					first = 2
				else:
					self.Total = int(splt[1])
					first = 1
		# Get number of clusters
		self.N = len(self.Centroids.keys())

#---------------Functions-----------------------------------------------------					

def pseudoAlign(t, q, l):
	# Performs local psuedo-alignment between query and target sequences
	# Returns weighted percent score, percent identical, and hit length
	cdef int score = 0
	cdef int p = 0
	cdef int idx
	cdef str kmer
	cdef str i
	cdef int hl = 0
	cdef int tl
	tl = len(t)
	for idx, kmer in enumerate(q):
		# Add 5 for hit, subtract 4 for miss
		if tl > idx:
			if t[idx] == kmer:
				# Match
				score += 5
				p += 1
				hl += 3
		elif tl > idx+1:
			if t[idx+1] == kmer:
				# Check for insertion
				score += 4
				idx += 1
				hl += 3
			elif tl > idx+2 and t[idx+2] == kmer:
				score += 3
				idx += 2
				hl += 3
		elif tl-idx >= 1:
			if t[idx-1] == kmer:
				# Check for deletions
				score += 4
				idx -= 1
				hl += 3
			elif tl-idx >= 2 and t[idx-2] == kmer:
				score += 3
				idx -= 2
				hl += 3
		else:
			# Miss
			score -= 4
	return score/(5*l), p/l, hl

def seed(t, q, end):
	# Returns first possible starting point for alignment
	cdef int idx
	cdef int ind
	for idx in range(len(q[:end])):
		for ind in range(len(t[:end])):
			if q[idx] == t[ind]:
				if q[idx+1] == t[ind+1]:
					if q[idx+2] == t[ind+2]:
						# Return indecies after 3 matches
						return idx, ind
	return -1,-1

def localAlign(alpha, total, k, cluster, query):
	# Aligns query to all targets in cluster
	cdef float score
	cdef float pid
	cdef float evalue = total/(4**k)
	cdef float e
	cdef list scores = []
	cdef list q
	for i in cluster:
		# Get 1/3 of average length as ending point to shorten loops
		end = floor((i.Total+query.Total)/6)
		minlen = query.Total - end
		e = 1.0
		t = i.Kmers
		q = query.Kmers
		while e > alpha:
			qstart, tstart = seed(t, q, end)
			if qstart >= 0:
				q = q[qstart:]
				t = t[tstart:]
				if len(q) > minlen:
					score, pid, hl = pseudoAlign(q, t, query.Total)
					if score > 0 and hl > 0:
						# Add length of uncounted nucleotides from last kmer
						hl = hl+k-1
						# Chance of getting random match of hit length with observed quality score of match
						e = evalue/(hl*score)
						if e <= alpha:
							if 0.0 < pid <= 1.0:
								scores.append([str(k), str(query.L), i.Name, str(hl), "{:.2%}".format(pid), str(e)])
								break
						else:
							# Index query to prevent endless loop
							q = q[1:]
					else:
						q = q[1:]
				else:
					# Break out of loop if sequences get too short
					break
			else:
				# Skip irresolvable sequences
				break
	return scores

def loadCluster(c, k, indir):
	# Reads in appropritate cluster file and converts entries to KMER class
	cdef list cluster = []
	cdef list fasta
	cdef str infile
	cdef list i
	infile = ("{}{}mer-Cluster{}.fa").format(indir, k, c)
	fasta,_,_ = ReadFasta(infile, stdout=False)
	for i in fasta:
		# Count kmers and add to list
		gene = KMER(k, i[0], i[1])
		# Erase seqs list to reduce memory load
		gene.Seqs = []
		cluster.append(gene)
	return cluster

def IdentifyGenes(outfile, alpha, model, fasta):
	# Assigns query to cluster and compares k-mers to all sequences in cluster
	scores = []
	ex = []
	query = KMER(model.K, fasta[0], fasta[1])
	if query.Total > 0:
		while not scores and len(ex) < 3:
			c = model.GetCluster(query, ex)
			cluster = loadCluster(c, model.K, model.Dir)
			scores = localAlign(alpha, model.Total, model.K, cluster, query)
			ex.append(c)
		AppendCSV(outfile, query.Name, scores)	
