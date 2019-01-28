#!/usr/bin/python2 -W ignore

from Bio import AlignIO 
import glob
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
import numpy as np
import argparse
import sys
import os

parser = argparse.ArgumentParser(usage = 'python2 filter_fasta_alns.py -o STR [-i STR -l INT -g INT -G INT -H INT -s INT -h]\n\nPython 2 required!!\n\nThis program filters aligned fastas for missing data and heterozygous sites. It first removes columns which contain a higher percentage of missing data (N or - characters) than the <max_gaps_per_site> threshold, if alignments are at least <min_length> long following this filter, it then removes sequences with a higher percentage of missing data than the <max_gaps_per_sequence> threshold and with a higher percentage of heterozygous sites than the <max_het_per_sequence> threshold. Finally, it outputs filtered alignments for all input files with at least <min_seqs> sequences remaining and produces a log file detailing the filtering in the output directory. All filters are optional, although sequences comprised completely of missing data are always removed. The output directory must be defined and the script must be run from within the input fasta directory\n\nAuthor: Owen G Osborne\n\n')
parser.add_argument("-i", "--input_extension", help="input file extension, e.g fa. Files must be in current working directory", default="fas")
parser.add_argument("-o", "--outdir", required=True, help="path to output folder (must exist)")
parser.add_argument("-l", "--min_length", help="Minimum alignment length", default=0, type=int)
parser.add_argument("-g", "--max_gaps_per_site", help="Maximum percentage missing data or gaps per site", default=100, type=int)
parser.add_argument("-G", "--max_gaps_per_sequence", help="Maximum percentage missing data or gaps per sequence", default=100, type=int)
parser.add_argument("-H", "--max_het_per_sequence", help="Maximum percentage of heterozygous sites per sequence", default=100, type=int)
parser.add_argument("-s", "--min_seqs", help="Minimum number of sequences to keep alignment", default=0, type=int)

args = parser.parse_args()

files = '{0}{1}'.format("*.",args.input_extension)
outdir = args.outdir
if not outdir.endswith(os.path.sep):
    outdir += os.path.sep
minLen = args.min_length
maxpcGap = args.max_gaps_per_sequence
minSeqs = args.min_seqs
maxpcGap_s = args.max_gaps_per_site
maxpcHet = args.max_het_per_sequence
numfiles = len(glob.glob(files))
logpath = '{0}{1}'.format(outdir,"fasta_alignment_filtering.log")

line1 = '{0} fasta files found, alignments which pass filters will be printed to {1}, and a log file will be printed to {2}\n'.format(numfiles, outdir, logpath)
print(line1)
log = open(logpath,'w')
log.write("filename,prefilter_len,postfilter_len,passed_len_filter,prefilter_nSeqs,nSeqs_removed_gap_filt,nSeqs_removed_het_filt,postfilter_nSeqs,passed_nSeq_filter" + '\n')

if len(glob.glob(files)) > 0:
	nPass = 0
	for fasta in glob.glob(files):#	get seqs
	
		fastaname = fasta.split('/')[-1]#	get fasta name without path
		aln = AlignIO.read(fasta, "fasta")#	extract alignment from fasta
		seqLen1 = aln.get_alignment_length()
		nSeq1 = len(aln) 
		filt1_aln = MultipleSeqAlignment([])
		aln_a = np.array([list(rec) for rec in aln], np.character)
		aln_r = range(seqLen1-1)
		delsites = list([])
		for i in aln_r:#	for each column in alignment
			site = list(aln_a[:,i])	
			nSeq = len(site) #	count seqs
			nGap_s = float(site.count('N') + site.count('n') + site.count('-')) #	count gaps and missing data
			pcGap_s = nGap_s / nSeq * 100
			if pcGap_s > maxpcGap_s:
				delsites.append(i)
			#	if proportion of seqs in the column with missing data is above threshold, delete column
		
		aln_a = np.delete(aln_a,delsites,1)
		 
		c = 0
		for current_seq in aln:  
			filt1_aln.add_sequence(current_seq.id, ''.join(map(str, list(aln_a[c,]))))
			c+= 1	
	
		length = filt1_aln.get_alignment_length()# if length of alignment after column-wise filter is above threshold, continue
		if length >= minLen: 
			filt2_aln = MultipleSeqAlignment([])
			for current_seq in filt1_aln:#	for each sequence 
				seq = str(current_seq.seq)
				nGap = float(seq.count('N') + seq.count('n') + seq.count('-'))
				pcGap = nGap / length * 100 
				if pcGap < maxpcGap: #	if proportion of missing data in seq is below threshold, print to filtered alignment
					filt2_aln.add_sequence(current_seq.id, str(current_seq.seq)) 
			filt3_aln = MultipleSeqAlignment([])
			for current_seq in filt2_aln:#	for each sequence 
				seq = str(current_seq.seq)
				nHet = float(seq.count('W') + seq.count('w') + seq.count('S') + seq.count('s') + seq.count('M') + seq.count('m') + seq.count(	'K') + seq.count('k') + seq.count('R') + seq.count('r') + seq.count('Y') + seq.count('y'))
				nATGC = float(seq.count('A') + seq.count('T') + seq.count('G') + seq.count('C') + seq.count('a') + seq.count('t') + seq.count	('g') + seq.count('c'))
				if nATGC+nHet > 0 :
					pcHet = nHet / (nHet+nATGC) * 100 
				else :
					pcHet = 0
				if pcHet < maxpcHet: #	if proportion of non gap or missing data bases that are hets in seq is below threshold, print to 	filtered alignment
					filt3_aln.add_sequence(current_seq.id, str(current_seq.seq))
			postFiltSeqs = len(filt3_aln)		
			if postFiltSeqs >= minSeqs: #	if number of remaining sequences is above threshold, print to output file
				outpath = '{0}{1}'.format(outdir,fastaname)
				f = open (outpath,'w')
				AlignIO.write(filt3_aln,f,"fasta")
				pasN = "Y"
				nPass+= 1
			else:
				pasN = "N"
			pasLen = "Y"
			remGap = (len(aln) - len(filt2_aln)) 
			remHet = (len(filt2_aln) - len(filt3_aln))
			
	
		else:
			pasLen = "N"
			remGap = "NA"
			remHet = "NA"
			postFiltSeqs = "NA"
			pasN = "NA"
	
		fastaname = str(fastaname)
		seqLen1 = str(seqLen1)
		length = str(length)
		pasLen = str(pasLen)
		nSeq1 = str(nSeq1)
		remGap = str(remGap)
		remHet = str(remHet)
		postFiltSeqs = str(postFiltSeqs)
		pasN = str(pasN)
		log.write(fastaname + "," + seqLen1 + "," + length + "," + pasLen + "," + nSeq1 + "," + remGap + "," + remHet + "," + postFiltSeqs + 	"," + pasN + "\n")
	log.close()
	line2 = 'All done! {0} of {1} alignments passed the filters. See {2} for more information\n'.format(nPass, numfiles, logpath)
	print (line2)
else:
	print ("\nNo fasta files found, check working directory contains fasta files and --input_extension is set correctly\n")
	sys.exit()


