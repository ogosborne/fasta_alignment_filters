#!/usr/bin/python2

from Bio import AlignIO 
import glob
from Bio.Align import MultipleSeqAlignment
from Bio.Align import AlignInfo
from Bio.SeqRecord import SeqRecord
import numpy as np
import argparse
import sys
import os
import Bio.Alphabet.IUPAC as alpha


parser = argparse.ArgumentParser(usage = 'python2 make_consensus_fastas.py -o STR [-i STR -I FLOAT -n STR -m -h]\n\nPython 2 required!!\n\nThis program creates consensus sequences from aligned fastas.  The output directory must be defined and the script must be run from within the input fasta directory\nRequires Bio, numpy, argparse \n\nAuthor: Owen G Osborne\n\n')
parser.add_argument("-i", "--input_extension", help="input file extension, e.g fa. Files must be in current working directory", default="fas")
parser.add_argument("-o", "--outdir", required=True, help="path to output folder (must exist)")
parser.add_argument("-I", "--identity", help="Minimum sequence identity to call consensus", default=1.0, type=float)
parser.add_argument("-n", "--name",  help="name of consensus sequence", default=1.0, type=str)
parser.add_argument("-m", "--include_missing_data", help="include missing data in calculation, otherwise it will be ignored", action='store_true')



args = parser.parse_args()
files = '{0}{1}'.format("*.",args.input_extension)
outdir = args.outdir
if not outdir.endswith(os.path.sep):
    outdir += os.path.sep
identity = args.identity
name = args.name
include_missing_data = args.include_missing_data
numfiles = len(glob.glob(files))
if (0 < identity <= 1): 
	line1 = '{0} fasta files found, consensus sequences will be printed to {1}\n'.format(numfiles, outdir)
	print(line1)
else:
	print ("\n--identity must be between 0 and 1!\n")
	sys.exit()

if len(glob.glob(files)) > 0:
	nPass = 0
	for fasta in glob.glob(files):#	get seqs
	
		fastaname = fasta.split('/')[-1]#	get fasta name without path
		aln = AlignIO.read(fasta, "fasta")#	extract alignment from fasta
#		nSeq1 = len(aln)
		if include_missing_data == False:
			gap_aln = MultipleSeqAlignment([])
			for current_seq in aln:#	for each sequence 
				seq = str(current_seq.seq)
				seq = str.replace(seq, "N","-") # change Ns to gaps to make them invisible to dumb_consensus
				gap_aln.add_sequence(current_seq.id, seq)
				aln = gap_aln
		summary_align = AlignInfo.SummaryInfo(aln)
		consensus_unamb = summary_align.dumb_consensus(threshold=identity, ambiguous="N",  consensus_alpha=alpha.IUPACUnambiguousDNA, require_multiple=0)
		consensus_unamb = str(consensus_unamb)
		consensus_unamb = str.replace(consensus_unamb, "n","N")

		outseq = MultipleSeqAlignment([])
		outseq.add_sequence(name, consensus_unamb)
		outpath = '{0}{1}{2}{3}{4}'.format(outdir,"cons_sim",identity,"_",fastaname)
		f = open (outpath,'w')
		AlignIO.write(outseq,f,"fasta")

		
		
else:
	print ("\nNo fasta files found, check working directory contains fasta files and --input_extension is set correctly\n")
	sys.exit()


