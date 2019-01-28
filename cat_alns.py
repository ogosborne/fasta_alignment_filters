#!/usr/bin/python2

from Bio import AlignIO 
import glob
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
import argparse
import sys


parser = argparse.ArgumentParser(usage = 'python2 cat_alns.py -o STR -t STR [-i STR -h]\n\nPython 2 only\n\nRequires: Bio, glob, argparse\n\nThis program concatenates multiple fasta alignments with missing sequences into a single fasta alignment. You must provide the taxa you want to keep in a comma-delimited list, if they are missing from input alignments they are replaced by "N"s. The script must be run from within the input fasta directory\n\nAuthor: Owen G Osborne\n\n')
parser.add_argument("-i", "--input_extension", help="input file extension, e.g fa. Files must be in current working directory", default="fas")
parser.add_argument("-o", "--outfile", required=True, help="path to output file")
parser.add_argument("-t", "--taxa",  help="comma delimited taxon names", required=True, type=str)

args = parser.parse_args()

files = '{0}{1}'.format("*.",args.input_extension)
outfile = args.outfile
taxa = args.taxa
taxa = taxa.split(",")
numfiles = len(glob.glob(files))
numtaxa = len(taxa)

line1 = '{0} fasta files and {1} taxa found, alignments will be concatenated and written to {2}\n'.format(numfiles, numtaxa, outfile)
print(line1)

if numfiles > 0:
	cataln = MultipleSeqAlignment([])
	for taxon in taxa:
		cataln.add_sequence(taxon, "") #	make alignment with all required taxa
	for fasta in glob.glob(files):	
		fastaname = fasta.split('/')[-1] #	get fasta name without path
		aln = AlignIO.read(fasta, "fasta") #	extract alignment from fasta
		seqLen = aln.get_alignment_length() 
		newaln = MultipleSeqAlignment([])
		seq = "X"
		for catrec in cataln: # for each taxon 
			
			catid = str(catrec.id) 
			for rec in aln:
				if str(rec.id) == catid: # find sequence in fasta alignment if it's there
					seq = str(rec.seq)
			if seq == "X":
				seq = ("N" * seqLen) # if not make a sequence of Ns of the correct length
			catseq = str(catrec.seq)+seq # concatenate the previously concatenate sequences and the new sequence
			newaln.add_sequence(catid,catseq) # add to a new alignment
			seq = "X"
		cataln = newaln # update the concatenated alignment
	
	f = open (outfile,'w')
	AlignIO.write(cataln,f,"fasta") # write the concatenated alignment to outfile
	print "Finished"
	sys.exit()

else:
	print "\nNo fasta files found, check working directory contains fasta files and --input_extension is set correctly\n"
	sys.exit()

