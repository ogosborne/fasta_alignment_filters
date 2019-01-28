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
import gffutils


parser = argparse.ArgumentParser(usage = 'python2 trim_alignments_to_CDS.py -o STR -c STR [-i STR -h]\n\nPython 2 required!!\n\nThis program trims aligned fastas to CDS defined by a gff3 file. The CDS gff file must only contain CDS (eg. by running grep "CDS" transdecoder.gff3) and contain no duplicate records. The script only outputs alignments for records in the positive orientation (TO DO, include reverse complementing of negative strand CDS alignments). The output directory must be defined and the script must be run from within the input fasta directory\nRequires Bio, numpy, gffutils,  argparse \n\nAuthor: Owen G Osborne\n\n')
parser.add_argument("-i", "--input_extension", help="fasta file input , e.g fa. Files must be in current working directory", default="fas")
parser.add_argument("-o", "--outdir", required=True, help="path to output folder (must exist)")
parser.add_argument("-c", "--cdsfile",  help="path to CDS only gff3 file", required=True, type=str)



args = parser.parse_args()
files = '{0}{1}'.format("*.",args.input_extension)
outdir = args.outdir
if not outdir.endswith(os.path.sep):
    outdir += os.path.sep
cdsfile = args.cdsfile

nfiles = len(glob.glob(files))
db = gffutils.create_db(cdsfile, dbfn=":memory:", id_spec=":seqid:", merge_strategy='error', verbose=True)
ncds = db.count_features_of_type("CDS")



line1 = "{0}{1}{2}{3}".format(ncds, " CDS found in -cdsfile, ", nfiles, " fasta files found in directory")
print line1


if len(glob.glob(files)) > 0:
	nPass = 0
	for fasta in glob.glob(files):#	get seqs
	
		fastaname = fasta.split('/')[-1]#	get fasta name without path
		alnname = fastaname.split('.')[0]
		aln = AlignIO.read(fasta, "fasta")#	extract alignment from fasta
		try:
			rec = db[alnname] # extract record from db
		except:
			rec = None # or if it isn't there, set as None
		if not rec is None: # if record exists, get strand, start and end
			strand = rec.strand 
			start = rec.start - 1
			end = rec.end

			if (end - start) % 3 == 0: # if cds ends in partial codon then trim it back
				end = end
			if (end - start) % 3 == 1:
				end = end - 1
			if (end - start) % 3 == 2:
				end = end - 2

			

			if strand == "+":  # if strand is "+", trim the alignment to the CDS start and end
				cdsaln = aln[:,start:end]
				
				alnlen = aln.get_alignment_length()
				alncdslen = cdsaln.get_alignment_length()
				rem3 = alncdslen % 3
				outpath = '{0}{1}'.format(outdir,fastaname)
				f = open (outpath,'w')
				AlignIO.write(cdsaln,f,"fasta")
				nPass+= 1

	line2 = "{0}{1}{2}".format(nPass, " CDS alignments printed to ", outdir)
	print line2
else:
	print ("\nNo fasta files found, check working directory contains fasta files and --input_extension is set correctly\n")
	sys.exit()

sys.exit()




