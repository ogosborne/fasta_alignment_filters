#!/usr/bin/python2

from Bio import AlignIO 
from Bio.Seq import Seq
import glob
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
import numpy as np
import argparse
import sys
import os
from Bio.Alphabet import generic_dna

# get options

parser = argparse.ArgumentParser(usage = 'python2 filter_CDS_alns_for_PAML.py -o STR -t STR [-i STR -l INT -s INT -e STR -v INT]\nThis script filters CDS alignments and outputs them in phylip format or fasta while maintaining CDS frame. \nAuthor: Owen G Osborne\n\n')
parser.add_argument("-i", "--input_extension", help="input file extension, e.g fa. Files must be in current working directory", default="fas")
parser.add_argument("-o", "--outdir", required=True, help="path to output folder (must exist)")
parser.add_argument("-l", "--min_length", help="Minimum alignment length in codons", default=0, type=int)
parser.add_argument("-s", "--min_unmasked_seqs_per_site", help="Minimum sequences per codon site", default=100, type=int)
parser.add_argument("-e", "--essential_taxa",  help="comma delimited list of taxon names which are required for each site to be included, these are also the taxa which must have SNPs between them if --variants_required is set", required=False, type=str, default=None)
parser.add_argument("-v", "--variants_required", help="Minimum SNPs required between essential taxa", default=0, type=int)
parser.add_argument("-t", "--taxa",  help="comma delimited taxon names", required=True, type=str)
parser.add_argument("-f", "--output_format",  help="output format, either fasta or phylip (default)", required=False, default="phylip", type=str)

args = parser.parse_args()
files = '{0}{1}'.format("*.",args.input_extension)
taxa = args.taxa
taxa = taxa.split(",")
numtaxa = len(taxa)
numfiles = len(glob.glob(files))
esstaxa = args.essential_taxa
outfmt = args.output_format
if not esstaxa == None:
	esstaxa = esstaxa.split(",")
	esstaxaind = list([])
	for esstaxon in esstaxa:
		esstaxaind.append(taxa.index(esstaxon))
else:
	esstaxaind = None

minLen = args.min_length

varReq = args.variants_required
minSeqs_s = args.min_unmasked_seqs_per_site
outdir = args.outdir
if not outdir.endswith(os.path.sep):
    outdir += os.path.sep
logpath = '{0}{1}'.format(outdir,"fasta_codon_filtering_for_paml.log")


#check sequence number 
nPass = 0
if numfiles < 1:
	print ("\nNo fasta files found, check working directory contains fasta files and --input_extension is set correctly\n")
	sys.exit()
else:
	if not esstaxaind == None:
		log = open(logpath,'w')
		log.write("filename,N_seqs,N_prefilter_codons,N_postfilter_codons,N_masked_seqs_premature_stop,N_codons_removed_missing_essential_taxa_"+"_".join(esstaxa)+",N_codons_removed_missing_data_threshold_"+str(minSeqs_s)+",N_aa_vars_essential_seqs,passed_len_filter_"+str(minLen)+",passed_min_var_filter_"+str(varReq)+"\n")
	else:
		log = open(logpath,'w')
		log.write("filename,N_seqs,N_prefilter_codons,N_postfilter_codons,N_masked_seqs_premature_stop,N_codons_removed_missing_essential_taxa,N_codons_removed_missing_data_threshold_"+str(minSeqs_s)+",N_aa_vars_essential_seqs,passed_len_filter_"+str(minLen)+",passed_min_var_filter_"+str(varReq)+"\n")
	minLen = minLen*3

# loop through all fastas

	for fasta in glob.glob(files):
		Npremstop = 0
		NMDEss = 0
		NMDThresh = 0
		passed_len = False
		Passed_var = True
		varcount = 0
		fastaname = fasta.split('/')[-1]
		fastabarename = fastaname.split('.')[0]
		aln = AlignIO.read(fasta, "fasta")
		seqLen = aln.get_alignment_length()
		seqlencodon = seqLen / 3
		if (float(seqLen)/3.0).is_integer() == False:
						print '{0}{1}'.format(fastaname, " does not have a multiple of 3 bases, skipping alignment\n")
		else:

# first make alignment with all sequences in -taxa

			fullaln = MultipleSeqAlignment([])
			seq = "X"
			for taxon in taxa: 
				for rec in aln:
					if str(rec.id) == taxon: # find sequence in fasta alignment if it's there
						seq = str(rec.seq)
				if seq == "X":
					seq = ("N" * seqLen) # if not make a sequence of Ns of the correct length

				fullaln.add_sequence(taxon,seq)
				seq = "X"

# convert partially missing codons to missing data. For each sequence, replace sequence with Ns if it has a premature stop codon

			aln_a = np.array([list(rec) for rec in fullaln], np.character)
			for i in range(0,seqLen,3):
				for j in range(numtaxa):
					if "N" in ''.join(aln_a[j,i:i+3]):
						aln_a[j,i] = "N"
						aln_a[j,i+1] = "N"
						aln_a[j,i+2] = "N"

# mask sequences with premature stop codons

			for j in range(numtaxa):
				isstop = False
				for i in range(0,seqLen,3):
					mycodon = ''.join(aln_a[j,i:i+3])
					myaa = Seq(mycodon, generic_dna).translate()
					if not i == (seqLen - 3):
						if myaa == '*':
							for x in range(seqLen):
								aln_a[j,x] = "N"
							isstop = True
					else:
						if myaa == '*':
							aln_a[j,i] = "N"
							aln_a[j,i+1] = "N"
							aln_a[j,i+2] = "N"


				if isstop == True:
						Npremstop+=1

# mask sites with missing data in essential seqs to missing data

			if not esstaxaind == None:
				essremsites = list([])
				for j in esstaxaind: 
					for i in range(seqLen):
						if aln_a[j,i] == "N":
							for x in range(numtaxa):
								aln_a[x,i] = "N"
								if i in list(range(0,seqLen,3)):
									essremsites.append(i)
				essremsites = list(set(essremsites))
				for n in essremsites:
					NMDEss+=1


# then filter by each site in the alignment for missing data, removing columns with less than minSeqs_s unmasked sequences (because of step 1, this can be by site rather than by codon, as if there is an N in one site of a codon there should also be Ns in all other sites of a codon)

			delsites = list([])
			for i in range(seqLen):
				site = list(aln_a[:,i])
				nMissing = float(site.count('N') + site.count('n') + site.count('-')) #	count gaps and missing data
				nNotMissing = numtaxa - nMissing
				if nNotMissing < minSeqs_s:
					delsites.append(i)
			aln_a = np.delete(aln_a,delsites,1)

# convert filtered alignment back to a MultipleSeqAlignment object

			filtaln = MultipleSeqAlignment([])
			c = 0
			for taxon in taxa:  
				filtaln.add_sequence(taxon, ''.join(map(str, list(aln_a[c,:]))))
				c+= 1

# if length of alignment after column-wise filter is above threshold, continue

			seqLen2 = filtaln.get_alignment_length()
			seqlencodonpost = seqLen2/3
			if seqLen2 >= minLen:
				passed_len = True
			else:
				passed_len = False


# make amino acid alignment of essential taxa

			if not esstaxaind == None:
				essaln= MultipleSeqAlignment([])
				for esstaxon in esstaxa:
					for rec in filtaln:
						if str(rec.id) == esstaxon:
							seq = str(rec.seq)
							essaln.add_sequence(esstaxon,seq)
				for rec in essaln:
					rec.seq = rec.seq.translate()

# count amino acid variants between essential taxa, if varcount is at least varReq, continue

				for i in range(essaln.get_alignment_length()):
					mysite = list(essaln[:,i])
					if not all(x==mysite[0] for x in mysite):
						varcount+=1
				if not varcount >= varReq:
					passed_var = False
				else:
					passed_var = True
			else:
				passed_var = True

# output filtered alignments for each sequence that passes the filters
			if passed_var == True and passed_len == True:
				if outfmt == "phylip":
					outpath = '{0}{1}{2}'.format(outdir,fastabarename,".phy")
					f = open (outpath,'w')
					AlignIO.write(filtaln,f,"phylip-sequential")
					nPass+=1
				elif outfmt == "fasta":
					outpath = '{0}{1}{2}'.format(outdir,fastabarename,".fas")
					f = open (outpath,'w')
					AlignIO.write(filtaln,f,"fasta")
					nPass+=1
				else:
					print "output format must be either fasta or phylip\n"
					sys.exit()

# write sequence info to log file
			NMDThresh = seqlencodon - seqlencodonpost - NMDEss
			log.write(fastaname + "," + str(numtaxa) + "," + str(seqlencodon) + "," + str(seqlencodonpost) + "," + str(Npremstop) + "," + str(NMDEss) + "," + str(NMDThresh) + "," +str(varcount)+ "," + str(passed_len) + "," + str(passed_var) + "\n")


print '{0}{1}{2}{3}{4}'.format("Finished, ", nPass, " alignments passed the filters. Filtered alignments were printed to ",outdir,"\n")
