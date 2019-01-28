# fasta_alignment_filters

This repository contains various python scripts for manipulating and filtering fasta alignments


Author: Owen G Osborne


# cat_alns.py
This script concatenates multiple fasta alignments (e.g. for different genes) into a single fasta alignment. Individual alignments can include missing samples and bases. You must provide the taxa you want to keep in a comma-delimited list, if they are missing from input alignments they are replaced by "N"s. Each individual alignments must include sequences with the exact same names as those in the taxa list. The script must be run from within the input fasta directory.

usage: python2 cat_alns.py -o STR -t STR [-i STR -h]
options:
  -h, --help            show this help message and exit
  -i INPUT_EXTENSION, --input_extension INPUT_EXTENSION
                        input file extension, e.g fa. Files must be in current
                        working directory
  -o OUTFILE, --outfile OUTFILE
                        path to output file
  -t TAXA, --taxa TAXA  comma delimited taxon names

example: cat_alns.py -o outfile.fasta -t species1,species2,species3 -i fas

# make_consensus_fastas.py
This script creates consensus sequences from aligned fastas. You can set the minumum identity to call a consensus base and decide whether to include missing data in the minimum identity calculation. 

usage: make_consensus_fastas.py -o STR [-i STR -I FLOAT -n STR -m -h]
options:
  -h, --help            show this help message and exit
  -i INPUT_EXTENSION, --input_extension INPUT_EXTENSION
                        input file extension, e.g fa. Files must be in current
                        working directory
  -o OUTDIR, --outdir OUTDIR
                        path to output folder (must exist)
  -I IDENTITY, --identity IDENTITY
                        Minimum sequence identity to call consensus
  -n NAME, --name NAME  name of consensus sequence
  -m, --include_missing_data
                        include missing data in calculation, otherwise it will
                        be ignored

# filter_fasta_alns.py

This script filters aligned fastas for missing data and heterozygous sites. It first removes columns which contain a higher percentage of missing data (N or - characters) than the <max_gaps_per_site> threshold, if alignments are at least <min_length> long following this filter, it then removes sequences with a higher percentage of missing data than the <max_gaps_per_sequence> threshold and with a higher percentage of heterozygous sites than the <max_het_per_sequence> threshold. Finally, it outputs filtered alignments for all input files with at least <min_seqs> sequences remaining and produces a log file detailing the filtering in the output directory. All filters are optional, although sequences comprised completely of missing data are always removed. The output directory must be defined and the script must be run from within the input fasta directory

usage: python2 filter_fasta_aln.py -o STR [-i STR -l INT -g INT -G INT -H INT -s INT -h]
options:
  -h, --help            show this help message and exit
  -i INPUT_EXTENSION, --input_extension INPUT_EXTENSION
                        input file extension, e.g fa. Files must be in current
                        working directory
  -o OUTDIR, --outdir OUTDIR
                        path to output folder (must exist)
  -l MIN_LENGTH, --min_length MIN_LENGTH
                        Minimum alignment length
  -g MAX_GAPS_PER_SITE, --max_gaps_per_site MAX_GAPS_PER_SITE
                        Maximum percentage missing data or gaps per site
  -G MAX_GAPS_PER_SEQUENCE, --max_gaps_per_sequence MAX_GAPS_PER_SEQUENCE
                        Maximum percentage missing data or gaps per sequence
  -H MAX_HET_PER_SEQUENCE, --max_het_per_sequence MAX_HET_PER_SEQUENCE
                        Maximum percentage of heterozygous sites per sequence
  -s MIN_SEQS, --min_seqs MIN_SEQS
                        Minimum number of sequences to keep alignment


# trim_alignments_to_CDS.py

This program trims aligned fastas to CDS defined by a gff3 file. The CDS gff file must only contain CDS (eg. by running grep "CDS" transdecoder.gff3 > transdecodercds.gff3) and contain no duplicate records. The script only outputs alignments for records in the positive orientation. The output directory must be defined and the script must be run from within the input fasta directory

usage: python2 trim_alignments_to_CDS.py -o STR -c STR [-i STR -h]
options:
-h, --help            show this help message and exit
-i INPUT_EXTENSION, --input_extension INPUT_EXTENSION
                      fasta file input , e.g fa. Files must be in current
                      working directory
-o OUTDIR, --outdir OUTDIR
                      path to output folder (must exist)
-c CDSFILE, --cdsfile CDSFILE
                      path to CDS only gff3 file

# filter_CDS_alns_for_PAML.py

This script filters CDS alignments and outputs them in phylip format (suitable for PAML) or fasta formate while maintaining CDS frame. Inputs should be aligned coding sequences. You start with a directory full of CDS alignments and provide the script with a list of taxa to include and a list of essential taxa (there must be no missing data in essential taxa for an alignment column to be output). It first makes an alignment of all individuals in the taxa list (as in cat_alns.py, whether they are in the input alignment or not - missing taxa are replaced by Ns). It then masks codons with any missing data (i.e. converts to NNN for only sequences with missing data), masks entire sequences that contain premature stop codons, masks codons that have missing data in essential taxa in the whole alignment. It then filters by each (codon) site in the alignment for missing data, removing columns with less than -s unmasked sequences. If -l, -e and -v are set, it then only outputs alignments which are longer than -l codons and have at least -e amino acid variants between the essential taxa (-e.)

usage: python2 filter_CDS_alns_for_PAML.py -o STR -t STR [-i STR -l INT -s INT -e STR -v INT]
options:
  -h, --help            show this help message and exit
  -i INPUT_EXTENSION, --input_extension INPUT_EXTENSION
                        input file extension, e.g fa. Files must be in current
                        working directory
  -o OUTDIR, --outdir OUTDIR
                        path to output folder (must exist)
  -l MIN_LENGTH, --min_length MIN_LENGTH
                        Minimum alignment length in codons
  -s MIN_UNMASKED_SEQS_PER_SITE, --min_unmasked_seqs_per_site MIN_UNMASKED_SEQS_PER_SITE
                        Minimum sequences per codon site
  -e ESSENTIAL_TAXA, --essential_taxa ESSENTIAL_TAXA
                        comma delimited list of taxon names which are required
                        for each site to be included, these are also the taxa
                        which must have SNPs between them if
                        --variants_required is set
  -v VARIANTS_REQUIRED, --variants_required VARIANTS_REQUIRED
                        Minimum SNPs required between essential taxa
  -t TAXA, --taxa TAXA  comma delimited taxon names
  -f OUTPUT_FORMAT, --output_format OUTPUT_FORMAT
                        output format, either fasta or phylip (default)

example: python2 filter_CDS_alns_for_PAML.py -o /output/directory/ -t species1,species2,species3,species4 -e species1,species2 -l 100 -s 3 -v 1
