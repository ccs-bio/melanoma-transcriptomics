#!/bin/sh

# ---------------------------------------------------------------------------------------------------------------------
#
#                                   		Center for Computational Science
#												http://www.ccs.miami.edu/
#                             			  			University of Miami
#
#   This software is a "University of Miami Work" under the terms of the United States Copyright Act.
#   Please cite the author(s) in any work or product based on this material.
#
#   OBJECTIVE:
#	The purpose of this program is to download a CDNA file from the Ensembl repository, and then prepare the necessary
#	files for obtaining the Parentgenes from the Pseudogenes.
#
#
#   NOTES:
#   This script should be run in the directory that will hold the CDNA data; said directory will also contain a
#	subdirectory called "blas" were all blast results will be deposited.
#	Do note that ParentGenes have a biotype of "protein_coding", and Pseudogenes have a biotype of ".*pseudogene" for
#	regular expression matching.  Please see the dependencies section below for the required libraries (if any).
#	Also note that this script needs to connect to the outside, so it cannot be run from a computational node.
#
#
#
#   DEPENDENCIES:
#
#       • wget — used to download the CDNA file from Ensembl.
#
#
#	AUTHOR:	Camilo Valdes (cvaldes3@miami.edu)
#			Computational Biology and Bioinformatics Group
#			Center for Computational Science (CCS), University of Miami
#
# ---------------------------------------------------------------------------------------------------------------------

echo ""
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Starting..."
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"


# --------------------------------------------------- File Download ---------------------------------------------------

echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Downloading CDNA file..."

DB_VERSION="76"
CDNA_FILE="Homo_sapiens.GRCh38.cdna.all.fa"
CDNA_FILE_COMPRESSED=$CDNA_FILE".gz"
ENSEMBL_CDNA_FILE_PATH="ftp://ftp.ensembl.org//pub/release-"$DB_VERSION"/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"

# Download the CDNA file using the WGET utility
wget $ENSEMBL_CDNA_FILE_PATH

echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Inflating CDNA file..."

# The file is compressed in ".gz" format, so we use gunzip to extract it
gunzip $CDNA_FILE_COMPRESSED

echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"


# -------------------------------------------------- Asset Extraction -------------------------------------------------


echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Extracting file assets..."
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"

# First we extract all the headers from the sequences.  We do this so that we can extract the information we need to
# create a mapping file of GeneID to Transcript ID.
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Extracting all fasta headers..."

FASTA_HEADERS_ALL="fasta_headers_for_CDNA_full.txt"

grep ">" $CDNA_FILE > $FASTA_HEADERS_ALL


# Once we have all the fasta headers, we pull out all the "protein_coding" as the Parentgenes are a subset.
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Extracting protein_coding headers..."

FASTA_HEADERS_PROTEIN_CODING="protein_coding-fasta_headers.txt"

grep "gene_biotype:protein_coding" Homo_sapiens.GRCh38.cdna.all.fa > $FASTA_HEADERS_PROTEIN_CODING


# We do the same for the pseudogenes and their RE pattern: ".*pseudogene"
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Extracting .*pseudogene headers..."

FASTA_HEADERS_PSEUDOGENES="pseudogenes-fasta_headers.txt"

grep "gene_biotype:.*pseudogene" Homo_sapiens.GRCh38.cdna.all.fa > $FASTA_HEADERS_PSEUDOGENES


# We also need a mapping file to identify ParentGene relationships after the BLAST step.  This mapping file
# maps (for pseudogenes) a TranscriptID to its encassing GeneID.
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Creating Pseudogene Mapping File..."

PSEUDOGENE_MAPPING_FILE="pseudogene_transcriptsIds-and-geneIds.txt"

cut -d' ' -f 1,4 $FASTA_HEADERS_PSEUDOGENES > $PSEUDOGENE_MAPPING_FILE

echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"



# Once we have the necessary fasta headers for the biotypes of interest, we can extract the GENE_ID's of said
# biotypes, but we only pick them up if they are in the standard 25 human chromosomes.

echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Obtaining GeneIDs for Protein Coding genes..."

GENE_ID_FILE_PROTEIN_CODING="gene_ids_for_protein_coding.txt"
grep "chromosome:GRCh38:[0-9]*:\|chromosome:GRCh38:MT:" $FASTA_HEADERS_PROTEIN_CODING | cut -d' ' -f 4  | sed 's/gene://' | sort -u > $GENE_ID_FILE_PROTEIN_CODING

NUMBER_OF_PROTEIN_CODING_GENE_IDS=`wc -l $GENE_ID_FILE_PROTEIN_CODING`
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] protein coding: "$NUMBER_OF_PROTEIN_CODING_GENE_IDS

echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"



echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Obtaining GeneIDs for Pseudogenes..."

GENE_ID_FILE_PSEUDOGENES="gene_ids_for_pseudogenes.txt"
grep "chromosome:GRCh38:[0-9]*:\|chromosome:GRCh38:MT:" $FASTA_HEADERS_PSEUDOGENES |cut -d' ' -f 4  | sed 's/gene://' | sort -u > $GENE_ID_FILE_PSEUDOGENES

NUMBER_OF_PSEUDOGENE_IDS=`wc -l $GENE_ID_FILE_PSEUDOGENES`
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] pseudogenes: "$NUMBER_OF_PSEUDOGENE_IDS

echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"


# -------------------------------------------------- FASTA Extraction -------------------------------------------------

echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Extracting FASTA Sequences..."
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"

PERL_APP_PATH='/nethome/cvaldes/data/references/ensembl/human/76/cdna/scripts'
BLAST_APP='/share/apps/blast+/2.2.28/bin'

CDNA_FASTA='/nethome/cvaldes/data/references/ensembl/human/76/cdna/Homo_sapiens.GRCh38.cdna.all.fa'


INPUT_FILE_WITH_PROTEIN_CODING_IDS=$GENE_ID_FILE_PROTEIN_CODING
INPUT_FILE_WITH_PSEUDOGENE_IDS=$GENE_ID_FILE_PSEUDOGENES

PREFIX_FOR_PROTEIN_CODING='protein_coding'
PREFIX_FOR_PSEUDOGENES='pseudogenes'

OUTPUT_DIR='/nethome/cvaldes/data/references/ensembl/human/'$DB_VERSION'/cdna'


# Run the script to extract the FASTA sequences of interest (these are used in the BLAST step to find 'Best_Hits')
$PERL_APP_PATH/fasta_extractor.pl -i $INPUT_FILE_WITH_PROTEIN_CODING_IDS -o $OUTPUT_DIR -c $CDNA_FASTA -p $PREFIX_FOR_PROTEIN_CODING
$PERL_APP_PATH/fasta_extractor.pl -i $INPUT_FILE_WITH_PSEUDOGENE_IDS -o $OUTPUT_DIR -c $CDNA_FASTA -p $PREFIX_FOR_PSEUDOGENES

echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"



# ------------------------------------------------------- BLAST -------------------------------------------------------
#
#	We BLAST against the ProteinCoding database.  That is, we create a BLAST database with the Protein Coding
#	sequences, and we blast the pseudogene sequences against it.
#

echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Creating BLAST databases..."

FASTA_FILE_FOR_PROTEIN_CODING=$OUTPUT_DIR'/'$PREFIX_FOR_PROTEIN_CODING'-cdna.fasta'

DB_TYPE='nucl'

# Make/Prep the BLAST DB we'll be using
$BLAST_APP/makeblastdb -in $FASTA_FILE_FOR_PROTEIN_CODING -dbtype $DB_TYPE -out $PREFIX_FOR_PROTEIN_CODING -parse_seqids

echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"



# Then run the actual BLAST QUERY (pseudogenes vs protein_coding)
# This will find the protein_coding GeneIds that have a high-level of sequence homology to the pseudogenes used as query,
# essentially finding the 'parent_genes' for the pseudogenes.

BLAST_DB_PROTEIN_CODING=$OUTPUT_DIR'/'$PREFIX_FOR_PROTEIN_CODING
BLAST_QUERY_FASTA_FILE=$OUTPUT_DIR'/'$PREFIX_FOR_PSEUDOGENES'-cdna.fasta'

BLAST_OUPUT_DIR=$OUTPUT_DIR'/blast'
mkdir -p $BLAST_OUPUT_DIR

# This output txt file will hold the blast results
OUTPUT_PSEUDOGENES_VS_PROTEIN_CODING=$BLAST_OUPUT_DIR'/pseudogenes_vs_proteinCoding.blast.out.txt'

BLAST_MAX_TARGETS='1'


echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Running BLAST..."

$BLAST_APP/blastn -db $BLAST_DB_PROTEIN_CODING -query $BLAST_QUERY_FASTA_FILE -outfmt "6 qseqid sseqid evalue bitscore stitle" -max_target_seqs $BLAST_MAX_TARGETS -out $OUTPUT_PSEUDOGENES_VS_PROTEIN_CODING

echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"



# --------------------------------------------- Parent Gene Identification --------------------------------------------

echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Identifying Parent Genes..."


SCRIPT_PATH='/nethome/cvaldes/data/references/ensembl/human/'$DB_VERSION'/cdna/scripts'

INPUT_BLAST_FILE=$OUTPUT_PSEUDOGENES_VS_PROTEIN_CODING

# A Pseudogene mapping file that maps a TranscriptID (pseudogene) to its encasing GeneID (pseudogene)
INPUT_MAPPING_FILE=$PSEUDOGENE_MAPPING_FILE

# Build & Analyze
$SCRIPT_PATH/parentGeneIdentifier.pl -i $INPUT_BLAST_FILE -m $INPUT_MAPPING_FILE -o $BLAST_OUPUT_DIR






echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Done."
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"
echo ""