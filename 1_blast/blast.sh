#!/bin/sh

PERL_APP_PATH='/projects/bioinf/data/cvaldes/tcga/skcm/examination'
BLAST_APP='/share/apps/blast+/2.2.28/bin'

CDNA_FASTA='/nethome/cvaldes/data/references/ensembl/human/72/cdna/Homo_sapiens.GRCh37.72.cdna.all.fa'

echo ""
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Starting..."

INPUT_FILE_WITH_PROTEIN_CODING_IDS='/nethome/cvaldes/data/references/ensembl/human/72/cdna/gene_ids_for_protein_coding.txt'
INPUT_FILE_WITH_PSEUDOGENE_IDS='/nethome/cvaldes/data/references/ensembl/human/72/cdna/gene_ids_for_pseudogenes.txt'

PREFIX_FOR_PROTEIN_CODING='protein_coding'
PREFIX_FOR_PSEUDOGENES='pseudogenes'

OUTPUT_DIR='/nethome/cvaldes/data/references/ensembl/human/72/cdna'

echo ""
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Extracting FASTA sequences..."

# Run the script to extract the FASTA sequences of interest (these are used in the BLAST step to find 'Best_Hits')
# $PERL_APP_PATH/fasta_extractor.pl -i $INPUT_FILE_WITH_PROTEIN_CODING_IDS -o $OUTPUT_DIR -c $CDNA_FASTA -p $PREFIX_FOR_PROTEIN_CODING
# $PERL_APP_PATH/fasta_extractor.pl -i $INPUT_FILE_WITH_PSEUDOGENE_IDS -o $OUTPUT_DIR -c $CDNA_FASTA -p $PREFIX_FOR_PSEUDOGENES


#------------------------------------------------------------------------------------------------------------
#
# Create the BLAST protein_coding database
#

FASTA_FILE_FOR_PROTEIN_CODING=$OUTPUT_DIR'/'$PREFIX_FOR_PROTEIN_CODING'-cdna.fasta'

DB_TYPE='nucl'

echo ""
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Creating BLAST database ("$PREFIX_FOR_PROTEIN_CODING")..."

# Make/Prep the BLAST DB we'll be using
$BLAST_APP/makeblastdb -in $FASTA_FILE_FOR_PROTEIN_CODING -dbtype $DB_TYPE -out $PREFIX_FOR_PROTEIN_CODING -parse_seqids


#------------------------------------------------------------------------------------------------------------
#
# Then run the actual BLAST QUERY (pseudogenes vs protein_coding)
# This will find the protein_coding GeneIds that have a high-level of sequence homology to the pseudogenes used as query,
# essentially finding the 'parent_genes' for the pseudogenes.

BLAST_DB_PROTEIN_CODING=$OUTPUT_DIR'/'$PREFIX_FOR_PROTEIN_CODING
BLAST_QUERY_FASTA_FILE=$OUTPUT_DIR'/'$PREFIX_FOR_PSEUDOGENES'-cdna.fasta'

BLAST_OUPUT_DIR=$OUTPUT_DIR'/blast'
mkdir -p $BLAST_OUPUT_DIR

OUTPUT_PSEUDOGENES_VS_PROTEIN_CODING=$BLAST_OUPUT_DIR'/pseudogenes_vs_proteinCoding.blast.out.txt'

BLAST_MAX_TARGETS='1'

echo ""
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Running BLAST..."

$BLAST_APP/blastn -db $BLAST_DB_PROTEIN_CODING -query $BLAST_QUERY_FASTA_FILE -outfmt "6 qseqid sseqid evalue bitscore stitle" -max_target_seqs $BLAST_MAX_TARGETS -out $OUTPUT_PSEUDOGENES_VS_PROTEIN_CODING

#------------------------------------------------------------------------------------------------------------
#
# Extract Unique Hits from the BLAST output
#

OUTPUT_SINGLE_HITS_PSEUDOGENE_VS_PROTEIN_CODING=$BLAST_OUPUT_DIR'/pseudogenes_vs_proteinCoding-single_hits.txt'

echo ""
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Extracting Unique hits..."

#sort -u $OUTPUT_PSEUDOGENES_VS_PROTEIN_CODING > $OUTPUT_SINGLE_HITS_PSEUDOGENE_VS_PROTEIN_CODING
