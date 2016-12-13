#!/bin/sh
#BSUB -J 0c2c67ce-c9ff-45ad-a11b-3abc15b7439b
#BSUB -o /path/to/logs/log_file.txt
#BSUB -W 168:00
#BSUB -N
#BSUB -u "cvaldes3@miami.edu"
#BSUB -x
#BSUB -n 16
#BSUB -R "span[ptile=16]"
#BSUB -M 33554432

# ---------------------------------------------------------------------------------------------------------
#
#											Main Alignment Run
#
#	The purpose of this script is to process (align, quantify) TCGA Melanoma (SKCM) samples to the human
#	genome.
#
#	Version: 2.0
#

echo ""
echo ""
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Starting Run..."


# ------------------------------------------ Sample Particulars -------------------------------------------

SAMPLE_BASE='0c2c67ce-c9ff-45ad-a11b-3abc15b7439b'
READ_BASE='UNCID_1286603.27bf4227-d8a9-4048-af7e-01c6d5d0e777'

# USER_DIR='/projects/bioinf/data/cvaldes'
USER_DIR='/scratch/cvaldes'

BASE_OUTPUT_DIR=$USER_DIR'/tcga/skcm'
mkdir -p $BASE_OUTPUT_DIR

# Main Output directory for Alignments, Quants, etc (all subsequent operations use the results from here)
ANALYSIS_BASE_DIR=$BASE_OUTPUT_DIR'/analysis_78'


# ------------------------------------------------ Tools --------------------------------------------------

# BAMTOOLS, used for converting BAM file to FASTQ
BAM2FASTQ_APP_PATH='/nethome/bioinfo/tools/bamtools/bam2fastq/1.1.0'

# SamTools, used for extracting statistics, parsing, and general processing of SAM/BAM output files
SAMTOOLS_PATH='/nethome/bioinfo/tools/samtools/0.1.19'

# TopHat, splice mapping tools used to align a RNA-Seq sample to the human genome
TOPHAT_APP_PATH='/nethome/bioinfo/tools/tophat/v2/2.0.12/'

# Cufflinks, used to Quantify & Assemble the mapped RNA-Seq reads
CUFFLINKS_APP_PATH='/nethome/bioinfo/tools/cufflinks/v2/2.2.1'

# Count extractor, a Perl script that generates the genomic-interval counts for a given gene assembly
COUNT_APP_PATH='/path/to/counting/script/perl'


# ------------------------------------------------ Tools --------------------------------------------------

# Number of Threads to dispatch
PROC_THREADS='32'


# ------------------------------------------ Genome Annotations -------------------------------------------

ANNOT_VERSION='72'

# Ensembl annotations file with known biotypes
GTF_FILE='/path/to/ensembl/human/'$ANNOT_VERSION'/gtf/Homo_sapiens-CLEAN.GRCh37.'$ANNOT_VERSION'.gtf'

# Mask file for ignoring highly abundant rRNAs
RRNA_MASK_GTF_FILE='/npath/to/references/ensembl/human/'$ANNOT_VERSION'/gtf/rRNAs.gtf'

# Genome index database(s) in the format for each algorithm
INDEX_FOR_BOWTIE2='/path/to/ref/annotations/ensembl/human/'$ANNOT_VERSION'/dna/bowtie2/ensembl_hs_'$ANNOT_VERSION'.dna'

# Human Genome FASTA File
GENOME_FASTA_FILE_DNA='/path/to/ref/annotations/ensembl/human/'$ANNOT_VERSION'/dna/bowtie2/ensembl_hs_'$ANNOT_VERSION'.dna.fa'


# -------------------------------------------- Miscellaneous ----------------------------------------------

DASHES=' ------------------------------------------------------------------'


# --------------------------------------------- BAM to FASTQ ----------------------------------------------

echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Inflating BAM to FASTQ..."
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"

OUTPUT_DIR_FOR_FASTQ=$BASE_OUTPUT_DIR'/fastq/'$SAMPLE_BASE
mkdir -p $OUTPUT_DIR_FOR_FASTQ

BAM_FILE='/path/to/tcga/skcm/raw/'$SAMPLE_BASE'/'$READ_BASE'.sorted_genome_alignments.bam'

FASTQ_BASE_NAME=$OUTPUT_DIR_FOR_FASTQ'/'$READ_BASE

# Template for the output files of the FASTQ conversion tool
FASTQ_BASE_NAME_MASK=$FASTQ_BASE_NAME'#'

# Run the conversion tool
$BAM2FASTQ_APP_PATH/bam2fastq --force --output $FASTQ_BASE_NAME_MASK $BAM_FILE

# FASTQ files of the BAM files are stored in the following two (2) files
READ_SAMPLE_1=$FASTQ_BASE_NAME'_1'
READ_SAMPLE_2=$FASTQ_BASE_NAME'_2'

echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Done."
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"


# ---------------------------------------------- Alignment ------------------------------------------------

echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Mapping sample using TopHat..."

OUTPUT_DIR_ALIGNMENTS_DNA=$ANALYSIS_BASE_DIR'/'$SAMPLE_BASE'/alignments/dna'
mkdir -p $OUTPUT_DIR_ALIGNMENTS_DNA

#
#	TopHat Parameters
#
TH_MIN_ANCHOR='6'

echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Aligning to DNA version..."

#
#	TopHat run
#
$TOPHAT_APP_PATH/tophat2 -p $PROC_THREADS --output-dir $OUTPUT_DIR_ALIGNMENTS_DNA \
										  --min-anchor $TH_MIN_ANCHOR \
										  --no-coverage-search \
										  --GTF $GTF_FILE \
										  $INDEX_FOR_BOWTIE2 \
										  $READ_SAMPLE_1 \
										  $READ_SAMPLE_2

echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Done."
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"

# ----------------------------------------------- SamTools ------------------------------------------------
#
#	The gimmick here is to split the Alignments into three (3) files: All, Unique, and Multi.
#	The Unique file contains alignment above a mapping-quality threshold (MAPQ=30 for 'unique'), while the
#	Multi file contains all the remaining alignments below the 'unique' threshold.
#
#	Note that TopHat does not report a uniform distribution of MAPQs, but rather, it bins alignments into classes of
#	MAPQ 255	= unique mapping
#	MAPQ 3 		= maps to 2 locations in the target
#	MAPQ 2 		= maps to 3 locations
#	MAPQ 1 		= maps to 4-9 locations
#	MAPQ 0 		= maps to 10 or more locations.
#

echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Samtools..."
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Extracting Alignment Quality Classes..."
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"

# Mapping quality threshold -- reads with this value (or above) will be filtered out
MAPQ_QUAL='1'

# The files that will point to the respective alignment classes, and that will be used in downstream analysis
ALIGNMENT_FILE_WITH_MULTI_HITS=$OUTPUT_DIR_ALIGNMENTS_DNA'/accepted_hits-multi.sorted.bam'
ALIGNMENT_FILE_WITH_MULTI_HITS_MAPQ_FILTERED=$OUTPUT_DIR_ALIGNMENTS_DNA'/accepted_hits-multi_mapq'$MAPQ_QUAL'.sorted.bam'
ALIGNMENT_FILE_WITH_UNIQUE_HITS=$OUTPUT_DIR_ALIGNMENTS_DNA'/accepted_hits-unique.sorted.bam'
ALIGNMENT_FILE_WITH_UNIQUE_HITS_SORTED_BY_NAME=$OUTPUT_DIR_ALIGNMENTS_DNA'/accepted_hits-unique.sorted_by_name.bam'
ALIGNMENT_FILE_WITH_UNMAPPED=$OUTPUT_DIR_ALIGNMENTS_DNA'unmapped.bam'

# ------------------------------------ All ------------------------------------

echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"$DASHES
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Samtools - MAPQ = All"
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"$DASHES

$SAMTOOLS_PATH/samtools sort $OUTPUT_DIR_ALIGNMENTS_DNA/accepted_hits.bam $OUTPUT_DIR_ALIGNMENTS_DNA/accepted_hits.sorted
$SAMTOOLS_PATH/samtools index $OUTPUT_DIR_ALIGNMENTS_DNA/accepted_hits.sorted.bam
$SAMTOOLS_PATH/samtools flagstat $OUTPUT_DIR_ALIGNMENTS_DNA/accepted_hits.sorted.bam > $OUTPUT_DIR_ALIGNMENTS_DNA/flagStats-all.txt

# ------------------------------------ Multi ----------------------------------

echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"$DASHES
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Samtools - MAPQ = Multi"
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"$DASHES

$SAMTOOLS_PATH/samtools view -b -f 256 $OUTPUT_DIR_ALIGNMENTS_DNA/accepted_hits.bam > $OUTPUT_DIR_ALIGNMENTS_DNA/accepted_hits-multi.bam
$SAMTOOLS_PATH/samtools sort $OUTPUT_DIR_ALIGNMENTS_DNA/accepted_hits-multi.bam $OUTPUT_DIR_ALIGNMENTS_DNA/accepted_hits-multi.sorted
$SAMTOOLS_PATH/samtools index $OUTPUT_DIR_ALIGNMENTS_DNA/accepted_hits-multi.sorted.bam
$SAMTOOLS_PATH/samtools flagstat $OUTPUT_DIR_ALIGNMENTS_DNA/accepted_hits-multi.sorted.bam > $OUTPUT_DIR_ALIGNMENTS_DNA/flagStats-multi.txt


# --------------------------------- Multi MAPQ --------------------------------

echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"$DASHES
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Samtools - MAPQ = Multi_MAPQ"
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"$DASHES

INPUT_FILE_WITH_MULTI_ALIGNMENTS=$OUTPUT_DIR_ALIGNMENTS_DNA'/accepted_hits-multi.sorted.bam'
OUTPUT_FILE_BASE_NAME='accepted_hits-multi_mapq'
OUTPUT_FILE_WITH_MULTI_ALIGNMENTS_FILTERED=$OUTPUT_DIR_ALIGNMENTS_DNA'/'$OUTPUT_FILE_BASE_NAME''$MAPQ_QUAL'.bam'
OUTPUT_FILE_WITH_MULTI_ALIGNMENTS_FILTERED_SORTED=$OUTPUT_DIR_ALIGNMENTS_DNA'/'$OUTPUT_FILE_BASE_NAME''$MAPQ_QUAL'.sorted'
OUTPUT_FILE_WITH_MULTI_ALIGNMENTS_FILTERED_SORTED_BAM=$OUTPUT_DIR_ALIGNMENTS_DNA'/'$OUTPUT_FILE_BASE_NAME''$MAPQ_QUAL'.sorted.bam'
OUTPUT_FILE_WITH_MULTI_FLAGSTATS=$OUTPUT_DIR_ALIGNMENTS_DNA'/flagstats-multi_mapq'$MAPQ_QUAL'.txt'

$SAMTOOLS_PATH/samtools view -b -q $MAPQ_QUAL $INPUT_FILE_WITH_MULTI_ALIGNMENTS > $OUTPUT_FILE_WITH_MULTI_ALIGNMENTS_FILTERED
$SAMTOOLS_PATH/samtools sort $OUTPUT_FILE_WITH_MULTI_ALIGNMENTS_FILTERED $OUTPUT_FILE_WITH_MULTI_ALIGNMENTS_FILTERED_SORTED
$SAMTOOLS_PATH/samtools index $OUTPUT_FILE_WITH_MULTI_ALIGNMENTS_FILTERED_SORTED_BAM
$SAMTOOLS_PATH/samtools flagstat $OUTPUT_FILE_WITH_MULTI_ALIGNMENTS_FILTERED_SORTED_BAM > $OUTPUT_FILE_WITH_MULTI_FLAGSTATS

# ----------------------------------- Unique ----------------------------------

echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"$DASHES
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Samtools - MAPQ = Unique"
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"$DASHES

$SAMTOOLS_PATH/samtools view -b -q 30 $OUTPUT_DIR_ALIGNMENTS_DNA/accepted_hits.bam > $OUTPUT_DIR_ALIGNMENTS_DNA/accepted_hits-unique.bam
# Sort by Position
$SAMTOOLS_PATH/samtools sort $OUTPUT_DIR_ALIGNMENTS_DNA/accepted_hits-unique.bam $OUTPUT_DIR_ALIGNMENTS_DNA/accepted_hits-unique.sorted
# Sort by Name
$SAMTOOLS_PATH/samtools sort -n $OUTPUT_DIR_ALIGNMENTS_DNA/accepted_hits-unique.bam $OUTPUT_DIR_ALIGNMENTS_DNA/accepted_hits-unique.sorted_by_name
# Index & Flagstat Position-based alignments
$SAMTOOLS_PATH/samtools index $OUTPUT_DIR_ALIGNMENTS_DNA/accepted_hits-unique.sorted.bam
$SAMTOOLS_PATH/samtools flagstat $OUTPUT_DIR_ALIGNMENTS_DNA/accepted_hits-unique.sorted.bam > $OUTPUT_DIR_ALIGNMENTS_DNA/flagStats-unique.txt


# Remove files that are no longer needed
rm -fR $OUTPUT_DIR_ALIGNMENTS_DNA/accepted_hits.bam
rm -fR $OUTPUT_DIR_ALIGNMENTS_DNA/accepted_hits.sorted.bam
rm -fR $OUTPUT_DIR_ALIGNMENTS_DNA/accepted_hits-multi.bam	# Remove the MULTI unsorted BAM file
rm -fR $OUTPUT_FILE_WITH_MULTI_ALIGNMENTS_FILTERED			# Remove the MAPQ_1 unsorted BAM file
rm -fR $OUTPUT_DIR_ALIGNMENTS_DNA/accepted_hits-unique.bam	# Remove the unique unsorted BAM file

echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Done."
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"


# ------------------------------------ FPKM Quantification & Assembly -------------------------------------

echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Cufflinks..."
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"

#
#	Cufflinks parameters
#
OVERHANG_TOLERANCE='8'			# The number of bp allowed to enter the intron of a transcript when determining if a read or another transcript is mappable to/compatible with it.
MIN_ISOFORM_FRACTION='0.00'		# Cufflinks filters out transcripts that it believes are very low abundance — 0.00 turns it off.
PRE_MRNA_FRACTION='0.00'		# Cufflinks uses this parameter to filter out alignments that lie within the intronic intervals implied by the spliced alignments.

OUTPUT_DIR_FOR_CUFFLINKS=$ANALYSIS_BASE_DIR'/'$SAMPLE_BASE'/quantification/dna/assemblies'
mkdir -p $OUTPUT_DIR_FOR_CUFFLINKS

ASSEMBLY_OUT_DIR_MULTI_DNA=$OUTPUT_DIR_FOR_CUFFLINKS'/multi'
ASSEMBLY_OUT_DIR_MULTI_DNA_MAPQ_FILTERED=$OUTPUT_DIR_FOR_CUFFLINKS'/multi_mapq'$MAPQ_QUAL
ASSEMBLY_OUT_DIR_UNIQUE_DNA=$OUTPUT_DIR_FOR_CUFFLINKS'/unique'

mkdir -p $ASSEMBLY_OUT_DIR_MULTI_DNA
mkdir -p $ASSEMBLY_OUT_DIR_MULTI_DNA_MAPQ_FILTERED
mkdir -p $ASSEMBLY_OUT_DIR_UNIQUE_DNA

# ------------------------------------ Cufflinks Multi ----------------------------------

echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"$DASHES
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Cufflinks Multi..."
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"$DASHES
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"

$CUFFLINKS_APP_PATH/cufflinks --no-update-check -p $PROC_THREADS \
												--no-faux-reads \
												--output-dir $ASSEMBLY_OUT_DIR_MULTI_DNA \
												--min-isoform-fraction $MIN_ISOFORM_FRACTION \
												--pre-mrna-fraction $PRE_MRNA_FRACTION \
												--overhang-tolerance $OVERHANG_TOLERANCE \
												--GTF-guide $GTF_FILE \
												--mask-file $RRNA_MASK_GTF_FILE \
												$ALIGNMENT_FILE_WITH_MULTI_HITS

# --------------------------------- Cufflinks Multi MAPQ --------------------------------

echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"$DASHES
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Cufflinks Multi_MAPQ..."
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"$DASHES
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"

$CUFFLINKS_APP_PATH/cufflinks --no-update-check -p $PROC_THREADS \
												--no-faux-reads \
												--output-dir $ASSEMBLY_OUT_DIR_MULTI_DNA_MAPQ_FILTERED \
												--min-isoform-fraction $MIN_ISOFORM_FRACTION \
												--pre-mrna-fraction $PRE_MRNA_FRACTION \
												--overhang-tolerance $OVERHANG_TOLERANCE \
												--GTF-guide $GTF_FILE \
												--mask-file $RRNA_MASK_GTF_FILE \
												$ALIGNMENT_FILE_WITH_MULTI_HITS_MAPQ_FILTERED

# -----------------------------------  Cufflinks Unique ---------------------------------

echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"$DASHES
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Cufflinks Unique..."
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"$DASHES
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"

$CUFFLINKS_APP_PATH/cufflinks --no-update-check -p $PROC_THREADS \
												--no-faux-reads \
												--output-dir $ASSEMBLY_OUT_DIR_UNIQUE_DNA \
												--min-isoform-fraction $MIN_ISOFORM_FRACTION \
												--pre-mrna-fraction $PRE_MRNA_FRACTION \
												--overhang-tolerance $OVERHANG_TOLERANCE \
												--GTF-guide $GTF_FILE \
												--mask-file $RRNA_MASK_GTF_FILE \
												$ALIGNMENT_FILE_WITH_UNIQUE_HITS


echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Done."
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"


# ------------------------------------------------ Counts -------------------------------------------------

echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Counts..."
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"

OUTPUT_DIR_FOR_COUNTS_BASE=$ANALYSIS_BASE_DIR'/'$SAMPLE_BASE'/counts/dna'
OUTPUT_DIR_FOR_COUNTS_UNIQUE=$OUTPUT_DIR_FOR_COUNTS_BASE'/unique_counts'
OUTPUT_DIR_FOR_COUNTS_MULTI=$OUTPUT_DIR_FOR_COUNTS_BASE'/multi_counts'
OUTPUT_DIR_FOR_COUNTS_MULTI_MAPQ_FILTERED=$OUTPUT_DIR_FOR_COUNTS_BASE'/multi_mapq'$MAPQ_QUAL'_counts'

mkdir -p $OUTPUT_DIR_FOR_COUNTS_BASE
mkdir -p $OUTPUT_DIR_FOR_COUNTS_UNIQUE
mkdir -p $OUTPUT_DIR_FOR_COUNTS_MULTI
mkdir -p $OUTPUT_DIR_FOR_COUNTS_MULTI_MAPQ_FILTERED

OUTPUT_DIR_FOR_EXONIC_COUNTS_BASE=$ANALYSIS_BASE_DIR'/'$SAMPLE_BASE'/counts_exonic/dna'
OUTPUT_DIR_FOR_EXONIC_COUNTS_UNIQUE=$OUTPUT_DIR_FOR_EXONIC_COUNTS_BASE'/unique_exonic_counts'

mkdir -p $OUTPUT_DIR_FOR_EXONIC_COUNTS_BASE
mkdir -p $OUTPUT_DIR_FOR_EXONIC_COUNTS_UNIQUE

#
#	Prefix for Intronic counts Script
#
PREFIX_FOR_COUNTS_MULTI='multi'
PREFIX_FOR_COUNTS_MULTI_MAPQ_FILTERED='multi_mapq'$MAPQ_QUAL
PREFIX_FOR_COUNTS_UNIQUE='unique'




# ------------------------------------- Counts Multi ------------------------------------

echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"$DASHES
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Counts Multi..."
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"$DASHES
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"

$COUNT_APP_PATH/countExtractor.pl -a $ALIGNMENT_FILE_WITH_MULTI_HITS \
								  -g $GTF_FILE \
								  -o $OUTPUT_DIR_FOR_COUNTS_MULTI \
								  -p $PREFIX_FOR_COUNTS_MULTI

# ---------------------------------- Counts Multi MAPQ ----------------------------------

echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"$DASHES
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Counts Multi_MAPQ..."
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"$DASHES
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"

$COUNT_APP_PATH/countExtractor.pl -a $ALIGNMENT_FILE_WITH_MULTI_HITS_MAPQ_FILTERED \
							      -g $GTF_FILE \
							      -o $OUTPUT_DIR_FOR_COUNTS_MULTI_MAPQ_FILTERED \
							      -p $PREFIX_FOR_COUNTS_MULTI_MAPQ_FILTERED

# ------------------------------------  Counts Unique -----------------------------------

echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"$DASHES
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Counts Unique..."
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"$DASHES
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"


echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Intronic Counts..."
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"

$COUNT_APP_PATH/countExtractor.pl -a $ALIGNMENT_FILE_WITH_UNIQUE_HITS \
							      -g $GTF_FILE \
							      -o $OUTPUT_DIR_FOR_COUNTS_UNIQUE \
							      -p $PREFIX_FOR_COUNTS_UNIQUE



# ------------------------------------  Exonic Counts Unique -----------------------------------
#
#	HTSeq Exonic Counts Unique
#

echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Exonic Counts..."
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"

OUTPUT_FILE_FOR_EXONIC_COUNTS_UNIQUE=$OUTPUT_DIR_FOR_EXONIC_COUNTS_UNIQUE'/unique_exonic_counts.txt'

htseq-count --quiet --format='bam' --order='name' $ALIGNMENT_FILE_WITH_UNIQUE_HITS_SORTED_BY_NAME $GTF_FILE > $OUTPUT_FILE_FOR_EXONIC_COUNTS_UNIQUE

echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Done."
echo ""



# ------------------------------------------------ CuffQuant -------------------------------------------------


echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"$DASHES
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] CuffQuant..."
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"$DASHES
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"

OUTPUT_DIR_FOR_CUFFQUANT=$ANALYSIS_BASE_DIR'/'$SAMPLE_BASE'/quantification/dna/cuffquant'
mkdir -p $OUTPUT_DIR_FOR_CUFFQUANT

CUFFQUANT_OUT_DIR_UNIQUE_DNA=$OUTPUT_DIR_FOR_CUFFQUANT'/unique'
mkdir -p $CUFFQUANT_OUT_DIR_UNIQUE_DNA

#
#	CuffQuant
#
$CUFFLINKS_APP_PATH/cuffquant   --no-update-check \
								-p $PROC_THREADS \
								-o $CUFFQUANT_OUT_DIR_UNIQUE_DNA \
								-M $RRNA_MASK_GTF_FILE \
								-b $GENOME_FASTA_FILE_DNA \
								$GTF_FILE \
								$ALIGNMENT_FILE_WITH_UNIQUE_HITS


echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Done."
echo ""


# ------------------------------------------------ Clean Up -------------------------------------------------

echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"$DASHES
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Cleaning up..."
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"$DASHES
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"

# Remove the FASTQ files
rm -fR $READ_SAMPLE_1
rm -fR $READ_SAMPLE_2

# Remove the Alignment files
rm -fR $ALIGNMENT_FILE_WITH_MULTI_HITS
rm -fR $ALIGNMENT_FILE_WITH_MULTI_HITS_MAPQ_FILTERED
# rm -fR $ALIGNMENT_FILE_WITH_UNIQUE_HITS
rm -fR $ALIGNMENT_FILE_WITH_UNMAPPED


echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Fin."
echo ""
echo ""