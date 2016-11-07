#!/share/opt/perl/5.18.1/bin/perl -w

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
#	The purpose of this program is to obtain pseudogene parent-gene information from two files: one is a
#   BLAST ouput formatted file, and the other is two (2) column file with pseudogene transcripts and gene_ids.
#
#
#	USAGE:
#
#		parentGeneIdentifier.pl -i pseudogenes_vs_proteinCoding.blast.out.txt -o /blast -m pseudogene_transcriptsIds-and-geneIds.txt
#
#
#   NOTES:
#   Please see the dependencies section below for the required libraries (if any).
#
#   DEPENDENCIES:
#
#       • Getopt::Long
#		• POSIX
#
#
#	AUTHOR:	Camilo Valdes (cvaldes3@miami.edu)
#			Computational Biology and Bioinformatics Group
#			Center for Computational Science (CCS), University of Miami
#
# ---------------------------------------------------------------------------------------------------------------------

# Perl Modules
use strict;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Benchmark;
use Cwd;
use POSIX qw/strftime/;

# ------------------------------------------------------- Main --------------------------------------------------------

my $bundle = "0.1";
my $build  = "A0001";

my $file;               # Input File with BLAST output format (-outfmt "6 qseqid sseqid evalue bitscore stitle")
my $outputDir;			# Output Directory
my $mappingFile;        # Two (2) column file mapping a pseudogene transcript id to its base (origin) gene_id

GetOptions ( 'i=s' => \$file, 'o:s' => \$outputDir, 'm=s' => \$mappingFile );

if( ! defined( $file) )
{ print( RED, "\n\nWARNING! Missing Files ", RESET ); print( ".. $!\n\n" ); exit(1); }

chomp($file);

my $currentWD = getcwd();       # Get a pointer to the Current Working Directory

# Default to the current working directory if none provided
if( !defined( $outputDir ) )
{ $outputDir = $currentWD; }

$| = 1;  # Flush STDOUT

print( "\nparentGeneIdentifier.pl, ($bundle, $build)\n" );
print( "===========================================================\n" );

# --------------------------------------------------- Output Files ----------------------------------------------------

my $fout1 = $outputDir . "/" . "parentGenes-GeneIDs.txt";
unless( open( OUTFILE, ">$fout1" ) ) { print "File $fout1 does not exist"; exit; }

print OUTFILE ( "Query GeneID\tSubject GeneID (Best Hit)\n" );

my $fout2 = $outputDir . "/" . "parentGenes-TranscriptIDs.txt";
unless( open( OUTFILE2, ">$fout2" ) ) { print "File $fout2 does not exist"; exit; }

print OUTFILE2 ( "Query TranscriptID\tSubject TranscriptID (Best Hit)\n" );

my $fout3 = $outputDir . "/" . "parentGenes-Biotypes.txt";
unless( open( OUTFILE3, ">$fout3" ) ) { print "File $fout3 does not exist"; exit; }

print OUTFILE3 ( "Subject GeneID\tSubject Biotype\n" );

my $fout4 = $outputDir . "/" . "parentGenes-All_Summary.txt";
unless( open( OUTFILE4, ">$fout4" ) ) { print "File $fout4 does not exist"; exit; }

print OUTFILE4 ( "Query GeneID\tQuery TranscriptID\tSubject TranscriptID (Best Hit)\tSubject GeneID (Best Hit)\tSubject Biotype\n" );

my $fout5 = $outputDir . "/" . "parentGenes-Evalues.txt";
unless( open( OUTFILE5, ">$fout5" ) ) { print "File $fout5 does not exist"; exit; }

print OUTFILE5 ( "Subject GeneID (Parent Gene)\tQuery GeneID (Pseudogene)\tE-value\n" );

# Output file for pseudogenes that have only 1 parent (1:1 ratios of parent and child)
my $fout6 = $outputDir . "/" . "parentGenes-Evalues-singletons.txt";
unless( open( OUTFILE6, ">$fout6" ) ) { print "File $fout6 does not exist"; exit; }

print OUTFILE6 ( "Subject GeneID (Parent Gene)\tQuery GeneID (Pseudogene)\tE-value\n" );



# ------------------------------------------------- Blast Hits Loader -------------------------------------------------

print( BOLD, GREEN, "\n--\n", RESET );
print( BOLD, "Loading BLAST Hits ... ", RESET );

unless( open( INFILE, $file ) )
{ die( "\n\nUnable to open Input File: $file.  Please verify.\n\n" ); }

# Maps a BLAST subject (mapping target) to its corresponding GeneID
# KEY: Subject TranscriptID
# VAL: Subject GeneID
my %subjects = ();

# Maps a subject GeneID to its corresponding Biotype
my %subjectBiotype = ();

# Maps a subject GeneID to its corresponding evalue (vs parentGene)
my %subjectEvalueHashAsTranscripts = ();

# Maps a subject GeneID to its corresponding bitscore (vs parentGene)
my %subjectBitScoreHash = ();

# Maps a BLAST query to it's subject
# But we only want the TOP HIT
# KEY: Query TranscriptID
# VAL: Subject TranscriptID
my %blastHits = ();

# Maps a QUERY transcriptID to a SUBJECT transcriptID
# KEY: Query TranscriptID
# VAL: Subject TranscriptID
my %queryToSubject = ();

my $numberOfLines = 0;

while( my $line = <INFILE> )
{
    chomp($line);
	$line   =~ s/\r//g;
    my @lineArray = split(/\t/,$line);

    my $queryTranscriptID   = $lineArray[ 0 ];
    my $subjectTranscriptID = $lineArray[ 1 ];
    my $subjectEvalue		= $lineArray[ 2 ];
    my $subjectBitScore     = $lineArray[ 3 ];
    my $subjectDescription  = $lineArray[ 4 ];

    my @arrayForSubjectDesc = split(/ /,$subjectDescription);
    my @geneIDTokenizer     = split(/:/,$arrayForSubjectDesc[ 2 ]);
    my @biotypeTokenizer    = split(/:/,$arrayForSubjectDesc[ 3 ]);

    my $subjectGeneID = $geneIDTokenizer[ 1 ];
    my $subjectGeneBiotype = $biotypeTokenizer[ 1 ];

    if( ! defined( $subjects{ $subjectTranscriptID } ) ) {
        $subjects{ $subjectTranscriptID } = $subjectGeneID;
    }

    if( ! defined( $blastHits{ $queryTranscriptID } ) ) {
        $blastHits{ $queryTranscriptID } = $subjectGeneID;
    }

    if( ! defined( $queryToSubject{ $queryTranscriptID } ) ) {
        $queryToSubject{ $queryTranscriptID } = $subjectTranscriptID;
    }

    if( ! defined( $subjectBiotype{ $subjectGeneID } ) ) {
        $subjectBiotype{ $subjectGeneID } = $subjectGeneBiotype;
    }

    $subjectEvalueHashAsTranscripts{ $subjectGeneID }{ $queryTranscriptID } = $subjectEvalue;


    $numberOfLines++;
}

close(INFILE);

print( CYAN, "\n Number of lines processed in BLAST file: ", RESET );
print( " " . &addCommas($numberOfLines) . "\n" );


# ------------------------------------------------ Mapping File Loader ------------------------------------------------

print( BOLD, GREEN, "\n--\n", RESET );
print( BOLD, "Loading Mapping File... ", RESET );

unless( open( INFILE2, $mappingFile ) )
{ die( "\n\nUnable to open Mapping File: $mappingFile.  Please verify.\n\n" ); }

# Mapping Hash for Transcript IDs to their Gene Ids
# KEY: TranscriptID
# VAL: GeneID
my %mappings = ();

my $numberOfLines2 = 0;

while( my $line = <INFILE2> )
{
    chomp($line);

    $line   =~ s/\r//g;
    $line   =~ s/>//g;

    my @lineArray = split(/ /,$line);

    my $transcriptID = $lineArray[ 0 ];
    my $geneID       = $lineArray[ 1 ];
    $geneID =~ s/gene://g;

    $mappings{ $transcriptID } = $geneID;

    $numberOfLines2++;
}

close(INFILE2);

print( CYAN, "\n Number of lines processed in mapping file: ", RESET );
print( " " . &addCommas($numberOfLines2) . "\n" );


# --------------------------------------------------- Associations ----------------------------------------------------
#
#   Loop through the Query Transcript ID's in the BLAST OUTPUT and create the proper mapping files of:
#
#       QueryGene       -- SubjectGene
#       QueryTranscript -- Subject Transcript
#

print( BOLD, GREEN, "\n--\n", RESET );
print( BOLD, "Generating Mappings... ", RESET );

for my $queryTranscriptId ( sort{ $a cmp $b || $a <=> $b } keys %blastHits )
{
    if( defined( $mappings{ $queryTranscriptId } ) ) {

        my $geneIdForQueryTranscript    = $mappings{ $queryTranscriptId }; # Query Transcript GeneID
        my $geneIdFromBlastSubjectHit   = $blastHits{ $queryTranscriptId }; # Subject GeneID from Query Transcript
        my $subjectTranscriptID         = $queryToSubject{ $queryTranscriptId }; # Subject TranscriptID
        my $subjectGeneBiotype          = $subjectBiotype{ $geneIdFromBlastSubjectHit }; # Subject Gene Biotype

        print OUTFILE ( "$geneIdForQueryTranscript\t$geneIdFromBlastSubjectHit\n" );
        print OUTFILE2 ( "$queryTranscriptId\t$subjectTranscriptID\n" );
        print OUTFILE3 ( "$geneIdFromBlastSubjectHit\t$subjectGeneBiotype\n" );
        print OUTFILE4 ( "$geneIdForQueryTranscript\t$queryTranscriptId\t$subjectTranscriptID\t$geneIdFromBlastSubjectHit\t$subjectGeneBiotype\n" );
    }
}

my $numberOfSubjects = scalar keys %subjects;
my $numberOfBlastHits = scalar keys %blastHits;
my $numberOfQueryToSubjectMappings = scalar keys %queryToSubject;
my $numberOfSubjectBiotypes = scalar keys %subjectBiotype;

print( "\n" );
print( " • Number of Unique TranscriptID Subjects: $numberOfSubjects\n" );
print( " • Number of Blast Hits: $numberOfBlastHits\n" );
print( " • Number of Query to Subject Transcript Mappings: $numberOfQueryToSubjectMappings\n" );
print( " • Number of Unique Subject Biotypes: $numberOfSubjectBiotypes\n" );

# ------------------------------------------------------ Evalues ------------------------------------------------------
#
#	Evalue Ratings
# 	We also want to answer the question of "for each parent gene, what are their children, and how good are their mappings
# 	Maps a ParentGeneID to a collection of (hashes) of child query genes

print( BOLD, GREEN, "\n--\n", RESET );
print( BOLD, "Rating Evalue scores... ", RESET );


for my $aParentGene ( sort{ $a cmp $b || $a <=> $b } keys %subjectEvalueHashAsTranscripts )
{
    my %childHits = %{ $subjectEvalueHashAsTranscripts{ $aParentGene } };

    my $numberOfChildHits += scalar keys %childHits;

    for my $aChildHit ( sort{ $a cmp $b || $a <=> $b } keys %childHits )
    {
        my $evalue = $childHits{ $aChildHit };

        my $geneIDForChildHit = $mappings{ $aChildHit };

        print OUTFILE5 ( "$aParentGene\t$geneIDForChildHit\t$evalue\n" );

        if ( $numberOfChildHits <= 1 ) {
            print OUTFILE6 ( "$aParentGene\t$geneIDForChildHit\t$evalue\n" );
        }
    }
}



close(OUTFILE);
close(OUTFILE2);
close(OUTFILE3);
close(OUTFILE4);
close(OUTFILE5);
close(OUTFILE6);

print("\n\n");

# Re-flush STDOUT
$| = 1;



# ------------------------------------------------ Functions & Methods ------------------------------------------------

#	addCommas
#	Simple method to format a number with commas
#
sub addCommas
{
    local($_) = shift;
    1 while s/^(-?\d+)(\d{3})/$1,$2/;
    return $_;
}


# Exit the program with a normal condition flag
exit( 0 );