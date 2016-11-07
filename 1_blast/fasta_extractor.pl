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
#	The purpose of this program is to generate a fasta file from a series of Ensembl IDs.
#
#
#   NOTES:
#   Please see the dependencies section below for the required libraries (if any).
#
#
#   DEPENDENCIES:
#
#       • Getopt::Long
#		• POSIX
#		• Bio::Perl
#		• Bio::DB::Sam
#		• Bio::SeqIO
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
use Bio::Perl;
use Bio::Seq;
use Bio::SeqIO;
use Bio::DB::Sam;

# ------------------------------------------------------- Main --------------------------------------------------------

my $bundle = "0.1";
my $build  = "A0001";

my $file;               # Input File with Ensembl Gene Ids
my $outputDir;			# Output Directory
my $cdna_file;          # File that holds all the Fasta records
my $prefix;             # Prefix name for output file

GetOptions ( 'i=s' => \$file, 'o:s' => \$outputDir, 'c:s' => \$cdna_file, 'p:s' => \$prefix );

if( ! defined( $file) )
{ print( RED, "\n\nWARNING! Missing Files ", RESET ); print( ".. $!\n\n" ); exit(1); }

chomp($file);

my $currentWD = getcwd();       # Get a pointer to the Current Working Directory

# Default to the current working directory if none provided
if( !defined( $outputDir ) )
{ $outputDir = $currentWD; }

if( !defined( $prefix ) )
{ $prefix = "fasta_file"; }

$| = 1;  # Flush STDOUT

print( "\nProgram Name, ($bundle, $build)\n" );
print( "===========================================================\n" );

# ------------------------------------------------------ Output File --------------------------------------------------

my $out1 = $outputDir . "/" . $prefix . "-cdna.fasta";

my %outfiles = ( 'out1' => Bio::SeqIO->new('-file' => ">$out1", '-format' => 'fasta') );

# --------------------------------------------------- GeneID File Loader ----------------------------------------------

print( BOLD, GREEN, "\n--\n", RESET );
print( BOLD, "Loading file(s) ... ", RESET );

unless( open( INFILE, $file ) )
{ die( "\n\nUnable to open Input File: $file.  Please verify.\n\n" ); }

# Hash to store GeneIDs — we are only interested in the keys
# KEY: GeneID
# VAL: Expression Meassure
my %geneIds = ();

my $numberOfLines = 0;

while( my $line = <INFILE> )
{
    chomp($line);
	$line   =~ s/\r//g;
    my @lineArray = split(/\t/,$line);

    my $ensemblGeneID       = $lineArray[ 0 ];

    $geneIds{ $ensemblGeneID } = 1;

    $numberOfLines++;
}

close(INFILE);

print( CYAN, "\n Number of lines processed in input file: ", RESET );
print( " " . &addCommas($numberOfLines) . "\n" );


# --------------------------------------------------- FASTA Loader ----------------------------------------------------

my $inseq = Bio::SeqIO->new('-file' => "<$cdna_file", '-format' => 'fasta');

while( my $seqin = $inseq->next_seq )
{
    my $fasta_header   = $seqin->desc;

    my $transcript_id  = $seqin->id;
    my $sequence       = $seqin->seq();

    # We'll parse the fasta header to extract the GeneID
    my @arrayForFastaHeader = split(/ /,$fasta_header);
    my @arrayForGeneToken   = split(/:/,$arrayForFastaHeader[ 2 ]);

    my $gene_id = $arrayForGeneToken[ 1 ];

    if( defined( $geneIds{ $gene_id } ) ) {
        $outfiles{'out1'}->write_seq($seqin);
    }

}


print("\n\n");

# Re-flush STDOUT
$| = 1;



# -------------------------------------------------- Functions & Methods ----------------------------------------------

#	addCommas
#	Simple method to format a number with commas
sub addCommas
{
    local($_) = shift;
    1 while s/^(-?\d+)(\d{3})/$1,$2/;
    return $_;
}


# Exit the program with a normal condition flag
exit( 0 );