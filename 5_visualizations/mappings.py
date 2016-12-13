#!/usr/bin/python
# coding: utf-8

# ---------------------------------------------------------------------------------------------------------------------
#
#                                   Center for Computational Science
#										http://www.ccs.miami.edu/
#                             			  University of Miami
#
#   This software is a "University of Miami Work" under the terms of the United States Copyright Act.
#   Please cite the author(s) in any work or product based on this material.
#
#   OBJECTIVE:
#	The purpose of this program is to map a Parental gene to a collection of Child Pseudogenes.
#
#
#   NOTES:
#   Please see the dependencies section below for the required libraries (if any).
#
#   DEPENDENCIES:
#
#       • NONE at this time.
#
#   The above libraries & modules are required. You can check the modules currently installed in your
#   system by running: python -c "help('modules')"
#
#   USAGE:
#   Run the program with the "--help" flag to see usage instructions.
#
#	AUTHOR:	Camilo Valdes (cvaldes3@med.miami.edu)
#			Computational Biology and Bioinformatics Group
#			Center for Computational Science (CCS), University of Miami
#
# ---------------------------------------------------------------------------------------------------------------------

# 	Python Modules
import os, sys
import argparse
import time
import math
import csv
import re

#------------------------------------------------------ Main ----------------------------------------------------------

sys.stdout.flush()

print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + "" )

# 	Pick up the command line arguments
parser = argparse.ArgumentParser()
parser.add_argument( "-m", "--mappingFile", required=True, type=str, help="Input File with Mappings" )
parser.add_argument( "-i", "--inputFile", required=True, type=str, help="Parental Genes of Interest")
parser.add_argument( "-o", "--out", required=False, type=str, help="Output Directory" )
args = parser.parse_args()

# 	Variables initialized from the command line arguments
filePathForMappingFile  = args.mappingFile
filePathForParentalGenesOfInterest = args.inputFile
output_directory        = args.out if args.out is not None else "./out"

#   Remove any whitespaces around the file paths
filePathForMappingFile.strip()
filePathForParentalGenesOfInterest.strip()

#   Maps a Parental gene to a nested array of pseudogenes
dictionaryForParentGeneMappings = {}

#------------------------------------------- Output Files & Directories -----------------------------------------------
#
# We'll check if the output directory exists — either the default (current) or the requested one
#
if not os.path.exists( output_directory ):
    print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + " Output directory does not exist.  Creating..." + "" )
    os.makedirs( output_directory )

outputFile = output_directory + "/mappings.txt"
writer = csv.writer( open( outputFile, "wb" ), delimiter='\t', lineterminator="\n" )
writer.writerow( [ 'Parentgene', 'Pseudogenes' ] )


# ------------------------------------------------- Mapping File ----------------------------------------------------

print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + "" )
print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + " Loading Mapping File..." + "" )

numberOfLinesInFile = 0

with open( filePathForMappingFile,'r' ) as INFILE:

    reader = csv.reader( INFILE, delimiter='\t' )

    # Input File Header
    inputFileHeader = next(reader, None)

    if inputFileHeader:
        next

    try:
        for row_line in reader:

            parentGeneID = row_line[ 1 ]
            pseudogeneID = row_line[ 0 ]

            #  ------------------ ParentGene Mappings ------------------
            if( parentGeneID in dictionaryForParentGeneMappings ):
                dictionaryForParentGeneMappings[ parentGeneID ].append(pseudogeneID)

            if( parentGeneID not in dictionaryForParentGeneMappings):
                dictionaryForParentGeneMappings[ parentGeneID ] = [ pseudogeneID ]

            numberOfLinesInFile += 1

    except csv.Error as e:
        sys.exit( "File %s, line %d: %s" % ( aFile, reader.line_num, e ) )

print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + " Number of Lines in file: " + '{:0,.0f}'.format(numberOfLinesInFile) + "" )
print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + "" )


# --------------------------------------------------- Input File ------------------------------------------------------
#
#   Parental genes of interest
#
print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + " Loading Parentals to extract..." + "" )

numberOfGeneIdsInFile_1 = 0

with open( filePathForParentalGenesOfInterest,'r' ) as INFILE_FOR_FILE_1:

    reader = csv.reader( INFILE_FOR_FILE_1, delimiter='\t' )

    try:
        for row_line in reader:
            
            newRowArray = []

            parentGeneID = row_line[ 0 ]

            arrayWithChildPseudogenes = dictionaryForParentGeneMappings[ parentGeneID ]
            
            newRowArray.append( parentGeneID )
            
            for childPseudogeneID in arrayWithChildPseudogenes:
                newRowArray.append( childPseudogeneID )

            writer.writerow( newRowArray )
            numberOfGeneIdsInFile_1 += 1

    except csv.Error as e:
        sys.exit( "File %s, line %d: %s" % ( aFile, reader.line_num, e ) )


print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + " Number of Pseudogenes (File 1): " + '{:0,.0f}'.format(numberOfGeneIdsInFile_1) + "" )
print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + "" )







#--------------------------------------------------- End of Line ------------------------------------------------------

print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + "" )
print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + "" )

sys.exit()
