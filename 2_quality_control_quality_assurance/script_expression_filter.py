#!/usr/bin/python
# coding: utf-8

#-----------------------------------------------------------------------------------------------------------
#
#                                   Center for Computational Science
#										http://www.ccs.miami.edu/
#                             			  University of Miami
#
#   This software is a "University of Miami Work" under the terms of the United States Copyright Act.
#   Please cite the author(s) in any work or product based on this material.
#
#   OBJECTIVE:
#	The purpose of this program is to filter a table file with expression results (fpkm or counts) based
#   on the (a) number of samples, and (b) the minimum expression value for a given gene in a collection
#   of samples.
#
#   For a gene to be present it has to have a value higher than a given minimum "T" in "N"  treatment
#   / tumor samples (as defined by the experimental condition).  A gene has to be at least present in
#   both experimental conditions in order for it to be called 'present'.  However, if a given gene has
#   an expression value of "0" in all the treatment/tumor samples, but has a value (above "T") in the
#   control samples, then the gene is kept.
#
#   NOTES:
#   Please see the dependencies section below for the required libraries (if any).
#
#   DEPENDENCIES:
#
#       • math, csv, re
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
#-----------------------------------------------------------------------------------------------------------

# 	Python Modules
import os, sys
import argparse
import time
import math
import csv
import re

#-----------------------------------------------------------------------------------------------------------
#
#													Main
#

sys.stdout.flush()

print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + " Python Starting" + "" )

# 	Pick up the command line arguments
parser = argparse.ArgumentParser()
parser.add_argument( "-i", "--input", required=True, type=str, help="Input File" )
parser.add_argument( "-e", "--interpretation", required=True, type=str, help="Experimental Interpretation" )
parser.add_argument( "-m", "--min_exp_val", required=True, type=float, help="Minimum Expression Value" )
parser.add_argument( "-p", "--presence", required=True, type=int, help="Percentage (%) of samples that a gene needs to have a value higher than 'm'." )
parser.add_argument( "-a", "--affix", required=False, type=str, help="Affix for Output files" )
parser.add_argument( "-o", "--out", required=False, type=str, help="Output Directory" )
parser.add_argument( "-c", "--ommit-last-column", action="store_true", default=False, dest="boolean_flag_for_last_column", required=False, help="BOOL Flag")
parser.add_argument( "--fpkm", action="store_true", default=False, dest="boolean_flag_for_fpkm", required=False, help="FPKM BOOL Flag")
parser.add_argument( "--counts", action="store_true", default=False, dest="boolean_flag_for_counts", required=False, help="Counts BOOL Flag")
args = parser.parse_args()

# 	Variables initialized from the command line arguments
filePathForInputFile    = args.input
expInterpretation       = args.interpretation
minExpressionValue      = args.min_exp_val
presenceValue           = args.presence
affixForOutputFile      = args.affix if args.affix is not None else "default_out"
output_directory        = args.out if args.out is not None else "./out"
ommitLastColumn         = args.boolean_flag_for_last_column
flagForUsingFPKM        = args.boolean_flag_for_fpkm
flagforUsingCounts      = args.boolean_flag_for_counts

# Remove any whitespaces around the file paths
filePathForInputFile.strip()

#-----------------------------------------------------------------------------------------------------------
#
#										Output Files & Directories
#

# We'll check if the output directory exists — either the default (current) or the requested one
if not os.path.exists( output_directory ):
    print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + " Output directory does not exist.  Creating..." + "" )
    os.makedirs( output_directory )

outputFile = output_directory + "/aggregated_expression_profile_for_" + affixForOutputFile + "_AllGenes.txt"

logFile = output_directory + "/filter_log.txt"
writerForLogFile = csv.writer( open( logFile, "wb" ), delimiter='\t' )

#-----------------------------------------------------------------------------------------------------------
#
#   Experimental Interpretation
#

# Get some initial stats on the number of samples we'll require to call a gene present
numberOfControlSamples   = expInterpretation.count('C')
numberOfTreatmentSamples = expInterpretation.count('T')

numberOfRequiredTreatmentSamples = (float( numberOfTreatmentSamples ) * float( presenceValue )) / 100.0
numberOfRequiredTreatmentSamples = math.ceil( numberOfRequiredTreatmentSamples )

numberOfRequiredControlSamples   = (float( numberOfControlSamples ) * float( presenceValue )) / 100.0
numberOfRequiredControlSamples   = math.ceil( numberOfRequiredControlSamples )

# And convert the Interpretation String into a list
listFromStringInterpretation = expInterpretation.split(',')

# This will hash an index as 'KEY' and type character as 'VALUE' so that we can cycle through the file columns
# and find the column's experimental condition
dictionaryForInterpretationIndices = {}

for (index, charType) in enumerate(listFromStringInterpretation):
    if( charType == "NONE" ):
        continue
    dictionaryForInterpretationIndices[ index ] = charType

# What type of expression metric are we using
weAreUsingString = "Using "

if(flagForUsingFPKM):
    weAreUsingString = weAreUsingString + " FPKM"

if(flagforUsingCounts):
    weAreUsingString = weAreUsingString + " Counts"

print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + " " + weAreUsingString )
print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + " Calculating Presence Requirements..." + "" )
print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + " Interpretation: " + expInterpretation )
print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + " No. Controls: " + '{:0,.0f}'.format(numberOfControlSamples) )
print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + " No. Treatments: " + '{:0,.0f}'.format(numberOfTreatmentSamples) )
print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + " Presence in Conditions: " + '{:0,}'.format(presenceValue) + "%" )
print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + " No. Control Samples Required: " + '{:0,.0f}'.format(numberOfRequiredControlSamples) )
print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + " No. Treatment Samples Required: " + '{:0,.0f}'.format(numberOfRequiredTreatmentSamples) )
print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + " Min Expression: " + '{:0,.3f}'.format(minExpressionValue) )

writerForLogFile.writerow( [ "", weAreUsingString ] )
writerForLogFile.writerow( [ " Calculating Presence Requirements...", "" ] )
writerForLogFile.writerow( [ " Interpretation", expInterpretation ] )
writerForLogFile.writerow( [ " No. Controls: ", '{:0,.0f}'.format(numberOfControlSamples) ] )
writerForLogFile.writerow( [ " No. Treatments: ", '{:0,.0f}'.format(numberOfTreatmentSamples) ] )
writerForLogFile.writerow( [ " Presence in Conditions: ", '{:0,}'.format(presenceValue), "%" ] )
writerForLogFile.writerow( [ " No. Control Samples Required: ", '{:0,.0f}'.format(numberOfRequiredControlSamples) ] )
writerForLogFile.writerow( [ " No. Treatment Samples Required: ", '{:0,.0f}'.format(numberOfRequiredTreatmentSamples) ] )
writerForLogFile.writerow( [ " Min Expression Value: ", '{:0,.3f}'.format(minExpressionValue) ] )


#-----------------------------------------------------------------------------------------------------------
#
#												File Loading
#

print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + " Loading Files..." + "" )
writerForLogFile.writerow( [ " Loading Files...", "" ] )

print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + filePathForInputFile + "" )
writerForLogFile.writerow( [ "", filePathForInputFile ] )


numberOfLinesInFile                             = 0
numberOfGenesPassingFilterOverall               = 0
numberOfGenesPassingFilterByControlCondition    = 0
numberOfGenesPassingFilterByTreatmentCondition  = 0

with open( filePathForInputFile,'r' ) as INFILE:

    reader = csv.reader( INFILE, delimiter='\t' )
    writer = csv.writer( open( outputFile, "wb" ), delimiter='\t', lineterminator="\n" )

    # Input File Header (THIS IS A LIST! Not a STRING with the contents of the line)
    inputFileHeader = next(reader, None)

    # if inputFileHeader:
#         if ( ommitLastColumn ):
#         inputFileHeader[ len(inputFileHeader) -1 ] = None # Remove the last column header


    writer.writerow( inputFileHeader )

    try:
        for row_line in reader:

            # Threshold is 'minExpressionValue'
            valuesGreaterThanThresholdForControlConditions   = 0
            valuesGreaterThanThresholdForTreatmentConditions = 0

            for( index, valueOfCell ) in enumerate( row_line ):

                if ( index == 0 ):  # Skip the GeneID column
                    continue

                # Not interested in the last column -- it contains information we'll recalculate (just to double-check)
                # if ( ommitLastColumn ):
#                     if( index == ( len( row_line ) - 1 ) ):
#                         continue

                # How many times does a gene in the 'T'/'Tumor'/'Treatment' condition above threshold
                if ( dictionaryForInterpretationIndices[ index ] == "T" ):
                    if (flagForUsingFPKM):
                        if ( float(valueOfCell) >= minExpressionValue ):
                            valuesGreaterThanThresholdForTreatmentConditions += 1

                    if (flagforUsingCounts):
                        if ( int(valueOfCell) >= minExpressionValue ):
                            valuesGreaterThanThresholdForTreatmentConditions += 1

                # How many times does a gene in the 'C'/'Control' condition above threshold
                if ( dictionaryForInterpretationIndices[ index ] == "C" ):
                    if (flagForUsingFPKM):
                        if ( float(valueOfCell) >= minExpressionValue ):
                            valuesGreaterThanThresholdForControlConditions += 1

                    if (flagforUsingCounts):
                        if ( int(valueOfCell) >= minExpressionValue ):
                            valuesGreaterThanThresholdForControlConditions += 1



            diagnosingString = ""

            # Genes that pass the threshold in BOTH conditions will be kept.  BUT, if a gene is 'present' in
            # the control samples, and has a value of zero ('0') in the treatments, it is also 'present.
            # Also keep track of the overall progress, as well as write-out the line that passed

            if ( ( ( valuesGreaterThanThresholdForTreatmentConditions >= numberOfRequiredTreatmentSamples ) or \
                    ( valuesGreaterThanThresholdForControlConditions >= numberOfRequiredControlSamples ) ) or \
                    ( ( valuesGreaterThanThresholdForControlConditions >= numberOfRequiredControlSamples ) and \
                    ( valuesGreaterThanThresholdForTreatmentConditions == 0 ) ) ):

                # row_line[ len(row_line) - 1 ] = None
                writer.writerow( row_line )
                numberOfGenesPassingFilterOverall += 1
                diagnosingString = "PASS — "

            # Statistics particular for Control samples
            if(valuesGreaterThanThresholdForControlConditions >= numberOfRequiredControlSamples):
                numberOfGenesPassingFilterByControlCondition += 1

            # Statistics particular for Treatment samples
            if(valuesGreaterThanThresholdForTreatmentConditions >= numberOfRequiredTreatmentSamples):
                numberOfGenesPassingFilterByTreatmentCondition += 1


            # print( diagnosingString + "Controls: " + str(valuesGreaterThanThresholdForControlConditions) + ", Treatments: " + str(valuesGreaterThanThresholdForTreatmentConditions) )

            numberOfLinesInFile += 1

    except csv.Error as e:
        sys.exit( "File %s, line %d: %s" % ( aFile, reader.line_num, e ) )


percentageOfGenesPassingFilterOverall               = ( float(numberOfGenesPassingFilterOverall) / float(numberOfLinesInFile) ) * 100.0
percentageOfGenesPassingFilterByControlCondition    = ( float(numberOfGenesPassingFilterByControlCondition) / float(numberOfLinesInFile) ) * 100.0
percentageOfGenesPassingFilterByTreatmentCondition  = ( float(numberOfGenesPassingFilterByTreatmentCondition) / float(numberOfLinesInFile) ) * 100.0

print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + " No. of Lines in file: " + '{:0,.0f}'.format(numberOfLinesInFile) + "" )
print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + " No. of Genes: " + '{:0,.0f}'.format(numberOfLinesInFile) + "" )
print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + " -------------------------------" )

print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + " No. of Genes Passing Filter Overall: " + '{:0,.0f}'.format(numberOfGenesPassingFilterOverall) + "" )
print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + " Percentage of Genes Passing Filter Overall: " + '{:0,.000f}'.format(percentageOfGenesPassingFilterOverall) + "%" )
print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + " -------------------------------" )

print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + " No. of Genes Passing Filter by Controls: " + '{:0,.0f}'.format(numberOfGenesPassingFilterByControlCondition) + "" )
print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + " Percentage of Genes Passing Filter by Controls: " + '{:0,.000f}'.format(percentageOfGenesPassingFilterByControlCondition) + "%" )
print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + " -------------------------------" )

print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + " No. of Genes Passing Filter by Treatments: " + '{:0,.0f}'.format(numberOfGenesPassingFilterByTreatmentCondition) + "" )
print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + " Percentage of Genes Passing Filter by Treatments: " + '{:0,.000f}'.format(percentageOfGenesPassingFilterByTreatmentCondition) + "%" )
print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + " -------------------------------" )

writerForLogFile.writerow( [ " No. of Lines in file: ", '{:0,.0f}'.format(numberOfLinesInFile) ] )
writerForLogFile.writerow( [ " No. of Genes: ", '{:0,.0f}'.format(numberOfLinesInFile) ] )

writerForLogFile.writerow( [ " No. of Genes Passing Filter Overall: ", '{:0,.0f}'.format(numberOfGenesPassingFilterOverall) ] )
writerForLogFile.writerow( [ " Percentage of Genes Passing Filter Overall: ", '{:0,.00f}'.format(percentageOfGenesPassingFilterOverall) ] )

writerForLogFile.writerow( [ " No. of Genes Passing Filter by Controls: ", '{:0,.0f}'.format(numberOfGenesPassingFilterByControlCondition) ] )
writerForLogFile.writerow( [ " Percentage of Genes Passing Filter by Controls: ", '{:0,.00f}'.format(percentageOfGenesPassingFilterByControlCondition) ] )

writerForLogFile.writerow( [ " No. of Genes Passing Filter by Treatments: ", '{:0,.0f}'.format(numberOfGenesPassingFilterByTreatmentCondition) ] )
writerForLogFile.writerow( [ " Percentage of Genes Passing Filter by Treatments: ", '{:0,.00f}'.format(percentageOfGenesPassingFilterByTreatmentCondition) ] )

#-----------------------------------------------------------------------------------------------------------
#
#                                               End of Line
#

print( "[ " + time.strftime('%d-%b-%Y %H:%M:%S',time.localtime()) + " ]" + " Done." + "\n" )
writerForLogFile.writerow( [ " Done.", "" ] )

sys.exit()
