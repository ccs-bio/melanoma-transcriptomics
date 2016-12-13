#
# R script to compute the differentially expressed genes using the
# Bioconductor package NOISeq
#

#---------------------------------------------------------------------------------------------
#
#                                         Main


library(NOISeq)

# Base Directory
setwd("/path/to/base/analysis/directory/noiseq")

print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Starting... ", sep="") )

algorithm = "noiseq"

pathForQuantResults = "/path/to/base/analysis/directory/results-103"

arrayForAlingmentTypes = c( 'counts' )
arrayForDataTypes = c( 'unique' )
arrayForFileType = c( 'Filtered' )
arrayForResultsTypes = c( 'results_pooled' )

#---------------------------------------------------------------------------------------------
#
#	Factors for Pooled Interpretation
#
factorsForPooledInterpretation = data.frame( Pooled = c("T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","C","C","C") )


#---------------------------------------------------------------------------------------------
#
#	Build & Analyze
#
for( alnType in arrayForAlingmentTypes ) {

	print( "" )
	print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Processing ", alnType, sep="") )

	for( dataType in arrayForDataTypes ) {

		for( fileType in arrayForFileType ) {	# O(3) anyone?

			print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Processing ", alnType, " - ", dataType, sep="") )

			inputFileName = paste( "aggregated_expression_profile_for_", alnType,"_AllGenes.txt", sep="" )
			inputFilePath = paste( pathForQuantResults, "/", alnType, "/", dataType, "/", fileType, "/", inputFileName, sep="" )

			countsTable = read.delim( inputFilePath, header = TRUE, row.names = 1 )
			countsTable$No..Samples = NULL # Delete the last column (contains frequency info)

			dataSetTable = NOISeq::readData( data = countsTable, factors = factorsForPooledInterpretation )

			param_k = 0.5					# Counts are thresholded to 0.5, while FPKMs are thresholded to 1.0
			param_norm = "tmm"				# Normalization
			param_factor = "Pooled"			# Experimental interpretation we're using
			param_replicates = "biological"	# Only used when we use the 'technical' replicates function

			if( alnType == 'fpkm' ) {
				param_k = 1.0
				param_norm = "n"
			}

			print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] with Params: K=", param_k, ", Norm=", param_norm, sep="") )

			# TECHNICAL replicates
			noiseqResults = noiseq( dataSetTable, k = param_k, norm = param_norm, factor = param_factor, replicates = param_replicates )

			# BIOLOGICAL replicates
			# noiseqResults = noiseqbio( dataSetTable, k = param_k, norm = param_norm, factor = param_factor )

			noiseqResultsFiltered = degenes( noiseqResults, q = 0.6, M = NULL )

			outputFile = paste( "", alnType, "/", "results_pooled", "/", dataType, "/", fileType, "/", "AllGenes_SigDiffExp-", dataType, ".txt", sep="" )

			write.table( noiseqResultsFiltered, file=outputFile, append=FALSE, sep="\t" )

		}

	} # End loop for data type

} # End loop for alignment type