#
# R script to compute the differentially expressed genes using the
# Bioconductor package NOISeq
#

#---------------------------------------------------------------------------------------------
#
#                                         Main


library(NOISeq)

# Base Directory
setwd("/Users/camilo/Documents/CCS/Enrico/Projects/Melanoma/Analysis_Full/differential_expression_103/_counts_exonic/noiseq")

print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Starting... ", sep="") )

algorithm = "noiseq"

pathForQuantResults = "/Users/camilo/Documents/CCS/Enrico/Projects/Melanoma/Analysis_Full/results-103"

arrayForAlingmentTypes = c( 'counts_exonic' )

arrayForDataTypes = c( 'unique' )
# arrayForFileType = c( 'Filtered_10_10', 'Filtered_66_100', 'Filtered_10_1000', 'Filtered_10_500' )
# arrayForFileType = c( 'Filtered_10_1000', 'Filtered_10_500' )
# arrayForFileType = c( 'Filtered_10_100' )
# arrayForFileType = c( 'Filtered_10_10', 'Filtered_10_100', 'Filtered_10_500', 'Filtered_10_1000', 'Filtered_66_100' )
arrayForFileType = c( 'Filtered_5_10', 'Filtered_20_10' )

arrayForResultsTypes = c( 'results_pooled' )

#---------------------------------------------------------------------------------------------
#
#	Factors for Pooled Interpretation
#
factorsForPooledInterpretation = data.frame( Pooled = c("T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","C","C") )


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

			noiseqResults = noiseq( dataSetTable, k = param_k, norm = param_norm, factor = param_factor, replicates = param_replicates )

			noiseqResultsFiltered = degenes( noiseqResults, q = 0.6, M = NULL )

			outputFile = paste( "", "counts", "/", "results_pooled", "/", dataType, "/", fileType, "/", "AllGenes_SigDiffExp-", dataType, ".txt", sep="" )

			write.table( noiseqResultsFiltered, file=outputFile, append=FALSE, sep="\t" )
			
			print( head( noiseqResultsFiltered, n=10 ) )
			
			
			# --------------------------------------------- Hierarhical Clustering ---------------------------------------------
			
			baseDirFigures = '/Users/camilo/Documents/CCS/Enrico/Projects/Melanoma/Analysis_Full/Figures_Numbers-Full_103/differential_expression/counts_exonic/noiseq'

			selectedGenesForNOISeq = noiseqResultsFiltered[noiseqResultsFiltered$prob > .80,]
			selectedGenesForNOISeq = selectedGenesForNOISeq[selectedGenesForNOISeq$M > 2,]
			
			print( head( selectedGenesForNOISeq, n=10 ) )
			
			figureFileNameForNOISeq = paste( fileType, ".png", sep="" )
			png( file=paste( baseDirFigures, "/", figureFileNameForNOISeq, sep="" ), width=1024, height=1024, units="px", res=120 )
			
			dataFrameForNOISeq = as.data.frame.matrix( countsTable[row.names(selectedGenesForNOISeq),] )
			
			heatmap( as.matrix( dataFrameForNOISeq[complete.cases(dataFrameForNOISeq),] ) )
			
			dev.off()

		}

	} # End loop for data type

} # End loop for alignment type