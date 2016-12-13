
# ---------------------------------------------------------------------------------------------------------------------
#
#
#                                   Center for Computational Science
#										http://www.ccs.miami.edu/
#                             			  University of Miami
#
#   This software is a "University of Miami Work" under the terms of the United States Copyright Act.
#   Please cite the author(s) in any work or product based on this material.
#
#   OBJECTIVE:
#	The purpose of this program is to fit a linear model onto the FC values for DEGs.
#
#
#
#	AUTHOR:	Camilo Valdes (cvaldes3@miami.edu)
#			Computational Biology and Bioinformatics Group
#			Center for Computational Science (CCS), University of Miami
#
# ---------------------------------------------------------------------------------------------------------------------

# 	Import any necessary libraries here

library(lars)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(data.table)
library(grid)
library(nnet)
library(MASS)
library(Rfit)


# 	Enable printing of large numbers. R defaults to "rounding" for display, and
# 	Cause the warnings to be printed as they occur.
options( width=220, digits=15, warn=1 )


# ------------------------------------------------- Project Setup -----------------------------------------------------

print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] ", sep="") )
print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] R Starting... ", sep="") )
print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] ", sep="") )


working_directory="/path/to/base/analysis/directory"
setwd( file.path( working_directory ) )

ggplot_output_dir="/path/to/base/analysis/directory/figures/lrm"

dir.create( file.path( working_directory, ggplot_output_dir ), showWarnings=FALSE, recursive=TRUE )



# ------------------------------------------------------ Main ---------------------------------------------------------

algorithms = c( 'counts_exonic' )

algorithm_names = c( 'DESeq', 'NOISeq', 'CuffDiff', 'Limma', 'GeneSpring' )


for( anAlgorithm in algorithms ) {

	fileWithFCsForDEGs		= paste( working_directory, "/merged_data/merged_data_", anAlgorithm, ".txt", sep="" )
	tableWithFCsForDEGs		= read.delim( fileWithFCsForDEGs, header=TRUE, row.names=1, sep="\t", check.names=FALSE )
	
	#	Convert the input table into a Data Frame, and change 'NA' to '0'
	dataFrameWithFCsForDEGs = as.data.frame(tableWithFCsForDEGs)
	dataFrameWithFCsForDEGs[is.na(dataFrameWithFCsForDEGs)] = 0
	
	colnames( dataFrameWithFCsForDEGs ) = algorithm_names
	
	print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Data Frame with Zeros", sep="") )
	print(head(dataFrameWithFCsForDEGs,n=20))
	print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] ", sep="") )
	
	#	Remove rows with all zeros
	dataFrameWithFCsForDEGs = dataFrameWithFCsForDEGs[apply( dataFrameWithFCsForDEGs[,-1], 1, function(x){ !all( as.numeric(x)==0 ) } ),]
	# rowSub = apply(dataFrameWithFCsForDEGs, 1, function(row) all(row !=0 ))
	#dataFrameWithFCsForDEGs = dataFrameWithFCsForDEGs[rowSub,]
	
	print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Data Frame with No Zeros Rows", sep="") )
	print(head(dataFrameWithFCsForDEGs,n=20))
	print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] ", sep="") )
	
	
	# -------------------------------------------------------------------------------------------------------------------------
	#
	#	Linear Model — lm()
	#
	
	# columnIndexWithResponseAlgorithm = 1	#	DESeq
	# columnIndexWithResponseAlgorithm = 2	#	NOISeq
	columnIndexWithResponseAlgorithm = 3	# 	CuffDiff
	# columnIndexWithResponseAlgorithm = 4	#	Limma
	# columnIndexWithResponseAlgorithm = 5	#	GeneSpring
	
	#
	#	Prepare the response vector
	#
	print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Response:", sep="") )
	responseDataFrame = dataFrameWithFCsForDEGs[ columnIndexWithResponseAlgorithm ]
	
	print( head( responseDataFrame, n=20) )
	print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] ", sep="") )
	
	#	Eliminate those rows with a Zero
	nonZeroRowSubset = apply(responseDataFrame, 1, function(row) all(row !=0 )) # Go through each row and determine if a value is zero
	responseDataFrame = responseDataFrame[nonZeroRowSubset,, drop=FALSE]		# drop=FALSE retains rownames! :)
	
	print( head( responseDataFrame, n=20) )
	
	print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Number of rows: ", nrow(responseDataFrame), sep="") )
	
	namesOfGenesInResponse = rownames(responseDataFrame)
	
	print( head( namesOfGenesInResponse, n=20) )
	
	
	#
	#	Prepare the Predictor Explanatory Matrix
	#
	print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] ", sep="") )
	print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Predictor:", sep="") )
	
	dataFrameWithJoinedResponse = dataFrameWithFCsForDEGs[namesOfGenesInResponse,]
	
		
	#
	#	Run the lm() function to get the linear model
	#
	
	print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] ", sep="") )
	print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] -------------------------------------------------------------", sep="") )
	print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Running lm()...", sep="") )
	print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] ", sep="") )
	
	
	# fittedDEGs = lm( CuffDiff ~ DESeq + NOISeq + Limma + GeneSpring, data=dataFrameWithFCsForDEGs )
	
	fittedDEGs = lm( CuffDiff ~ DESeq + NOISeq + Limma + GeneSpring, data=dataFrameWithJoinedResponse )
	# fittedDEGs = lm( NOISeq ~ DESeq + CuffDiff + Limma + GeneSpring, data=dataFrameWithJoinedResponse )
	# fittedDEGs = lm( DESeq ~ NOISeq + CuffDiff + Limma + GeneSpring, data=dataFrameWithJoinedResponse )
	# fittedDEGs = lm( Limma ~ DESeq + NOISeq + CuffDiff + GeneSpring, data=dataFrameWithJoinedResponse )
	# fittedDEGs = lm( GeneSpring ~ DESeq + NOISeq + CuffDiff + Limma, data=dataFrameWithJoinedResponse )
	
	print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] -------------------------------------------------------------", sep="") )
	print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Summary:", sep="") )
	
	print(summary(fittedDEGs))	
	
	
	#
	#	After inspecting the Coefficients, we run again the regression (lm) with only the SIGNIFICANT coefficients
	#	Significant Coefficients are those that have a "Signif. Code" less than 0.05 (or three *** stars) in the Summary output.
	#
	fittedDEGsWithSignificantCoefficients = lm( CuffDiff ~ DESeq + Limma + GeneSpring, data=dataFrameWithJoinedResponse )	#	CuffDiff
	# fittedDEGsWithSignificantCoefficients = lm( NOISeq ~ DESeq + CuffDiff + Limma, data=dataFrameWithJoinedResponse )		# 	NOISeq
	# fittedDEGsWithSignificantCoefficients = lm( DESeq ~ NOISeq + CuffDiff + Limma, data=dataFrameWithJoinedResponse )		#	DESeq
	# fittedDEGsWithSignificantCoefficients = lm( GeneSpring ~ DESeq + NOISeq + CuffDiff + Limma, data=dataFrameWithJoinedResponse )	#	GeneSpring
	
	
	print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] -------------------------------------------------------------", sep="") )
	print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Significant Summary:", sep="") )
	
	print(summary(fittedDEGsWithSignificantCoefficients))
	
	
	print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] ", sep="") )
	print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] ", sep="") )
	print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] ", sep="") )
	print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] -------------------------------------------------------------", sep="") )
	print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Predicted Values:", sep="") )
	
	#
	#	Predicted Fitted
	#
	predictedGeneExpression = predict.lm( fittedDEGsWithSignificantCoefficients, interval="prediction" )
	
	
	OutputFileForFittedGeneExpression = paste( working_directory, "/linear_models/", anAlgorithm, "/round_robin_update/Fitted_Expression-Response_Col-", columnIndexWithResponseAlgorithm, ".txt", sep="" )
	write.table( predictedGeneExpression, file=OutputFileForFittedGeneExpression, append=FALSE, row.names=TRUE, col.names=NA, quote=FALSE, sep="\t" )
	
	
	print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] -------------------------------------------------------------", sep="") )
	print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Diagnostic Plots:", sep="") )
	
	plot(fittedDEGs, which=1:6)
	termplot(fittedDEGs, data=dataFrameWithJoinedResponse)
	
	#	Q-Q Plot diagnostic
	qqplotLM.stdres = rstandard(fittedDEGs)
	qqnorm(qqplotLM.stdres, ylab="Standardized Residuals", xlab="Normal Scores", main="Q-Q Plot")
	qqline(qqplotLM.stdres)
	
}





# ------------------------------------------------------ END ---------------------------------------------------------

print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] ", sep="") )
print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Done.", sep="") )
print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] ", sep="") )
