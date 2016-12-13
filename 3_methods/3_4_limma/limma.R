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
#	The purpose of this program is to compute the differentially expressed genes using the limma package.  It takes as
#	input FPKM values calculated with the Cufflinks package.
#
#
#
#	AUTHOR:	Camilo Valdes (cvaldes3@med.miami.edu)
#			Computational Biology and Bioinformatics Group
#			Center for Computational Science (CCS), University of Miami
#
# ---------------------------------------------------------------------------------------------------------------------

# 	Import any necessary libraries here

library(Biobase)
library(edgeR)
library(reshape2)
library(ggplot2)
library(limma)

options( digits = 15 ) # Enable printing of large numbers. R defaults to "rounding" for display.


# ----------------------------------------------------- Setup ---------------------------------------------------------

print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] R Starting... ", sep="") )

working_directory="/path/to/base/analysis/directory"

setwd( file.path( working_directory ) )


# ------------------------------------------------------ Main ---------------------------------------------------------

# The Design Matrix to use for the Interpretation of C v T
design_pool = model.matrix(~0+c( rep('test',101), rep('control',2) ) )

# Re-head the table to display what is what
colnames(design_pool) = c('test','control')

# The file 'list' we wish to process
files = c( 'AllGenes' )

# Read-alignment Types for which counts have been extracted
types = c( 'unique' )

# Filtered or Raw data
data_types = c( 'All', 'Filtered_10_1' )


for( aFile in files ) {

	for( aType in types ) {

		for( dataType in data_types )
		{
			
			baseDirFigures = '/path/to/base/analysis/directory/differential_expression/fpkm_cuffdiff/limma'
			
			
			# ------------------------------------------------- Read Data In ---------------------------------------------------
			
			print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] ", dataType, sep="") )
			
			
			fileToProcess = paste( "results-103/fpkm_cuffdiff/", aType, "/", dataType, "/aggregated_expression_profile_for_fpkm_cuffdiff_", aFile, ".txt", sep="" )

			fpkmTable = read.delim( fileToProcess, header=TRUE, row.names=1)
			
			
			
			# --------------------------------------------- BoxPlot 1 ------------------------------------------------
			
			outputFileForBoxPlot_1 	= paste(  baseDirFigures, "/boxplot-", dataType, ".pdf", sep="" )
			
			expressionTableMelted_1 = melt( as.matrix(fpkmTable), id.vars = "GeneID", measure.vars = c(2:104))
			
			print( head( expressionTableMelted_1, n=5 ) )
			
			# Y-axis display limits
			yAxis_lower = -20
			yAxis_upper = 15
			
			boxplot_1 =   ggplot( expressionTableMelted_1, aes(x=Var2, y=log2(value), fill=Var2 ) ) +
							xlab("Samples") +
							ylab( paste( "CuffDiff FPKM"," (log2)", sep="" ) ) +
							geom_boxplot( outlier.colour = "red" ) +
							ggtitle( paste( "CuffDiff FPKM", dataType, " - limma", sep="" ) ) +
							guides(fill=FALSE) + ylim( yAxis_lower, yAxis_upper ) +
							theme( text = element_text(size=5), axis.text.x = element_text(angle=270, vjust=1) )

			ggsave( outputFileForBoxPlot_1 )
			
			
			
			# ------------------------------------------------ Classic Analysis ------------------------------------------------
			
			print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Classic", sep="") )
			
			
			# Log-Transform the FPKMs (keeping parity with GeneSpring)
			fpkmTable_log = log2(fpkmTable)
			
			print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Log Table", sep="") )
			
			print( head( fpkmTable_log, n=10 ) )
			

			# Genes could have an FPKM of zero (not detected), so we clean the log(0) result of '-INF' to
			# a format that limma can understand, namely 'NA'
			is.na(fpkmTable_log) <- sapply(fpkmTable_log, is.infinite)

			# Creates appropriate contrast matrix to perform the comparisons
			cm = makeContrasts(test-control, levels=design_pool)

			# Fits a linear model for each gene based on the given series of samples
			fit = lmFit( fpkmTable_log, design_pool )

			# Compute a contrast matrix with all desired comparisons
			# and estimated coefficients and standard errors for a given set of contrasts.
			fit2 = contrasts.fit( fit,cm )

			# Computes moderated t-statistics and log-odds of differential expression by
			# empirical Bayes shrinkage of the standard errors towards a common value.
			fit3 = eBayes( fit2 )

			# Generates a list of top differentially expressed genes sorted by B-values ('sort.by=B')
			# for each of the comparison groups ('coef=1')
			# The summary table contains the following information: logFC is the log2-fold change, the
			# the AveExpr is the average expression value accross all samples, the moderated t-statistic (t)
			# is the logFC to its standard error, the P.Value is the associated p-value, the adj.P.Value
			# is the p-value adjusted for multiple testing and the B-value (B) is the log-odds that a
			# gene is differentially expressed (the-higher-the-better)
			#
			# Note: Usually one wants to base gene selection on the adjusted P-value rather than the
			# t- or B-values.
			#
			topTable( fit3, coef=1, adjust="fdr", sort.by="B", number=10 )
			
			
			
			# ---------------------------------------------------- Output ------------------------------------------------------
			
			outputFile_classic 	= paste( "differential_expression_103/_fpkm-cuffdiff/limma/results_pooled/", aType, "/", dataType, "/", aFile, "_SigDiffExp-", aType, ".txt", sep="" )
			
			topTableForClassic	= topTable( fit3, coef=1, adjust="fdr", sort.by="P", number=50000)
			
			write.table( topTableForClassic, file=outputFile_classic, row.names=T, col.names=NA, sep="\t")

			
			
			# --------------------------------------------- Hierarhical Clustering ---------------------------------------------
			#
			# 	Classic Analysis
			
			print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Heatmap...", sep="") )
			
			# Select the genes below a given threshold (pvalue < 0.05 && fc > 2)
			selectedGenesForClassic = topTableForClassic[topTableForClassic$adj.P.Val < 0.05,]
			selectedGenesForClassic = selectedGenesForClassic[selectedGenesForClassic$logFC > 2,]
			selectedGenesForClassic_2 = selectedGenesForClassic[, c("logFC", "adj.P.Val")]
			
			print( head( selectedGenesForClassic_2, n=10) )
			
			figureFileNameForClassic = paste( dataType, ".png", sep="" )
			
			png( file=paste( baseDirFigures, "/", figureFileNameForClassic, sep="" ), width=1024, height=1024, units="px", res=120 )
			
			dataFrameForClassic = as.data.frame.matrix( fpkmTable_log[row.names(selectedGenesForClassic),] )
			
			print( head( dataFrameForClassic, n=10) )
			
			heatmap( as.matrix( dataFrameForClassic[complete.cases(dataFrameForClassic),] ) )
			
			dev.off()
			
			
			

		} # end data type (filtered vs raw)

	} # end type to process (multi & unique)

} # end file names to process
