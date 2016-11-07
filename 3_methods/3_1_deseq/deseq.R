
#-----------------------------------------------------------------------------------------------------------
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
#	The purpose of this program is to use the Bioconductor package "DESeq" to analyze a count-based RNA-Seq
#	dataset.
#
#	AUTHOR:	Camilo Valdes (cvaldes3@med.miami.edu)
#			Computational Biology and Bioinformatics Group
#			Center for Computational Science (CCS), University of Miami
#
#-----------------------------------------------------------------------------------------------------------

# 	Import any necessary libraries here

library(DESeq)

options( digits = 15 ) # Enable printing of large numbers. R defaults to "rounding" for display.

#-----------------------------------------------------------------------------------------------------------
#
#                                         			Main


# Working Directory
working_directory="/Users/camilo/Documents/CCS/Enrico/Projects/Melanoma/Analysis_Full"

setwd( file.path( working_directory ) )

# Experimental Grouping and Interpretation for Pooled Testing
conds_pool = c("T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","C","C")

# The file 'list' we wish to process
files = c( 'AllGenes' )

# Read-alignment Types for which counts have been extracted
types = c( 'unique' )

# File types
# fileTypes = c( 'Filtered_10_10', 'Filtered_66_100', 'Filtered_10_1000', 'Filtered_10_500' )
# fileTypes = c( 'Filtered_10_100' )
fileTypes = c( 'Filtered_5_10', 'Filtered_20_10' )

#---------------------------------------------------------------------------------------------
#
#	Pool Interpretation — Treat all TCGA samples as replicates
#

print( "Testing Pooled Interpretation..." )

for( aFile in files ) {

	for( aType in types ) {

		for( fileType in fileTypes ) {

			fileToProcess = paste( "results-103/counts_exonic/", aType, "/", fileType, "/aggregated_expression_profile_for_counts_exonic_", aFile, ".txt", sep="" )

			print(fileToProcess)

			countsTable = read.delim( fileToProcess, header=TRUE, row.names=1, sep="\t", check.names=FALSE )

			cds	= newCountDataSet( countsTable, conds_pool )
			cds	= estimateSizeFactors(cds)
			cds	= estimateDispersions(cds, method="per-condition", sharingMode="fit-only", fitType="parametric" )  # If the parametric test fails, use the fitType="local" parameter

			results = nbinomTest( cds, "T", "C" )
			pvalues = results$pval
			adjustedpvalues = p.adjust( pvalues, method="fdr" )
			
			
			
			# ---------------------------------------------------- Output ------------------------------------------------------
			
			outputDir = paste( "differential_expression_103/_counts_exonic/deseq/counts/results_pooled/", aType, "/", fileType, sep="" )
			
			dir.create( file.path( working_directory, outputDir ), showWarnings=FALSE, recursive=TRUE )

			outputFile = paste(  outputDir, "/", aFile, "_SigDiffExp-", aType, ".txt", sep="" )

			write.table( results, file=outputFile, append=FALSE, sep="\t" )
			
			
			
			# --------------------------------------------- Hierarhical Clustering ---------------------------------------------
			
			baseDirFigures = '/Users/camilo/Documents/CCS/Enrico/Projects/Melanoma/Analysis_Full/Figures_Numbers-Full_103/differential_expression/counts_exonic/deseq'

			selectedGenesForDeseq = results[results$padj < 0.05,]
			selectedGenesForDeseq = selectedGenesForDeseq[selectedGenesForDeseq$foldChange > 2,]
			
			# Deseq clobbers the row names of the dataframe, so we have to rename them
			rownames(selectedGenesForDeseq) = selectedGenesForDeseq[,1]			

			figureFileNameForDeseq = paste( fileType, ".png", sep="" )
			png( file=paste( baseDirFigures, "/", figureFileNameForDeseq, sep="" ), width=1024, height=1024, units="px", res=120 )
			
			dataFrameForDeseq = as.data.frame.matrix( countsTable[row.names(selectedGenesForDeseq),] )
			
			heatmap( as.matrix( dataFrameForDeseq[complete.cases(dataFrameForDeseq),] ) )
			
			dev.off()


		}
	}
}

