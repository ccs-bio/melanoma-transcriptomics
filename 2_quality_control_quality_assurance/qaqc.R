
# -----------------------------------------------------------------------------------------------------------
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
#	The purpose of this program is to generate some quality assurance and quality control (QA/QC) 
#	visualizations to use as an idicator of the quality of the data.  The script generates three (3)
#	plots:
#			• Boxplot — Displays a samples experssion value distribution.
#			• Heatmap — Expression values create a hierarchical clustering of Genes & Samples.
#			• Correlation Matrix — Expression values are used to display sample associations.
#			• PCA — Displays a PCA plot for analyzing variation and clustering
#	
#
#	AUTHOR:	Camilo Valdes (cvaldes3@med.miami.edu)
#			Computational Biology and Bioinformatics Group
#			Center for Computational Science (CCS), University of Miami
#
# -----------------------------------------------------------------------------------------------------------

# 	Import any necessary libraries here

library(corrgram)
library(gplots)
library(reshape2)
library(ggplot2)

options( digits = 15 ) # Enable printing of large numbers. R defaults to "rounding" for display.


# ------------------------------------------------- Main ----------------------------------------------------


working_directory="/Users/camilo/Documents/CCS/Enrico/Projects/Melanoma/Analysis_Full"

setwd( file.path( working_directory ) )

# Directory in which the results are stored at
resultsDirectory = "/results-103"

inputDirectory = paste( working_directory, resultsDirectory, sep="" )

# Visualization friendly "sample names"
col_headers = c( "GeneID", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", "58", "59", "60", "61", "62", "63", "64", "65", "66", "67", "68", "69", "70", "71", "72", "73", "74", "75", "76", "77", "78", "79", "80", "81", "82", "83", "84", "85", "86", "87", "88", "89", "90", "91", "92", "93", "94", "95", "96", "97", "98", "99", "100", "101", "C1", "C2" )

# The different quantification units we'll be plotting
quantTypes = c( 'fpkm', 'counts_exonic', 'counts_intronic' )
# quantTypes = c( 'counts_exonic', 'counts_intronic' )
#quantTypes = c( 'fpkm' )

# Read-alignment Types
mappingTypes = c( 'unique' )

# File types, e.g., All or Filtered at some level
fileTypes = c( 'All' )
# fileTypes = c( 'Filtered_10_10', 'Filtered_10_100', 'Filtered_10_500', 'Filtered_10_1000', 'Filtered_66_100' )
# fileTypes = c( 'Filtered_10_0.5' )



#----------------------------------------------- Build & Run ------------------------------------------------


print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] R Starting... ", sep="") )

for( quantType in quantTypes ) {
	
	for( mappingType in mappingTypes ) {
		
		for( fileType in fileTypes ) {
			
			# ------------------------------------------ Input ---------------------------------------------
			
			fileToProcess = paste( inputDirectory, "/", quantType, "/", mappingType, "/", fileType, "/aggregated_expression_profile_for_", quantType, "_AllGenes.txt", sep="" )

			expressionTable = read.delim( fileToProcess, sep="\t", header=TRUE, check.names = FALSE )
			
			# Rename Column Headers for Display
			names(expressionTable) = col_headers
			
			
			# ------------------------------------------ Output --------------------------------------------
			
			outputDir = paste( "Figures_Numbers-Full_103/QA_QC/", quantType, "/", mappingTypes, "/", fileTypes, sep="" )
			
			dir.create( file.path( working_directory, outputDir ), showWarnings=FALSE, recursive=TRUE )

			# Boxplot output file (as PDF)
			outputFileForBoxPlot 	= paste(  outputDir, "/boxplot.pdf", sep="" )
			
			# HeatMap output file (as PNG, pdf is too big)
			outputFileForHeatMap 		= paste(  outputDir, "/heatmap.png", sep="" )
			outputFileForHeatMapLog2 	= paste(  outputDir, "/log2-heatmap.png", sep="" ) 
			
			# Correlation Matrix ouput file (PNG)
			outputFileForCorrMatrix 	= paste(  outputDir, "/correlationMatrix.png", sep="" )
			outputFileForCorrMatrixLog2	= paste(  outputDir, "/log2-correlationMatrix.png", sep="" )
			
			# PCA plot output files (PNG)
			outputFileForPCA1	= paste(  outputDir, "/pca_components.png", sep="" )
			outputFileForPCA2	= paste(  outputDir, "/pca.png", sep="" )
			
			
			# ------------------------------------------ BoxPlot -------------------------------------------
			
			print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Boxplot... ", sep="") )
			
			# Melt the data for box_plot format
			expressionTableMelted = melt( expressionTable, id.vars = "GeneID", measure.vars = c(2:104))
			
			#
			#	Boxplot
			#
			# boxplot =   ggplot( expressionTableMelted, aes(x=variable, y=log2(value), fill=variable ) ) +
# 							xlab("Samples") +
# 							ylab( paste( quantType," (log2)", sep="" ) ) +
# 							geom_boxplot( outlier.colour = "red" ) +
# 							ggtitle( paste( quantType, " - ", fileType, sep="" ) ) +
# 							guides(fill=FALSE) +
# 							theme( text = element_text(size=5), axis.text.x = element_text(angle=270, vjust=1) )
#
# 			ggsave( outputFileForBoxPlot )
			
			#
			#	FPKM requires an extra BoxPlot that is Zoomed-In on the Y-Axis
			#
			if( quantType == 'fpkm' ) {
				
				print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Boxplot (FPKM zoomed)... ", sep="") )
				
				# Y-axis display limits
				yAxis_lower = -20
				yAxis_upper = 15
				
				# boxplotFPKM =   ggplot( expressionTableMelted, aes(x=variable, y=log2(value), fill=variable ) ) +
# 								xlab("Samples") +
# 								ylab( paste( quantType," (log2)", sep="" ) ) +
# 								geom_boxplot( outlier.colour = "red" ) +
# 								ggtitle( paste( quantType, " - ", fileType, sep="" ) ) +
# 								guides(fill=FALSE) + ylim( yAxis_lower, yAxis_upper ) +
# 								theme( text = element_text(size=5), axis.text.x = element_text(angle=270, vjust=1) )
#
# 				ggsave( paste(  outputDir, "/boxplot_zoomed.pdf", sep="" ) )
			}
			
			
			
			# ------------------------------------------ HeatMap -------------------------------------------
			
			print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] HeatMap... ", sep="") )
			
			# The heatmap() function does not take Tables, so we need to convert the input table to a Data Frame
			dataFrameForHeatMap = as.data.frame.matrix( expressionTable[,2:104] )
			
			# RowNames have to be renamed
			rownames(dataFrameForHeatMap) = expressionTable[,1]
			
			
			# ----------------------------- Heatmap -------------------------------
			
			# Open the graphics device to save to
			# png( file=outputFileForHeatMap, width=2056, height=2056, units="px", res=120 )
#
# 			heatmap.2( 	as.matrix( dataFrameForHeatMap[complete.cases(dataFrameForHeatMap),] ),
# 						col=redgreen(75),
# 						scale="row",
# 						key=TRUE,
# 						keysize=.5,
# 			          	density.info="none",
# 						trace="none",
# 						labRow=NA,
#						main=paste( quantType, " - ", fileType, sep="" ) )
#
#
# 			dev.off()	# Close the graphics device
			
			
			# --------------------------- Heatmap Log2 ----------------------------
			
			print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] HeatMap Log2()... ", sep="") )
			
			# Zeros give us problems when we Log2, so we replace them with a minimum
			dataFrameForHeatMapLog2 = dataFrameForHeatMap
			
			# Replace Zero values with the Minimum value in each Sample (column)
			# Why not do it across Samples — The minimum value of a gene (row)???
			dataFrameForHeatMapLog2 = apply(dataFrameForHeatMapLog2,2,function(x){x[x==0] = min(x[x!=0]);x})
			
			dataFrameForHeatMapLog2 = log2( dataFrameForHeatMapLog2 )
			
			# Open the graphics device to save to
			# png( file=outputFileForCorrMatrixLog2, width=2056, height=2056, units="px", res=120 )
#
# 			heatmap.2( 	as.matrix( dataFrameForHeatMapLog2 ),
# 						col=redgreen(75),
# 						scale="row",
# 						key=TRUE,
# 						keysize=1.5,
# 			          	density.info="none",
# 						trace="none",
# 						cexCol=0.9,
# 						labRow=NA,
#						main=paste( quantType, " - ", fileType, sep="" ) )
#
# 			dev.off()	# Close the graphics device
			
			
			
			# ----------------------------------- Correlation Matrix ---------------------------------------
			
			print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Correlation Plot... ", sep="") )
			
			# Need to reformat the input data structs to calculate the correlation matrix
			expressionDataFrameForCorrMatrix = as.data.frame.matrix( expressionTable[,2:104] )
			
			rownames(expressionDataFrameForCorrMatrix) = expressionTable[,1]
			
			correlationMatrix = cor( expressionDataFrameForCorrMatrix, use="complete.obs", method="pearson" )
			
			# vector of colors (purple to orange)
			purple_orange = c("#8073AC","#B2ABD2","white","#E08214","#FDB863")
			
			
			# ----------------------------- Correlation Plot -------------------------------
			
			# Open the graphics device to save to
			# png( file=outputFileForCorrMatrix, width=2056, height=2056, units="px", res=120 )
#
# 			corrgram( 	correlationMatrix,
# 						order=TRUE,
# 						lower.panel=panel.shade,
# 						upper.panel=panel.shade,
# 						text.panel=panel.txt,
# 						main=paste( quantType, " - ", fileType, sep="" ),
# 					  	col.regions = colorRampPalette( purple_orange ) )
#
# 			dev.off()	# Close the graphics device
			
			
			# --------------------------- Correlation Plot Log2 ----------------------------
		
			print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Correlation Plot Log2()... ", sep="") )
			
			# Open the graphics device to save to
			# png( file=outputFileForCorrMatrixLog2, width=2056, height=2056, units="px", res=120 )
#
# 			corrgram( 	log2(correlationMatrix),
# 						order=TRUE,
# 						lower.panel=panel.shade,
# 						upper.panel=panel.shade,
# 						text.panel=panel.txt,
# 						main=paste( quantType, " - ", fileType, sep="" ),
# 					  	col.regions = colorRampPalette( purple_orange ) )
#
# 			dev.off()	# Close the graphics device
			
			
			
			# ----------------------------------- PCA ----------------------------------------
			#
			# 	WARNING: 'prcomp' cannot deal with "NA", so it COULD be removed using 'na.omit()' 
			#	but this function will remove an entire row if any of its cells contain an 'NA'.  
			#	'complete.cases' is another option as it allows partial selection by using part
			#	of the dataframe.
			#
			
			dataFrameForPCA = as.data.frame.matrix( expressionTable[,2:104] )
			
			# RowNames have to be renamed
			rownames(dataFrameForHeatMap) = expressionTable[,1]

			# Replace Zero values with the Minimum value in each Sample (column)
			dataFrameForPCAwithMin = apply(dataFrameForPCA,2,function(x){x[x==0] = min(x[x!=0]);x})
			
			# Log2() Transform for comparinsons
			dataFrameForPCAwithMinLog2 = log2( dataFrameForPCAwithMin )
			
			# Need to transpose the input matrix so that we cluster on Samples (columns)
			transposedDataFrameForPCA 	  = t( dataFrameForPCAwithMin )
			
			transposedDataFrameForPCAlog2 = t( dataFrameForPCAwithMinLog2 )
			
			
			# --------------------------- PCA Plot 1 ----------------------------
			# 	Plot of the variances (y-axis) associated with the PCs (x-axis).
			# 	The Figure below is useful to decide how many PCs to retain for further analysis.
			
			print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] PCA Plot 1 ... ", sep="") )
			
			pca1 = prcomp( na.omit(transposedDataFrameForPCA), scale=FALSE )
			
			# Open the graphics device to save to
			png( file=outputFileForPCA1, width=2056, height=2056, units="px", res=120 )


			plot( 	pca1, 
					main=paste( quantType, " - ", fileType, " Refinement - PCs", sep="" ), 
					type="l" )
					

			dev.off() # Close the graphics device
			
			
			# --------------------------- PCA Plot 2 ----------------------------
			
			print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] PCA Plot 2 ... ", sep="") )
			
			pca2 = prcomp( na.omit(transposedDataFrameForPCAlog2), scale=FALSE )
			
			# Open the graphics device to save to
			png( file=outputFileForPCA2, width=2056, height=2056, units="px", res=120 )
			
			
			plot( 	pca2$x[,1:2], 
					main=paste(  quantType, " - ", fileType, " Refinement - PCs", sep="" ), 
					pch=19, 
					col="blue", 
					cex=.9 )
					
			text(	pca2$x[,1], 
					pca2$x[,2], 
					labels=c( "T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10", "T11", "T12", "T13", "T14", "T15", "T16", "T17", "T18", "T19", "T20", "T21", "T22", "T23", "T24", "T25", "T26", "T27", "T28", "T29", "T30", "T31", "T32", "T33", "T34", "T35", "T36", "T37", "T38", "T39", "T40", "T41", "T42", "T43", "T44", "T45", "T46", "T47", "T48", "T49", "T50", "T51", "T52", "T53", "T54", "T55", "T56", "T57", "T58", "T59", "T60", "T61", "T62", "T63", "T64", "T65", "T66", "T67", "T68", "T69", "T70", "T71", "T72", "T73", "T74", "T75", "T76", "T77", "T78", "T79", "T80", "T81", "T82", "T83", "T84", "T85", "T86", "T87", "T88", "T89", "T90", "T91", "T92", "T93", "T94", "T95", "T96", "T97", "T98", "T99", "T100", "T101", "C1", "C2" ), 
					pos=3 )
			
			
			dev.off() # Close the graphics device
			
			
			
			
			
			

		} # End loop for file type (All, etc.)
		
	} # End loop for mapping type (fpkm, counts_exonic, etc.)
	
} # End loop for quality type (unique, etc.)

