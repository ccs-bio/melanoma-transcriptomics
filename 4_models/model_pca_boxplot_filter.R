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
#	The purpose of this program is to create scatter plots for comparing Raw and Permuted averages
#
#
#
#	AUTHOR:	Camilo Valdes (cvaldes3@miami.edu)
#			Computational Biology and Bioinformatics Group
#			Center for Computational Science (CCS), University of Miami
#
# ---------------------------------------------------------------------------------------------------------------------

# 	Import any necessary libraries here

library(reshape2)
library(RColorBrewer)
library(ggplot2)
require(grid)
library(plyr)

# 	Enable printing of large numbers. R defaults to "rounding" for display, and
# 	Cause the warnings to be printed as they occur.
options( digits=15, warn=1 )


# ------------------------------------------------- Project Setup -----------------------------------------------------

working_directory="/path/to/pca_coefficients"

ggplot_output_dir="/path/to/figures/directory/for/PCA_Reduction/coefficients"

setwd( file.path( working_directory ) )


# ------------------------------------------------ Custom Functions ---------------------------------------------------
#
#

#
#	Append data to an Empty Data.frame
#	Missing values are replaced with NA
#
cbind.all <- function (...) 
{
    nm <- list(...)
    nm <- lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow))
    do.call(cbind, lapply(nm, function(x) rbind(x, matrix(, n - 
        nrow(x), ncol(x)))))
}


# ------------------------------------------------------ Main ---------------------------------------------------------

print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] ", sep="") )
print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] R Starting... ", sep="") )
print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] ", sep="") )


algorithms = c( 'counts_exonic_deseq', 'counts_exonic_noiseq', 'counts_intronic_deseq', 'counts_intronic_noiseq', 'cuffdiff', 'fpkm_limma', 'fpkm_genespring')

globalDataFrame = data.frame()

#
#	Plots for individual Algorithms
#
for( anAlgorithm in algorithms ) {
	
	
	# --------------------------------------------------- Data Prep ------------------------------------------------------
	#
	#
	
	rawFile 		= paste( working_directory, "/raw_averages/coefficients/pca_coefficients-", anAlgorithm, ".txt", sep="" )
	rawTable		 = read.delim( rawFile, header=TRUE, sep="\t", row.names=1, check.names=FALSE )
	rawTableSorted 		= rawTable[order(row.names(rawTable)),]
	
	
	#	Structs for Plotting
	dataFrameWithPC1 = data.frame( seq(from=1,to=length(rawTableSorted[,"PC1"]),by=1), rawTableSorted[,"PC1"] )
	colnames(dataFrameWithPC1) = c("x", "y")
	rownames(dataFrameWithPC1) = rownames(rawTableSorted)
	maximumValue = max( rawTableSorted[,"PC1"] )
	
	#	Standard Deviation
	standardDeviationForPC1 = sapply(rawTableSorted, sd)
	
	print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Standard Deviation: ", standardDeviationForPC1, sep="") )
		
	
	# Add the interesting column to the global dataframe
	globalDataFrame = cbind.all( globalDataFrame, rawTableSorted[,"PC1"] )
	
		
	# ----------------------------------------------------- Plots --------------------------------------------------------
	#
	#
	
	#	Color Ramps
	colorRampForPositives = colorRampPalette(brewer.pal( n=8, "Blues" ))(1000)
	colorRampForNegatives = colorRampPalette(brewer.pal( n=8, "Oranges" ))(1000)
	
	#
	#	Bar Plot
	#
	ggplot( dataFrameWithPC1, aes( x=x, y=y, fill=y ) ) + 
			scale_fill_gradient2( low=rev(colorRampForNegatives), high=colorRampForPositives, limits=c(-maximumValue, maximumValue) ) +
			geom_bar(stat="identity", position="identity") +
			theme( legend.title=element_blank() ) +
			ylim( -maximumValue, maximumValue ) +
			ggtitle( paste( "", anAlgorithm, sep="" ) ) +
			xlab( "Coefficients") + ylab("Coefficient Values")

	ggsave( paste( ggplot_output_dir, "/bar_plots/bar_plot-raw-", anAlgorithm, ".png" ,sep="" ) )
		
	
	#
	#	Histogram
	#
	ggplot( dataFrameWithPC1, aes( x=y ) ) +
			geom_histogram( aes(y = ..density..) ) +
			geom_density( colour="red" ) +
			ggtitle( paste( "", anAlgorithm, sep="" ) ) +
			xlab( "Coefficient") + ylab("Coefficient Frequency")
	
	ggsave( paste( ggplot_output_dir, "/histograms/histogram-raw-", anAlgorithm, ".png" ,sep="" ) )
	
	
	
	
	# ----------------------------------------------- BoxPlot Extraction -------------------------------------------------
	#
	#	Extraction based on the BoxPlot of Genes with interesting (i.e., significant) coefficients.
	#
	print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] ", anAlgorithm, " Data Frame...", sep="") )
	print(head(dataFrameWithPC1,n=10))
	
	#
	#	Quick, single view, Box Plot
	#
	boxplot( dataFrameWithPC1$y, main=paste("", anAlgorithm, sep="") )
	
	#	Get the boxplot stats
	boxplotStats = boxplot( dataFrameWithPC1$y, plot=FALSE )$stats
	
	statsUpperHinge = boxplotStats[4,]
	statsLowerHinge = boxplotStats[2,]
	
	print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Stats...", sep="") )
	print(boxplotStats)
	print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Upper Hinge: ", statsUpperHinge, sep="") )
	print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Lower Hinge: ", statsLowerHinge, sep="") )
	
	print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"]", sep="") )
	
	genesGreaterThanUpperHinge = dataFrameWithPC1[dataFrameWithPC1$y > statsUpperHinge,]
	genesLowerThanLowerHinge   = dataFrameWithPC1[dataFrameWithPC1$y < statsLowerHinge,]
	
	print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Upper Hinge Genes...", sep="") )
	print(nrow(genesGreaterThanUpperHinge))
	print(head(genesGreaterThanUpperHinge,n=10))
	
	print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Lower Hinge Genes...", sep="") )
	print(nrow(genesLowerThanLowerHinge))
	print(head(genesLowerThanLowerHinge,n=10))
	
	
	outputFileforUpperHinge = paste( working_directory, "/raw_averages/boxplot/hinges/upper_hinge-", anAlgorithm, ".txt", sep="" )
	write.table( genesGreaterThanUpperHinge, file=outputFileforUpperHinge, append=FALSE, row.names=TRUE, col.names=NA, quote=FALSE, sep="\t", na="0" )
	
	outputFileforLowerHinge = paste( working_directory, "/raw_averages/boxplot/hinges/lower_hinge-", anAlgorithm, ".txt", sep="" )
	write.table( genesLowerThanLowerHinge, file=outputFileforLowerHinge, append=FALSE, row.names=TRUE, col.names=NA, quote=FALSE, sep="\t", na="0" )
	
	
	# ---------------------------------------------------- Outliers ------------------------------------------------------
	#
	#	Grab the outliers (their values)
	outliersFromBoxPlot = boxplot( dataFrameWithPC1$y, plot=FALSE )$out
	
	#	Extract the outliers (Gene labels) from the original data frame
	outlierGeneLabels = dataFrameWithPC1[dataFrameWithPC1$y %in% outliersFromBoxPlot,]
	
	
	print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Outlier Values...", sep="") )
	print(nrow(outlierGeneLabels))
	print(head(outliersFromBoxPlot,n=10))
	
	print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Outlier Labels...", sep="") )
	print(nrow(outlierGeneLabels))
	print(head(outlierGeneLabels,n=10))
	
	
	outputFileforOutliers = paste( working_directory, "/raw_averages/boxplot/outliers/outliers-", anAlgorithm, ".txt", sep="" )
	write.table( outlierGeneLabels, file=outputFileforOutliers, append=FALSE, row.names=TRUE, col.names=NA, quote=FALSE, sep="\t", na="0" )
}

colnames(globalDataFrame) = algorithms
print(head(globalDataFrame,n=10))



# ------------------------------------------------- Global BoxPlot ----------------------------------------------------
#
#

meltedDataFrameForGlobalDataFrame = melt( globalDataFrame )

print(head(meltedDataFrameForGlobalDataFrame,n=10))

ggplot( meltedDataFrameForGlobalDataFrame, aes( x=Var2, y=value, fill=Var2 ) ) +
		geom_boxplot() + 
		theme( legend.position="none", axis.text.x=element_text(angle=90, hjust=1) ) +
		ggtitle( "Coefficients" ) +
		xlab( "Algorithms")

ggsave( paste( ggplot_output_dir, "/box_plots/boxplot-raw-", anAlgorithm, ".png" ,sep="" ) )




# ------------------------------------------------- Global SD Plot ----------------------------------------------------
#
#



ggplot( meltedDataFrameForGlobalDataFrame, aes( x=Var2, y=value ) ) +
	geom_point() +
	stat_summary( fun.y=mean, colour="red", geom="point", shape=18, size=3, show_guide=FALSE ) +
	stat_summary( fun.data=mean_cl_normal, colour="red", geom="errorbar", mult=1 ) +
	theme( legend.position="none", axis.text.x=element_text( angle=90, hjust=1 ) ) +
	ggtitle( "Coefficients" ) +
	xlab( "Algorithms")

ggsave( paste( ggplot_output_dir, "/standard_deviations/sd-raw-", anAlgorithm, ".png" ,sep="" ) )









# ------------------------------------------------------ END ---------------------------------------------------------

print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] ", sep="") )
print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Done.", sep="") )
print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] ", sep="") )
