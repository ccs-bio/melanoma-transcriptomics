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
#	The purpose of this program is to create PCA plots for the DEG results from the algorithms
#
#
#
#	AUTHOR:	Camilo Valdes (cvaldes3@miami.edu)
#			Computational Biology and Bioinformatics Group
#			Center for Computational Science (CCS), University of Miami
#
# ---------------------------------------------------------------------------------------------------------------------

# 	Import any necessary libraries here

library(KernSmooth)
library(reshape2)
library(ggplot2)
require(grid)

# 	Enable printing of large numbers. R defaults to "rounding" for display, and
# 	Cause the warnings to be printed as they occur.
options( digits=15, warn=1 )


# ------------------------------------------------- Project Setup -----------------------------------------------------

working_directory="/path/to/base/directory/for/PCA_Reduction_103"
setwd( file.path( working_directory ) )

ggplot_output_dir="/path/to/figures/for/PCA_Reduction/raw_averages"


# ------------------------------------------------------ Main ---------------------------------------------------------


print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] ", sep="") )
print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] R Starting... ", sep="") )
print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] ", sep="") )




#
#	Path to the files with the DEG results
#

print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Loading DEG files...", sep="") )

degFileCountsExonicForDeseq = paste( working_directory, "/pca_data_lists/", "counts_exonic/deseq", "/deg_with_expression.txt", sep="" )
degFileCountsExonicForNoiseq = paste( working_directory, "/pca_data_lists/", "counts_exonic/noiseq", "/deg_with_expression.txt", sep="" )

degFileCountsIntronicForDeseq = paste( working_directory, "/pca_data_lists/", "counts_intronic/deseq", "/deg_with_expression.txt", sep="" )
degFileCountsIntronicForNoiseq = paste( working_directory, "/pca_data_lists/", "counts_intronic/noiseq", "/deg_with_expression.txt", sep="" )

degFilePoissonQuartileForCuffdiff = paste( working_directory, "/pca_data_lists/", "cuffdiff_poisson_quartile/cuffdiff", "/deg_with_expression.txt", sep="" )

degFileFpkmForLimma = paste( working_directory, "/pca_data_lists/", "fpkm/limma", "/deg_with_expression.txt", sep="" )
degFileFpkmForGenespring = paste( working_directory, "/pca_data_lists/", "fpkm/genespring", "/deg_with_expression.txt", sep="" )



#
#	Tables with DEG results
#

print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Creating DEG tables...", sep="") )

degTableCountsExonicForDeseq 	= read.delim( degFileCountsExonicForDeseq, header=TRUE, row.names=1, sep="\t", check.names=FALSE )
degTableCountsExonicForNoiseq 	= read.delim( degFileCountsExonicForNoiseq, header=TRUE, row.names=1, sep="\t", check.names=FALSE )

degTableCountsIntronicForDeseq 	= read.delim( degFileCountsIntronicForDeseq, header=TRUE, row.names=1, sep="\t", check.names=FALSE )
degTableCountsIntronicForNoiseq = read.delim( degFileCountsIntronicForNoiseq, header=TRUE, row.names=1, sep="\t", check.names=FALSE )

degTablePoissonQuartileForCuffdiff = read.delim( degFilePoissonQuartileForCuffdiff, header=TRUE, row.names=1, sep="\t", check.names=FALSE )

degTableFpkmForLimma		= read.delim( degFileFpkmForLimma, header=TRUE, row.names=1, sep="\t", check.names=FALSE )
degTableFpkmForGenespring	= read.delim( degFileFpkmForGenespring, header=TRUE, row.names=1, sep="\t", check.names=FALSE )



# --------------------------------------------------- Log2 Transformation -----------------------------------------------------
#
#	Use the log2() function to transform the data into a Normal distribution
#
#	Careful with zeros ("0"), change them to 1 and then Log!
#

degTableCountsExonicForDeseq[ degTableCountsExonicForDeseq == 0 ] = 1
degTableCountsExonicForNoiseq[ degTableCountsExonicForNoiseq == 0 ] = 1

degTableCountsIntronicForDeseq[ degTableCountsIntronicForDeseq == 0 ] = 1
degTableCountsIntronicForNoiseq[ degTableCountsIntronicForNoiseq == 0 ] = 1

degTablePoissonQuartileForCuffdiff[ degTablePoissonQuartileForCuffdiff == 0 ] = 1

degTableFpkmForLimma[ degTableFpkmForLimma == 0 ] = 1
degTableFpkmForGenespring[ degTableFpkmForGenespring == 0 ] = 1



degTableCountsExonicForDeseq		= log2(degTableCountsExonicForDeseq)
degTableCountsExonicForNoiseq		= log2(degTableCountsExonicForNoiseq)

degTableCountsIntronicForDeseq		= log2(degTableCountsIntronicForDeseq)
degTableCountsIntronicForNoiseq		= log2(degTableCountsIntronicForNoiseq)

degTablePoissonQuartileForCuffdiff	= log2(degTablePoissonQuartileForCuffdiff)

degTableFpkmForLimma				= log2(degTableFpkmForLimma)
degTableFpkmForGenespring			= log2(degTableFpkmForGenespring)


# ------------------------------------------------------ Raw Averages ---------------------------------------------------------


# 	Counts Exonic
geneAverageInTumorsForCountsExonicInDeseq 	= apply( degTableCountsExonicForDeseq[,1:101], 1, mean )
geneAverageInControlsForCountsExonicInDeseq = apply( degTableCountsExonicForDeseq[,102:103], 1, mean )

geneAverageInTumorsForCountsExonicInNoiseq 	 = apply( degTableCountsExonicForNoiseq[,1:101], 1, mean )
geneAverageInControlsForCountsExonicInNoiseq = apply( degTableCountsExonicForNoiseq[,102:103], 1, mean )


# 	Counts Intronic
geneAverageInTumorsForCountsIntronicInDeseq 	= apply( degTableCountsIntronicForDeseq[,1:101], 1, mean )
geneAverageInControlsForCountsIntronicInDeseq 	= apply( degTableCountsIntronicForDeseq[,102:103], 1, mean )

geneAverageInTumorsForCountsIntronicInNoiseq 	= apply( degTableCountsIntronicForNoiseq[,1:101], 1, mean )
geneAverageInControlsForCountsIntronicInNoiseq	= apply( degTableCountsIntronicForNoiseq[,102:103], 1, mean )


#	Cuffdiff
geneAverageInTumorsForPoissonQuartileInCuffdiff 	= apply( degTablePoissonQuartileForCuffdiff[,1:101], 1, mean )
geneAverageInControlsForPoissonQuartileInCuffdiff 	= apply( degTablePoissonQuartileForCuffdiff[,102:103], 1, mean )


#	FPKM
geneAverageInTumorsForFPKMsInLimma 		= apply( degTableFpkmForLimma[,1:101], 1, mean )
geneAverageInControlsForFPKMsInLimma 	= apply( degTableFpkmForLimma[,102:103], 1, mean )

geneAverageInTumorsForFPKMsInGenespring 	= apply( degTableFpkmForGenespring[,1:101], 1, mean )
geneAverageInControlsForFPKMsInGenespring	= apply( degTableFpkmForGenespring[,102:103], 1, mean )


# ------------------------------------------------------ Merge Averages ---------------------------------------------------------
#
#	Merge individual Tumor and Control averages into a data.frame type we can use with PCA
#


#	Counts Exonic
dataFrameWithAverageTumorsAndControlsForCountsExonicInDeseq = merge( geneAverageInTumorsForCountsExonicInDeseq, geneAverageInControlsForCountsExonicInDeseq, by="row.names", all.x=TRUE, all.y=TRUE )
rownames(dataFrameWithAverageTumorsAndControlsForCountsExonicInDeseq) = dataFrameWithAverageTumorsAndControlsForCountsExonicInDeseq[,1]
dataFrameWithAverageTumorsAndControlsForCountsExonicInDeseq[,1] = NULL

dataFrameWithAverageTumorsAndControlsForCountsExonicInNoiseq = merge( geneAverageInTumorsForCountsExonicInNoiseq, geneAverageInControlsForCountsExonicInNoiseq, by="row.names", all.x=TRUE, all.y=TRUE )
rownames(dataFrameWithAverageTumorsAndControlsForCountsExonicInNoiseq) = dataFrameWithAverageTumorsAndControlsForCountsExonicInNoiseq[,1]
dataFrameWithAverageTumorsAndControlsForCountsExonicInNoiseq[,1] = NULL


#	Counts Intronic
dataFrameWithAverageTumorsAndControlsForCountsIntronicInDeseq = merge( geneAverageInTumorsForCountsIntronicInDeseq, geneAverageInControlsForCountsIntronicInDeseq, by="row.names", all.x=TRUE, all.y=TRUE )
rownames(dataFrameWithAverageTumorsAndControlsForCountsIntronicInDeseq) = dataFrameWithAverageTumorsAndControlsForCountsIntronicInDeseq[,1]
dataFrameWithAverageTumorsAndControlsForCountsIntronicInDeseq[,1] = NULL

dataFrameWithAverageTumorsAndControlsForCountsIntronicInNoiseq = merge( geneAverageInTumorsForCountsIntronicInNoiseq, geneAverageInControlsForCountsIntronicInNoiseq, by="row.names", all.x=TRUE, all.y=TRUE )
rownames(dataFrameWithAverageTumorsAndControlsForCountsIntronicInNoiseq) = dataFrameWithAverageTumorsAndControlsForCountsIntronicInNoiseq[,1]
dataFrameWithAverageTumorsAndControlsForCountsIntronicInNoiseq[,1] = NULL


#	CuffDiff
dataFrameWithAverageTumorsAndControlsForPoissonQuartileInCuffDiff = merge( geneAverageInTumorsForPoissonQuartileInCuffdiff, geneAverageInControlsForPoissonQuartileInCuffdiff, by="row.names", all.x=TRUE, all.y=TRUE )
rownames(dataFrameWithAverageTumorsAndControlsForPoissonQuartileInCuffDiff) = dataFrameWithAverageTumorsAndControlsForPoissonQuartileInCuffDiff[,1]
dataFrameWithAverageTumorsAndControlsForPoissonQuartileInCuffDiff[,1] = NULL


#	FPKMs
dataFrameWithAverageTumorsAndControlsForFPKMsInLimma = merge( geneAverageInTumorsForFPKMsInLimma, geneAverageInControlsForFPKMsInLimma, by="row.names", all.x=TRUE, all.y=TRUE )
rownames(dataFrameWithAverageTumorsAndControlsForFPKMsInLimma) = dataFrameWithAverageTumorsAndControlsForFPKMsInLimma[,1]
dataFrameWithAverageTumorsAndControlsForFPKMsInLimma[,1] = NULL

dataFrameWithAverageTumorsAndControlsForFPKMsInGeneSpring = merge( geneAverageInTumorsForFPKMsInGenespring, geneAverageInControlsForFPKMsInGenespring, by="row.names", all.x=TRUE, all.y=TRUE )
rownames(dataFrameWithAverageTumorsAndControlsForFPKMsInGeneSpring) = dataFrameWithAverageTumorsAndControlsForFPKMsInGeneSpring[,1]
dataFrameWithAverageTumorsAndControlsForFPKMsInGeneSpring[,1] = NULL



#
#	Save averages for Comparinson
#
write.table( dataFrameWithAverageTumorsAndControlsForCountsExonicInDeseq, file=paste( working_directory, "/pca_coefficients/raw_averages/averages/averages-counts_exonic_deseq.txt", sep="" ), append=FALSE, row.names=TRUE, col.names=NA, quote=FALSE, sep="\t" )
write.table( dataFrameWithAverageTumorsAndControlsForCountsExonicInNoiseq, file=paste( working_directory, "/pca_coefficients/raw_averages/averages/averages-counts_exonic_noiseq.txt", sep="" ), append=FALSE, row.names=TRUE, col.names=NA, quote=FALSE, sep="\t" )

write.table( dataFrameWithAverageTumorsAndControlsForCountsIntronicInDeseq, file=paste( working_directory, "/pca_coefficients/raw_averages/averages/averages-counts_intronic_deseq.txt", sep="" ), append=FALSE, row.names=TRUE, col.names=NA, quote=FALSE, sep="\t" )
write.table( dataFrameWithAverageTumorsAndControlsForCountsIntronicInNoiseq, file=paste( working_directory, "/pca_coefficients/raw_averages/averages/averages-counts_intronic_noiseq.txt", sep="" ), append=FALSE, row.names=TRUE, col.names=NA, quote=FALSE, sep="\t" )

write.table( dataFrameWithAverageTumorsAndControlsForPoissonQuartileInCuffDiff, file=paste( working_directory, "/pca_coefficients/raw_averages/averages/averages-cuffdiff.txt", sep="" ), append=FALSE, row.names=TRUE, col.names=NA, quote=FALSE, sep="\t" )

write.table( dataFrameWithAverageTumorsAndControlsForFPKMsInLimma, file=paste( working_directory, "/pca_coefficients/raw_averages/averages/averages-fpkm_limma.txt", sep="" ), append=FALSE, row.names=TRUE, col.names=NA, quote=FALSE, sep="\t" )
write.table( dataFrameWithAverageTumorsAndControlsForFPKMsInGeneSpring, file=paste( working_directory, "/pca_coefficients/raw_averages/averages/averages-fpkm_genespring.txt", sep="" ), append=FALSE, row.names=TRUE, col.names=NA, quote=FALSE, sep="\t" )




# ------------------------------------------------------------ PCA --------------------------------------------------------------


#	Counts Exonic
pcaForGenesInCountsExonicDeseq = prcomp( na.omit( t( dataFrameWithAverageTumorsAndControlsForCountsExonicInDeseq ) ), scale=FALSE )
print( summary( pcaForGenesInCountsExonicDeseq ) )

pcaForGenesInCountsExonicNoiseq = prcomp( na.omit( t( dataFrameWithAverageTumorsAndControlsForCountsExonicInNoiseq ) ), scale=FALSE )
print( summary( pcaForGenesInCountsExonicNoiseq ) )


#	Counts Intronic
pcaForGenesInCountsIntronicDeseq = prcomp( na.omit( t( dataFrameWithAverageTumorsAndControlsForCountsIntronicInDeseq ) ), scale=FALSE )
print( summary( pcaForGenesInCountsIntronicDeseq ) )

pcaForGenesInCountsIntronicNoiseq = prcomp( na.omit( t( dataFrameWithAverageTumorsAndControlsForCountsIntronicInNoiseq ) ), scale=FALSE )
print( summary( pcaForGenesInCountsIntronicNoiseq ) )


#	CuffDiff
pcaForGenesInPoissonQuartileForCuffDiff = prcomp( na.omit( t( dataFrameWithAverageTumorsAndControlsForPoissonQuartileInCuffDiff ) ), scale=FALSE )
print( summary( pcaForGenesInPoissonQuartileForCuffDiff ) )


#	FPKMs
pcaForGenesInFPKMsForLimma = prcomp( na.omit( t( dataFrameWithAverageTumorsAndControlsForFPKMsInLimma ) ), scale=FALSE )
print( summary( pcaForGenesInFPKMsForLimma ) )

pcaForGenesInFPKMsForGenespring = prcomp( na.omit( t( dataFrameWithAverageTumorsAndControlsForFPKMsInGeneSpring ) ), scale=FALSE )
print( summary( pcaForGenesInFPKMsForGenespring ) )



# ------------------------------------------------------ PCA Coefficients -------------------------------------------------------


#	Counts Exonic
coefficientsFromPCAForCountsExonicInDeseq = as.data.frame(pcaForGenesInCountsExonicDeseq$rotation)
print( head( coefficientsFromPCAForCountsExonicInDeseq, n=10 ) )
write.table( coefficientsFromPCAForCountsExonicInDeseq, file=paste( working_directory, "/pca_coefficients/raw_averages/coefficients/pca_coefficients-counts_exonic_deseq.txt", sep="" ), append=FALSE, row.names=TRUE, col.names=NA, quote=FALSE, sep="\t" )

coefficientsFromPCAForCountsExonicInNoiseq = as.data.frame(pcaForGenesInCountsExonicNoiseq$rotation)
print( head( coefficientsFromPCAForCountsExonicInNoiseq, n=10 ) )
write.table( coefficientsFromPCAForCountsExonicInNoiseq, file=paste( working_directory, "/pca_coefficients/raw_averages/coefficients/pca_coefficients-counts_exonic_noiseq.txt", sep="" ), append=FALSE, row.names=TRUE, col.names=NA, quote=FALSE, sep="\t" )


#	Counts Intronic
coefficientsFromPCAForCountsIntronicInDeseq = as.data.frame(pcaForGenesInCountsIntronicDeseq$rotation)
print( head( coefficientsFromPCAForCountsIntronicInDeseq, n=10 ) )
write.table( coefficientsFromPCAForCountsIntronicInDeseq, file=paste( working_directory, "/pca_coefficients/raw_averages/coefficients/pca_coefficients-counts_intronic_deseq.txt", sep="" ), append=FALSE, row.names=TRUE, col.names=NA, quote=FALSE, sep="\t" )

coefficientsFromPCAForCountsIntronicInNoiseq = as.data.frame(pcaForGenesInCountsIntronicNoiseq$rotation)
print( head( coefficientsFromPCAForCountsIntronicInNoiseq, n=10 ) )
write.table( coefficientsFromPCAForCountsIntronicInNoiseq, file=paste( working_directory, "/pca_coefficients/raw_averages/coefficients/pca_coefficients-counts_intronic_noiseq.txt", sep="" ), append=FALSE, row.names=TRUE, col.names=NA, quote=FALSE, sep="\t" )


#	CuffDiff
coefficientsFromPCAForPoissonQuartileInCuffDiff = as.data.frame(pcaForGenesInPoissonQuartileForCuffDiff$rotation)
print( head( coefficientsFromPCAForPoissonQuartileInCuffDiff, n=10 ) )
write.table( coefficientsFromPCAForPoissonQuartileInCuffDiff, file=paste( working_directory, "/pca_coefficients/raw_averages/coefficients/pca_coefficients-cuffdiff.txt", sep="" ), append=FALSE, row.names=TRUE, col.names=NA, quote=FALSE, sep="\t" )


#	FPKMs
coefficientsFromPCAForFPKMsInLimma = as.data.frame(pcaForGenesInFPKMsForLimma$rotation)
print( head( coefficientsFromPCAForFPKMsInLimma, n=10 ) )
write.table( coefficientsFromPCAForFPKMsInLimma, file=paste( working_directory, "/pca_coefficients/raw_averages/coefficients/pca_coefficients-fpkm_limma.txt", sep="" ), append=FALSE, row.names=TRUE, col.names=NA, quote=FALSE, sep="\t" )

coefficientsFromPCAForFPKMsInGenespring = as.data.frame(pcaForGenesInFPKMsForGenespring$rotation)
print( head( coefficientsFromPCAForFPKMsInGenespring, n=10 ) )
write.table( coefficientsFromPCAForFPKMsInGenespring, file=paste( working_directory, "/pca_coefficients/raw_averages/coefficients/pca_coefficients-fpkm_genespring.txt", sep="" ), append=FALSE, row.names=TRUE, col.names=NA, quote=FALSE, sep="\t" )



# --------------------------------------------------------- PCA Scores -----------------------------------------------------------


#	Counts Exonic
scoresFromPCAforCountsExonicInDeseq = as.data.frame(pcaForGenesInCountsExonicDeseq$x)
print( head( scoresFromPCAforCountsExonicInDeseq, n=10 ) )
write.table( scoresFromPCAforCountsExonicInDeseq, file=paste( working_directory, "/pca_scores/raw_averages/pca_scores-counts_exonic_deseq.txt", sep="" ), append=FALSE, row.names=TRUE, col.names=NA, quote=FALSE, sep="\t" )

scoresFromPCAForCountsExonicInNoiseq = as.data.frame(pcaForGenesInCountsExonicNoiseq$x)
print( head( scoresFromPCAForCountsExonicInNoiseq, n=10 ) )
write.table( scoresFromPCAForCountsExonicInNoiseq, file=paste( working_directory, "/pca_scores/raw_averages/pca_scores-counts_exonic_noiseq.txt", sep="" ), append=FALSE, row.names=TRUE, col.names=NA, quote=FALSE, sep="\t" )


#	Counts Intronic
scoresFromPCAForCountsIntronicInDeseq = as.data.frame(pcaForGenesInCountsIntronicDeseq$x)
print( head( scoresFromPCAForCountsIntronicInDeseq, n=10 ) )
write.table( scoresFromPCAForCountsIntronicInDeseq, file=paste( working_directory, "/pca_scores/raw_averages/pca_scores-counts_intronic_deseq.txt", sep="" ), append=FALSE, row.names=TRUE, col.names=NA, quote=FALSE, sep="\t" )

scoresFromPCAForCountsIntronicInNoiseq = as.data.frame(pcaForGenesInCountsIntronicNoiseq$x)
print( head( scoresFromPCAForCountsIntronicInNoiseq, n=10 ) )
write.table( scoresFromPCAForCountsIntronicInNoiseq, file=paste( working_directory, "/pca_scores/raw_averages/pca_scores-counts_intronic_noiseq.txt", sep="" ), append=FALSE, row.names=TRUE, col.names=NA, quote=FALSE, sep="\t" )



#	CuffDiff
scoresFromPCAForPoissonQuartileInCuffDiff = as.data.frame(pcaForGenesInPoissonQuartileForCuffDiff$x)
print( head( scoresFromPCAForPoissonQuartileInCuffDiff, n=10 ) )
write.table( scoresFromPCAForPoissonQuartileInCuffDiff, file=paste( working_directory, "/pca_scores/raw_averages/pca_scores-cuffdiff.txt", sep="" ), append=FALSE, row.names=TRUE, col.names=NA, quote=FALSE, sep="\t" )



#	FPKMs
scoresFromPCAForFPKMsInLimma = as.data.frame(pcaForGenesInFPKMsForLimma$x)
print( head( scoresFromPCAForFPKMsInLimma, n=10 ) )
write.table( scoresFromPCAForFPKMsInLimma, file=paste( working_directory, "/pca_scores/raw_averages/pca_scores-fpkm_limma.txt", sep="" ), append=FALSE, row.names=TRUE, col.names=NA, quote=FALSE, sep="\t" )

scoresFromPCAForFPKMsInGenespring = as.data.frame(pcaForGenesInFPKMsForGenespring$x)
print( head( scoresFromPCAForFPKMsInGenespring, n=10 ) )
write.table( scoresFromPCAForFPKMsInGenespring, file=paste( working_directory, "/pca_scores/raw_averages/pca_scores-fpkm_genespring.txt", sep="" ), append=FALSE, row.names=TRUE, col.names=NA, quote=FALSE, sep="\t" )




# ----------------------------------- PCA PLot: Principal Components -----------------------------------------


#	Counts Exonic
png( paste( ggplot_output_dir, "/pca_genes-Counts_Exonic-DESeq.png", sep="" ), width=2056, height=2056, units="px", res=120 )
plot( pcaForGenesInCountsExonicDeseq, type="l", main="PCA Counts Exonic - DESeq - Tumors & Controls - Raw Averages")
dev.off()

png( paste( ggplot_output_dir, "/pca_genes-Counts_Exonic-Noiseq.png", sep="" ), width=2056, height=2056, units="px", res=120 )
plot( pcaForGenesInCountsExonicNoiseq, type="l", main="PCA Counts Exonic - NOISeq - Tumors & Controls - Raw Averages")
dev.off()


#	Counts Intronic
png( paste( ggplot_output_dir, "/pca_genes-Counts_Intronic-DESeq.png", sep="" ), width=2056, height=2056, units="px", res=120 )
plot( pcaForGenesInCountsIntronicDeseq, type="l", main="PCA Counts Intronic - DESeq - Tumors & Controls - Raw Averages")
dev.off()

png( paste( ggplot_output_dir, "/pca_genes-Counts_Intronic-Noiseq.png", sep="" ), width=2056, height=2056, units="px", res=120 )
plot( pcaForGenesInCountsIntronicNoiseq, type="l", main="PCA Counts Intronic - NOISeq - Tumors & Controls - Raw Averages")
dev.off()



#	CuffDiff
png( paste( ggplot_output_dir, "/pca_genes-Poisson_Quartile-CuffDiff.png", sep="" ), width=2056, height=2056, units="px", res=120 )
plot( pcaForGenesInPoissonQuartileForCuffDiff, type="l", main="PCA Poisson Quartile - CuffDiff - Tumors & Controls - Raw Averages")
dev.off()



#	Counts Intronic
png( paste( ggplot_output_dir, "/pca_genes-FPKMs_Limma.png", sep="" ), width=2056, height=2056, units="px", res=120 )
plot( pcaForGenesInFPKMsForLimma, type="l", main="PCA FPKMs - Limma - Tumors & Controls - Raw Averages")
dev.off()

png( paste( ggplot_output_dir, "/pca_genes-FPKMs_GeneSpring.png", sep="" ), width=2056, height=2056, units="px", res=120 )
plot( pcaForGenesInFPKMsForGenespring, type="l", main="PCA FPKMs - GeneSpring - Tumors & Controls - Raw Averages")
dev.off()






#
#	End
#
print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] ", sep="") )
print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Done.", sep="") )
print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] ", sep="") )
