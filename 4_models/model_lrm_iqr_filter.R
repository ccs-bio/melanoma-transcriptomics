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
library(grid)


# 	Enable printing of large numbers. R defaults to "rounding" for display, and
# 	Cause the warnings to be printed as they occur.
options( width=220, digits=15, warn=1 )


# ----------------------------------------------------- Setup ---------------------------------------------------------

print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] R Starting... ", sep="") )

working_directory="/path/to/data/with/linear_models/IQR_2"
setwd( file.path( working_directory ) )


fileToProcess = paste( working_directory, "/data-averages_for_IQR.txt", sep="" )

expressionTable = read.csv( fileToProcess, sep="\t", header=TRUE, row.names=1, check.names=FALSE )

dataFrameForIQR = as.data.frame( expressionTable )

print( head(dataFrameForIQR, n=10))


# ----------------------------------------------- Coefficients & IQR -------------------------------------------------
#
#	Clean up the coefficients data structure and use the IQR strategy to extract the significant genes
#

# columnToInspect = "Fitted.with.Significant.Coefficients"
columnToInspect = "Fitted.with.Significant.Coefficients.No.Limma"

#	Constant that will be used to filter the iqr ranges
constantForIQR = 2


dataFrameForBarPlot = data.frame( seq(from=1,to=length(dataFrameForIQR[,c(columnToInspect)]),by=1), dataFrameForIQR[,c(columnToInspect)] )
colnames(dataFrameForBarPlot) = c("x", "y")
rownames(dataFrameForBarPlot) = rownames(dataFrameForIQR)

print(head(dataFrameForBarPlot,n=5))


rangeIQR = IQR( dataFrameForBarPlot$y )

IQRtimesConstant = rangeIQR * constantForIQR

print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] IQR: ", rangeIQR, sep="") )
print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] IQR * ", constantForIQR, ": ", IQRtimesConstant, sep="") )
print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"]", sep="") )


#
#	Extract the new outliers
#
genesGreaterThanIQR 	= dataFrameForBarPlot[dataFrameForBarPlot$y > IQRtimesConstant,]
genesLowerThanIQR   	= dataFrameForBarPlot[dataFrameForBarPlot$y < -IQRtimesConstant,]

numberGreaterThanIQR 	= nrow(genesGreaterThanIQR)
numberLessThanIQR    	= nrow(genesLowerThanIQR)
TotalOutliersByIQR 		= numberGreaterThanIQR + numberLessThanIQR


print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Higher: ", numberGreaterThanIQR, sep="") )
print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Lower: ", numberLessThanIQR, sep="") )
print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Total: ", TotalOutliersByIQR, sep="") )


appendedDataFrameFromIQR = rbind(genesGreaterThanIQR, genesLowerThanIQR)

outputFileforIQROutliers = paste( working_directory, "/iqr_threshold_", constantForIQR,"/outliers-", columnToInspect, ".", constantForIQR, "-SAMMSON.txt", sep="" )
write.table( appendedDataFrameFromIQR, file=outputFileforIQROutliers, append=FALSE, row.names=TRUE, col.names=NA, quote=FALSE, sep="\t", na="0" )


# ----------------------------------------------- BoxPlot -------------------------------------------------

dataFrameForIQRMelted = melt( as.matrix(dataFrameForIQR), id.vars="ROW_HEADER", measure.vars=c(1:2))

print(head(dataFrameForIQRMelted,n=10))

IQRboxPlot = ggplot(dataFrameForIQRMelted, aes(x=Var2, y=value)) + 
				geom_boxplot() + 
				geom_hline( yintercept=IQRtimesConstant, colour="red", linetype="dashed", size=0.3 ) +
				geom_hline( yintercept=-1.8703421, colour="blue", size=0.3 ) +
				geom_hline( yintercept=-1.1526562, colour="green", size=0.3 ) +
				geom_hline( yintercept=-IQRtimesConstant, colour="red", linetype="dashed", size=0.3 ) +
				ylab("LRM Coefficients") +
				xlab("Condition")
						
IQRboxPlot
