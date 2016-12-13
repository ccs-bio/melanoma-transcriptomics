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
#	The purpose of this program is to create a scatter plot of Parental and child Pseudogene associations.
#
#
#
#	AUTHOR:	Camilo Valdes (cvaldes3@miami.edu)
#			Computational Biology and Bioinformatics Group
#			Center for Computational Science (CCS), University of Miami
#
# ---------------------------------------------------------------------------------------------------------------------

#
# 	Import any necessary system or 3rd-party libraries & frameworks here
#
library(ggplot2)
library(RColorBrewer)
library(reshape2)

#
#	Import any R files here.  These contain custom classes, functions, etc.
#


# 	Set some options to alter the default behavior of the REPL.
options( width=512, digits=15, warn=1, echo=TRUE )


#
#	Custom Colors for Color Blind folks
#
# 	The palette with black:
cbbPalette_black = c( "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7" )

# 	The palette with grey:
cbPalette_grey = c( "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7" )





# ------------------------------------------------- Project Setup -----------------------------------------------------

print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] ", sep="") )
print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] R Starting... ", sep="") )
print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] ", sep="") )


working_directory="/path/to/Paper/Figures/Scatter_Plot"
setwd( file.path( working_directory ) )

figures_base_dir=paste(working_directory, "/figures" , sep="")
dir.create( file.path( figures_base_dir ), showWarnings=FALSE, recursive=TRUE )


# ------------------------------------------------------ Functions ---------------------------------------------------------

#'	Returns a fitted linear-model object for calculating R^2
#'  @param df	A Data Frame with the columns to calculate the R^2 from
#'  @returns	m A linear model object
lm_eqn = function( df )
{
    y <- df[[1]]
    x <- df[[2]]

    m = lm( y ~ x, df )
    return(m)
}

# ------------------------------------------------------ Main ---------------------------------------------------------

#	Read-in the Parental-Pseudogene associations
fileWithAssociations	= paste( working_directory, "/data/all_parental_pseudogenes.txt", sep="" )
fileWithAssociationsNumberCol = max(count.fields(fileWithAssociations, sep="\t"))
associationsDataTable	= read.table( 	fileWithAssociations, header=FALSE, check.names=FALSE, sep="\t", 
										fill=TRUE, col.names=1:fileWithAssociationsNumberCol, na.strings=c(""," ", "NA"),
										stringsAsFactors=FALSE )
																				
#	Read-in the Parental FCs
parentalGeneFCs 		 = paste( working_directory, "/data/parentGenes_deg_out.txt", sep="" )
parentalGeneFCsDataTable = read.table(parentalGeneFCs, header=TRUE, check.names=FALSE, sep="\t", stringsAsFactors=FALSE)
parentalGeneFCsDataTable[is.na(parentalGeneFCsDataTable)] = 0

#	Read-in the Pseudogene FCs
pseudogeneFCs			 = paste( working_directory, "/data/pseudogene_deg_out.txt", sep="" )
pseudogeneFCsDataTable	 = read.table(pseudogeneFCs, header=TRUE, check.names=FALSE, sep="\t", stringsAsFactors=FALSE)
pseudogeneFCsDataTable[is.na(pseudogeneFCsDataTable)] = 0


#	Get the MAX FCs for each Parental & Pseudogene — these are the FC values we'll plot
parentalGeneFCsDataTable$max = apply(parentalGeneFCsDataTable[, 2:6], 1, function(x) x[which.max(abs(x))] )
pseudogeneFCsDataTable$max   = apply(pseudogeneFCsDataTable[, 2:6], 1, function(x) x[which.max(abs(x))] )


columnNames = c( "Parent.ID", paste("pseudo_", seq(1:9), sep="") )
colnames(associationsDataTable) = columnNames

#	Melt the associations table so we have it in long-format
associationsDataTable[is.na(associationsDataTable)] = 0
associationsDataTableMelted = melt(associationsDataTable, id.vars="Parent.ID")
associationsDataTableMelted = subset(associationsDataTableMelted, associationsDataTableMelted[,"value"] != 0)
associationsDataTableMelted$Row.ID = seq(1:nrow(associationsDataTableMelted))

colnames(associationsDataTableMelted) = c("Parent.ID", "Pseudo.Child.Number", "Pseudogene.ID", "Row.ID")

#	Grab the corresponding Parental & Pseudogene FCs
associationsDataTableMelted$Parental.FC = parentalGeneFCsDataTable$max[match(associationsDataTableMelted$Parent.ID, parentalGeneFCsDataTable$GeneID)]
associationsDataTableMelted$Pseudogene.FC = pseudogeneFCsDataTable$max[match(associationsDataTableMelted$Pseudogene.ID, pseudogeneFCsDataTable$GeneID)]

#	Remove any extreme outliers - Pseudogenes are the problem
thresholdForOutliers = 10
associationsDataTableMelted = associationsDataTableMelted[associationsDataTableMelted$Parental.FC != 0,]
associationsDataTableMelted = associationsDataTableMelted[associationsDataTableMelted$Pseudogene.FC != 0,]

#	Parental Threshold
associationsDataTableMelted = associationsDataTableMelted[associationsDataTableMelted$Parental.FC > -thresholdForOutliers,]
associationsDataTableMelted = associationsDataTableMelted[associationsDataTableMelted$Parental.FC < thresholdForOutliers,]

#	Pseudogene Threshold
associationsDataTableMelted = associationsDataTableMelted[associationsDataTableMelted$Pseudogene.FC > -thresholdForOutliers,]
associationsDataTableMelted = associationsDataTableMelted[associationsDataTableMelted$Pseudogene.FC < thresholdForOutliers,]

numberOfAssociations = nrow(associationsDataTableMelted)
numberOfAssociations



# -------------------------------------------------- Scatter Plot -----------------------------------------------------

# 	Compute the R2 for the two controls
correlationCoeff = lm_eqn( associationsDataTableMelted[ c( "Parental.FC", "Pseudogene.FC" ) ] )
r2_for_comparo  = format( summary(correlationCoeff)$r.squared, digits=2)

r2_label = paste("R^2 == ", r2_for_comparo, sep="")

#	Position of the label is Top-Left (min.x and max.y)
r2_label_y = max(associationsDataTableMelted$Parental.FC)
r2_label_x = min(associationsDataTableMelted$Pseudogene.FC)


scatterPlot = ggplot(associationsDataTableMelted, aes(x=Pseudogene.FC, y=Parental.FC)) +
				geom_point(shape=1) +
				geom_smooth(method=lm) +
				geom_hline( yintercept=0, colour="black", size=0.6 ) +
				geom_vline( xintercept=0, colour="black", size=0.6 ) +
				ggtitle( paste( "5 Method Parental Gene and Pseudogene Associations", sep="" ) ) +
				xlab( "Pseudogene Expression (log2 FC)") + 
				ylab( "Parental Expression (log2 FC)" ) # +
# 				annotate("text", x=r2_label_x, y=r2_label_y, vjust=3, hjust=-0.3, label=r2_label, parse=TRUE, color='blue')

scatterPlot

outputFileForScatterPlot = paste(figures_base_dir, "/", "scatter_plot-", numberOfAssociations, "-@", thresholdForOutliers, ".png", sep="")
ggsave(outputFileForScatterPlot)


# ------------------------------------------------------ END ---------------------------------------------------------

print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] ", sep="") )
print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Done.", sep="") )
print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] ", sep="") )

