
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
#	The purpose of this program is to create a figure of Parental and child Pseudogene lanes.
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


working_directory="/path/to/base/directory/for/Figures/Parental_Pseudo_Lanes"
setwd( file.path( working_directory ) )

figures_base_dir=paste(working_directory, "/figures" , sep="")
png_output_dir=paste( figures_base_dir, "/png", sep="" )

dir.create( file.path( png_output_dir ), showWarnings=FALSE, recursive=TRUE )



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

fileWithAssociations	= paste( working_directory, "/data/_plot_data/associations.txt", sep="" )
associationsDataTable	= read.delim( fileWithAssociations, header=FALSE, check.names=FALSE, sep="\t" )

fileWithPseudogenesFromLRM	= paste( working_directory, "/data/_plot_data/pseudogenes_from_lrm.txt", sep="" )
pseudogenesFromLRMDataTable	= read.delim( fileWithPseudogenesFromLRM, header=FALSE, check.names=FALSE, sep="\t" )
pseudogenesFromLRMDataTable$IsInLRM = 1
pseudogenesFromLRMDataTable

fileWithPseudogeneFCs	= paste( working_directory, "/data/_plot_data/pseudogene_fcs.txt", sep="" )
pseudogeneFCsDataTable	= read.delim( fileWithPseudogeneFCs, header=FALSE, check.names=FALSE, sep="\t" )
colnames(pseudogeneFCsDataTable) = c("Pseudogene.Child", "Child.FC", "Child.Name")


#
#	Melt the data into Long-Format
#
namesOfColumnVariables = c("V1", "V2", "V3")

associationsMelted = melt( associationsDataTable, id.vars=namesOfColumnVariables)
colnames(associationsMelted) = c("GeneID", "Gene.Name", "Parental.FC", "Pseudogene.Child.Number", "Pseudogene.Child")
associationsMelted$Row.ID = seq(1:nrow(associationsMelted))
associationsMelted$Pseudogene.Child.FC = pseudogeneFCsDataTable$Child.FC[match(associationsMelted$Pseudogene.Child, pseudogeneFCsDataTable$Pseudogene.Child)]
associationsMelted$Pseudogene.Child.FC[is.na(associationsMelted$Pseudogene.Child.FC)] = 0
associationsMelted$Pseudogene.Child.Name = pseudogeneFCsDataTable$Child.Name[match(associationsMelted$Pseudogene.Child, pseudogeneFCsDataTable$Pseudogene.Child)]


#	Flag those values in the LRM model
associationsMelted$isInLRM = 0
associationsMelted$isInLRM[associationsMelted$Pseudogene.Child %in% pseudogenesFromLRMDataTable$V1] = 1


#	Remove any extreme outliers
thresholdForOutliers = 40
thresholdForFC = 2
geneToIsolate = "HMGA1"
associationsMelted = associationsMelted[associationsMelted$Pseudogene.Child.FC > -thresholdForOutliers,]
associationsMelted = associationsMelted[associationsMelted$Pseudogene.Child.FC < thresholdForOutliers,]

# associationsMelted = associationsMelted[associationsMelted$Parental.FC > -thresholdForOutliers,]
# associationsMelted = associationsMelted[associationsMelted$Parental.FC < thresholdForOutliers,]

# associationsMelted = associationsMelted[associationsMelted$Pseudogene.Child.FC >= thresholdForFC,]
# associationsMelted = associationsMelted[associationsMelted$Pseudogene.Child.FC <= thresholdForFC,]

associationsMelted = associationsMelted[associationsMelted$Gene.Name == geneToIsolate,]

#
#	Remove DUPLICATES
#
# associationsMelted = subset(associationsMelted,!duplicated(associationsMelted$Gene.Name))


#	Keep only those associations in the LRM model
# associationsMelted = subset(associationsMelted, associationsMelted[,"isInLRM"] != 0)
associationsMelted = subset(associationsMelted, associationsMelted[,"Pseudogene.Child.FC"] != 0)
associationsMelted = associationsMelted[order(associationsMelted["GeneID"]),]

#	Recompute the Row.Ids so we can have a nice ordered plot
associationsMelted$Row.ID = seq(1:nrow(associationsMelted))

#	Create the demarcation markers for the Parent Lanes
spacerSequence = seq(1,nrow(associationsMelted) + 1, 0.5)
spacerSequence = spacerSequence[seq(2, length(spacerSequence), 2)]
associationsMelted$Spacer.locs = spacerSequence

head(associationsMelted,n=100)
nrow(associationsMelted)
as.vector(associationsMelted$Gene.Name)


# --------------------------------------------------- Error Plot ------------------------------------------------------


errorPlot = ggplot(associationsMelted, aes(x=Row.ID, y=Pseudogene.Child.FC, ymax=Pseudogene.Child.FC, ymin=0)) +
				geom_pointrange(aes(group=Gene.Name, colour=Pseudogene.Child.FC>0), size=1) +
				coord_flip() +
				geom_hline( yintercept=0, colour="blue" ) +
				geom_hline( yintercept=2, colour="blue", linetype="dashed", size=0.5 ) +
				geom_hline( yintercept=-2, colour="blue", linetype="dashed", size=0.5 ) +
				# geom_vline( aes(xintercept=Spacer.locs), colour="gray65", size=0.2, linetype="dashed" ) +
				ggtitle( paste( "Parental Gene and Pseudogene Associations", sep="" ) ) +
				ylab( "Pseudogene Expression (log2)" ) +
				# theme(	legend.position="top", axis.text.y=element_text(size=6), axis.ticks=element_blank(),
# 						panel.grid.major = element_blank(),
# 						panel.grid.minor = element_blank()) +
				theme(	legend.position="top", axis.text.y=element_text(size=2)) +
				scale_colour_manual(name="Pseudogene Regulation", labels=c("Down Regulated", "Up Regulated"), values=c("#41A317", "red")) +
				scale_x_discrete(name="Pseudogene Labels", labels=associationsMelted$Pseudogene.Child.Name) # +
# 				ylim(-20,9)
				# scale_x_discrete(name=geneToIsolate, labels=element_blank()) +

errorPlot

output_file_error_plot = paste(png_output_dir, "/", "parental_lane_plot.png", sep="")
ggsave(output_file_error_plot)


q()

# -------------------------------------------------- Scatter Plot -----------------------------------------------------

# 	Compute the R2 for the two controls
correlationCoeff = lm_eqn( associationsMelted[ c( "Parental.FC", "Pseudogene.Child.FC" ) ] )
r2_for_comparo  = format( summary(correlationCoeff)$r.squared, digits=2)

r2_label = paste("R^2 == ", r2_for_comparo, sep="")

#	Position of the label is Top-Left (min.x and max.y)
r2_label_y = max(associationsMelted$Parental.FC)
r2_label_x = min(associationsMelted$Pseudogene.Child.FC)


scatterPlot = ggplot(associationsMelted, aes(x=Pseudogene.Child.FC, y=Parental.FC)) +
				geom_point(shape=1) +
				geom_smooth(method=lm) +
				geom_hline( yintercept=0, colour="black", size=0.6 ) +
				geom_vline( xintercept=0, colour="black", size=0.6 ) +
				ggtitle( paste( "LRM Parental Gene and All Pseudogene Associations", sep="" ) ) +
				xlab( "Pseudogene Expression (log2 FC)") + 
				ylab( "Parental Expression (log2 FC)" ) +
				annotate("text", x=r2_label_x, y=r2_label_y, vjust=3, hjust=-0.3, label=r2_label, parse=TRUE, color='blue')

scatterPlot

outputFileForScatterPlot = paste(png_output_dir, "/", "scatter_plot.png", sep="")
ggsave(outputFileForScatterPlot)


# ------------------------------------------------------ END ---------------------------------------------------------

print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] ", sep="") )
print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Done.", sep="") )
print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] ", sep="") )
