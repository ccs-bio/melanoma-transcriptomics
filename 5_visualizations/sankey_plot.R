
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
#	The purpose of this program is to...
#
#
#
#	AUTHOR:	Camilo Valdes (cvaldes3@miami.edu)
#			Computational Biology and Bioinformatics Group
#			Center for Computational Science (CCS), University of Miami
#
# ---------------------------------------------------------------------------------------------------------------------

# 	Import any necessary libraries here

library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(arcdiagram)
require(rCharts)
require(grid)
require(igraph)


# 	Enable printing of large numbers. R defaults to "rounding" for display, and
# 	Cause the warnings to be printed as they occur.
options( digits=15, warn=1 )


#
#	Custom Colors for Color Blind folks
#
# 	The palette with black:
cbPalette_black = c( "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7" )

# 	The palette with grey:
cbPalette_grey = c( "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7" )



# ------------------------------------------------- Project Setup -----------------------------------------------------

print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] ", sep="") )
print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] R Starting... ", sep="") )
print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] ", sep="") )


working_directory="/path/to/linear_model-vs-PCA/Sankey_Plots"
setwd( file.path( working_directory ) )

figures_base_dir=paste( working_directory, "/figures", sep="" )
html_output_dir=paste( figures_base_dir, "/html", sep="" )
dir.create( file.path( html_output_dir ), showWarnings=FALSE, recursive=TRUE )


# ------------------------------------------------------ Main ---------------------------------------------------------

inputData = paste( working_directory, "/data/sankey_input.txt",sep="" )

inputTableWithEdges = read.delim( inputData, header=FALSE, sep="\t", check.names=FALSE )

dataFrameWithEdges = as.data.frame(inputTableWithEdges)
colnames(dataFrameWithEdges) = c("source","target","value")

print( dataFrameWithEdges )

dataFrameWithEdges$source = as.character(dataFrameWithEdges$source)
dataFrameWithEdges$target = as.character(dataFrameWithEdges$target)


sankeyPlot = rCharts$new()
sankeyPlot$setLib('/path/to/R/rCharts_d3_sankey-gh-pages')
sankeyPlot$setTemplate(script = "/path/to/R/rCharts_d3_sankey-gh-pages/layouts/chart.html")

sankeyPlot$set(
				  data = dataFrameWithEdges,
				  nodeWidth = 20,
				  nodePadding = 25,
				  layout = 32,
				  width = 1260,
				  height = 656
)


#
#	BUG!
#	
#	In order to create & display the image we have to build the call as below.!!
#
sankeyPlot
sankeyPlot$save( paste( html_output_dir, "/sankey.html", sep="" ), standalone=TRUE )


# ------------------------------------------------------ Graph ---------------------------------------------------------
#
#	Graph View
#
graphRepresentation = graph.data.frame( d=dataFrameWithEdges, directed=TRUE )

edges = get.edgelist(graphRepresentation)

#	Inspect the edges and verify
edges

# plot( graphRepresentation, layout=layout.reingold.tilford( graphRepresentation, root=1 ) )


# ------------------------------------------------------ Arcs ---------------------------------------------------------
#
#
set.seed(120)

#
#
#
arcplot( 	edges, 
			show.nodes=TRUE,
			sorted=TRUE,
			labels=paste("node",1:21,sep="-"), 
			lwd.arcs=4*runif(10,.5,2), 
			col.arcs=hsv(runif(9,0.6,0.8),alpha=0.4), 
			pch.nodes=21, 
			cex.nodes=runif(10,1,3), col.nodes="gray80", bg.nodes="gray90", lwd.nodes=2 )





# ------------------------------------------------------ END ---------------------------------------------------------

print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] ", sep="") )
print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Done.", sep="") )
print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] ", sep="") )
