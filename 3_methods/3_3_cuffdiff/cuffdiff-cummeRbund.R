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
#	The purpose of this program is to run the CummeRbund analysis workflow a set of Cuffdiff results.
#
#
#
#	AUTHOR:	Camilo Valdes (cvaldes3@med.miami.edu)
#			Computational Biology and Bioinformatics Group
#			Center for Computational Science (CCS), University of Miami
#
# ---------------------------------------------------------------------------------------------------------------------

# 	Import any necessary libraries here

library(cummeRbund)

options( digits = 15 ) # Enable printing of large numbers. R defaults to "rounding" for display.


# ----------------------------------------------------- Setup ---------------------------------------------------------

print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] R Starting... ", sep="") )

working_directory="/path/to/base/analysis/directory/cuffdiff"

setwd( file.path( working_directory ) )



# ------------------------------------------------------ Main ---------------------------------------------------------

#	Annotations file
gtf_file='/path/to/ensembl/72/gtf/Homo_sapiens-CLEAN.GRCh37.72.gtf'


# 	Dispersion methods tested
# dispertion_methods = c( 'blind', 'per-condition', 'poisson', 'pooled' )
# dispertion_methods = c( 'per-condition', 'poisson', 'pooled' )
dispertion_methods = c( 'poisson' )


# 	Normalization methods tested
# normalization_methods = c( 'classic-fpkm', 'geometric', 'quartile' )
normalization_methods = c( 'quartile' )


# 	Base directory with cuffdiff results
cuffdiff_results_base = '/path/to/base/results/directory/Melanoma/cuffdiff'


# Read-alignment Types for which counts have been extracted
mapping_types = c( 'unique' )



#
#	Build & Run
#
for( dispersionMethod in dispertion_methods ) {
	
	for( normalizationMethod in normalization_methods) {
		
		for( mappingType in mapping_types ) {
			
			
			print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] ", dispersionMethod, " - ", normalizationMethod, sep="") )
		
		
			# ------------------------------------------ Loading ------------------------------------------
		
			results_name = paste( "results_", dispersionMethod, "_", normalizationMethod, sep="" )
		
			cuffdiff_results_dir = paste( cuffdiff_results_base, "/", results_name, sep="" )
		
			# 	The CuffDiff results
			# cuffdiff_results = readCufflinks( dir=cuffdiff_results_dir, rebuild=TRUE )
			cuffdiff_results = readCufflinks( dir=cuffdiff_results_dir )

			cuffdiff_results
		
		
		
			# ---------------------------------- Differential Expression ----------------------------------

			print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Differential Expression...", sep="") )

			sigGeneIds = getSig( cuffdiff_results, alpha=0.05, level="genes" )

			print( head( sigGeneIds, n=10 ) )
		
			print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] No. Diff. Exp. Genes: ", length( sigGeneIds ), sep="") )
			
				
		
			# ---------------------------------------- DEG Output ----------------------------------------
		
			print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Writing DEG Table...", sep="") )
			
			# 	Get the Full Output table (pvalues, Fold Changes, etc.)
			#	This is the table that we'll write-out to disk
			tableWithDifferentialData = diffData( genes( cuffdiff_results ) )
			
		
			outputDir = paste( mappingType, "/", results_name, sep="" )
		
			dir.create( file.path( working_directory, outputDir ), showWarnings=FALSE, recursive=TRUE )


			outputFile = paste(  outputDir, "/AllGenes_SigDiffExp-", mappingType, ".txt", sep="" )

			write.table( tableWithDifferentialData, file=outputFile, append=FALSE, sep="\t" )
			
			
			
			# --------------------------------- Per-Replicate FPKM Table Output ---------------------------
			
			
			print( paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Writing Replicate FPKM Table...", sep="") )
			
			#
			#	Replicate Matrix
			#
			tableWithReplicateFPKMs = repFpkmMatrix( genes( cuffdiff_results ) )
			
			print( head( tableWithReplicateFPKMs, n=10 ) )
			
			
			
			outputFile = paste(  outputDir, "/aggregated_expression_profile_for_fpkm_AllGenes.txt", sep="" )

			write.table( tableWithReplicateFPKMs, file=outputFile, append=FALSE, quote=FALSE, sep="\t" )
			
			
			
			
		}
		
	}
	
}



