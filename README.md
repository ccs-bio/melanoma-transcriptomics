# Melanoma Transcriptomics

The contents of this repository are the source code scripts for the paper "**Ensemble Model Approach Targeting Pseudogenes in the Melanoma Transcriptome**".  The software contained here covers the methods mentioned in the paper, along with miscellaneous programs, scripts and notes.


## How To Get Started

- [Download the Code](https://github.com/ccs-bio/melanoma-transcriptomics) and follow the standard instructions on how to install the framework/library.

## Communication

- If you **found a bug**, open an issue and _please provide detailed steps to reliably reproduce it_.
- If you have **feature request**, open an issue.
- If you **would like to contribute**, please submit a pull request.

## Requirements

The basic requirements for the software available here are:

- [R, 3.3.1 or above](https://www.r-project.org/)
- [Python 2.7.x](https://www.python.org/)
- [Perl 5.18.x](https://www.perl.org/)
- [GeneSpring GX](http://www.genomics.agilent.com/en/Microarray-Data-Analysis-Software/GeneSpring-GX/?cid=AG-PT-130&tabId=AG-PR-1061)
- [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)

You will also need a reasonably modern computer and operating system that can support the above.  Note that some packages have their own dependencies, and you should refer to their respective user manuals and installation instructions for details.

### BioPerl
[BioPerl](http://bioperl.org/) is required for parsing the results from BLAST.

### Bowtie2

[Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) is required for the alignment step.  See the Bowtie2 [Getting Started](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#getting-started-with-bowtie-2-lambda-phage-example) guide for instructions on installing and running Bowtie2.

### Samtools

[Samtools](http://www.htslib.org/) is required for parsing the alignments and generating miscellaneous alignment statistics.  See the [Samtools Documentation](http://www.htslib.org/doc/) for instructions on installing and running Samtools.

### R

The latest version of [R](http://www.r-project.org/) is required for calculating the bacterial p-values, and visualizing the read-hit distribution.  You will also need the following R libraries & packages from [CRAN](http://cran.r-project.org/):

- [ggplot2](http://ggplot2.org/)
- [reshape2](http://cran.r-project.org/web/packages/reshape2/index.html)
- [plyr](http://cran.r-project.org/web/packages/plyr/index.html)
- [DESeq](http://bioconductor.org/packages/release/bioc/html/DESeq.html)
- [limma](https://bioconductor.org/packages/release/bioc/html/limma.html)
- [NOISeq](https://www.bioconductor.org/packages/release/bioc/html/NOISeq.html)



## Contact

Contact [Camilo Valdes](mailto:cvaldes3@miami.edu) for pull requests, bug reports, good jokes and coffee recipes.

### Maintainers

- [Camilo Valdes](mailto:cvaldes3@miami.edu)


## License

The software in this repository is available under the [GNU General Public License, Version 3](https://github.com/ccs-bio/melanoma-transcriptomics/blob/master/LICENSE).  See the [LICENSE](https://github.com/ccs-bio/melanoma-transcriptomics/blob/master/LICENSE) file for more information.
