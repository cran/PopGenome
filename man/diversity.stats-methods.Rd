 \name{diversity.stats-methods}

\docType{methods}

\alias{diversity.stats,GENOME-method}
\alias{diversity.stats-methods}


\title{Diversities}

\description{
 A generic function to calculate some nucleotide & haplotype diversities.
}

\usage{

\S4method{diversity.stats}{GENOME}(object,new.populations=FALSE,subsites=FALSE,pi=FALSE)

}

 \arguments{

  \item{object}{An object of class \code{"GENOME"}}
  \item{new.populations}{list of populations. default=\code{FALSE}}
  \item{subsites}{ 

                  \code{"transitions"}: SNPs which are transitions. \cr
                  \code{"transversions"}: SNPs which are transversions. \cr
                  \code{"syn"}: synonymous sites. \cr
                  \code{"nonsyn"}:  nonsynonymous sites. \cr
                  \code{"exon"}: SNPs in exon regions. \cr
                  \code{"intron"}: SNPs in intron regions. \cr
                  \code{"coding"}: SNPs in coding regions (CDS). \cr
                  \code{"utr"}:  SNPs in UTR regions. \cr
                  \code{"gene"}: SNPs in genes. \cr
		  \code{"intergenic"} : SNPs in intergenic regions. \cr
  } 

  \item{pi}{Nei's pi calculation}

}


\value{

Returned value is an modified object of class \code{"GENOME"} \cr  \cr
--------------------------------------------------------- \cr
Following slots will be modified in the \code{"GENOME"}  object \cr
--------------------------------------------------------- \cr

\tabular{rlll}{
            \tab Slot                	    \tab Reference   \tab Description\cr
          1.\tab \code{nuc.diversity.within}	\tab [1]         \tab Nucleotide diversity (within the population) \cr
          2.\tab \code{Pi}	                \tab [2]         \tab Diversity from Nei  (within the population)    \cr 
          3.\tab \code{hap.diversity.within}	\tab [1]         \tab Haplotype diversity (within the population)\cr
          }

}

\references{
[1] Hudson, R. R., M. Slatkin, and W.P. Maddison (1992). \emph{Estimating of levels of gene flow from DNA sequence data.} Gentics 13(2),583-589 \cr \cr
[2] Nei, M. (1987). \emph{Molecular Evolutionary Genetics.} Columbia Univ. Press, New York. \cr
}

\examples{

# GENOME.class <- readData("\home\Alignments")
# GENOME.class
# GENOME.class <- diversity.stats(GENOME.class)
# GENOME.class <- diversity.stats(GENOME.class,list(1:4,5:10))
# GENOME.class <- diversity.stats(GENOME.class,
# list(c("seq1","seq5","seq3"),c("seq2","seq8")))
# show the result:
# GENOME.class@nuc.diversity.within

}


\keyword{methods}
