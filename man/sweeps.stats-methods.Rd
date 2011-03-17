\name{sweeps.stats-methods}
\docType{methods}
\alias{sweeps.stats,GENOME-method}
\alias{sweeps.stats-methods}
\alias{get.sweeps,GENOME-method}
\alias{get.sweeps-methods}


\title{Selective Sweeps}


\description{
 This module contains some statistics to detect selective sweeps.
}

\usage{

\S4method{sweeps.stats}{GENOME}(object,new.populations=FALSE,subsites=FALSE)
\S4method{get.sweeps}{GENOME}(object)

}

\arguments{

  \item{object}{an object of class \code{"GENOME"}}
  \item{new.populations}{list of populations. default:\code{FALSE}}
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
          default:\code{FALSE}}
  
}



\value{


returned value is an modified object of class \code{"GENOME"} \cr \cr
-------------------------------------------------------------- \cr
Following Slots will be modified in the \code{"GENOME"} object \cr
-------------------------------------------------------------- \cr

   \item{CL}{Composite Likelihood of SNPs}
   \item{CLmax}{Max. Composite Likeliood of SNPs}

}

\references{
Nielson, R. (2005). \emph{Genomic scans for selective sweeps using SNP data}
Genome Res. 2005 15: 1566-1575
}

\examples{

# GENOME.class <- readData("\home\Alignments")
# GENOME.class
# GENOME.class <- sweeps.stats(GENOME.class)
# GENOME.class <- sweeps.stats(GENOME.class,list(1:7,8:12))
# show the result:
# sweeps.values <- get.sweeps(GENOME.class)
# sweeps.values[[1]] # first population !
# GENOME.class@CLmax
 

}



\keyword{methods}
