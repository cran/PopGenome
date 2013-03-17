\name{MKT-methods}
\docType{methods}
\alias{MKT,GENOME-method}
\alias{MKT-methods}
\alias{get.MKT,GENOME-method}
\alias{get.MKT-methods}

\title{McDonald-Kreitman Test (McDonald \& Kreitman 1991)}


\description{
 This generic function calculates an approximate version of the McDonald-Kreitman Test.
}

\usage{

\S4method{MKT}{GENOME}(object,new.populations=FALSE)
\S4method{get.MKT}{GENOME}(object)

}

\arguments{

  \item{object}{an object of class \code{"GENOME"}}
  \item{new.populations}{list of populations. default:\code{FALSE}}
}

\details{

This function assumes that in population genetic analysis the probability of two variants \cr
in one codon is very small. Due to this only single nucleotide polymorphisms (SNPs) are examined. \cr 
When there was no gff-file specified, an alignment in the right reading frame is expected. \cr
The outgroup has to be defined as a population !

}



\value{


  Returned value is an modified object of class "GENOME" \cr \cr
  --------------------------------------------------------- \cr
  Following slots will be modified in the "GENOME" object \cr
  --------------------------------------------------------- \cr

  \item{MKT}{a matrix which includes following values:}

   \tabular{rll}{
            \tab Columns               	                 \tab Description\cr
          1.\tab \code{P_nonsyn}                         \tab nonsynonymous sites\cr
          2.\tab \code{P_syn}         	                 \tab synonymous sites\cr
          3.\tab \code{D_nonsyn}                         \tab fixed nonsynonymous sites\cr      
          4.\tab \code{D_syn}         	                 \tab fixed synonymous sites\cr
	  5.\tab \code{neutrality.index}           	 \tab $(P_{nonsyn}/P_{syn})/(D_{nonsyn}/D_{syn})$\cr
          6.\tab \code{alpha}                            \tab 1-neutrality.index\cr        
        }
  
  
}

\references{ 

McDonald, J. H.; Kreitman, M. (1991). \emph{Adaptive protein evolution at the Adh locus in Drosophila.} Nature 351 (6328): 652-654 \cr \cr 

}


\examples{

# GENOME.class <- readData("\home\Alignments")
# GENOME.class
# GENOME.class <- MKT(GENOME.class)
# GENOME.class <- MKT(GENOME.class,list(1:7,8:12))
# show the result:
# get.MKT(GENOME.class)
}



\keyword{methods}
