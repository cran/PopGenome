\name{get.codons-methods}
\docType{methods}
\alias{get.codons,GENOME-method}
\alias{get.codons-methods}


\title{Detail informations about codon-substitutions}


\description{
This generic function returns some informations about the codon-changes resulting 
from the observed SNPs.
}


\usage{
\S4method{get.codons}{GENOME}(object, regionID, reading.start.pos, ref.chr)
}

\arguments{
  \item{object}{an object of class \code{"GENOME"}}
  \item{regionID}{what region should be analyzed ?}
  \item{reading.start.pos}{start position of the reading frame}
  \item{ref.chr}{reference chromosome as a FASTA-file}
}


\value{
The function \code{get.codons} returns a data.frame with following informations, \cr

\item{1}{Position of the SNPs} \cr
\item{2}{Mayor Codon} \cr
\item{3}{Minor Codon} \cr  
\item{4}{Mayor Protein}  \cr  
\item{5}{Minor Protein} \cr
\item{6}{synonymous (TRUE/FALSE)} \cr
\item{7}{Polarity of the Mayor Protein} \cr
\item{8}{Polarity of the Minor Protein} \cr  

}

\examples{
# GENOME.class <- readSNP("MPI",CHR=1)
# GENOME.class.split <- splitting.data(GENOME.class,positions=list(10000:20000),type=2)
# get.codons(GENOME.class.split,1,1000,"chromosome1.fas")
}



\keyword{methods}
