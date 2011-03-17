\name{set.synnonsyn-methods}
\docType{methods}
\alias{set.synnonsyn,GENOME-method}
\alias{set.synnonsyn-methods}


\title{Sets synonymous positions for SNP data}


\description{
 This generic function verifies which of the observed biallelic positions are synonymous. 
}

\usage{

\S4method{set.synnonsyn}{GENOME}(object,ref.chr)



}

\arguments{

  \item{object}{an object of class \code{"GENOME"}}
  \item{ref.chr}{the reference chromosome in FASTA format}
}



\value{


returned value is an modified object of class \code{"GENOME"} \cr 
  
  
}

\examples{

# GENOME.class <- readData("\home\Alignments",gffpath="GFF.Folder")
# GENOME.class 
# GENOME.class <- set.synnonsyn(GENOME.class,ref.chr="chromosome1.fas")
 

}



\keyword{methods}
