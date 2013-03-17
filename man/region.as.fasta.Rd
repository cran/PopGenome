\name{region.as.fasta-methods}
\docType{methods}

\alias{region.as.fasta,GENOME-method}
\alias{region.as.fasta-methods}


\title{Extract a region as a FASTA file}


\description{
This generic function writes a FASTA file of the observed biallelic positions to the current workspace.
}

\usage{

\S4method{region.as.fasta}{GENOME}(object,region.id=FALSE,filename=FALSE,type=1,ref.chr=FALSE)

}

\arguments{

  \item{object}{an object of class \code{"GENOME"}}
  \item{region.id}{region of the genome}
  \item{filename}{name of the FASTA file}
  \item{type}{1: extract SNPs 2: extract all nucleotides}
  \item{ref.chr}{reference sequence}
}



\examples{

#GENOME.class       <- readSNP("Arabidopsis",CHR=1)
#GENOME.class.split <- splitting.data(GENOME.class, positions=list(100:200))
#region.as.fasta(GENOME.class.split,1,"my_fasta_file.fas")

}



\keyword{methods}

