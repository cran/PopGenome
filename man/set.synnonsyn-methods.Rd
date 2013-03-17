\name{set.synnonsyn-methods}
\docType{methods}
\alias{set.synnonsyn,GENOME-method}
\alias{set.synnonsyn-methods}


\title{Set synonymous positions for SNP data}


\description{
 This generic function classifies the observed biallelic positions read from SNP data files into synonymous and non-synonymous SNPs.\cr
}

\usage{

\S4method{set.synnonsyn}{GENOME}(object,ref.chr)



}

\arguments{

  \item{object}{an object of class \code{"GENOME"}}
  \item{ref.chr}{the reference chromosome in FASTA format}
}



\value{


The return value is a modified object of class \code{"GENOME"} \cr 
  
  
}

\note{
The data has to be read in with a corresponding GFF/GTF file (including the CDS fields); 
otherwise a correct classification is not possible.
}

\examples{

# GENOME.class <- readData("\home\Alignments",gffpath="GFF.Folder")
# GENOME.class 
# GENOME.class <- set.synnonsyn(GENOME.class,ref.chr="chromosome1.fas")
 

}



\keyword{methods}
