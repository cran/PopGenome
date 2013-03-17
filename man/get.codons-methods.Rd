\name{get.codons-methods}
\docType{methods}
\alias{get.codons,GENOME-method}
\alias{get.codons-methods}


\title{Detailed information about the nature of codon changes}


\description{
This generic function returns some information about the codon changes resulting 
from the observed SNPs.
}


\usage{
\S4method{get.codons}{GENOME}(object, regionID, reading.start.pos, ref.chr, 
SNP.DATA=FALSE, reverse.strand=FALSE)
}

\arguments{
  \item{object}{an object of class \code{"GENOME"}}
  \item{regionID}{what region/alignment should be analyzed ?}
  \item{reading.start.pos}{start position of the reading frame (nucleotide position of the reference seq.)}
  \item{ref.chr}{reference chromosome as a FASTA file}
  \item{SNP.DATA}{in case of SNP data. see details !}
  \item{reverse.strand}{in case the coding region is on a reverse strand}
}
\details{
In case of calculations based on alignments, the parameters \cr
\code{reading.start.pos} and \code{ref.chr} are not required. \cr
The alignment data has to be read in with the readData function \cr
and the corresponding GFF-files to verify synonymous and nonsynonymous \cr
sites. In case of SNP-data the GFF-file is not neccessary, but the start position \cr
of the reading frame is needed. When splitting the data into coding regions \cr
(\code{splitting.data}), the corresponding information is in the slot \code{GENOME.class.split@region.names}
}

\value{
The function \code{get.codons} returns a data.frame with the following information \cr

\item{1}{Position of the SNPs} 
\item{2}{Major Codon} 
\item{3}{Minor Codon}   
\item{4}{Major amino acid}    
\item{5}{Minor amino acid} 
\item{6}{synonymous (TRUE/FALSE)} 
\item{7}{Polarity of the major amino acid} 
\item{8}{Polarity of the minor amino acid}   

}

\examples{
# GENOME.class <- readSNP("MPI",CHR=1)
# GENOME.class.split <- splitting.data(GENOME.class,positions=list(10000:20000),type=2)
# get.codons(GENOME.class.split,1,1000,"chromosome1.fas",SNP.DATA=TRUE)
# When splitting into coding regions (splitting.data(,subsites="coding")) 
# the reading.start.position is in split@region.names
# GENOME.class <- readData("Alignments", gffpath="GFF_folder")
# get.codons(GENOME.class.split,2)
}



\keyword{methods}
