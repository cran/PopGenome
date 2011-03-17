\name{readVCF}

\docType{methods}
\alias{readVCF}

\title{Reading SNP data from the 1000 Genome-project}


\description{
 This function reads tabixed VCF-files from the 1000 Genome-project. 
}

\usage{

readVCF(filename, numcols, tid, frompos, topos,
        samplenames=NA, gffpath = FALSE )

}

\arguments{
 
  \item{filename}{the corresponding VCF-file}
  \item{numcols}{number of SNPs should be read in as a chunk}
  \item{tid}{which chromosome ? (character)}
  \item{frompos}{start of the region}
  \item{topos}{end of the region}
  \item{samplenames}{a vector of individuals}
  \item{gffpath}{the corresponding GFF-file}
}

\details{

The ff-package we use to store the SNP informations is limited by \cr
individuals * (number of SNPs) <= .Machine$integer.max \cr
Use the function \code{vcf_handle <-.Call("VCF_open", filename)} \cr 
to open a VCF-file and \code{.Call("VCF_getSampleNames",vcf_handle)} \cr 
to get the individual names. \cr

}

\value{


  The function creates an object of class "GENOME" \cr \cr
  --------------------------------------------------------- \cr
  Following Slots will be filled in the "GENOME" object \cr
  --------------------------------------------------------- \cr

            \tabular{rll}{
            \tab Slot                	        \tab Description\cr
          1.\tab \code{n.sites}                 \tab total number of sites\cr
          2.\tab \code{n.biallelic.sites}       \tab number of biallelic sites \cr
          3.\tab \code{region.data}            	\tab some detail data informations\cr
          4.\tab \code{region.names}            \tab names of each region \cr 
        }


}


\examples{

# GENOME.class <- readVCF("...\chr1.vcf.gz", 1000, "1", 1, 100000)
# GENOME.class
# GENOME.class@region.names
# GENOME.class <- neutrality.stats(GENOME.class,FAST=TRUE)
# show the result:
# get.sum.data(GENOME.class)
# GENOME.class@region.data

}



\keyword{methods}
