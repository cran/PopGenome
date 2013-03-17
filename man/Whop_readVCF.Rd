\name{Whop_readVCF}

\docType{methods}
\alias{Whop_readVCF}

\title{Reading tabixed VCF files (An interface to WhopGenome)}


\description{
 This function provides an interface to the WhopGenome package which 
 is specialized on reading tabixed VCF files.

}

\usage{

Whop_readVCF(v, numcols, tid, frompos, topos,
        samplenames=NA, gffpath = FALSE, include.unknown=FALSE)

}

\arguments{
 
  \item{v}{An vcf_handle returned from \code{vcf_open()}}
  \item{numcols}{number of SNPs should be read in as a chunk}
  \item{tid}{which chromosome ? (character)}
  \item{frompos}{start of the region}
  \item{topos}{end of the region}
  \item{samplenames}{a vector of individuals}
  \item{gffpath}{the corresponding GFF-file}
  \item{include.unknown}{including unknown positions}
}

\details{
WhopGenome is required ! \code{require(WhopGenome)}
WhopGenome provides some powerful filter meachanisms which can be applied to 
the VCF reading process. The filter rules can be set via WhopGenome functions.
\code{Whop_readVCF} expects a \code{vcf_handle} returned from \code{vcf_open}.
The Whop_readVCF function expects a tabixed VCF with a diploid GT-field. \cr
In case of haploid data, the GT-field have to be transformed to a pseudo 
diploid field like (0 -> 0|0). An alternative is readData(..., format="VCFhap") 
which can read non-tabixed haploid VCFs directly. \cr 
The ff-package we use is limited by
individuals * (number of SNPs) <= .Machine$integer.max \cr
In case of really huge data sets the bigmemory package will be used. (slower) \cr
See also readData(..., format="VCF") ! \cr
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

# require(WhopGenome)
# vcf_handle   <- vcf_open("chr2.vcf.gz")
# GENOME.class <- Whop_readVCF(vcf_handle, 1000, "2", 1, 100000)
# GENOME.class
# GENOME.class@region.names

}



\keyword{methods}
