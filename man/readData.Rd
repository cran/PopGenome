\name{readData}

\docType{methods}
\alias{readData}


\alias{get.sum.data,GENOME-method}
\alias{get.sum.data-methods}

\title{Reading alignments and calculating summary data}


\description{
 This function reads alignments in FASTA format and calculates some summary data. 
}

\usage{

readData(path,populations=FALSE,outgroup=FALSE,include.unknown=TRUE,
gffpath=FALSE)

\S4method{get.sum.data}{GENOME}(object)

}

\arguments{
 
  \item{object}{object of class \code{"GENOME"}}
  \item{path}{the basepath of the alignments}
  \item{outgroup}{vector of outgroup sequences}
  \item{include.unknown}{if unknown positions should be considered. default:\code{TRUE}}
  \item{populations}{list of populations.default:\code{FALSE}}
  \item{gffpath}{the basepath of the corresponding gff-files. default:\code{FALSE}}
}

\details{
When there is no gff-file specified, an alignment in the right reading frame is expected. Otherwise the examination of synonymous and nonsynonymous positions is useless. \cr \cr
When there is no population defined the whole alignment is considered.
}

\value{


  The function creates an object of class "GENOME" \cr \cr
  --------------------------------------------------------- \cr
  Following Slots will be filled in the "GENOME" object \cr
  --------------------------------------------------------- \cr

            \tabular{rll}{
            \tab Slot                	        \tab Description\cr
          1.\tab \code{n.sites}                        \tab total number of sites\cr
          2.\tab \code{n.biallelic.sites}          	\tab number of biallelic sites \cr
          3.\tab \code{n.gaps}                         \tab number of sites with gaps\cr      
          4.\tab \code{n.unknowns}          	        \tab number of sites with unknown nucleotides\cr
	  5.\tab \code{n.valid.sites}            	\tab number of valid sites\cr
          6.\tab \code{n.polyallelic.sites}            \tab number of sites with >2 nucleotides \cr      
          7.\tab \code{trans.transv.ratio}          	\tab transition/transversion ratio of biallelic sites\cr
	  8.\tab \code{region.data}            	\tab some detail data informations\cr
        }


}


\examples{

# GENOME.class <- readData("...\Alignments")
# GENOEM.class <- readData("...\Alignments", include.unknown=F)
# GENOME.class
# show the result:
# get.sum.data(GENOME.class)
# GENOME.class@region.data

}



\keyword{methods}
