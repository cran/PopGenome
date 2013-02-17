\name{readData}

\docType{methods}
\alias{readData}


\alias{get.sum.data,GENOME-method}
\alias{get.sum.data-methods}

\title{Reading alignments and calculating summary data}


\description{
 This function reads alignments/SNP-data in several formats and calculates some summary data. 
}


\usage{

readData(path,populations=FALSE,outgroup=FALSE,include.unknown=FALSE,
         gffpath=FALSE,format="fasta",parallized=FALSE,
         progress_bar_switch=TRUE, FAST=FALSE,big.data=FALSE,
         SNP.DATA=FALSE
        )

\S4method{get.sum.data}{GENOME}(object)

}

\arguments{
 
  \item{object}{object of class \code{"GENOME"}}
  \item{path}{the basepath of the alignments}
  \item{outgroup}{vector of outgroup sequences}
  \item{include.unknown}{if unknown positions should be considered.}
  \item{populations}{list of populations.default:\code{FALSE}}
  \item{gffpath}{the basepath of the corresponding gff-files. default:\code{FALSE}}
  \item{format}{Data formats. \code{"fasta"} is default. See detail !}
  \item{parallized}{parallele processing. see detail !}
  \item{progress_bar_switch}{progress_bar}
  \item{FAST}{Fast computation. See detail !}
  \item{big.data}{using the ff-package}
  \item{SNP.DATA}{important for reference positions, should be TRUE, if you use SNP-data in alignment format}

}

\details{
If there is no gff-file specified, an alignment in the right reading frame is expected. Otherwise the examination of synonymous and nonsynonymous positions is useless. \cr \cr

format: \cr \cr 

\code{"fasta"},\code{"nexus"},\code{"phylip"}, \cr
\code{"MAF"},\code{"MEGA"},\code{"HapMap"},\code{"VCF"}, 
\code{"VCFhap"} (haploid), \cr
\code{"RData"} \cr \cr

parallized: \cr \cr

- only works on UNIX, because of the multicore package. \cr
- will speed up calculation if you use a huge amount of alignments \cr \cr

FAST:   \cr \cr

- fast computation of biallelic matrix, biallelic sites, transversions/transitions \cr
  and biallelic substitutions                    


big.data: \cr \cr

- using the ff-package \cr \cr
- ff mechanism for biallelic.matrix and gff/gtf information \cr 
- is done automatically for readVCF or readSNP \cr
- Note! should switch to TRUE, if you use big chunks \cr 
  and you want to concatenate them in the PopGenome framework \cr
  (for example: sliding window of the whole data).\cr

SNP.DATA: \cr \cr

- should be switched to \code{TRUE}, if you use SNP-data in alignment format. \cr \cr

}

\value{


  The function creates an object of class "GENOME" \cr \cr
  --------------------------------------------------------- \cr
  Following Slots will be filled in the "GENOME" object \cr
  --------------------------------------------------------- \cr

            \tabular{rll}{
            \tab Slot                	                \tab Description\cr
          1.\tab \code{n.sites}                         \tab total number of sites\cr
          2.\tab \code{n.biallelic.sites}          	\tab number of biallelic sites \cr
          3.\tab \code{n.gaps}                          \tab number of sites with gaps\cr      
          4.\tab \code{n.unknowns}          	        \tab number of sites with unknown nucleotides\cr
	  5.\tab \code{n.valid.sites}            	\tab number of valid sites\cr
          6.\tab \code{n.polyallelic.sites}             \tab number of sites with >2 nucleotides \cr      
          7.\tab \code{trans.transv.ratio}          	\tab transition/transversion ratio of biallelic sites\cr
          8.\tab \code{region.names}          	        \tab names of each region\cr
	  9.\tab \code{region.data}            	        \tab some detail data informations\cr
          
        }


}


\examples{

# GENOME.class <- readData("...\Alignments", FAST=TRUE)
# GENOME.class@region.names
# GENOME.class <- readData("...\Alignments", big.data=TRUE)
# object.size(GENOME.class)
# GENOME.class <- readData("...\Alignments",gffpath="...\Alignments_GFF")
# GENOME.class
# show the result:
# get.sum.data(GENOME.class)
# GENOME.class@region.data

}



\keyword{methods}
