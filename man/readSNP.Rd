\name{readSNP}

\docType{methods}
\alias{readSNP}

\title{Reading SNP data from the 1001 Genome-project}


\description{
 This function reads SNP data (quality_variant) from the 1001 Genome-project. 
}

\usage{

readSNP(folder,populations=FALSE,outgroup=FALSE,gffpath=FALSE,CHR=FALSE,
        ref.chr=FALSE,snp.window.size=FALSE,parallized=FALSE)

}

\arguments{
 
  \item{folder}{the basepath of the variant_calls}
  \item{outgroup}{vector of outgroup sequences}
  \item{populations}{list of populations}
  \item{gffpath}{the corresponding gff-file}
  \item{CHR}{which chromosome ?, default: all chromosomes}
  \item{ref.chr}{reference chromosome (to verfiy synonymous positions)}
  \item{snp.window.size}{scan SNP chunks}
  \item{parallized}{multicore computation}
}

\details{
The ff-package we use to store the SNP informations is limited by \cr
individuals * (number of SNPs) <= .Machine$integer.max \cr
The text-files containing the SNP-information of each individuell have to be stored in one folder. \cr
The slots \code{transitions}, \code{biallelic.sites} and \code{biallelic.substitutions} of the class \cr
\code{"regions.data"} will be filled. \cr
Synonymous and non-synonymous positions are not considered. Thus, the McDonald-Kreitman test will not work correctly. \cr
Use the function \code{set.synnonsyn} to set the syn/nonsyn sites.

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
          3.\tab \code{region.data}            	\tab some detail data informations\cr
          4.\tab \code{region.names}            \tab names of each region \cr 
        }


}


\examples{

# GENOME.class <- readSNP("...\SNPData")
# GENOME.class <- readSNP("...\SNPData", CHR="1")
# GENOME.class <- readSNP("...\SNPData", CHR="1", gffpath="Gff_file.gff")
# GENOME.class
# GENOME.class <- neutrality.stats(GENOME.class,FAST=TRUE)
# show the result:
# get.sum.data(GENOME.class)
# GENOME.class@region.data

}



\keyword{methods}
