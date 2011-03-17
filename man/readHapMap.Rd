\name{readHapMap}

\docType{methods}
\alias{readHapMap}

\title{Reading SNP data from the HapMap consortium}


\description{
 This function reads HapMap data. 
}

\usage{

readHapMap(folder,hap_gffpath,populations=FALSE,outgroup=FALSE)

}

\arguments{
 
  \item{folder}{the basepath of the variant_calls}
  \item{hap_gffpath}{the basepath of the corresponding gff-files. Note! This gff_file does not contain information about subsites. see    details!}
  \item{populations}{list of populations}
  \item{outgroup}{vector of outgroup sequences}
  
}

\details{
PopGenome reads the gff file distributed on the HapMap plattform only to verify the reference positions of the chromosomes.
In the next release this function will also handle GFF/GTF files to get informations about subsites like exons and so on.
The folder includes the files of different individuals of one chromosome. This enables also FST calculations of the HapMap data.
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
        }


}


\examples{

# GENOME.class <- readHapMap("...\HapMapData")
# GENOME.class
# show the result:
# get.sum.data(GENOME.class)
# GENOME.class@region.data

}



\keyword{methods}
