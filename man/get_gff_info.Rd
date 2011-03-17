\name{get_gff_info}

\docType{methods}
\alias{get_gff_info}


\title{Annotation info}


\description{
This function extract annotation informations from a gtf/gff file.
}

\usage{

get_gff_info(object=FALSE,gff.file,chr,position)

}

\arguments{

  \item{object}{object of class \code{GENOME}} 
  \item{gff.file}{basepath of the gtf/gff file}
  \item{chr}{the chromosome}
  \item{position}{reference position or region id (when \code{object} is specified)}

}  

\details{
coming soon !
}


\examples{

# get_gff_info("Arabidopsis.gff",chr=1,200202)
# get_gff_info(GENOME.class,"Arabidopsis.gff",chr=1,position=3)
}



\keyword{methods}
