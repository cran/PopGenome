\name{concatenate_to_whole_genome}

\docType{methods}
\alias{concatenate_to_whole_genome}

\title{Concatenate}


\description{
 This function concatenates the regions/chunks. 
}

\usage{

concatenate_to_whole_genome(obj, n.chunks)

}

\arguments{
 
  \item{obj}{object of class \code{GENOME}}
  \item{n.chunks}{number of regions/chunks}
  
}


\value{
  The function creates an object of class "GENOME"  
}


\examples{

# GENOME.class <- readData("Three_Alignments/")
# WHOLE        <- concatenate_to_whole_genome(GENOME.class,3)
# slide        <- sliding.window.transform(WHOLE,10,10)
# slide        <- neutrality.stats(slide)
# slide@Tajima.D

}

\keyword{methods}
