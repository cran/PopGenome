\name{popGetBial-methods}

\docType{methods}

\alias{popGetBial,GENOME-method}
\alias{popGetBial-methods}


\title{Get the biallelic matrix}

\description{
This function returns the biallelic matrix of a specific region.
}

\usage{

\S4method{popGetBial}{GENOME}(object,bialmatNr)

}

 \arguments{

  \item{object}{An object of class \code{"GENOME"}}
  \item{bialmatNr}{id of the region}

}


\value{

Biallelic matrix \cr
rows:    names of each individuell \cr
columns: biallelic sites \cr

}


\examples{

# GENOME.class <- readData("\home\Alignments")
# popGetBial(GENOME.class,7) # biallelic matrix of the 7th alignment
}


\keyword{methods}
