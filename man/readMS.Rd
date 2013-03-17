\name{readMS}

\docType{methods}
\alias{readMS}

\title{Reading data from MS and MSMS output}


\description{

This function reads data produced from the coalescent simulation programs \cr
MS (Hudson, 2002 ) and MSMS (Greg, 2010)

}

\usage{

readMS(file)

}

\arguments{
  \item{file}{the basepath of the MS/MSMS output}
}

\value{

An object of class \code{GENOME}

}

\references{ 

Hudson, R. R. (2002). \emph{Generating samples under a Wright-Fisher neutral model of genetic variation.}
Bioinformatics 18: 337-338 \cr \cr 

Gregory Ewing and Joachim Hermisson, \emph{MSMS: A Coalescent Simulation Program Including Recombination, Demographic Structure, and Selection at a Single Locus.} Bioinformatics 2010, doi: 10.1093/bioinformatics/btq322 

}

\examples{

# GENOME.class <- readMS("ms.output.txt")
# GENOME.class@region.names

}

\keyword{methods}
