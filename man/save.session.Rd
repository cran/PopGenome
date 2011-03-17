\name{save.session}

\docType{methods}
\alias{save.session}


\title{Saving a PopGenome-session}


\description{
This function saves the PopGenome-session in the current workspace.
Can be loaded again with load.session()
}

\usage{

save.session(object)

}

\arguments{

  \item{object}{object of class \code{"GENOME"}}

}  

\details{
Saving R and ff-objects created by the ff-package.
}


\examples{

# GENOME.class <- readData("...\Alignments")
# save.session(GENOME.class)
# load.session()
}



\keyword{methods}
