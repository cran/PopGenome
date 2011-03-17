\name{MS}
\docType{methods}
\alias{MS}

\title{Coalescent Simulation}


\description{
 This function uses Hudson's MS to compare simulated data with the observed data. 
}

\usage{

MS(GENO,niter=10,thetaID="user",params=FALSE,detail=FALSE,
neutrality=FALSE,linkage=FALSE,F_ST=FALSE)

}

\arguments{

  \item{GENO}{an object of class \code{"GENOME"}}
  \item{niter}{number of samples per loci}
  \item{thetaID}{\code{"Tajima"},\code{"Watterson"} or \code{"user"}. default:\code{"user"}}
  \item{neutrality}{Calculate neutrality tests. default=\code{FALSE}}
  \item{linkage}{Calculate linkage disequilibrium. default=\code{FALSE}}
  \item{F_ST}{Calculate fixation index. default=\code{FALSE}}
  \item{params}{an object of class \code{"test.params"}. see \code{?test.params}}
  \item{detail}{detail statistics. Note:slower! default=\code{FALSE}}

}


\details{
You can choose different mutation rate estimators to generate simulation data.
When \code{thetaID="user"}, you have to define the theta values in an object of class \code{"test.params"}. 
The \code{"test.params"} class can also be used to specify some additional parameter like migration 
and/or recombination rates... (\code{?test.params}).
}


\value{


  The function creates an object of class \code{"cs.stats"} 

}


\references{ 
Hudson, R. R. (2002). Generating samples under a Wright-Fisher neutral model of genetic variation.
Bioinformatics 18: 337-338
}




\examples{

# GENOME.class <- readData("...\Alignments")
# GENOME.class <- neutrality.stats(GENOME.class,list(1:6))
# MS.class <- MS(GENOME.class,thetaID="Tajima",neutrality=TRUE)
# MS.class
# MS.class@obs.val

}

\note{
The executable file ms from Hudson have to be stored in the current workspace.
}
\keyword{methods}
