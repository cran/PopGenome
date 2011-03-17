\name{mult.linkage.stats-methods}
\docType{methods}
\alias{mult.linkage.stats,GENOME-method}
\alias{mult.linkage.stats-methods}


\title{Multilocus linkage statistics}


\description{
 This generic function calculates the linkage disequilibrium between regions.
}

\usage{

\S4method{mult.linkage.stats}{GENOME}(object,new.populations="list",lower.bound=0,upper.bound=1)

}

\arguments{

  \item{object}{an object of class \code{"GENOME"}}
  \item{new.populations}{list of populations. default:\code{FALSE}}
  \item{lower.bound}{sites with minor-allele-frequency=lower.bound are considered}
  \item{upper.bound}{sites with minor-allele-frequency=upper.bound are considered}
}



\value{


returned value is an modified object of class \code{"GENOME"} \cr \cr
-------------------------------------------------------------- \cr
Following Slots will be modified in the \code{"GENOME"} object \cr
-------------------------------------------------------------- \cr

   \item{mult.Linkage}{Some values about linkage for each pair of site (R2...)}
  
}

\examples{

# GENOME.class <- readData("...\Alignments")
# GENOME.class
# GENOME.class <- mult.linkage.stats(GENOME.class)
# show the result:
# GENOME.class@mult.Linkage
 

}



\keyword{methods}
