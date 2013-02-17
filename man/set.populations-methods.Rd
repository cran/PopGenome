\name{set.populations-methods}
\docType{methods}

\alias{set.populations,GENOME-method}
\alias{set.populations-methods}


\title{Define populations}


\description{
This generic function defines the populations. \cr
The advantage of this function is, that you don't have to specify the populations for each calculation.
The populations are set for each statistic module.  
}

\usage{

\S4method{set.populations}{GENOME}(object,new.populations=FALSE, diploid=FALSE)

}

\arguments{

  \item{object}{an object of class \code{"GENOME"}}
  \item{new.populations}{list of populations. default:\code{FALSE}}
  \item{diploid}{if diploid data is present}
}



\examples{

# GENOME.class <- readData("\home\Alignments")
# pop.1 <- c("seq1","seq2")
# pop.2 <- c("seq3","seq4","seq1")
# GENOME.class <- set.populations(GENOME.class,list(pop.1,pop.2))
# GENOME.class <- neutrality.stats(GENOME.class)
}



\keyword{methods}

