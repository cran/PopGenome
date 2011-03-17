 \name{F_ST.stats.2-methods}

\docType{methods}

\alias{F_ST.stats.2,GENOME-method}
\alias{F_ST.stats.2-methods}


\title{Fixation Index (2)}

\description{
 A generic function to calculate some FST measurenments.
}

\usage{

\S4method{F_ST.stats.2}{GENOME}(object,new.populations="list",snn=TRUE)

}

 \arguments{

  \item{object}{An object of class \code{"GENOME"}}
  \item{new.populations}{list of populations. default=\code{FALSE}}
  \item{snn}{Snn statistic fron Hudson}

}


\value{

Returned value is an modified object of class \code{"GENOME"} \cr  \cr
--------------------------------------------------------- \cr
Following slots will be modified in the \code{"GENOME"}  object \cr
--------------------------------------------------------- \cr

\tabular{rlll}{
            \tab Slot                	\tab Reference   \tab Description\cr
          1.\tab \code{Hudson.Snn}	\tab [1]         \tab Snn statistic from Hudson (2000) \cr
          }

}

\references{
[1] Hudson, R. R. (2000).\emph{A new statistic for detecting genetic differentiation.} Genetics 155: 2011-2014.
}

\examples{

# GENOME.class <- readData("\home\Alignments")
# GENOME.class
# GENOME.class <- F_ST.stats.2(GENOME.class)
# GENOME.class <- F_ST.stats.2(GENOME.class,list(1:4,5:10))
# GENOME.class <- F_ST.stats.2(GENOME.class,
# list(c("seq1","seq5","seq3"),c("seq2","seq8")))
# show the result:
# GENOME.class@Hudson.Snn

}


\keyword{methods}
