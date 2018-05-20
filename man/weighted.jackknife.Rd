\name{weighted.jackknife}
\docType{methods}

\alias{weighted.jackknife,GENOME-method}
\alias{weighted.jackknife-methods}

\title{Weighted Jackknife}


\description{
 Implements the weighted drop-one jackknife procedure.   
}

\usage{

\S4method{weighted.jackknife}{GENOME}(object,
do.D=TRUE, 
do.BDF=TRUE, 
per.region=FALSE, 
block.size=1
	)
}


\arguments{

  \item{object}{an object of class \code{"GENOME"}}
  \item{do.D}{default:\code{TRUE}}
  \item{do.BDF}{default:\code{TRUE}}
  \item{per.region}{\code{FALSE} drop a genomic region, \code{TRUE} drop single SNPs}
  \item{block.size}{in case \code{per.region=TRUE} the number of SNPs to be dropped}
}

\value{
The function creates a transformed object of class \code{"GENOME"}.\\
The following slots are filled: \code{BDF.z},\code{BDF.pval},\code{BDF.SE},\code{D.z},\code{D.pval}
and \code{D.SE}.
}

\note{
The \code{introgression.stats} module has to be performed prior to the anaylses.
}


\examples{
# GENOME.class      <- readVCF("...", ...)
# GENOME.class      <- set.populations()
# GENOME.class.slide <- sliding.window.transform(...)
# GENOME.class.slide <- introgression.stats(...)
# GENOME.class.slide <- weighted.jackknife(GENOME.class.slide, per.region=FALSE)
# GENOME.class.slide@BDF.z
# GENOME.class.slide@D.z
}



\keyword{methods}
