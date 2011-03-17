\name{detail.stats-methods}
\docType{methods}
\alias{detail.stats,GENOME-method}
\alias{detail.stats-methods}
\alias{get.detail,GENOME-method}
\alias{get.detail-methods}


\title{Several statistics}


\description{
 This generic function calculates some mixed statistics.
}

\usage{

\S4method{detail.stats}{GENOME}(object,new.populations=FALSE,subsites=FALSE,biallelic.structure=TRUE,mismatch.distribution=TRUE)
\S4method{get.detail}{GENOME}(object)



}

\arguments{

  \item{object}{an object of class \code{"GENOME"}}
  \item{new.populations}{list of populations. default:\code{FALSE}}
  \item{subsites}{\code{"syn"} for synonymous sites and \code{"nonsyn"} for nonsynonymous sites. default:\code{FALSE}}
  \item{biallelic.structure}{fixed and shared polymorphisms (stored in \code{GENOME.class@region.stats}). default:\code{TRUE}}
  \item{mismatch.distribution}{statistics based on mismatch distribution. default:\code{TRUE}}
  
}



\value{


returned value is an modified object of class \code{"GENOME"} \cr \cr
-------------------------------------------------------------- \cr
Following Slots will be modified in the \code{"GENOME"} object \cr
-------------------------------------------------------------- \cr

   \item{MDSD}{}
   \item{MDG1}{}
   \item{MDG2}{}
   \item{region.stats}{the slot  \code{biallelic.structure} will be filled}
   
  
  
}

\examples{

# GENOME.class <- readData("\home\Alignments")
# GENOME.class
# GENOME.class <- detail.stats(GENOME.class)
# GENOME.class <- detail.stats(GENOME.class,list(1:7,8:12))
# show the result:
# mismatch.values <- get.detail(GENOME.class)
# GENOME.class@region.stats@biallelic.structure
# GENOME.class@region.stats@biallelic.structure[[1]]

}



\keyword{methods}
