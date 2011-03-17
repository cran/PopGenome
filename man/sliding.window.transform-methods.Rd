\name{sliding.window.transform-methods}
\docType{methods}

\alias{sliding.window.transform,GENOME-method}
\alias{sliding.window.transform-methods}

\title{Sliding Window Transformation.}


\description{
 This generic function transforms the existing object of class \code{"GENOME"} in another object of class \code{"GENOME"}, so that 
 sliding window calculations are possible.
}

\usage{

\S4method{sliding.window.transform}{GENOME}(object,width=7,jump=5,type=1)


}

\arguments{

  \item{object}{an object of class \code{"GENOME"}}
  \item{width}{window size. default:\code{7}}
  \item{jump}{jump size. default:\code{5}}
  \item{type}{\code{1} scan only biallelic positions, \code{2} scan the whole alignment. default:\code{1}}

}

\value{

The function creates an object of class \code{"GENOME"} 

}

\note{
In this version of PopGenome Sliding Window applications will scan alignments separately.
What can be the file set to be analyzed depends on how much workspace is available. 
In a future version of PopGenome this problem is circumvented by moving data swapped.
}


\examples{

# GENOME.class  <- readData("...\Alignments")
# slide.GENOME.class <- sliding.window.transform(slide.GENOME.class)
# slide.GENOME.class <- neutrality.stats(slide.GENOME.class)
# values <- get.neutrality(slide.GENOME.class)

}



\keyword{methods}
