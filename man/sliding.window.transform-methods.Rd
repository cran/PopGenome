\name{sliding.window.transform-methods}
\docType{methods}

\alias{sliding.window.transform,GENOME-method}
\alias{sliding.window.transform-methods}

\title{Sliding Window Transformation.}


\description{
 This generic function transforms the existing object of class \code{"GENOME"} in another object of class \code{"GENOME"}, so that 
 sliding window calculations become possible without loosing the method-capacity of the PopGenome environment.   
 
}

\usage{

\S4method{sliding.window.transform}{GENOME}(object,width=7,jump=5,type=1,whole.data=TRUE)


}

\arguments{

  \item{object}{an object of class \code{"GENOME"}}
  \item{width}{window size. default:\code{7}}
  \item{jump}{jump size. default:\code{5}}
  \item{type}{\code{1} scan only biallelic positions (SNPs), \code{2} scan the genome. default:\code{1}}
  \item{whole.data}{scanning the whole data by concatenating the regions. If FALSE, the regions are scaned seperately}

}

\value{

The function creates an transformed object of class \code{"GENOME"} 

}

\note{

You should not use the big.data option in the readData function, if you want to scan regions seperately.

 
}


\examples{

# GENOME.class       <- readData("...\Alignments")
# slide.GENOME.class <- sliding.window.transform(GENOME.class)
# slide.GENOME.class <- sliding.window.transform(GENOME.class,100,100)
# slide.GENOME.class <- neutrality.stats(slide.GENOME.class)
# slide.GENOME.class@region.names
# values             <- get.neutrality(slide.GENOME.class)

}



\keyword{methods}
