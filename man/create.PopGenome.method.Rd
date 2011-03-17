\name{create.PopGenome.method}

\docType{methods}
\alias{create.PopGenome.method}


\title{Integration of own functions in the PopGenome-framework}


\description{
This function generates a skeleton for a PopGenome function.
Ensures an effortless integration of new methods in the PopGenome-framework.  
}

\usage{

create.PopGenome.method(function.name,population.specific=TRUE)

}

\arguments{

  \item{function.name}{name of your function}
  \item{population.specific}{\code{TRUE}:function returns one value per population.\code{FALSE}:function returns one value for all populations. for example FST measurements}

}  

\details{
This mechanism enables you to use your own functions in the PopGenome environment. Sliding window analysis and 
calculation of subsites will work well.\cr
Please look at the generated function, which documents where to place your own function in detail. 
}


\examples{

# GENOME.class <- readData(".../Alignments")
# create.PopGenome.method("myFunction")
# edit myFunction.R
# source("myFunction")
# value <- myFunction(test)
# value

}



\keyword{methods}
