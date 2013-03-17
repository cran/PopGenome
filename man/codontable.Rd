\name{codontable}

\docType{methods}
\alias{codontable}


\title{Prints the Codontable which is used in the PopGenome framework}


\description{
This functions prints the Nucleotide Triplets as numerical values and 
the correspnding protein character strings.
}

\usage{

codontable()

}

\arguments{
no arguments
}  

\details{
The returned value is a list including two matrices. \cr
The first matrix contains the proteins and the second matrix the \cr
corresponding nucleotide triplets. In the PopGenome Vignette you can see how \cr
to manipulate these tables to use your own defined genetic code.
}


\examples{

# table <- codontable()
# table$Proteins 
# table$Triplets
}



\keyword{methods}
