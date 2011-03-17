\name{PopGenome}
\alias{PopGenome}
\alias{PopGenome}
\docType{package}
\title{

 PopGenome



}
\description{

R-package for Population genetic analysis

}
\details{
\tabular{ll}{
Package: \tab PopGenome\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-03-17\cr
Depends: \tab R (>= 1.8.0)\cr
License: \tab What license is it under?\cr
LazyLoad: \tab yes\cr
Packaged: \tab Tue Sep 20 15:33:01 2011; bapfe\cr
Built: \tab R 2.8.1; $x86_{64}$-unknown-linux-gnu; 2011-09-20 15:33:04; unix\cr
}

Index:
\preformatted{
F_ST.stats              Fixation index
MKT                     McDonald & Kreitman test
MS                      Coalescent simulation
detail.stats            Several statistics
linkage.stats           Linkage disequilibrium
neutrality.stats        Neutrality statistics
readData                Reading alignments and calculating summary data
sliding.window.transform
                        Sliding window transformation.
test.params             Set parameter for coalescent simulation.
}

}

\author{
Bastian Pfeifer
Maintainer: Bastian Pfeifer <Bastian.Pfeifer@uni-duesseldorf.de>

}

\keyword{ package }
\seealso{

?readData
\link{readData}
}
\examples{

# GENOME.class <- readData("...\Alignments")
# GENOME.class <- neutrality.stats(GENOME.class)
# values       <- get.neutrality(GENOME.class)
# GENOME.class <- F_ST.stats(GENOME.class,list(1:5,6:10))
# values       <- get.F_ST(GENOME.class)

}
