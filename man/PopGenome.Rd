\name{PopGenome}
\alias{PopGenome}
\alias{PopGenome}
\docType{package}
\title{

 PopGenome



}
\description{

R-package for Population genetic & genomic analysis

}
\details{

Index:
\preformatted{
F_ST.stats               Fixation index
diversity.stats          Diversities
MKT                      McDonald & Kreitman test
MS                       Coalescent simulation
detail.stats             Several statistics
linkage.stats            Linkage disequilibrium
neutrality.stats         Neutrality statistics
readData                 Reading alignments and calculating summary data
readSNP                  Read data from the 1001 genome project
readVCF                  Read data from the 1000 genome project
readHapMap               Read HapMap data
sliding.window.transform Sliding window transformation.
splitting.data           Splitting data in subsites
test.params              Set parameter for coalescent simulation.

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
