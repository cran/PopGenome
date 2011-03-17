\name{test.params-class}
\docType{class}
\alias{test.params}
\alias{test.params-class}

\title{

Set parameter for Coalescent Simulation.

}

\description{

This object can be passed to the function MS after having set parameter 
values. This class eases the process of passing on all necessary values
to the MS function.
 
}

\arguments{
 
  \item{theta}{mutation parameter theta (4Nmu), where N is the diplod population size and mu the mutation rate per locus. It needs to be provided as vector of length n.regions}  
  \item{seeds}{specify 3 random number seeds. a vector of length 3 with positive values is expected}
  \item{fixedSegsites}{usually the number of segregating sites varies in each iteration. Please provide a  single numeric value if the number of segregating sites needs to be fixed.}
  \item{recombination}{provide a vector of format: c(p, nsites), p = cross over parameter rate, nsites is the number of sites between recombination occurs}
  \item{geneConv}{in addition to recombination intra-locus non-cross-over exchange gene conversion can be included in simulation, expected format is c(f, gamma),  f denote the ratio, g/r, where r is the probability per generation  of crossing-over between adjacent sites. (see Wiuf and Hein 2000), gamma is the mean conversion tract length}
  \item{growth}{population size is measured by $N(t) = N0 exp^{alpha*t}$. provide alpha as integer value. negative values indicate that population was larger in the past than present}
  \item{migration}{specify the migration rate between popualtions. Please provide a single numeric value.}
  \item{demography}{vector of length 3 or 4 with first value denoted as 'type' \cr
					valid 'types' for vectors of length 3 are as following: \cr
					- \code{1} to set a growth rate change alpha at a certain time t:\cr 
					  		c(1, t, alpha)\cr
\cr
					- \code{2}	set all subpop to size $x * N_{0}$ and growth rate to zero:\cr 
					  		c(2, t, x)\cr
\cr
					- \code{3}      set all elements of migration matrix to $x/(npop-1)$:\cr 
					 		c(3, t, x)					\cr
\cr
							valid 'types' for vector of length 4 with the following values:\cr	
					- \code{4}      set growth rate of subpop i to alpha at time z: \cr
					 		 c(4, t, i, alpha)\cr
\cr
					- \code{5}      set subpop i size to $x * N_{0}$ at time t and growth rate to zero:\cr 
					  		c(5, t, i, x)\cr
\cr
					- \code{6}      split subpopulation i into subpopulation i and a new subpopulation,\cr 
					  		labeled npop + 1. Each ancestral lineage in subpopulation i is randomly \cr 
					  		assigned to subpopulation i with probability p and subpopulation  \cr
					 		npop + 1 with probability 1 - p. The size of subpopulation npop + 1 is\cr  
					  		set to $N_{0}$. Migration rates to and from the new subpopulation are assumed\cr  
					  		to be zero and the growth rate of the new subpopulation is set to zero: \cr
					 		c(6, t, i, p)\cr
\cr
					- \code{7}      move all lineages in subpopulation i to subpopulation j at time t.\cr 
					  		Migration rates from subpopulation i are set to zero: \cr
					  		c(7, t, i, j)\cr
  }

}


\author{
Bastian Pfeifer
}

\seealso{
\link{MS}
}
\examples{

# params           <- new("test.params")
# params@theta     <- rep(5,n.regions)
# params@migration <- 3

}

\keyword{classes}

