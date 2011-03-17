\name{GENOME-class}
\docType{class}
\alias{GENOME-class}

\alias{neutrality.stats}
\alias{linkage.stats}
\alias{sliding.window.transform}
\alias{F_ST.stats}
\alias{detail.stats}
\alias{MKT}
\alias{set.populations}
\alias{show.slots}
\alias{get.status}
\alias{get.neutrality}
\alias{get.F_ST}
\alias{get.diversity}
\alias{get.linkage}
\alias{get.detail}
\alias{get.MKT}
\alias{get.sum.data}

\alias{F_ST.stats.2,GENOME-method}
\alias{getBayes,GENOME-method}

\alias{getMS,GENOME-method}

\alias{popFSTN,GENOME-method}
\alias{show,GENOME-method}
\alias{usage,GENOME-method}

\title{Class "GENOME"}
\description{A class where all values are stored}

\section{Objects from the Class}{

coming soon ...

}
\section{Slots}{
	 \describe{
    \item{\code{basepath}:}{Object of class \code{"character"} ~~ }
    \item{\code{project}:}{Object of class \code{"character"} ~~ }
    \item{\code{populations}:}{Object of class \code{"list"} ~~ }
    \item{\code{poppairs}:}{Object of class \code{"vector"} ~~ }
    \item{\code{outgroup}:}{Object of class \code{"vector"} ~~ }
    \item{\code{region.names}:}{Object of class \code{"character"} ~~ }
    \item{\code{genelength}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{n.sites}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{n.biallelic.sites}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{n.gaps}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{n.unknowns}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{n.valid.sites}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{n.polyallelic.sites}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{trans.transv.ratio}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{Pop_Neutrality}:}{Object of class \code{"list"} ~~ }
    \item{\code{Pop_FSTN}:}{Object of class \code{"list"} ~~ }
    \item{\code{Pop_FSTH}:}{Object of class \code{"list"} ~~ }
    \item{\code{Pop_Linkage}:}{Object of class \code{"list"} ~~ }
    \item{\code{Pop_Slide}:}{Object of class \code{"list"} ~~ }
    \item{\code{Pop_MK}:}{Object of class \code{"list"} ~~ }
    \item{\code{Pop_Detail}:}{Object of class \code{"list"} ~~ }
    \item{\code{FSTNLISTE}:}{Object of class \code{"list"} ~~ }
    \item{\code{nucleotide.F_ST}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{nucleotide.F_ST2}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{nuc.diversity.between}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{nuc.diversity.within}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{nuc.F_ST.pairwise}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{nuc.F_ST.vs.all}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{n.haplotypes}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{hap.diversity.within}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{hap.diversity.between}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{Pi}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{PIA_nei}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{haplotype.counts}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{haplotype.F_ST}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{hap.F_ST.pairwise}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{Nei.G_ST.pairwise}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{hap.F_ST.vs.all}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{Nei.G_ST}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{Hudson.G_ST}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{Hudson.H_ST}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{Hudson.K_ST}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{Hudson.Snn}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{hap.pair.F_ST}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{MKT}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{Tajima.D}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{SLIDE}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{Fay.Wu.H}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{Zeng.E}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{theta_Tajima}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{theta_Watterson}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{theta_Fu.Li}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{theta_Achaz.Watterson}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{theta_Achaz.Tajima}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{theta_Fay.Wu}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{theta_Zeng}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{Fu.Li.F}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{Fu.Li.D}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{Yach}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{n.segregating.sites}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{Rozas.R_2}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{Fu.F_S}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{Strobeck.S}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{Kelly.Z_nS}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{Rozas.ZZ}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{Rozas.ZA}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{Wall.B}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{Wall.Q}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{MDSD}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{MDG1}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{MDG2}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{genes}:}{Object of class \code{"list"} ~~ }
    \item{\code{region.data}:}{Object of class \code{"region.data"} ~~ }
    \item{\code{region.stats}:}{Object of class \code{"region.stats"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{detail.stats}{\code{signature(object = "GENOME")}: Several Statistics}
    \item{F_ST.stats.2}{\code{signature(object = "GENOME")}: ... }
    \item{F_ST.stats}{\code{signature(object = "GENOME")}: Fixation Index }
    \item{getBayes}{\code{signature(object = "GENOME")}: ... }
    \item{get.detail}{\code{signature(object = "GENOME")}: ... }
    \item{get.diversity}{\code{signature(object = "GENOME")}: ... }
    \item{get.F_ST}{\code{signature(object = "GENOME")}: ... }
    \item{get.linkage}{\code{signature(object = "GENOME")}: ... }
    \item{get.MKT}{\code{signature(object = "GENOME")}: ... }
    \item{getMS}{\code{signature(object = "GENOME")}: ... }
    \item{get.neutrality}{\code{signature(object = "GENOME")}: ... }
    \item{get.status}{\code{signature(object = "GENOME")}: ... }
    \item{get.sum.data}{\code{signature(object = "GENOME")}: ... }
    \item{linkage.stats}{\code{signature(object = "GENOME")}: Linkage Disequilibrium }
    \item{MKT}{\code{signature(object = "GENOME")}: MKT Test}
    \item{neutrality.stats}{\code{signature(object = "GENOME")}: Neutrality Statistics }
    \item{popFSTN}{\code{signature(object = "GENOME")}: ... }
    \item{set.populations}{\code{signature(object = "GENOME")}: Define the populations}
    \item{show}{\code{signature(object = "GENOME")}: ... }
    \item{show.slots}{\code{signature(object = "GENOME")}: ... }
    \item{sliding.window.transform}{\code{signature(object = "GENOME")}: Sliding Window }
    \item{usage}{\code{signature(object = "GENOME")}: ... }
	 }
}
\references{ ~put references to the literature/web site here ~ }
\author{Bastian Pfeifer}
\note{ test test}



\seealso{
	bla bla
        
}
\examples{
showClass("GENOME")
}
\keyword{classes}
