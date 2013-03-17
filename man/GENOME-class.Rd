\name{GENOME-class}
\docType{class}
\alias{GENOME-class}

\alias{neutrality.stats}
\alias{linkage.stats}
\alias{calc.R2}
\alias{sliding.window.transform}
\alias{F_ST.stats}
\alias{F_ST.stats.2}
\alias{detail.stats}
\alias{diversity.stats}
\alias{MKT}
\alias{recomb.stats}
\alias{sweeps.stats}
\alias{Achaz.stats}
\alias{set.populations}
\alias{set.outgroup}
\alias{get.individuals}
\alias{show.slots}
\alias{get.status}
\alias{get.neutrality}
\alias{get.F_ST}
\alias{get.diversity}
\alias{get.linkage}
\alias{get.detail}
\alias{get.MKT}
\alias{get.sum.data}
\alias{get.recomb}
\alias{get.sweeps}
\alias{set.ref.positions}
\alias{set.synnonsyn}
\alias{splitting.data}
\alias{mult.linkage.stats}
\alias{getBayes}
\alias{get.biallelic.matrix}
\alias{region.as.fasta}
\alias{get.codons}
\alias{PG_plot.biallelic.matrix}


\alias{getMS,GENOME-method}

\alias{popFSTN,GENOME-method}
\alias{show,GENOME-method}
\alias{usage,GENOME-method}

\title{Class "GENOME"}

\description{A class where all values are stored}

\section{Slots}{
    \describe{
    \item{\code{BIG.BIAL}:}{Biallelic matrix as an ff-object}
    \item{\code{SLIDE.POS}:}{Positions of biallelic sites (Sliding window mode)}
    \item{\code{big.data}:}{ff-package ?}
    \item{\code{gff.info}:}{Gff information ?}
    \item{\code{snp.data}:}{SNP data ? }
    \item{\code{basepath}:}{The basepath of the data}
    \item{\code{project}:}{----}
    \item{\code{populations}:}{Populations definded before reading data}
    \item{\code{poppairs}:}{---}
    \item{\code{outgroup}:}{A vector of outgroup sequences}
    \item{\code{region.names}:}{Names/identifier of each region}
    \item{\code{genelength}:}{Number of regions }
    \item{\code{keep.start.pos}:}{position of start for sliding window}
    \item{\code{n.sites}:}{Total number of sites calculated}
    \item{\code{n.sites2}:}{Total number of sites calculated}
    \item{\code{n.biallelic.sites}:}{Number of biallelic sites (SNPs) }
    \item{\code{n.gaps}:}{Number of gaps observed in the data}
    \item{\code{n.unknowns}:}{Number of unknown.positions }
    \item{\code{n.valid.sites}:}{Sites without gaps}
    \item{\code{n.polyallelic.sites}:}{Sites with more than two variants}
    \item{\code{trans.transv.ratio}:}{Transition-transversion ratio}
    \item{\code{Coding.region}:}{Number of nucleotides in the CDS region }
    \item{\code{UTR.region}:}{Number of nucleotides in the UTR region }
    \item{\code{Intron.region}:}{Number of nucleotides in Intron region }
    \item{\code{Exon.region}:}{Number of nucleotides in Exon region }
    \item{\code{Gene.region}:}{Number of nucleotides in Gene region }
    \item{\code{Pop_Neutrality}:}{Populations defined in the neutrality module}
    \item{\code{Pop_FSTN}:}{Populations defined in the FST (nucleotide) module}
    \item{\code{Pop_FSTH}:}{Populations defined in the FST (haplotype) module}
    \item{\code{Pop_Linkage}:}{Populations defined in the Linkage module}
    \item{\code{Pop_Slide}:}{---}
    \item{\code{Pop_MK}:}{Populations defined in the MK module}
    \item{\code{Pop_Detail}:}{Populations defined in the Detail module}
    \item{\code{Pop_Recomb}:}{Populations defined in the Recombination module}
    \item{\code{Pop_Sweeps}:}{Populations defined in the Selective sweeps module}
    \item{\code{FSTNLISTE}:}{---}
    \item{\code{nucleotide.F_ST}:}{Nucleotide FST }
    \item{\code{nucleotide.F_ST2}:}{---}
    \item{\code{nuc.diversity.between}:}{Nucleotide diversity between the populations}
    \item{\code{nuc.diversity.within}:}{Nucleotide diversity within the populations}
    \item{\code{nuc.F_ST.pairwise}:}{For each pair of populations}
    \item{\code{nuc.F_ST.vs.all}:}{FST for one population vs the rest}
    \item{\code{n.haplotypes}:}{---}
    \item{\code{hap.diversity.within}:}{Haplotype diversity withing the popualtions}
    \item{\code{hap.diversity.between}:}{Haplotype diversity between the populations}
    \item{\code{Pi}:}{Pi from Nei}
    \item{\code{PIA_nei}:}{Pi between the poulations}
    \item{\code{haplotype.counts}:}{Counts of the haplotypes observed}
    \item{\code{haplotype.F_ST}:}{Haplotype FST}
    \item{\code{hap.F_ST.pairwise}:}{Haplotype diversity of each pair of populations}
    \item{\code{Nei.G_ST.pairwise}:}{Haplotype diversity of each pair of populations}
    \item{\code{hap.F_ST.vs.all}:}{FST for one population vs the rest}
    \item{\code{Nei.G_ST}:}{GST from Nei}
    \item{\code{Hudson.G_ST}:}{GST from Hudson}
    \item{\code{Hudson.H_ST}:}{HST from Hudson}
    \item{\code{Hudson.K_ST}:}{KST from Hudson}
    \item{\code{Hudson.Snn}:}{Snn from Hudson}
    \item{\code{Phi_ST}:}{Fixation index from Excoffier}
    \item{\code{hap.pair.F_ST}:}{---}
    \item{\code{MKT}:}{Mcdonald-Kreitman values}
    \item{\code{Tajima.D}:}{Tajima D statistic}
    \item{\code{SLIDE}:}{---}
    \item{\code{Fay.Wu.H}:}{}
    \item{\code{Zeng.E}:}{}
    \item{\code{theta_Tajima}:}{}
    \item{\code{theta_Watterson}:}{}
    \item{\code{theta_Fu.Li}:}{}
    \item{\code{theta_Achaz.Watterson}:}{}
    \item{\code{theta_Achaz.Tajima}:}{}
    \item{\code{theta_Fay.Wu}:}{}
    \item{\code{theta_Zeng}:}{}
    \item{\code{Fu.Li.F}:}{}
    \item{\code{Fu.Li.D}:}{}
    \item{\code{Yach}:}{}
    \item{\code{n.segregating.sites}:}{total number of segregating sites}
    \item{\code{Rozas.R_2}:}{}
    \item{\code{Fu.F_S}:}{}
    \item{\code{Strobeck.S}:}{}
    \item{\code{Kelly.Z_nS}:}{}
    \item{\code{Rozas.ZZ}:}{}
    \item{\code{Rozas.ZA}:}{}
    \item{\code{Wall.B}:}{}
    \item{\code{Wall.Q}:}{}
    \item{\code{mult.Linkage}:}{Linkage disequilibrium between regions}
    \item{\code{RM}:}{Minimum number of recombination events (Hudson)}
    \item{\code{CL}:}{Composite Likelihood of SNPs (Nielsen et. al)}
    \item{\code{CLmax}:}{Max. Composite Likelihood of SNPs (Nielsen et.al)}
    \item{\code{CLR}:}{Composite Likelihood ratio test (Nielsen et. al)}
    \item{\code{MDSD}:}{}
    \item{\code{MDG1}:}{}
    \item{\code{MDG2}:}{}
    \item{\code{genes}:}{}
    \item{\code{region.data}:}{Detail information about the data}
    \item{\code{region.stats}:}{Detail (site-specific) statistics}
  }
}
\section{Methods}{
  \describe{
    \item{detail.stats}{Several statistics}
    \item{diversity.stats}{Haplotype and Nucleotide Diversities}
    \item{F_ST.stats.2}{Snn Hudson }
    \item{F_ST.stats}{Fixation index }
    \item{getBayes}{Get the input for BayeScanR}
    \item{get.detail}{Get the results from the Detail module}
    \item{get.codons}{Get informations about the codon-changes}
    \item{get.diversity}{Get diversities from the FST module}
    \item{get.F_ST}{Get FST values from the FST module}
    \item{get.linkage}{Get the values from the Linkage module}
    \item{get.MKT}{Mcdonald-Kreitman values}
    \item{getMS}{---}
    \item{get.neutrality}{Get the statistic-values from the Neutrality module}
    \item{get.status}{Status of calculation}
    \item{get.sum.data}{Get some data observed from the alignments}
    \item{linkage.stats}{Linkage disequilibrium}
    \item{calc.R2}{Linkage disequilibrium}
    \item{mult.linkage.stats}{Linkage disequilibrium between regions }
    \item{recomb.stats}{Recombination statistics }
    \item{sweeps.stats}{Selective sweeps}
    \item{Achaz.stats}{Achaz statistics}
    \item{get.recomb}{Get the values from the Recombination module}
    \item{get.sweeps}{Get the values frome the Selective Sweep module}
    \item{set.ref.positions}{Set the SNP positions}
    \item{set.synnonsyn}{Verifies synonymous positions}
    \item{splitting.data}{Splits the data in subsites}
    \item{MKT}{MKT Test}
    \item{neutrality.stats}{Neutrality Statistics }
    \item{popFSTN}{Internal function}
    \item{get.biallelic.matrix}{Prints the biallelic.matrix}
    \item{set.populations}{Defines the populations}
    \item{set.outgroup}{Defines the outgroup}
    \item{get.individuals}{Individual names}
    \item{region.as.fasta}{Extract the region as a fasta file}
    \item{show}{---}
    \item{show.slots}{Show slots of the class \code{GENOME}}
    \item{sliding.window.transform}{Sliding window analysis}
    \item{usage}{---}
    \item{PG_plot.biallelic.matrix}{Plots the Biallelic Matrix}
  }
}

\references{see documentation for each module}
\author{Bastian Pfeifer}

\examples{
#GENOME.class <- readData("Alignments")
#GENOME.class@n.sites
#GENOME.class@region.names
}

\keyword{classes}
