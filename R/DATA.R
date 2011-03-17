setClass("region.data", representation(

populations      = "list",
populations2     = "list",
popmissing       = "list", 
outgroup         = "list",
outgroup2        = "list",

CodingSNPS       = "list",
UTRSNPS          = "list",
IntronSNPS       = "list",

transitions      = "list",  # matrix_sv  transition war eine 1
biallelic.matrix = "list",  # matrix_pol
biallelic.sites  = "list", # matrix_pos
matrix_codonpos  = "list", # codonpos. of biallelics
synonymous       = "list", # synnonsyn
matrix_freq      = "list",
n.singletons     = "list", # unic
trans.transv.ratio     = "list",
n.valid.sites        = "list", # algsites
n.biallelic.sites       = "list",
polyallelic.sites = "list", # mhitbp
n.nucleotides    = "list", # sum_sam
biallelic.compositions  = "list", # TCGA
biallelic.substitutions = "list", # subst
minor.alleles    = "list", # mutations
codons           = "list",
sites.with.gaps  = "list", # gaps
sites.with.unknowns = "list"
 
))

#### SHOW ######
setMethod("show", "region.data",
 function(object){
 # print(summary(object))
 cat("-----\n")
 cat("SLOTS:\n")
 cat("-----\n")
 out <- data.frame(Slots=c("populations","outgroup","transitions","biallelic.matrix","n.singletons",
                   "biallelic.sites","n.nucleotides","biallelic.compositions","synonymous","biallelic.substitutions","polyallelic.sites","sites.with.gaps","sites.with.unknowns","minor.alleles","codons"),         
               

                Description=c("samples of each population","samples of outgroup",
               "biallelic site transitions","biallelic matrix","number of singletons",
               "position of biallelic sites","number of nucleotides per sequence","nucleotides per sequence (biallelic)",
                "synonymous biallelic sites","biallelic substitutions","sites with >2 nucleotides","sites with gap positions","sites with unknown positions","minor alleles","codons of biallelic substitutions") )
 print(out)
 cat("\n---------------\n")
 print("These are the Slots (class region.data)")
  
})

#setMethod("summary","DATA",function(object){
#cat("Some basic informations: \n")
#cat("---------------------------- \n")
#cat("Path of Gene: ",object@genename,"\n")
#cat("unic:       number of single mutations:",object@unic,"\n")
#cat("totalmis:   number of missing positions",object@totalmis,"\n")
#cat("bial_sites: number of biallelic positions",object@bial_sites,"\n")
#})
