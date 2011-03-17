
readData <- function(path,populations=FALSE,outgroup=FALSE,include.unknown=TRUE,gffpath=FALSE){

methods <- "DATA"


########################
 if (.Platform$OS.type == "unix") {
 path_C_code      <- file.path(.path.package("PopGenome"),"libs","PopGenome.so")
 }else{
 ## muss je nach 32 64 Win geändert werden !
 path_C_code      <- file.path(.path.package("PopGenome"),"libs","PopGenome.dll")
 }
 dyn.load(path_C_code,PACKAGE="PopGenome")
########################
  
 
 npops    <- length(populations)
 popnames <- paste("pop",1:npops)
 
 # Get the alignments
 liste    <- list.files(path,full.names=TRUE)
 liste2   <- list.files(path)
 liste3   <- gsub(".fas","",liste2)
 

 gff_objects  <- vector("list",length(liste))
 ## GET THE GFF FILES -----------------------------
 if(gffpath[1]!=FALSE){ 
 
    
    #gff_objects  <- lapply(gff_objects,function(x){return(FALSE)})
    gff_liste    <- list.files(gffpath,full.names=TRUE)
    gff_liste2   <- list.files(gffpath)
    gff_liste3   <- gsub(".gff","",gff_liste2)             
    treffer      <- match(gff_liste3,liste3)
    
    # gff_objects  <- lapply(gff_liste,function(x){read.table(x,fill=TRUE)})
    for(xx in 1:length(treffer)){
      gff_objects[[treffer[xx]]] <- read.table(gff_liste[xx],fill=TRUE)
    }
    
 }else{}
 

 # ---------------------------------------------



  
# --- Uli
# liste    <- list.files(path,full.names=T)				        # @uw: modified to exclude subdirectories
# liste    <- liste[ as.matrix(file.info(liste)$isdir)==F ]		        # subdirectory called RData/ is created to store the codematrices 
#-- Uli



################# Get the poppairs names ##################################
if(npops>1){
 if(outgroup[1]!=FALSE){
  poppairs <- choose(npops+1,2) # Outgroup is included !!
  pairs    <- combn(1:(npops+1),2)
 }else{
  poppairs <- choose(npops,2)   # Outgroup is not included !!
  pairs    <- combn(1:(npops),2)
 } 
###########################################################################

 
#### --- Names of population pairs --- ####################################### 
 nn <- paste("pop",pairs[1,1],"/pop",pairs[2,1],sep="")
 if(dim(pairs)[2]>1){ # more than 2 Populations
  for(xx in 2:dim(pairs)[2]){
    m  <- paste("pop",pairs[1,xx],"/pop",pairs[2,xx],sep="")
    nn <- c(nn,m)
  } 
 }#END if
}# End npops > 1
else{poppairs <- 1;nn <- "pop1"} 
##### ------------------------------ ######################################### 


 
 ### --- ------------------------------ #### 
 sizeliste              <- length(liste)
 genome                 <- new("GENOME")  
 genome@basepath        <- file.path(path)  
 genome@project         <- file.path(path) 
 genome@genelength      <- sizeliste
 #populationsX           <- as.matrix(populations)
 #rownames(populationsX) <- popnames
 #colnames(populationsX) <- "Number of Samples"
 if(!is.list(populations)){genome@populations <- list(NULL)}else{genome@populations <- populations}
 genome@poppairs        <- nn
 genome@outgroup        <- outgroup
 genome@region.names    <- liste2
 DATABOOL               <- is.element("DATA",methods)

 nsites                 <- integer(sizeliste)

 # if(ALL_METHODS){
 # genelist  <- vector("list",length(liste))      # list of objects (class GEN )
 # datalist  <- vector("list",length(liste))    # list of objects (class DATA)
 # } 

 # INIT
 # do this calculation every time
 region.data  <- new("region.data")
 region.stats <- new("region.stats")
 init         <- vector("list",sizeliste)
 init2        <- numeric(sizeliste)
 init3        <- rep(NaN,sizeliste)

      # region.data init
      populationsX      <- init
      populations2      <- init
      popmissing        <- init 
      outgroupX         <- init
      outgroup2         <- init

      CodingSNPS        <- init
      UTRSNPS           <- init
      IntronSNPS        <- init

      transitions         <- init  # matrix_sv  transition war eine 1
      biallelic.matrix    <- init  # matrix_pol
      biallelic.sites     <- init # matrix_pos
      matrix_codonpos     <- init # codonpos. of biallelics
      synonymous          <- init # synnonsyn
      matrix_freq         <- init
      n.singletons        <- init # unic
      polyallelic.sites   <- init # mhitbp
      n.nucleotides       <- init # sum_sam
      biallelic.compositions  <- init # TCGA
      biallelic.substitutions <- init # subst
      minor.alleles       <- init # mutations
      codons              <- init
      sites.with.gaps     <- init # gaps
      sites.with.unknowns <- init
      
      # GENOME data  init
      n.valid.sites       <- init2
      n.gaps              <- init2
      n.unknowns          <- init2
      n.polyallelic.sites <- init2
      n.biallelic.sites   <- init2
      trans.transv.ratio  <- init3


## PROGRESS #########################
 progr <- progressBar()
#####################################

 # -----------------------------#
 #if(!is.list(populations)){
 #allseq <- T}else{allseq <- F}
 # -----------------------------#

 for(xx in 1:sizeliste){ # 
    
     
    gen <- try(PopGenread(liste[xx]),silent=TRUE)  
     
      # gen <-  .Call("readdna",liste[xx])
      # print("######## -- ALIGNMENT -- ##########################")
      # print(liste[xx])
      # print("###################################################")
    
    if(is.matrix(gen)){
    
         #if(allseq){
         #  populations <- list(1:dim(gen)[1]) <----- schon in popgen
         #}
        
    nsites[xx] <- dim(gen)[2] 
    
    #------------------------------------------------------#
    result    <- popgen(gen,Populations=populations,outgroup=outgroup,methods=methods,include.unknown=include.unknown,gff=gff_objects[[xx]])
    #------------------------------------------------------#
    rm(gen) # delete gen
    
    }else{result <- NA;next}

# PROGRESS #######################################################
    progr <- progressBar(xx,sizeliste, progr)
    if(xx==ceiling(sizeliste/2)){gc()}
###################################################################
 
    if(is.list(result)){
      
      # fill region.data
      populationsX[[xx]]        <- result$populations
      populations2[[xx]]        <- result$populations2
      popmissing[[xx]]          <- result$popmissing
      outgroupX[[xx]]           <- result$outgroup
      outgroup2[[xx]]           <- result$outgroup2

      datt                      <- result$data.sum
      CodingSNPS[[xx]]          <- datt$CodingSNPS
      UTRSNPS[[xx]]             <- datt$UTRSNPs
      IntronSNPS[[xx]]          <- datt$IntronSNPs

      transitions[[xx]]         <- datt$transitions       # matrix_sv  transition war eine 1
      biallelic.matrix[[xx]]    <- datt$biallelic.matrix  # matrix_pol
      biallelic.sites[[xx]]     <- datt$biallelic.sites   # matrix_pos
      matrix_codonpos[[xx]]     <- datt$matrix_codonpos   # codonpos. of biallelics
      synonymous[[xx]]          <- datt$synonymous        # synnonsyn
      matrix_freq[[xx]]         <- datt$matrix_freq
      n.singletons[[xx]]        <- datt$n.singletons      # unic
      polyallelic.sites[[xx]]   <- datt$polyallelic.sites # mhitbp
      n.nucleotides[[xx]]       <- datt$n.nucleotides     # sum_sam
      biallelic.compositions[[xx]]  <- datt$biallelic.compositions # TCGA
      biallelic.substitutions[[xx]] <- datt$biallelic.substitutions# subst
      minor.alleles[[xx]]       <- datt$minor.alleles # mutations
      codons[[xx]]              <- datt$codons
      sites.with.gaps[[xx]]     <- datt$sites.with.gaps
      sites.with.unknowns[[xx]] <- datt$sites.with.unknowns
      
      # fill Genome data
      n.valid.sites[xx]       <- datt$n.valid.sites
      n.gaps[xx]              <- length(datt$sites.with.gaps)
      n.unknowns[xx]          <- length(datt$sites.with.unknowns)
      n.polyallelic.sites[xx] <- length(datt$polyallelic.sites)
      n.biallelic.sites[xx]   <- length(datt$biallelic.sites)
      trans.transv.ratio[xx]  <- datt$trans.transv.ratio
    
    }# else{warnings("No biallelic position !")}

        
 }# End of For  
 
      # region.data
      region.data@populations      <- populationsX
      region.data@populations2     <- populations2
      region.data@popmissing       <- popmissing
      region.data@outgroup         <- outgroupX
      region.data@outgroup2        <- outgroup2

      region.data@CodingSNPS       <- CodingSNPS
      region.data@UTRSNPS          <- UTRSNPS
      region.data@IntronSNPS       <- IntronSNPS

      region.data@transitions      <- transitions  # matrix_sv  transition war eine 1
      region.data@biallelic.matrix <- biallelic.matrix# matrix_pol
      region.data@biallelic.sites  <- biallelic.sites # matrix_pos
      region.data@matrix_codonpos  <- matrix_codonpos # codonpos. of biallelics
      region.data@synonymous       <- synonymous# synnonsyn
      region.data@matrix_freq      <- matrix_freq
      region.data@n.singletons     <- n.singletons # unic
      region.data@polyallelic.sites <- polyallelic.sites # mhitbp
      region.data@n.nucleotides    <- n.nucleotides # sum_sam
      region.data@biallelic.compositions  <- biallelic.compositions  # TCGA
      region.data@biallelic.substitutions <- biallelic.substitutions # subst
      region.data@minor.alleles    <- minor.alleles # mutations
      region.data@codons           <- codons
      region.data@sites.with.gaps  <- sites.with.gaps # gaps
      region.data@sites.with.unknowns <- sites.with.unknowns

      region.stats@nucleotide.diversity   <- init
      region.stats@haplotype.diversity    <- init
      region.stats@haplotype.counts       <- init       # sfreqh
      region.stats@minor.allele.freqs     <- init       # JFD
      region.stats@biallelic.structure    <- init       # SXX
      region.stats@linkage.disequilibrium <- init

      # GENOME data
      names(nsites)              <- liste2
      genome@n.sites             <- nsites
      genome@n.valid.sites       <- n.valid.sites    
      genome@n.gaps              <- n.gaps     
      genome@n.unknowns          <- n.unknowns   
      genome@n.polyallelic.sites <- n.polyallelic.sites
      genome@n.biallelic.sites   <- n.biallelic.sites
      genome@trans.transv.ratio  <- trans.transv.ratio


genome@region.data  <- region.data
genome@region.stats <- region.stats

genome@Pop_Neutrality$calculated <- FALSE
genome@Pop_FSTN$calculated       <- FALSE
genome@Pop_FSTH$calculated       <- FALSE
genome@Pop_MK$calculated         <- FALSE
genome@Pop_Linkage$calculated    <- FALSE
genome@Pop_Slide$calculated      <- FALSE
genome@Pop_Detail$calculated     <- FALSE

return(genome)

}# End of Function

