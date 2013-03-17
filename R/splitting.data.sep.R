setGeneric("splitting.data.sep", function(object,subsites=FALSE, positions=FALSE, type=1) standardGeneric("splitting.data.sep"))
 setMethod("splitting.data.sep", "GENOME",

 function(object,subsites,positions,type){


if(subsites[1]!=FALSE & !object@gff.info){
stop("No GFF/GTF informations !")
}

# What information do we need from the original object ?!
if(subsites=="coding"){
N.subsites <- sapply(object@region.data@Coding.matrix,function(x){
              if(length(dim(x)[1])==0){return(0)}else{return(dim(x)[1])}
              })
}
if(subsites=="gene"){
N.subsites <- sapply(object@region.data@Gene.matrix,function(x){
              if(length(dim(x)[1])==0){return(0)}else{return(dim(x)[1])}
              })
}
if(subsites=="exon"){
N.subsites <- sapply(object@region.data@Exon.matrix,function(x){
              if(length(dim(x)[1])==0){return(0)}else{return(dim(x)[1])}
              })
}
if(subsites=="utr"){
N.subsites <- sapply(object@region.data@UTR.matrix,function(x){
              if(length(dim(x)[1])==0){return(0)}else{return(dim(x)[1])}
              })
}
if(subsites=="intron"){
N.subsites <- sapply(object@region.data@Intron.matrix,function(x){return(dim(x)[1])
              if(length(dim(x)[1])==0){return(0)}else{return(dim(x)[1])}
              })
}

# Check the blocks of the Gene.matrix
# bial.pos.ids <- matrix(,length(N.subsites),2)
# start <- 1
# for(xx in 1:length(N.subsites)){
#    bial.pos.ids[xx,] <- c(start,N.subsites[xx])
#   start <- N.subsites[xx] + 1
# }
###########################################

N.Bial.Sites <- object@n.biallelic.sites

# Check the blocks of the Gene matrix
gene.pos.ids <- matrix(,length(N.subsites),2)
start        <- 1
ADDX         <- 0

for(xx in 1:length(N.subsites)){
   ADDX              <- ADDX + N.subsites[xx]
   gene.pos.ids[xx,] <- c(start,ADDX)
   start             <- ADDX + 1
}
###########################################
REGION.NAMES <- object@region.names


# ------------------------------------------------------------

if(object@genelength>1){
# object@snp.data = TRUE # to manipulate the concatenation
object          <- concatenate_to_whole_genome(object,object@genelength)
}


genomeobj               <-  new("GENOME") 
ddatt                   <-  new("region.data")
XXX                     <-  object@region.data 


genomeobj@BIG.BIAL[[1]] <- object@region.data@biallelic.matrix[[1]]
 
if(!is.list(positions)){

 if(subsites=="coding"){
    Gene.matrix <- object@region.data@Coding.matrix[[1]][,]
    LENGTH      <- GLOBAL.GFF$Coding
 }
 if(subsites=="exon"){
    Gene.matrix <- object@region.data@Exon.matrix[[1]][,]
    LENGTH      <- GLOBAL.GFF$Exon
 }
 if(subsites=="utr"){
    Gene.matrix <- object@region.data@UTR.matrix[[1]][,]
    LENGTH      <- GLOBAL.GFF$UTR
 }
 if(subsites=="intron"){
    Gene.matrix <- object@region.data@Intron.matrix[[1]][,]
    LENGTH      <- GLOBAL.GFF$Intron
 }
 if(subsites=="gene"){
    Gene.matrix <- object@region.data@Gene.matrix[[1]][,]
    LENGTH      <- GLOBAL.GFF$Gene
 }
 
 if(length(dim(Gene.matrix)[1])==0 || dim(Gene.matrix)[1]==0){
  cat("\n")
  stop(paste("No ",subsites, "defined !"))
 }
 
}
 
 # init -------------------------------------------------------
  if(!is.list(positions)){
     init             <- vector("list",dim(Gene.matrix)[1])
     region.names     <- character(dim(Gene.matrix)[1])
     n.sites          <- numeric(dim(Gene.matrix)[1])
  }else{
     init             <- vector("list",length(positions))    
     region.names     <- character(length(positions))  
     n.sites          <- numeric(length(positions))
  }

  SLIDE.POS        <- init            
  biallelic.matrix <- init
  biallelic.sites  <- init
  outgroup         <- init
  populations      <- init
  popmissing       <- init
  transitions      <- init
  synonymous       <- init
 
  #-------------------------------------------------------------

    if(!is.list(positions) & !object@snp.data){
    
     ### Progress
     progr <- progressBar()
     ###
  
     count <- 1
     for(zz in 1:dim(Gene.matrix)[1]){
   
        bialpos                       <- Gene.matrix[zz,1]:Gene.matrix[zz,2]
        bialpos                       <- is.element(XXX@biallelic.sites[[1]],bialpos)
        bialpos                       <- which(bialpos)
   
        if(zz==(gene.pos.ids[count,2]+1)){count <- count + 1} 
        CHUNK.name                    <- REGION.NAMES[count]
        region.names[zz]              <- paste(CHUNK.name,": ",Gene.matrix[zz,1],"-",Gene.matrix[zz,2])
        n.sites[zz]                   <- (Gene.matrix[zz,2]-Gene.matrix[zz,1]) + 1
        SLIDE.POS[[zz]]               <- bialpos
        outgroup[[zz]]                <- XXX@outgroup[[1]]
        populations[[zz]]             <- XXX@populations[[1]]
        synonymous[[zz]]              <- XXX@synonymous[[1]][bialpos]
        transitions[[zz]]             <- XXX@transitions[[1]][bialpos]
        biallelic.sites[[zz]]         <- XXX@biallelic.sites[[1]][bialpos]  

     ## Progress
     progr <- progressBar(zz,dim(Gene.matrix)[1], progr)
     ####
       
     }
    }
 
    if(!is.list(positions) & object@snp.data){

     ### Progress
     progr <- progressBar()
     ###
     
     # Check the blocks of the baillelic.matrix
     bial.pos.ids  <- matrix(,length(N.Bial.Sites),2)
     start         <- 1
     ADDX          <- 0

     for(xx in 1:length(N.Bial.Sites)){
      ADDX               <- ADDX + N.Bial.Sites[xx]
      bial.pos.ids [xx,] <- c(start,ADDX)
      start              <- ADDX + 1
     }
     count <- 1
     
###########################################

     bial.sites <- XXX@biallelic.sites[[1]][bial.pos.ids[count,1]:bial.pos.ids[count,2]] # first chunk
     JUMP <- 0
     for(zz in 1:dim(Gene.matrix)[1]){
   
        #print(JUMP)
        #print(bial.pos.ids)
        #print(N.Bial.Sites)
        #print(gene.pos.ids)

        if(zz==(gene.pos.ids[count,2]+1)){ # end of chunk
           JUMP       <- JUMP    + N.Bial.Sites[count]  
           count      <- count + 1
           bial.sites <- XXX@biallelic.sites[[1]][bial.pos.ids[count,1]:bial.pos.ids[count,2]] # next chunk        
        } 

         #print("WAT")
         CHUNK.name                    <- REGION.NAMES[count]
         region.names[zz]              <- paste(CHUNK.name,":",Gene.matrix[zz,1],"-",Gene.matrix[zz,2])
         n.sites[zz]                   <- Gene.matrix[zz,2] - Gene.matrix[zz,1] + 1

         #bialpos                       <- Gene.matrix[zz,1]:Gene.matrix[zz,2]
         #bialpos                       <- match(bialpos,XXX@biallelic.sites2[[1]])
         #bialpos                       <- bialpos[!is.na(bialpos)]
          
         #print(bial.sites) 
         #print(Gene.matrix[zz,1])
         #print(Gene.matrix[zz,2])
         bialpos                       <- .Call("find_windowC",bial.sites, Gene.matrix[zz,1],Gene.matrix[zz,2],1)
         # bialpos                       <- XXX@biallelic.sites2[[1]][bialpos[1]:bialpos[2]] 
         if(length(bialpos)==0){next} 
         bialpos                       <- bialpos[1]:bialpos[2]     

         #start.region                  <- XXX@biallelic.sites[[1]][bialpos[1]]
         #end.region                    <- XXX@biallelic.sites[[1]][bialpos[length(bialpos)]]
         
        #region.names[zz]             <- paste(Gene.matrix[zz,1],"-",Gene.matrix[zz,2])
       

        
        SLIDE.POS[[zz]]               <- bialpos + JUMP
        outgroup[[zz]]                <- XXX@outgroup[[1]]
        populations[[zz]]             <- XXX@populations[[1]]
        synonymous[[zz]]              <- XXX@synonymous[[1]][bialpos + JUMP]
        transitions[[zz]]             <- XXX@transitions[[1]][bialpos + JUMP]
        biallelic.sites[[zz]]         <- XXX@biallelic.sites[[1]][bialpos + JUMP]  
        
     ## Progress
     progr <- progressBar(zz,dim(Gene.matrix)[1], progr)
     ####
       
     }
    }
   


   if(is.list(positions) & type==1) {

     ### Progress
     progr <- progressBar()
     ###

      for(zz in 1:length(positions)){

        bialpos                       <- positions[[zz]]
        region.names[zz]              <- paste(positions[[zz]][1],"-",positions[[zz]][length(positions[[zz]])])
        SLIDE.POS[[zz]]               <- bialpos
        outgroup[[zz]]                <- XXX@outgroup[[1]]
        populations[[zz]]             <- XXX@populations[[1]]
        synonymous[[zz]]              <- XXX@synonymous[[1]][bialpos]
        transitions[[zz]]             <- XXX@transitions[[1]][bialpos]
        biallelic.sites[[zz]]         <- XXX@biallelic.sites[[1]][bialpos]  

     ## Progress
     progr <- progressBar(zz,length(positions), progr)
     ####

      }
   }

   if(is.list(positions) & type==2) {


     ### Progress
     progr <- progressBar()
     ###
   
     for(zz in 1:length(positions)){
     
        bialpos                       <- positions[[zz]]
        bialpos                       <- is.element(XXX@biallelic.sites[[1]],bialpos)
        bialpos                       <- which(bialpos)
        
        region.names[zz]              <- paste(positions[[zz]][1],"-",positions[[zz]][length(positions[[zz]])])
        SLIDE.POS[[zz]]               <- bialpos
        outgroup[[zz]]                <- XXX@outgroup[[1]]
        populations[[zz]]             <- XXX@populations[[1]]
        synonymous[[zz]]              <- XXX@synonymous[[1]][bialpos]
        transitions[[zz]]             <- XXX@transitions[[1]][bialpos]
        biallelic.sites[[zz]]         <- XXX@biallelic.sites[[1]][bialpos] 
  
     ## Progress
     progr <- progressBar(zz,length(positions), progr)
     ####
     
     }
   } 



ddatt@UTR.matrix       <- object@region.data@UTR.matrix
ddatt@Exon.matrix      <- object@region.data@Exon.matrix
ddatt@Gene.matrix      <- object@region.data@Gene.matrix
ddatt@Intron.matrix    <- object@region.data@Intron.matrix
ddatt@Coding.matrix    <- object@region.data@Coding.matrix
ddatt@biallelic.matrix <- biallelic.matrix
ddatt@biallelic.sites  <- biallelic.sites
ddatt@populations      <- populations
ddatt@outgroup         <- outgroup
ddatt@popmissing       <- popmissing
ddatt@transitions      <- transitions
ddatt@synonymous       <- synonymous

genomeobj@SLIDE.POS                 <- SLIDE.POS
genomeobj@populations               <- object@populations
genomeobj@region.names              <- region.names
genomeobj@n.sites                   <- n.sites
genomeobj@genelength                <- length(SLIDE.POS)
genomeobj@region.data               <- ddatt
genomeobj@Pop_Neutrality$calculated <- FALSE
genomeobj@Pop_FSTN$calculated       <- FALSE
genomeobj@Pop_FSTH$calculated       <- FALSE
genomeobj@Pop_MK$calculated         <- FALSE
genomeobj@Pop_Linkage$calculated    <- FALSE
genomeobj@Pop_Recomb$calculated     <- FALSE
genomeobj@Pop_Slide$calculated      <- TRUE
genomeobj@Pop_Detail$calculated     <- FALSE
genomeobj@big.data                  <- object@big.data
genomeobj@snp.data                  <- object@snp.data

return(genomeobj)

})
  
