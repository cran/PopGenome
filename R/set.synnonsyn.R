setGeneric("set.synnonsyn", function(object,ref.chr=FALSE) standardGeneric("set.synnonsyn"))
 setMethod("set.synnonsyn", "GENOME",
 function(object,ref.chr){


if(ref.chr[1]==FALSE){
stop("Please verify the reference sequence")
}

for (xyz in 1:length(ref.chr)){
# erstmal nur für ein Chunk

Coding.matrix            <- object@region.data@Coding.matrix2[[xyz]][,] # weil ff object, 2 because (fitting GFF)
biallelic.sites2         <- object@region.data@biallelic.sites2[[xyz]]
biallelic.sites          <- object@region.data@biallelic.sites[[xyz]]
START                    <- object@region.data@reading.frame[[xyz]][,] # weil ff object
CodingSNPS               <- object@region.data@CodingSNPS[[xyz]]

START                    <- START[,1] + START[,2]
Coding.matrix            <- Coding.matrix

#print(length(START)) ist ok :)
#print(dim(Coding.matrix)) ist ok :)

# ---
# unique
# wenn regionen doppelt
#ids                     <- !duplicated(Coding.matrix[,1])
#START                   <- START[ids]
#Coding.matrix           <- Coding.matrix[ids,,drop=FALSE]


# Create Region and save size of region 
SIZE  <<- numeric(dim(Coding.matrix)[xyz])
count <<- 1
erg  <- apply(Coding.matrix,1,function(xx){
 region        <- xx[1]:xx[2] 
 SIZE[count]   <<- length(region)
 count         <<- count + 1 
 return(region)

})

# what are the real positions ?
 erg         <- unlist(erg)
 bial.pos    <- .Call("my_match_C",erg,biallelic.sites2)
 bial.pos[bial.pos==-1] <- NaN
 #return(bial.pos)
 bial.pos    <- biallelic.sites[bial.pos]

# Create Start Vector
 count <<- 1
 vec <- sapply(START,function(x){
      gg    <- rep(x,SIZE[count])       
      count <<- count + 1
 return(gg)
 })

 START.vec   <- unlist(vec)
 cod.pos     <- (bial.pos - START.vec)%%3
 

 # Schmeisse NaNs raus
 cod.pos     <- cod.pos[!is.na(bial.pos)]
 bial.pos    <- bial.pos[!is.na(bial.pos)]
 
 ids         <- !duplicated(bial.pos)
 cod.pos     <- cod.pos[ids]
 bial.pos    <- bial.pos[ids]


# Create the codons
# bial.pos and cod.pos

codons <- matrix(,length(cod.pos),3)
for (xx in 1:length(cod.pos)){
    if(cod.pos[xx]==0){codons[xx,]=c(bial.pos[xx],bial.pos[xx]+1,bial.pos[xx]+2);next}
    if(cod.pos[xx]==1){codons[xx,]=c(bial.pos[xx]-1,bial.pos[xx],bial.pos[xx]+1);next}
    if(cod.pos[xx]==2){codons[xx,]=c(bial.pos[xx]-2,bial.pos[xx]-1,bial.pos[xx]);next}
}

## Reading the reference chromosome
file.info <- .Call("get_dim_fasta",ref.chr[xyz])
CHR       <- .Call("get_ind_fasta",ref.chr,1,file.info[[1]][2])


# Create codons with nucleotides
Nuc.codons    <- CHR[codons]
Nuc.codons    <- matrix(Nuc.codons,ncol=3)

ALT           <- Nuc.codons
REF           <- Nuc.codons
Subst         <- object@region.data@biallelic.substitutions[[xyz]]
minor         <- Subst[1,CodingSNPS]
major         <- Subst[2,CodingSNPS]

for(xx in 1: dim(Nuc.codons)[1]){
 if(cod.pos[xx]==0){REF[xx,1] <- minor[xx];ALT[xx,1]<-major[xx];next}
 if(cod.pos[xx]==1){REF[xx,2] <- minor[xx];ALT[xx,2]<-major[xx];next}
 if(cod.pos[xx]==2){REF[xx,3] <- minor[xx];ALT[xx,3]<-major[xx];next}
} 

# Coding Codons ...
ALT <- codonise64(ALT)
REF <- codonise64(REF)

CC  <- codontable()

ALT <- CC$Protein[1,ALT]
REF <- CC$Protein[1,REF]

CHECK <- cbind(ALT,REF)

erg <- apply(CHECK,1,function(x){return(length(unique(x)))})
erg[erg==2] <- FALSE #nonsyn
erg[erg==1] <- TRUE  #syn

# Change object of class GENOME
change <- object@region.data
change@synonymous[[xyz]][CodingSNPS] <- erg
object@region.data <- change 

}# End Iteration over chunks or chromosomes

return(object)

})
