setGeneric("region.as.fasta", function(object, region.id= FALSE, filename = FALSE) standardGeneric("region.as.fasta"))
setMethod("region.as.fasta", "GENOME",
function(object,region.id, filename){

if(region.id[1]==FALSE){
 stop("Define a region !")
}


bial    <- popGetBial(object,region.id)
subst   <- object@region.data@biallelic.substitutions[[region.id]]
minor   <- subst[1,]
mayor   <- subst[2,]

for(xx in 1:dim(bial)[2]){

  vek       <- bial[,xx]
  ids.minor <- vek==1
  ids.mayor <- vek==0

  bial[ids.minor,xx] <- minor[xx]
  bial[ids.mayor,xx] <- mayor[xx]

}

 

 number     <- c(1,1,1,1,2,2,3,3,4,4,5,5,5,6)
 nuc        <- c("T","t","U","u","C","c","G","g","A","a","N","n","?","-")


 bial <- apply(bial,1,function(x){return(nuc[match(x,number)])})
 bial <- t(bial)

write.dna(bial,file=filename,colsep="",format="fasta")

return(bial)

})
