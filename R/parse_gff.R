parse_gff <- function(tab){



Coding_ids <- which(tab[,3]=="CDS")
Coding     <- tab[Coding_ids,4:5,drop=FALSE]
readframe  <- tab[Coding_ids,8]

if(length(readframe)>0){
 for (xx in 1: dim(Coding)[1]){

     if(readframe[xx]==1){Coding[xx,1] <- Coding[xx,1] + 1}
     if(readframe[xx]==2){Coding[xx,1] <- Coding[xx,1] + 2}
 }
}

Intron_ids <- which(tab[,3]=="intron")
Intron     <- tab[Intron_ids,4:5,drop=FALSE]

UTR_ids    <- which(tab[,3]=="UTR")
UTR        <- tab[UTR_ids,4:5,drop=FALSE]



return(list(Coding=Coding,Intron=Intron,UTR=UTR))

}

