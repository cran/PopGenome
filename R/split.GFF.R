split.GFF <- function(GFF,pos1,pos2){

# find the subset
GFF.pos   <- GFF[,4]
regionX   <- .Call("find_windowC",GFF.pos,pos1,pos2,1)
GFF.pos   <- GFF[,5]
regionY   <- .Call("find_windowC",GFF.pos,pos1,pos2,1)

region    <- c(regionX,regionY)

if(length(region)!=0){
 region   <- sort(region)
 start    <- region[1]
 end      <- region[length(region)]
}else{return(NULL)}

return(GFF[start:end,,drop=FALSE])
}
