
# check each population. calculate the segragating sites (positions)

get_segsites <- function(matrix_pol,populations){

if(missing(populations)){

populations      <- vector("list",1)
populations[[1]] <- 1:dim(matrix_pol)[1]

}

npops    <- length(populations)
Segpos   <- vector("list",npops)

for(xx in 1:npops){
   popmatrix <- matrix_pol[populations[[xx]],,drop=FALSE]
   erg <- apply(popmatrix,2,function(x){
          check   <- unique(x)
          nongaps <- !is.na(check)
          check   <- check[nongaps]
          if(length(check)==2){return(TRUE)}else{return(FALSE)}
          })
   Segpos[[xx]] <- which(erg)
   
}
return(Segpos)
}
