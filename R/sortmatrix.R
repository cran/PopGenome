########################################################
################### sortmatrix #########################
########################################################
# A Matlab Function

sortmatrix <- function (matrix){

 m <- dim(matrix)[1]
 n <- dim(matrix)[2]

if(m==1 || n==1){
  return(sort(matrix))

}else{
     
     rownames(matrix) <- 1:m   
     
     
     for(xx in 1:n){
         mat <- unique(matrix[,xx])
         if(length(mat)>1){
          break;
         }
     }
     
     sortmatrix <- sort(matrix[,xx])
     sortmatrix <- matrix[names(sortmatrix),]
}

return(sortmatrix)
}