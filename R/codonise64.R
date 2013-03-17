########################################################
##################### codonise64 #######################
########################################################

codonise64 <- function(s){


m<-dim(s)[2]
n<-dim(s)[1]

if((m%%3)>0){
 stop("length of coding sequence cannot divide by 3!") 
}

t <- matrix(,n,m/3)
p <- seq(1,m,3)
q <- 1:(m/3)

s[s>=5] <- 100 # if gaps or unknown in the codon

t[ ,q] <- (s[ ,p]-1)*16 + (s[ ,p+1]-1)*4 + (s[ ,p+2]-1)+1
t[t>64]<-65

return(t)
}
