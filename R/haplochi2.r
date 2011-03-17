haplochi2 <- function(populations,nh,sfreqh){ # Chi2Test

npops <- length(populations)

# Get sample size of each population
n <- vector(,npops)
for(xx in 1:npops){
   n[xx] <- length(populations[[xx]])
}

# Number of haplotypes
 K              <- dim(sfreqh)[2]
 Sumnij              <- colSums(sfreqh)
 p              <- Sumnij/K
 
 return(list(K=K,Sumnij=Sumnij,p=p))
 }