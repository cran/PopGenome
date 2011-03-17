complike <- function(bial,populations){

freqs    <- jointfreqdist(bial,populations,outgroup=FALSE)

for(xx in 1:length(freqs)){
 
  CL      <-  sum(log(freqs[[xx]][freqs[[xx]]>0]))
  
  Counts  <- table(freqs[[xx]])
  p       <- as.numeric(names(Counts))
  val     <- Counts^p
  CLmax   <- sum(log(val[val>0]))  
}

return(list(CL=CL,CLmax=CLmax))

}
