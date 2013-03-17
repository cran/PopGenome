readMS <- function(file){

out     <- read.ms.output(file.ms.output=file)
gametes <- out$gametes

dir.create("SwapMS")

for(xx in 1:length(gametes)){

 d <- gametes[[xx]]
 d <- list(matrix=d,positions=NaN)
 samplename <- paste("ms_sample_",xx,".RD",sep="")
 save(d,file= file.path ("SwapMS",samplename) )

}

test <- readData("SwapMS", SNP.DATA=TRUE, FAST=TRUE, format="RData")

unlink("SwapMS",recursive=TRUE)

return(test)

}
