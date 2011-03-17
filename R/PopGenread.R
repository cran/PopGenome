PopGenread <- function(filepath,format="fasta") {

  gen           <-  .Call("readdna",filepath,PACKAGE="PopGenome")
  rownames(gen) <-  gsub(" ","",rownames(gen))
  return(gen)

}
