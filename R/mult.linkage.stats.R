setGeneric("mult.linkage.stats", function(object,new.populations="list",lower.bound=0,upper.bound=1) standardGeneric("mult.linkage.stats"))
 setMethod("mult.linkage.stats", "GENOME",

 function(object,new.populations,lower.bound,upper.bound){



# Check the populations ----------------------------------------

 if(!missing(new.populations)){
   NEWPOP <- TRUE
   npops <- length(new.populations)
   populations <- vector("list",npops)
  }else{
   NEWPOP <- FALSE
  } 

     if(!NEWPOP){

	npops       <- length(object@region.data@populations[[1]])
     	populations <- object@region.data@populations[[1]] # if there is no new population !!!!

     }else{
      	
	 bial  <- popGetBial(object,1) # !!!!!!!!!!!!!!!!
	 for(yy in 1:npops){
           if(is.character(new.populations[[yy]])){
              populations[[yy]] <- match(new.populations[[yy]],rownames(bial))
              naids             <- which(!is.na(populations[[yy]]))
              populations[[yy]] <- populations[[yy]][naids]
           }else{
              populations[[yy]] <- new.populations[[yy]]
              ids               <- which(populations[[yy]]>dim(bial)[1])
              if(length(ids)>0){ populations[[yy]] <- populations[[yy]][-ids]}   
           }   
       }
       
     } 


# End of checking populations ---------------------------------

 # define an evironment
 multGLOBAL <- new.env()

 pairs      <- combn(length(object@region.names),2)
 n.pairs    <- dim(pairs)[2]

 multGLOBAL$res        <- vector("list",n.pairs)
 multGLOBAL$iter       <- 1

#### NAMES ----------------------------------------
pp <- pairs
nn <- paste("",pp[1,1],"/",pp[2,1],sep="")
if(dim(pp)[2]>1){ # more than 2 sites
 for(yy in 2:dim(pp)[2]){
    m <- paste("",pp[1,yy],"/",pp[2,yy],sep="")
    nn <- c(nn,m)
 }
}#END if
#### ---------


## PROGRESS #########################
 progr <- progressBar()
#####################################

exx <- apply(pairs,2,function(xx){
		
	bial1 <- popGetBial(object,xx[1])
	bial2 <- popGetBial(object,xx[2])

        if(length(bial1)==0){multGLOBAL$iter <- multGLOBAL$iter + 1;return(0)}
	if(length(bial2)==0){multGLOBAL$iter <- multGLOBAL$iter + 1;return(0)}

	# Get the frequencies
        freq1 <- jointfreqdist(bial1,populations) 
        freq2 <- jointfreqdist(bial2,populations) 
 
	freq1 <- freq1$jfd        
	freq2 <- freq2$jfd

        subsites1 <- (freq1 >= lower.bound) & (freq1 <=upper.bound)
	subsites2 <- (freq2 >= lower.bound) & (freq2 <=upper.bound)
       	

        bial1 <- bial1[,subsites1,drop=FALSE]
        bial2 <- bial2[,subsites2,drop=FALSE]

 	multGLOBAL$res[[multGLOBAL$iter]]   <- pair_linkdisequ(bial1,bial2,populations)
     
 # PROGRESS #######################################################
    progr <- progressBar(multGLOBAL$iter,n.pairs, progr)
 ################################################################### 

        multGLOBAL$iter          <- multGLOBAL$iter + 1

   


 })

res           <- as.matrix(multGLOBAL$res)
rownames(res) <- nn
colnames(res) <- "region-pairwise Linkage"
object@mult.Linkage <- res

return(object)

})
