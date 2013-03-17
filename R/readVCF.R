
readVCF <- function( filename, numcols, tid, frompos, topos, samplenames=NA, gffpath = FALSE, include.unknown=FALSE, approx=TRUE)
{

frompos <- as.integer(frompos)
topos   <- as.integer(topos)

if(file.exists("SNPRObjects")){unlink("SNPRObjects",recursive=TRUE)}
if(file.exists("GFFRObjects")){unlink("GFFRObjects",recursive=TRUE)}


GFF <- FALSE
### Prepare GFF INFO # ----------------

if(gffpath[1]!=FALSE){

GFF <- TRUE

dir.create("GFFRObjects")

cat("\n")
cat("GFF information ...")
cat("\n")
  
   tid2 <- tid
   if(length(grep("Chr",tid2))==1 | length(grep("chr",tid2))==1){
   tid2 <- substr(tid2,4,nchar(tid2))
   }
   if(nchar(tid2)==1){CHR <- c(tid2,"z")}else{CHR <- unlist(strsplit(tid2,""))}
   my_pos      <- .Call("find_lines_GFF_Human", gffpath, CHR)
   gff.table   <- read.table(gffpath,sep="\t",colClasses=c(rep("character",3),rep("integer",2),rep("character",2),"character","NULL"),
                             skip = my_pos[1], nrows = my_pos[2] - my_pos[1]
                            )


# Sort the GFF- file -------
POS        <- gff.table[,4]
names(POS) <- 1:length(POS)
POS        <- sort(POS)
ids        <- as.numeric(names(POS))
gff.table  <- gff.table[ids,,drop=FALSE]

## vorletzte Zeile integer draus wegen GFF3 "." -> 0
vorl       <- gff.table[,8]
punkte     <- which(vorl==".")
if(length(punkte)>0){
  vorl[punkte]  <- 0
  gff.table[,8] <- as.integer(vorl) 
}else{
  gff.table[,8] <- as.integer(gff.table[,8])
}
# -------------------------

gffpath <- "GFFRObjects"

}






##-------------------------------------


	outdir <- "SNPRObjects"

        #
	resvec = dir.create( path = outdir , recursive=T, showWarnings=F )
	oldcurdir = setwd(outdir)	# will produce an error if that does not work
	setwd(oldcurdir)
	
	#
	#
	v <- .Call("VCF_open", filename )
	if( is.null( v ) )
	{
		# setwd(curdir)
		return(FALSE);
	}
	
	#
	#
	alltids = .Call("VCF_getContigNames",v)
	print("Available ContigIdentifiers (parameter tid):")
	print(alltids)
	if( ! (tid %in% alltids) )
	{
		stop("Parameter <tid> invalid : not a chromosome/contig Identifier found in the VCF file!")
	}
	
	#
	#
	sn <- .Call("VCF_getSampleNames",v)
        # sn <- sn[1:(length(sn)-1)] # FIXME	
         
	#
	#
	if( ! any( is.na( samplenames ) ) )
	{
		schnittmengesamples = samplenames[ samplenames %in% sn ]
	}
	else
	{
		schnittmengesamples = sn
	}
	
	#
	stopifnot( length( schnittmengesamples ) > 0 )
	
	#print( schnittmengesamples )
	.Call("VCF_selectSamples",v,schnittmengesamples )
	
	#
	#
	regset <- .Call("VCF_setRegion",v,tid,frompos,topos)
	if( regset == FALSE )
	{
		stop("Region could not be set!");
	}

	#16050408,17059916,		17066020
	#
#	mm <- matrix( nrow=12,ncol=10,data="-", dimnames=list(sl,rep("x",10)) )
	mi <- matrix( nrow=length(schnittmengesamples),ncol=numcols,data=as.integer(0), dimnames=list(schnittmengesamples,rep("x",numcols)) )
	if(approx==FALSE){
	# save names for diploid data
	dottwo    <- paste(schnittmengesamples,".2", sep="") 
  	diplNAMES <- as.vector(rbind(schnittmengesamples,dottwo)) 
	}

	#
	#
	numusedcols=numcols
	filenum=0
	fileprefix=paste(sep="","chr:","_",tid,"_",frompos,"-",topos,"_")

	#
	while( numusedcols == numcols )
	{
 
		#
		if(approx){
		.Call("VCF_readIntoCodeMatrix",v,mi)
		}else{
		.Call("VCF_readIntoCodeMatrixdiploid",v,mi)
		}
		#
		setcols     <- as.integer( colnames(mi) ) > 0
		numusedcols <- sum(setcols)

		#if(numusedcols==0){
                #   .Call("VCF_close",v)
                #   stop("no SNPs !")
                #}
		#

		cn <- as.integer( colnames(mi)[1:numusedcols] )
		#cat( paste(sep="","used = ", numusedcols,":\n" ) )
		#print( cn )
		
		#
		fullfilename = paste( sep="", outdir, "/", filenum,":","chr",tid,"_",cn[1],"-",cn[numusedcols],"_",".RData" )
		filenum      = filenum + 1
		#print( fullfilename )
		
		#
		
		#
		#print(cn[1])
		#print(cn[numusedcols])
		#print(mi)
		if(approx){
		o_b_j_sub  <- list(matrix=mi,reference=NaN, positions=cn)
                }else{ # create diploid matrix
                diplmi <- as.character(mi)
                diplmi <- strsplit(diplmi,split="")
                diplmi <- as.numeric(unlist(diplmi))
                diplmi <- matrix(diplmi, length(diplNAMES), dim(mi)[2]) 
		#diplmi     <- mi
		# modify matrix
		#diplmi     <- matrix(as.character(diplmi),nrow=dim(diplmi)[1],ncol=dim(diplmi)[2]) #apply(diplmi,2,as.character)
		#diplmi     <- apply(diplmi,2,function(x){as.numeric(sapply(strsplit(x,split=""),rbind))})
                rownames(diplmi) <- diplNAMES		
		o_b_j_sub  <- list(matrix=diplmi,reference=NaN, positions=cn)
                }

		save( file = fullfilename , o_b_j_sub  )

                if(GFF){

                  SUB_GFF      <- split.GFF(gff.table,cn[1],cn[numusedcols])
                  
                  if(length(SUB_GFF)!=0){
                   o_b_j_sub    <- SUB_GFF
                  
                   fullfilename <- paste( sep="","GFFRObjects", "/", (filenum-1),":","chr",tid,"_",cn[1],"-",cn[numusedcols],"_",".RData" )
                   save( file = fullfilename , o_b_j_sub  )
                  }

                }


		#
	}#..while( 
	
       # .Call("VCF_close",v)
	# print(fileprefix)


res         <- readData(outdir,populations = FALSE, outgroup=FALSE, include.unknown=include.unknown, gffpath=gffpath,format="RData",
                parallized=FALSE,progress_bar_switch=TRUE,
                FAST=TRUE,big.data=TRUE,SNP.DATA=TRUE)


if(res@genelength>1){

 res              <- concatenate_to_whole_genome(res,res@genelength)

}
 
 res@region.names   <- paste(tid,":",frompos,"-",topos,sep=" ")
 res@n.sites        <- as.numeric((topos - frompos + 1))
 res@keep.start.pos <- frompos
 res@gff.info       <- GFF

# Delete
unlink("SNPRObjects", recursive=TRUE)
if(GFF){
unlink("GFFRObjects", recursive=TRUE)
}

return( res )
        
}
