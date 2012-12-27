
readVCF <- function( filename, numcols, tid, frompos, topos, samplenames=NA, gffpath = FALSE )
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
  
   if(nchar(tid)==1){CHR <- c(tid,"z")}else{CHR <- tid}
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

	#
	#
	numusedcols=numcols
	filenum=0
	fileprefix=paste(sep="","chr:","_",tid,"_",frompos,"-",topos,"_")

	#
	while( numusedcols == numcols )
	{
		#
		.Call("VCF_readIntoCodeMatrix",v,mi)
		
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
		o_b_j_sub  <- list(matrix=mi,reference=NaN, positions=cn)
		save( file = fullfilename , o_b_j_sub  )

                if(GFF){

                  SUB_GFF      <- split.GFF(gff.table,cn[1],cn[numusedcols])
                  o_b_j_sub    <- SUB_GFF
                  fullfilename <- paste( sep="","GFFRObjects", "/", (filenum-1),":","chr",tid,"_",cn[1],"-",cn[numusedcols],"_",".RData" )
                  save( file = fullfilename , o_b_j_sub  )

                }


		#
	}#..while( 
	
       # .Call("VCF_close",v)
	# print(fileprefix)


res         <- readData(outdir,populations = FALSE, outgroup=FALSE, include.unknown=FALSE, gffpath=gffpath,format="RData",
                parallized=FALSE,progress_bar_switch=TRUE,
                FAST=TRUE,big.data=TRUE,SNP.DATA=TRUE)


if(res@genelength>1){

 res <- concatenate_to_whole_genome(res,res@genelength)

}

# Delete
unlink("SNPRObjects", recursive=TRUE)
if(GFF){unlink("GFFRObjects", recursive=TRUE)}

return( res )
        
}
