get_gff_info <- function(object=FALSE, gff.file, chr, position){


if(is(object)=="logical"){
 chr    <- as.character(chr)
 region <- .Call("find_lines_GFF",gff.file,chr)
 start  <- region[1]
 end    <- region[2]
 info   <- .Call("get_gff_info_C",start,end,gff.file,position)
return(info)
}

if(is(object)=="GENOME"){

 chr    <- as.character(chr)
 region <- .Call("find_lines_GFF",gff.file,chr)
 start  <- region[1]
 end    <- region[2]

RET.INFO <- vector("list",length(position))

for(xx in 1:length(position)){

 region <- position[xx]
 region <- object@region.data@biallelic.sites[[region]]
 info   <- sapply(region,function(x){return(.Call("get_gff_info_C",start,end,gff.file,x))})
 names(info)    <- region
 RET.INFO[[xx]] <- info

}

RET.INFO           <- as.matrix(RET.INFO)
rownames(RET.INFO) <- position

return(RET.INFO)
}


}
