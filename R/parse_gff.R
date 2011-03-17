parse_gff <- function(tab,transform=FALSE,SNP.DATA){



#if(length(tab$seqname)!=0){ # if is a list

#SEQNAME <- tab$seqname
#SOURCE  <- tab$source
#FEATURE <- tab$feature
#START   <- tab$start
#END     <- tab$end
#SCORE   <- tab$score 
#STRAND  <- tab$strand
#FRAME   <- tab$frame
#ATTR    <- tab$attributes
#tab     <- cbind(SEQNAME,SOURCE,FEATURE,START,END,SCORE,STRAND,FRAME,ATTR)

#}

if(transform){return(tab)}

Coding_ids  <- which(tab[,3]=="CDS")
Coding      <- tab[Coding_ids,4:5,drop=FALSE]
#Coding.info <- tab[Coding_ids,9]

if(!SNP.DATA){
readframe   <- tab[Coding_ids,8]

if(length(readframe)>0){

 for (xx in 1: dim(Coding)[1]){

     if(readframe[xx]==1){Coding[xx,1] <- as.numeric(Coding[xx,1]) + 1}
     if(readframe[xx]==2){Coding[xx,1] <- as.numeric(Coding[xx,1]) + 2}
 }
}
}# END SNP


Exon_ids   <- which(tab[,3]=="exon")
Exon       <- tab[Exon_ids,4:5,drop=FALSE]
#Exon.info  <- tab[Exon_ids,9]

Gene_ids   <- which(tab[,3]=="gene")
Gene       <- tab[Gene_ids,4:5,drop=FALSE]
#Gene.info  <- tab[Gene_ids,9]

Intron_ids  <- which(tab[,3]=="intron")
Intron      <- tab[Intron_ids,4:5,drop=FALSE]
#Intron.info <- tab[Intron_ids,9]

UTR_ids    <- which(tab[,3]=="UTR")
UTR        <- tab[UTR_ids,4:5,drop=FALSE]
#UTR.info   <- tab[UTR_ids,9]

### GLOBAL important for SNP DATA
READ <- NULL

if(length(GLOBAL.GFF$GFF)!=0){

 CHECK                    <- GLOBAL.GFF$GFF
 GLOBAL.GFF$Coding        <- CHECK[Coding_ids]
 GLOBAL.GFF$Exon          <- CHECK[Exon_ids]
 GLOBAL.GFF$Gene          <- CHECK[Gene_ids]
 GLOBAL.GFF$Intron        <- CHECK[Intron_ids]
 GLOBAL.GFF$UTR           <- CHECK[UTR_ids]
 #GLOBAL.GFF$GFF          <- NULL
 READ                     <- cbind(GLOBAL.GFF$start[Coding_ids],tab[Coding_ids,8])
 if(dim(READ)[1]>0){
  READ                    <- ff(READ,dim=dim(READ))
 }
}

return(list(Coding=Coding,Intron=Intron,UTR=UTR,Exon=Exon,Gene=Gene,reading.frame= READ))
#,Coding.info=Coding.info,Intron.info=Intron.info,
# Gene.info=Gene.info,UTR.info=UTR.info))

}

# Ausm Internet

getAttributeField <- function (x, field, attrsep = ";") {
     s = strsplit(x, split = attrsep, fixed = TRUE)
     sapply(s, function(atts) {
         a = strsplit(atts, split = "=", fixed = TRUE)
         m = match(field, sapply(a, "[", 1))
         if (!is.na(m)) {
             rv = a[[m]][2]
         }
         else {
             rv = as.character(NA)
         }
         return(rv)
     })
}



gffRead <- function(gffFile, nrows = -1) {
     #cat("Reading ", gffFile, ": ", sep="")
     gff = read.table(gffFile, sep="\t", as.is=TRUE, quote="",
     header=FALSE, comment.char="#", nrows = nrows,
     colClasses=c("character", "character", "character", "integer",  
    "integer",
    "character", "character", "character", "character")
)
     colnames(gff) = c("seqname", "source", "feature", "start", "end",
             "score", "strand", "frame", "attributes")
     #cat("found", nrow(gff), "rows with classes:",
         #paste(sapply(gff, class), collapse=", "), "\n")
     stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))
     return(gff)
}

#Now you can do stuff like

#gff <- gffRead(gfffile)
#gff$Name <- getAttributeField(gff$attributes, "Name")
#gff$ID <- getAttributeField(gff$attributes, "ID")
