# ------------------------------------------------------------
# get Codon Informations 
# ------------------------------------------------------------
setGeneric("get.codons", function(object, regionID, reading.start.pos, ref.chr) standardGeneric("get.codons"))
 setMethod("get.codons", "GENOME",
 function(object, regionID, reading.start.pos, ref.chr){


CodonInfo <- codontable()
Triplets  <- CodonInfo$Triplets
Protein   <- CodonInfo$Protein

# Generate polarity Codon Vector
POLARITY <- character(64)
POLARITY[c(1:4,11,17:24,33:44,49:52)] <- "nonpolar"
POLARITY[c(5:10,13:14,53:58,61:62)]   <- "polar"
POLARITY[c(12,15:16)]                 <- "stop codon"
POLARITY[c(25:32,59:60,63:64)]        <- "basic"
POLARITY[45:48]                       <- "acidic"
# --------------------------------


# Generate Codon Tables 
bial.pos     <- object@region.data@biallelic.sites[[regionID]]
minor        <- object@region.data@biallelic.substitutions[[regionID]][1,]
major	     <- object@region.data@biallelic.substitutions[[regionID]][2,]


cod.pos      <- (bial.pos - reading.start.pos) %%3
# Generate Codon-Positions # Rows are SNP-positions
codons <- matrix(,length(cod.pos),3)
for (yy in 1:length(cod.pos)){
    if(cod.pos[yy]==0){codons[yy,]=c(bial.pos[yy],bial.pos[yy]+1,bial.pos[yy]+2);next}
    if(cod.pos[yy]==1){codons[yy,]=c(bial.pos[yy]-1,bial.pos[yy],bial.pos[yy]+1);next}
    if(cod.pos[yy]==2){codons[yy,]=c(bial.pos[yy]-2,bial.pos[yy]-1,bial.pos[yy]);next}
}

## Reading the reference chromosome
file.info <- .Call("get_dim_fasta",ref.chr)
CHR       <- .Call("get_ind_fasta",ref.chr,1,file.info[[1]][2])

# Create codons with nucleotides
Nuc.codons    <- CHR[codons]
Nuc.codons    <- matrix(Nuc.codons,ncol=3)

ALT           <- Nuc.codons
REF           <- Nuc.codons

for(yy in 1: dim(Nuc.codons)[1]){
 if(cod.pos[yy]==0){REF[yy,1] <- minor[yy];ALT[yy,1] <- major[yy];next}
 if(cod.pos[yy]==1){REF[yy,2] <- minor[yy];ALT[yy,2] <- major[yy];next}
 if(cod.pos[yy]==2){REF[yy,3] <- minor[yy];ALT[yy,3] <- major[yy];next}
} 

# Character Codons
MATCHNUCS <- c("T","C","G","A","N")
NucALT    <- matrix(MATCHNUCS[ALT], ncol=3) 
NucREF    <- matrix(MATCHNUCS[REF], ncol=3) 

NucALT <- apply(NucALT,1,function(x){paste(x,collapse="")})
NucREF <- apply(NucREF,1,function(x){paste(x,collapse="")})

# Coding Codons ...
ALT <- codonise64(ALT)
REF <- codonise64(REF)

# Polarity
ALTPOL <- POLARITY[ALT]
REFPOL <- POLARITY[REF]

# Proteins
ALT <- Protein[1,ALT]
REF <- Protein[1,REF]

# ALT and REF Proteins
CHECK <- cbind(ALT,REF)

# Check syn/nonsyn
erg         <- apply(CHECK,1,function(x){return(length(unique(x)))})
erg[erg==2] <- FALSE #nonsyn
erg[erg==1] <- TRUE  #syn


# Collect information and store in data.frame 
XX    <- c("Position","Codons (major)","Codons (minor)","Protein (major)","Protein (minor)","synonymous", "Polarity (major)","Polarity (minor)")
DATA  <- data.frame(bial.pos,NucALT, NucREF, ALT, REF, as.logical(erg), ALTPOL, REFPOL)
colnames(DATA) <- XX

return(DATA)

})
