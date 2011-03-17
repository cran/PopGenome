setClass("region.stats", representation(

######## Pops for each Test Module

Pop_Neutrality  =  "list",
Pop_FSTN        =  "list",
Pop_FSTH        =  "list",
Pop_Linkage     =  "list",
Pop_MK          =  "list",


genename               = "list",
nucleotide.diversity   = "list",
haplotype.diversity    = "list",
haplotype.counts       = "list" ,      # sfreqh
minor.allele.freqs     = "list",       # JFD
biallelic.structure    = "list",       # SXX
linkage.disequilibrium = "list"       # Link

))


#"THETA","S (THETA)","T (THETA)","SA (THETA)","TA (THETA)","FL (THETA)","L (THETA)","FW (THETA)"
#"FREQUENCY (THETA) VALUES","Watterson Frequency","Tajima Frequency","Achaz, Watterson Frequency","Achaz, Tajima Frequency","Fu & Li  
# Frequency","Zeng Frequency","Fay & Wu Frequency",


#### SHOW ######
setMethod("show", "region.stats",
 function(object){
cat("-----\n")
 cat("SLOTS:\n")
 cat("-----\n")
 out <- data.frame(Slots=c("nucleotide.diversity","haplotype.diversity","haplotype.counts", "minor.allele.freqs","linkage.disequilibrium","biallelic.structure"),         
               

                Description=c(
               " nucleotide diversity","haplotype diversity","haplotype distribution","minor allele frequencies","linkage disequilibrium","shared and fixed polymorphisms"),

               Module=c("FST","FST","FST","Neutrality","Linkage","Detail"))
 
              


 print(out)
 cat("\n---------------\n")
 print("These are the Slots (class region.data)") 
})




