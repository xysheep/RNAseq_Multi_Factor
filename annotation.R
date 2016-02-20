setwd("C:/Users/xyang88/OneDrive/MyProject/RNAseq_HARA")
# set the above path to whatever you want to store your results



library(GenomicAlignments)
library(rtracklayer)
fls<-list.files("C:/Users/xyang88/Downloads/tophat_all_HA_RA_new",pattern="bam$",full.names=T)
# Set the above to path to where you store your bam files

bamfls<-BamFileList(fls)
flag<-scanBamFlag(isSecondaryAlignment = FALSE, isProperPair = TRUE)
param<-ScanBamParam(flag=flag)


library(TxDb.Hsapiens.UCSC.hg19.knownGene)
exByGn <- exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, "gene")

gnCnt<-summarizeOverlaps(exByGn,bamfls,mode="Union",ignore.strand=TRUE,single.end=FALSE,param=param)

save(gnCnt,file='Cnt.RData')

# A file named "gnCnt.RData" will be generated.
