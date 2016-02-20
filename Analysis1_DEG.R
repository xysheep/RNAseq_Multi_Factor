
load('Cnt.RData')
d = assay(gnCnt)
d = d[rowSums(d)>1,]
# Combine technical replication separately
# if no technical replication, no need to run the following 5 lines.
shortd = d[,1:12]
for (i in 1:12){
  shortd[,i] = rowSums(d[,((i-1)*4+1):(i*4)],na.rm=T)
}
d=shortd

coldata <- DataFrame(read.csv('design.csv',header=T))# prepare a design table
#############################################################Below
# You need to change the following code!!!!!!!!!!!!!!!!!!!!!!!!!!!

# Subtypes1 vs. subtypes2
dds1 <- DESeqDataSetFromMatrix(d, colData=coldata, design = ~ subtype + donor)
dds1<-DESeq(dds1)
l <- list(results(dds1,contrast=c('subtype','Sub1','Sub2')))

# Responder vs. non-responder
dds2 <- DESeqDataSetFromMatrix(d, colData=coldata, design = ~ subtype + responder)
dds2<-DESeq(dds2)
l[length(l)+1] <- results(dds2,contrast=c('responder','yes','no'))

# Responder vs. non-responder on sub1 
dds_rnor <- DESeqDataSetFromMatrix(d[,coldata$subtype=='Sub1'], 
                                   colData=coldata[coldata$subtype=='Sub1',], design = ~ responder)
dds_rnor<-DESeq(dds_rnor)
l[length(l)+1] <- results(dds_rnor)
####################################################################Above

q = lapply(l,function(x) data.frame(x['padj']))

exp.fc = lapply(l,function(x) data.frame(x['log2FoldChange']))




fdr = 0.1
score = list_score(exp.fc,q);

h_path_msig<-lapply(score,analysis_depath,hdb,fdr,hdbinfo)
c2_path_msig<-lapply(score,analysis_depath,c2db,fdr,c2dbinfo)
c7_path_msig<-lapply(score,analysis_depath,c7db,fdr,c7dbinfo)
kegg_path_msig<-lapply(score,analysis_depath,kegg.gs,fdr,kegg.info)
c3tft_path_msig<-lapply(score,analysis_depath,c3tftdb,fdr,c3tftdbinfo)


norm_cnts <- log10(counts(dds1,normalized=T))
norm_cnts[norm_cnts<0] = 0




