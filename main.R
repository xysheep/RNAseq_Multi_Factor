#setwd('C:/Users/xyang88/OneDrive/MyProject/RNAseq_CATTILL')
setwd('C:/Users/xysheep/OneDrive/MyProject/RNAseq_share')
source('func_RNAseq_report.R')
source('makedb.R')

# Change the following code to whatever you like
# Below: where to put the output. Remember the "/" at last
outdir<-"./output_xingyu/" 
# Below: Prefix of output files names. 
# Below: This example has three elements, if you have more comparison, you can add more
outnames<-c('Sub1 vs. Sub2','Responder vs. non-responder','Responder vs. non-responder on sub1 ')


source('Analysis1_DEG.R')
fdr = 0.1
source('Analysis2_Overlap.R')
source('output.R')
source('make_preranked.R')


colnames(norm_cnts) = rownames(coldata)
for (i in 1:length(q)){
  pdf(paste(outdir,outnames[i],'.pdf',sep=''))
  used_gene_idx1 <- (!is.na(q[[i]]) & q[[i]]<fdr)
  draw_samplecluster(used_gene_idx1,norm_cnts,outnames[i])
  draw_heatmap(used_gene_idx1,norm_cnts,outnames[i])
  dev.off()
}

