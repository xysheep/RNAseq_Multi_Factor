library(gage)
library(preprocessCore)
library(GenomicAlignments)
library(DESeq2)
library(Rtsne)
library(colorRamps)
library("RColorBrewer")
library(org.Hs.eg.db)
library(annotate)
library(pathview)
library(rtracklayer)


draw_samplecluster<-function(used_gene_idx,norm_cnts,title){
  deg_norm_cnts<-norm_cnts[used_gene_idx,]
  deg_norm_cnts<-data.frame(deg_norm_cnts)
  
  ## PCA
  pca_cnts<-prcomp(deg_norm_cnts,na.action=na.omit)
  pca<-pca_cnts$rotation
  text_label<-c(names(deg_norm_cnts))
  plot(pca[,1],pca[,2],cex=3,pch=22,main=title,xlab='PC1',ylab='PC2')
  text(pca[,1],pca [,2]-0.025,text_label,cex=0.75)
  
  ## Rtsne
  names(deg_norm_cnts)<-text_label
  set.seed(9464)
  tsne_results1<-Rtsne(t(deg_norm_cnts),perplexity=2)
  tsne_cnts1<-tsne_results1$Y
  plot(tsne_cnts1[,1],tsne_cnts1[,2],cex=3,pch=22,main='tSNE',xlab='tSNE1',ylab='tSNE2')
  text(tsne_cnts1[,1],tsne_cnts1[,2]-0.025,text_label,cex=0.75)
  
}

draw_heatmap<-function(used_gene_idx,norm_cnts,title){
  
  hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
  hmcol <- blue2yellow(100)
  temp_cnts<-norm_cnts[used_gene_idx,]
  #select <- order(apply(temp_cnts,1,median),decreasing=TRUE)[1:nrow(temp_cnts)]
  heatmap(temp_cnts#[select,]
          ,Colv=NA,#Rowv=NA,#
          scale="none",na.rm=T,col=hmcol,
          labCol=gsub('ub','',gsub('onor','',colnames(norm_cnts))),#c(paste('HA',1:4,sep=''),paste(rep(c('RA','RActl'),c(4,4)),rep(1:4,2),sep='')),
          labRow='',margins=c(5,5),xlab='Sample',ylab='Genes',main=title)
}

analysis_depath<-function(score,db,fdr,dbinfo){
  #gage
  fc.kegg.p <- gage(score, gsets = db, ref = NULL, samp = NULL)#subject to change
  sel <- fc.kegg.p$greater[, "q.val"] < fdr & !is.na(fc.kegg.p$greater[, "q.val"])
  path.ids <- cbind(rownames(fc.kegg.p$greater)[sel],fc.kegg.p$greater[sel, "q.val"],rep('up',sum(sel)))
  sel.l <- fc.kegg.p$less[, "q.val"] < fdr &  !is.na(fc.kegg.p$less[,"q.val"])
  path.ids.l <- cbind(rownames(fc.kegg.p$less)[sel.l],fc.kegg.p$less[sel.l, "q.val"],rep('down',sum(sel.l)))
  path.ids.all <- rbind(path.ids, path.ids.l)
  data.frame(cbind(path.ids.all,unlist(dbinfo[path.ids.all[,1]])))
}


analysis_deg<-function(cnts,subtype,donor){
  #DEG

  coldata=DataFrame(subtype=factor(subtype),donor=factor(donor))
  dds <- DESeqDataSetFromMatrix(cnts, colData=coldata, design = ~ subtype + donor)
  dds<-DESeq(dds)
  dds
}

analysis_deg_onefactor<-function(cnts,responder){
  #DEG
  
  coldata=DataFrame(responder=responder)
  dds <- DESeqDataSetFromMatrix(cnts, colData=coldata, design = ~ responder)
  dds<-DESeq(dds)
  dds
}

output_deg_path<-function(updegids,downdegids,db,depath,str){
  uppathdegs<-find_deg_inpath(db,depath,updegids)  
  downpathdegs<-find_deg_inpath(db,depath,downdegids)
  outresults=cbind(depath,uppathdegs[[1]],uppathdegs[[2]],downpathdegs[[1]],downpathdegs[[2]])
  if (nrow(depath)>0){
    names(outresults)=c('pathway name','q-value','up/down','information','up DEGs id','up DEGs symb','down DEGs id','down DEGs symb')
    }
  write.table(outresults,paste(outdir,
                               paste('D. pathways',outnames[i],'(',str,').csv'),
                               sep=''),
              row.names=F,sep=',')
}
find_deg_inpath<-function(db,paths,gids){
  l=lapply(db[as.character(paths[,1])],intersect,gids)
  s1=sapply(l,paste,collapse = ' ')
  s2=sapply(lapply(l,(function (x) as.character(id2symbols[x,]))),paste,collapse = ' ')
  list(s1,s2)
}
cal_score<-function(fc,qv){
  s=1-qv;
  s[is.na(s)]=0
  s[!is.na(fc) & fc<0] = -s[!is.na(fc) & fc<0];
  s
}
list_score<-function(exp.fc,q){
  s = q
  for(i in 1:length(q)){
    s[[i]] =  cal_score(exp.fc[[i]],q[[i]])
  }
  s
}