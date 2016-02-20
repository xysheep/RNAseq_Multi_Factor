   
dir.create(outdir,showWarnings=FALSE)
dir.create('./KEGGview/',showWarnings=FALSE)
dir.create(paste(outdir,'KEGG',sep=''),showWarnings=FALSE)
id2symbols<-data.frame(getSYMBOL(rownames(d),data='org.Hs.eg'))

for (i in 1:length(q)){
  index.deg<- q[[i]]<fdr & !is.na(q[[i]])
  index.sort<-order(q[[i]][index.deg])
  # output lists of DEGs
  outresults=cbind(rownames(id2symbols),id2symbols,q[[i]],exp.fc[[i]])
  names(outresults)<-c('EntrezID','Symbol','q-value','log2_fold_change')
  write.table(outresults[index.deg,][index.sort,],
              paste(outdir,
                    paste('DEGs',outnames[i],'(entrez gene id).csv'),
                    sep=''),
              row.names=F,sep=',')
  # output DE pathways 
  foldchange=as.vector(exp.fc[[i]])
  foldchange[is.na(foldchange)]=0
  updegids=rownames(id2symbols)[index.deg & foldchange>0]
  downdegids=rownames(id2symbols)[index.deg & foldchange<0]
  
  output_deg_path(updegids,downdegids,hdb,h_path_msig[[i]],'Msigdb Hallmark')
  output_deg_path(updegids,downdegids,c2db,c2_path_msig[[i]],'Msigdb C2')
  output_deg_path(updegids,downdegids,c7db,c7_path_msig[[i]],'Msigdb C7')
  output_deg_path(updegids,downdegids,kegg.gs,kegg_path_msig[[i]],'KEGG')
  output_deg_path(updegids,downdegids,c3tftdb,c3tft_path_msig[[i]],'c3tft')
  
  # output KEGG visualization
  if (length(kegg_path_msig[[i]][[1]])>0){
    sapply(kegg_path_msig[[i]][[1]],function(pid) pathview(gene.data=score[[i]],
                                                           pathway.id = pid,species='hsa',
                                                           out.suffix=outnames[i],
                                                           kegg.dir='./KEGGview/'))
    out_k_f=paste(kegg_path_msig[[i]][[1]],outnames[i],'png',sep='.')
    tmpoutdir = paste(outdir,'KEGG',outnames[i],sep='/')
    dir.create(tmpoutdir,showWarnings=FALSE,recursive = T)
    file.copy(out_k_f,tmpoutdir)
    file.remove(out_k_f)
  }
}
