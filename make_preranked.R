outdir_glist<-paste(outdir,"preranked_genelist/up/",sep='')
dir.create(outdir_glist,showWarnings=FALSE,recursive = T)
for (i in 1:length(q)){
  d=cbind(id2symbols,data.frame(score[[i]]))
  write.table(d,paste(outdir_glist,
                      paste('Ranked_Gene_List',outnames[i],'.rnk'),
                      sep=''),row.names=F,col.names=F,sep='\t',quote=F)
}

outdir_glist<-paste(outdir,"preranked_genelist/down/",sep='')
dir.create(outdir_glist,showWarnings=FALSE,recursive = T)
for (i in 1:length(q)){
  d=cbind(id2symbols,data.frame(-score[[i]]))
  write.table(d,paste(outdir_glist,
                      paste('Ranked_Gene_List',outnames[i],'.rnk'),
                      sep=''),row.names=F,col.names=F,sep='\t',quote=F)
}