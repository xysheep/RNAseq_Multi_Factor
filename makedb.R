makedb<-function(dbfile){
  x<-scan(dbfile,what="",sep="\n",quiet=T)
  y <- strsplit(x, "\t")
  names(y)<-sapply(y,`[[`,1)
  y<-lapply(y,`[`,c(-1,-2))
}
makedbinfo<-function(dbfile){
  x<-scan(dbfile,what="",sep="\n",quiet=T)
  y <- strsplit(x, "\t")
  names(y)<-sapply(y,`[[`,1)
  y<-lapply(y,`[`,2)
}

hdb = makedb('GENESETS/h.all.v5.1.entrez.gmt')
hdbinfo = makedbinfo('GENESETS/h.all.v5.1.entrez.gmt')

c1db = makedb('GENESETS/c1.all.v5.0.entrez.gmt')
c1dbinfo = makedbinfo('GENESETS/c1.all.v5.0.entrez.gmt')

c2db = makedb('GENESETS/c2.all.v5.0.entrez.gmt')
c2dbinfo = makedbinfo('GENESETS/c2.all.v5.0.entrez.gmt')

c7db = makedb('GENESETS/c7.all.v5.0.entrez.gmt')
c7dbinfo = makedbinfo('GENESETS/c7.all.v5.0.entrez.gmt')

data(kegg.gs)
kegg.info<-as.list(substr(names(kegg.gs),10,100))
names(kegg.gs)<-substr(names(kegg.gs),1,8)
names(kegg.info)<-names(kegg.gs)

c3tftdb = makedb('GENESETS/c3.tft.v5.0.entrez.gmt')
c3tftdbinfo = makedbinfo('GENESETS/c3.tft.v5.0.entrez.gmt')

c3mirdb = makedb('GENESETS/c3.mir.v5.0.entrez.gmt')
c3mirdbinfo = makedbinfo('GENESETS/c3.mir.v5.0.entrez.gmt')




