

#???ð?
library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

inputFile="all.txt"         #?????ļ?
gmtFile="c2.cp.kegg.v7.4.symbols.gmt"      #???????ļ?


#??ȡ?ļ?,?????????ļ?????????
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
logFC=as.vector(rt[,2])
names(logFC)=as.vector(rt[,1])
logFC=sort(logFC, decreasing=T)

#???????????ļ?
gmt=read.gmt(gmtFile)

#????????
kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$pvalue<0.05,]
write.table(kkTab,file="GSEA.result.txt",sep="\t",quote=F,row.names = F)

#????ʵ???鸻????ͼ??
termNum=5      #չʾͨ·????Ŀ
kkUp=kkTab[kkTab$NES>0,]
if(nrow(kkUp)>=termNum){
	showTerm=row.names(kkUp)[1:termNum]
	gseaplot=gseaplot2(kk, showTerm, base_size=8, title="Enriched in RSA")
	pdf(file="GSEA.RSA.pdf", width=7, height=5.5)
	print(gseaplot)
	dev.off()
}

#?????????鸻????ͼ??
termNum=5      #չʾͨ·????Ŀ
kkDown=kkTab[kkTab$NES<0,]
if(nrow(kkDown)>=termNum){
	showTerm=row.names(kkDown)[1:termNum]
	gseaplot=gseaplot2(kk, showTerm, base_size=8, title="Enriched in Control")
	pdf(file="GSEA.con.pdf", width=7, height=5.5)
	print(gseaplot)
	dev.off()
}

