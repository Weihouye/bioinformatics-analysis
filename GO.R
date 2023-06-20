
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library(GOplot)

pvalueFilter=0.05       #pֵ????????
qvalueFilter=1      #????????pֵ????????

#??????ɫ
colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}

#setwd("C:\\biowolf\\neuralDiagnostic\\09.GO")       #???ù???Ŀ¼
rt=read.table("diff.txt", header=T, sep="\t", check.names=F)     #??ȡ?????ļ?

#????????ת??Ϊ????id
colnames(rt)[1]="Gene"
genes=as.vector(rt[,1])
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        #ȥ??????idΪNA?Ļ???
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#GO????????
kk=enrichGO(gene=gene,OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable =T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
#?????????????Ľ???
write.table(GO,file="GO.txt",sep="\t",quote=F,row.names = F)

#??????ʾTerm??Ŀ
showNum=5
if(nrow(GO)<30){
	showNum=nrow(GO)
}

library(ggplot2)
pdf(file="barplot.pdf", width=10, height=8)
bar <- barplot(kk, drop = TRUE, showCategory =showNum, split="ONTOLOGY", color = colorSel) + 
  facet_grid(ONTOLOGY~., scale='free') +
  theme(legend.position = "right",
        legend.key.size = unit(0.7, "cm"),
        plot.margin = unit(c(1, 1, 1, 6), "cm"))
print(bar)
dev.off()


		
#????ͼ
pdf(file="bubble.pdf", width=10, height=8)
bub=dotplot(kk,showCategory = showNum, orderBy = "GeneRatio",split="ONTOLOGY", color = colorSel) +
  facet_grid(ONTOLOGY~., scale='free')+
  theme(legend.position = "right",
        legend.key.size = unit(0.7, "cm"),
        plot.margin = unit(c(1, 1, 1, 6), "cm"))
print(bub)
dev.off()


#??ȡGO??Ϣ
go=data.frame(Category=GO$ONTOLOGY, ID=GO$ID, Term=GO$Description, Genes = gsub("/", ", ", GO$geneID), adj_pval = GO$p.adjust)
#??ȡ??????logFC
genelist <- data.frame(ID=rt$Gene, logFC=rt$logFC)
row.names(genelist)=genelist[,1]
#????Ȧͼ????
circ <- circle_dat(go, genelist)
termNum =8       #??ʾGO??Ŀ
termNum=ifelse(nrow(go)<termNum,nrow(go),termNum)
geneNum=200      #?޶???????Ŀ
geneNum=ifelse(nrow(genelist)<geneNum, nrow(genelist), geneNum)
#????Ȧͼ
chord <- chord_dat(circ, genelist[1:geneNum,], go$Term[1:termNum])
pdf(file="GOcircos.pdf", width=12, height=12)
GOChord(chord, 
        space = 0.001,           #????֮???ļ???
        gene.order = 'logFC',    #????logFCֵ?Ի???????
        gene.space = 0.25,       #????????ԲȦ?????Ծ???
        gene.size = 5,           #????????????С 
        border.size = 0.1,       #??????ϸ
        process.label = 6)       #GO??????С
dev.off()

#????ͼ
pdf(file="GOcluster.pdf",width=12, height=10)
GOCluster(circ, 
          go$Term[1:termNum], 
          lfc.space = 0.2,        #logFC????֮???Ŀ?϶??С
          lfc.width = 1,          #logFC??ԲȦ????
          term.space = 0.2,       #logFC??GO֮????϶?Ĵ?С
          term.width = 1)         #GOԲȦ?Ŀ???
dev.off()          




