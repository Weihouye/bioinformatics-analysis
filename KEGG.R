

library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library("GOplot")

pvalueFilter=0.05       #pֵ????????
qvalueFilter=1    #????????pֵ????????

#??????ɫ
colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}
	
 
rt=read.table("diff.txt", header=T, sep="\t", check.names=F)      #??ȡ?????ļ?

#????????ת??Ϊ????id
colnames(rt)[1]="Gene"
genes=as.vector(rt[,1])
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=cbind(rt,entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        #ȥ??????idΪNA?Ļ???
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#kegg????????
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
#去除后缀
kk@result$Description<- gsub("- Homo sapiens \\(human\\)", "",kk@result$Description)
KEGG=as.data.frame(kk)

KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$Gene[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)


showNum=20
if(nrow(KEGG)<showNum){
	showNum=nrow(KEGG)
}

#??״ͼ
pdf(file="KEGG_barplot.pdf", width=10, height=8)
barplot(kk, drop = TRUE, showCategory = showNum, color = colorSel)+
  theme(legend.position = "right",
        legend.key.size = unit(0.7, "cm"),
        plot.margin = unit(c(1, 1, 1, 6), "cm"))
dev.off()

#????ͼ
pdf(file="KEGG_bubble.pdf", width=10, height=8)
dotplot(kk, showCategory = showNum, orderBy = "GeneRatio",color = colorSel)+
  theme(legend.position = "right",
        legend.key.size = unit(0.7, "cm"),
        plot.margin = unit(c(1, 1, 1, 6), "cm"))
dev.off()


#??ȡͨ·??Ϣ
kegg=data.frame(Category="ALL", ID = KEGG$ID, Term=KEGG$Description, Genes = gsub("/", ", ", KEGG$geneID), adj_pval = KEGG$p.adjust)
#??ȡ?????Ĳ?????Ϣ
genelist <- data.frame(ID=rt$Gene, logFC=rt$logFC)
row.names(genelist)=genelist[,1]
#????Ȧͼ????
circ <- circle_dat(kegg, genelist)
termNum =8       #??ʾͨ·????Ŀ
termNum=ifelse(nrow(kegg)<termNum,nrow(kegg),termNum)
geneNum=200      #??ʾ????????Ŀ
geneNum=ifelse(nrow(genelist)<geneNum, nrow(genelist), geneNum)
#????ͨ·??Ȧͼ
chord <- chord_dat(circ, genelist[1:geneNum,], kegg$Term[1:termNum])
pdf(file="KEGG_circos.pdf", width=13, height=12)
GOChord(chord, 
        space = 0.01,           
        gene.order = 'logFC',    
        gene.space = 0.2,       
        gene.size = 3,           
        border.size = 0.01,       
        process.label = 9)       
dev.off()

#ͨ·?ľ???ͼ
pdf(file="KEGG_cluster.pdf",width=12, height=10)
GOCluster(circ, 
          kegg$Term[1:termNum], 
          lfc.space = 0.2,        
          lfc.width = 1,          
          term.space = 0.2,       
          term.width = 1)         
dev.off()          




