
library(limma)
library(pheatmap)

inputFile="03.GSE165004_已处理完成.txt"       
logFCfilter=0.585            
adj.P.Val.Filter=0.05       



rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
sampleName1=c()
files=dir()
files=grep("S1.txt$", files, value=T)
for(file in files){
    rt=read.table(file, header=F, sep="\t", check.names=F)      
    geneNames=as.vector(rt[,1])      
    uniqGene=unique(geneNames)       
    sampleName1=c(sampleName1, uniqGene)
}


sampleName2=c()
files=dir()
files=grep("s2.txt$", files, value=T)
for(file in files){
    rt=read.table(file, header=F, sep="\t", check.names=F)      
    geneNames=as.vector(rt[,1])     
    uniqGene=unique(geneNames)       
    sampleName2=c(sampleName2, uniqGene)
conData=data[,sampleName1]
treatData=data[,sampleName2]
data=cbind(conData,treatData)
conNum=ncol(conData)
treatNum=ncol(treatData)

Type=c(rep("con",conNum),rep("treat",treatNum))
design <- model.matrix(~0+factor(Type))
colnames(design) <- c("con","treat")
fit <- lmFit(data,design)
cont.matrix<-makeContrasts(treat-con,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

allDiff=topTable(fit2,adjust='fdr',number=200000)
allDiffOut=rbind(id=colnames(allDiff),allDiff)
write.table(allDiffOut, file="all.txt", sep="\t", quote=F, col.names=F)

outData=rbind(id=paste0(colnames(data),"_",Type),data)
write.table(outData, file="normalize.txt", sep="\t", quote=F, col.names=F)

diffSig=allDiff[with(allDiff, (abs(logFC)>logFCfilter & P.Value < adj.P.Val.Filter )), ]
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut, file="diff.txt", sep="\t", quote=F, col.names=F)




diffGeneExp=data[row.names(diffSig),]
diffGeneExpOut=rbind(id=paste0(colnames(diffGeneExp),"_",Type), diffGeneExp)
write.table(diffGeneExpOut, file="diffGeneExp.txt", sep="\t", quote=F, col.names=F)




geneNum=50
diffSig=diffSig[order(as.numeric(as.vector(diffSig$logFC))),]
diffGeneName=as.vector(rownames(diffSig))
diffLength=length(diffGeneName)
hmGene=c()
if(diffLength>(2*geneNum)){
    hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
}else{
    hmGene=diffGeneName
}
hmExp=data[hmGene,]
Type=c(rep("Con",conNum),rep("Treat",treatNum))
names(Type)=colnames(data)
Type=as.data.frame(Type)
pdf(file="heatmap.pdf", width=10, height=8)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 8,
         fontsize_row=7,
         fontsize_col=8)
dev.off()


library(dplyr)
library(ggplot2)
library(ggrepel)

rtVol = read.table("all.txt", header=T, sep="\t", check.names=F)

Sig=ifelse((rtVol$P.Value<adj.P.Val.Filter) & (abs(rtVol$logFC)>logFCfilter), ifelse(rtVol$logFC>logFCfilter,"Up","Down"), "Not")


rtVol = mutate(rtVol, Sig=Sig)
p = ggplot(rtVol, aes(logFC, -log10(adj.P.Val)))+
  geom_point(aes(col=Sig))+
  scale_color_manual(values=c("green", "black","red"))+
  labs(title = " ")+
  theme(plot.title = element_text(size = 16, hjust = 0.5， face = "bold"))

p1=p+geom_label_repel(data=filter(rtVol, ((rtVol$adj.P.Val<adj.P.Val.Filter) & (abs(rtVol$logFC)>logFCfilter))),
                      box.padding=0.1, point.padding=0.1, min.segment.length=0.05,
                      size=1.8, aes(label=id)) + theme_bw()

pdf(file="vol.pdf", width=7, height=6.1)
print(p1)
dev.off()

