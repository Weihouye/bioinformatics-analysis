

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


library(limma)              #???ð?
expFile="4GeneExp.txt"     #?????????ļ?
diffFile="diff.txt"         #???????????ļ?
setwd("C:\\biowolf\\neuralDiagnostic\\14.geneScore")      #???ù???Ŀ¼

#??ȡ?????ļ????????????ļ?????
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#??ȡ???????????ļ?
diffRT=read.table(diffFile, header=T, sep="\t", check.names=F, row.names=1)
diffRT=diffRT[row.names(data),]

#?????��?
dataUp=data[diffRT[,"logFC"]>0,]
dataDown=data[diffRT[,"logFC"]<0,]
dataUp2=t(apply(dataUp,1,function(x)ifelse(x>median(x),1,0)))
dataDown2=t(apply(dataDown,1,function(x)ifelse(x>median(x),0,1)))

#?????????��ֵĽ???
outTab=rbind(dataUp2, dataDown2)
outTab=rbind(id=colnames(outTab), outTab)
write.table(outTab, file="geneScore.txt", sep="\t", quote=F, col.names=F)


######Vi