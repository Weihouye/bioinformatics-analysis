
##加载R包
library(GEOquery)
library(Biobase)
library(limma)  
library(dplyr)
library(tibble)
library(affy)
library(dplyr)
dir.create("sampFile")
samPath = "./sampFile"  
list.files()
untar(list.files()[i], exdir = samPath)
dat = list.files(samPath, pattern = "txt")dat
dat=read.maimages(files= dat,path = "./sampFile",green.only = TRUE,source="agilent")
ep=dat$E
qx <- as.numeric(quantile(ep, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) { ep[which(ep <= 0)] <- NaN
ep <- log2(ep+1)
print("log2 transform finished")}else{print("log2 transform not needed")}
par(cex = 0.9)
cols <- rainbow(length(ep) * 1.2)
boxplot(ep, col = cols, xlab = "Sample", ylab = "Log intensity")
RG = backgroundCorrect(dat, "normexp", normexp.method = "rma", offset=50) 
E = normalizeBetweenArrays(RG, method="quantile")
plotMA(dat,array = 1)
plotMA(RG,array = 1)
plotMA(E,array = 1)
plotDensities(dat)
plotDensities(RG)
plotDensities(E)
ep=as.data.frame( E )
ep=ep[,6:length(colnames(ep))]
par(cex = 0.9)
cols <- rainbow(length(ep) * 1.2)
boxplot(ep, col = cols, xlab = "Sample", ylab = "Log intensity")
expr <- as.data.frame( E )
expr$ID <- E[["Genes"]]$ProbeName
GPL=read.table("GPL16699-15607.txt",
                  header = TRUE,fill = T,sep = "\t",
                  comment.char = "#",
                  stringsAsFactors = FALSE,
                  quote = "")
colnames(GPL)
annotation <- as.data.frame(GPL) %>%
  dplyr::select("NAMES","GENE_SYMBOL")
colnames(annotation) 

colnames(expr)
colnames(annotation)[1]="ID"

expr2 <- merge(annotation,expr, by='ID')
expr2=expr2[,-c(1,3,4,5,6,7)]
colnames(expr2)
m=colnames(expr2)[2:length(colnames(expr2))]
m= substr(m, 1, 10) 
colnames(expr2)[2:length(colnames(expr2))]=m  #反过来修改列名
head(expr2)
colnames(expr2)
exprSet2 <- aggregate(x = expr2[,2:ncol(expr2)],
                      by = list(expr2$GENE),
                      FUN = mean)

exprSet6=data.frame(annotation[match(exprSet2$Group.1,annotation$GENE),],exprSet2)  
head(exprSet6)
exprSet6=exprSet6[,-c(1)]
colnames(exprSet6)[2:length(colnames(exprSet6))]=
head(exprSet6)
write.table(exprSet6, file="GSE165004.txt", sep="\t", quote=F, col.names=T)





