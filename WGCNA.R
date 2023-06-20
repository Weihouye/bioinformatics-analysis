
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install(c("GO.db", "preprocessCore", "impute","limma"))

#install.packages(c("matrixStats", "Hmisc", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival")) 
#install.packages("WGCNA")


#???ð?
library(limma)
library(WGCNA)
expFile="normalize.txt"      #?????????ļ?
setwd("C:\\biowolf\\wgcnaDiagnostic\\11.WGCNA")     #???ù???Ŀ¼

#??ȡ?????ļ????????????ļ?????
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
#data=data[apply(data,1,sd)>0.7,]     #ɾ??????С?Ļ???

#??ȡ????????ʵ????????Ʒ??Ϣ
Type=gsub("(.*)\\_(.*)", "\\2", colnames(data))
conCount=length(Type[Type=="con"])
treatCount=length(Type[Type=="treat"])
datExpr0=t(data)

###????ȱʧֵ
gsg = goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

###??Ʒ????
sampleTree = hclust(dist(datExpr0), method = "average")
pdf(file = "1_sample_cluster.pdf", width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
###??????
abline(h = 20000, col = "red")
dev.off()

###ɾ???????????µ???Ʒ
clust = cutreeStatic(sampleTree, cutHeight = 20000, minSize = 10)
table(clust)
keepSamples = (clust==1)
datExpr0 = datExpr0[keepSamples, ]


###׼???ٴ?????
traitData=data.frame(Con=c(rep(1,conCount),rep(0,treatCount)),
                     Treat=c(rep(0,conCount),rep(1,treatCount)))
row.names(traitData)=colnames(data)
fpkmSamples = rownames(datExpr0)
traitSamples =rownames(traitData)
sameSample=intersect(fpkmSamples,traitSamples)
datExpr0=datExpr0[sameSample,]
datTraits=traitData[sameSample,]


###??Ʒ????,?õ???Ʒ????????ͼ
sampleTree2 = hclust(dist(datExpr0), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE)
pdf(file="2_sample_heatmap.pdf",width=12,height=12)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

###powerֵɢ??ͼ
enableWGCNAThreads()   #???̹߳???
powers = c(1:20)       #??ָ????Χ1:20
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
pdf(file="3_scale_independence.pdf",width=9,height=5)
par(mfrow = c(1,2))
cex1 = 0.9
###????ָ????powerֵɢ??ͼ
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red") #?????޸?
###ƽ??��ͨ????powerֵɢ??ͼ
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

###?ڽӾ???ת??
sft #?鿴????powerֵ
softPower =sft$powerEstimate #????powerֵ
adjacency = adjacency(datExpr0, power = softPower)
softPower

#?޳߶??????????˹?ϵ
pdf(file="3_softConnectivity.pdf", width=9, height=5)
k <- softConnectivity(datE=datExpr0,power=softPower)
#sizeGrWindow(10, 5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k,main="Check Scale free topology\n")
dev.off()

###TOM????
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

###????????
geneTree = hclust(as.dist(dissTOM), method = "average");
pdf(file="4_gene_clustering.pdf",width=12,height=9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()


###??̬????ģ??ʶ??
minModuleSize = 60      #ģ????????Ŀ
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
pdf(file="5_Dynamic_Tree.pdf",width=8,height=6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()


###??ģ?????о???,?ҳ?????ģ??????
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average")
pdf(file="6_Clustering_module.pdf",width=7,height=6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.25  #???и߶ȿ??޸?
abline(h=MEDissThres, col = "red")
dev.off()


###????ģ???ϲ?
merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
pdf(file="7_merged_dynamic.pdf", width = 9, height = 6)
plotDendroAndColors(geneTree, mergedColors,"Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()
moduleColors = mergedColors
table(moduleColors)
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs


###ģ??????״??????ͼ
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
pdf(file="8_Module_trait.pdf", width=6, height=5.5)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(5, 10, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


###????MM??GSֵ
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
traitNames=names(datTraits)
geneTraitSignificance = as.data.frame(cor(datExpr0, datTraits, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep="")

###????ģ????Ҫ?Ե?ͼ??
y=datTraits[,1]
GS1=as.numeric(cor(y, datExpr0, use="p"))
GeneSignificance=abs(GS1)
ModuleSignificance=tapply(GeneSignificance, mergedColors, mean, na.rm=T)
pdf(file="9_GeneSignificance.pdf", width=11, height=7)
plotModuleSignificance(GeneSignificance, mergedColors)
dev.off()

###??��??????״??ģ??ɢ??ͼ
trait="Treat"
traitColumn=match(trait,traitNames)  
for (module in modNames){
	column = match(module, modNames)
	moduleGenes = moduleColors==module
	if (nrow(geneModuleMembership[moduleGenes,]) > 1){
	    outPdf=paste("10_", trait, "_", module,".pdf",sep="")
	    pdf(file=outPdf,width=7,height=7)
	    par(mfrow = c(1,1))
	    verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
	                         abs(geneTraitSignificance[moduleGenes, traitColumn]),
	                         xlab = paste("Module Membership in", module, "module"),
	                         ylab = paste("Gene significance for ",trait),
	                         main = paste("Module membership vs. gene significance\n"),
	                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
	    abline(v=0.8,h=0.5,col="red")
	    dev.off()
	}
}


###????GS_MM????
probes = colnames(datExpr0)
geneInfo0 = data.frame(probes= probes,
                       moduleColor = moduleColors)
for (Tra in 1:ncol(geneTraitSignificance))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneTraitSignificance[,Tra],
                         GSPvalue[, Tra])
  names(geneInfo0) = c(oldNames,names(geneTraitSignificance)[Tra],
                       names(GSPvalue)[Tra])
}

for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[,mod],
                         MMPvalue[, mod])
  names(geneInfo0) = c(oldNames,names(geneModuleMembership)[mod],
                       names(MMPvalue)[mod])
}
geneOrder =order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]
write.table(geneInfo, file = "GS_MM.xls",sep="\t",row.names=F)


###????ÿ??ģ???Ļ???
for (mod in 1:nrow(table(moduleColors)))
{  
  modules = names(table(moduleColors))[mod]
  probes = colnames(datExpr0)
  inModule = (moduleColors == modules)
  modGenes = probes[inModule]
  write.table(modGenes, file =paste0("module_",modules,".txt"),sep="\t",row.names=F,col.names=F,quote=F)
}

###????ÿ??ģ???ĺ??Ļ???
geneSigFilter=0.5         #??????Ҫ?ԵĹ???????
moduleSigFilter=0.8       #??????ģ???????ԵĹ???????
datMM=cbind(geneModuleMembership, geneTraitSignificance)
datMM=datMM[abs(datMM[,ncol(datMM)])>geneSigFilter,]
for(mmi in colnames(datMM)[1:(ncol(datMM)-2)]){
	dataMM2=datMM[abs(datMM[,mmi])>moduleSigFilter,]
	write.table(row.names(dataMM2), file =paste0("hubGenes_",mmi,".txt"),sep="\t",row.names=F,col.names=F,quote=F)
}


######Vi