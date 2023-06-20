


library(randomForest)
 
inputFile="GeneExp.txt"       



data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))


rf=randomForest(as.factor(group)~., data=data, ntree=1000)
pdf(file="forest.pdf", width=6, height=6)
plot(rf, main=" ", lwd=2)
dev.off()


optionTrees=which.min(rf$err.rate[,1])
optionTrees
rf2=randomForest(as.factor(group)~., data=data, ntree=optionTrees)

#查看基因的重要性
importance=importance(x=rf)

#绘制基因的重要性图
pdf(file="geneImportance.pdf", width=6.2, height=5.8)
varImpPlot(rf2, main="")
dev.off()

#挑选疾病特征基因
rfGenes=importance[order(importance[,"MeanDecreaseGini"], decreasing = TRUE),]
rfGenes=names(rfGenes[rfGenes>1])    
#rfGenes=names(rfGenes[1:30])         
write.table(rfGenes, file="rfGenes.txt", sep="\t", quote=F, col.names=F, row.names=F)




