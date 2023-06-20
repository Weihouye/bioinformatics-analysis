
library(neuralnet)
library(NeuralNetTools)
inputFile="geneScore.txt"      
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
data=as.data.frame(t(data))
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))
data$con=ifelse(group=="con", 1, 0)
data$treat=ifelse(group=="treat", 1, 0)


fit=neuralnet(con+treat~., data, hidden=5)
fit$result.matrix
fit$weight
#plot(fit)

pdf(file="neuralnet.pdf", width=9, height=7)
plotnet(fit)
dev.off()


net.predict=compute(fit, data)$net.result
net.prediction=c("con", "treat")[apply(net.predict, 1, which.max)]
predict.table=table(group, net.prediction)
predict.table
conAccuracy=predict.table[1,1]/(predict.table[1,1]+predict.table[1,2])
treatAccuracy=predict.table[2,2]/(predict.table[2,1]+predict.table[2,2])
paste0("Con accuracy: ", sprintf("%.3f", conAccuracy))
paste0("Treat accuracy: ", sprintf("%.3f", treatAccuracy))


colnames(net.predict)=c("con", "treat")
outTab=rbind(id=colnames(net.predict), net.predict)
write.table(outTab, file="neural.predict.txt", sep="\t", quote=F, col.names=F)


