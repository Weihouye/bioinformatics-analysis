


library(e1071)
library(kernlab)
library(caret)

inputFile="vennGeneExp.txt"        



data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))

#SVM-RFE
Profile=rfe(x=data,
            y=as.numeric(as.factor(group)),
            sizes = c(2,4,6,8, seq(10,40,by=3)),
            rfeControl = rfeControl(functions = caretFuncs, method = "cv"),
            methods="svmRadial")


pdf(file="SVM-RFE_Accuracy.pdf", width=6, height=5.5)
par(las=1)
x = Profile$results$Variables
y = Profile$results$RMSE
plot(x, y, xlab="Number of Features", ylab="5 X CV Accuracy", col="darkgreen")
lines(x, y, col="darkgreen")
wmin=which.min(y)
wmin.x=x[wmin]
wmin.y=y[wmin]
dev.off()


pdf(file="SVM-RFE_Error.pdf", width=6, height=5.5)
par(las=1)
x = Profile$results$Variables
y = 1-Profile$results$RMSE
plot(x, y, xlab="Number of Features", ylab="5 X CV Error", col="darkgreen")
lines(x, y, col="darkgreen")
wmin=which.min(y)
wmin.x=x[wmin]
wmin.y=y[wmin]
dev.off()



featureGenes=Profile$optVariables
write.table(file="SVM-RFE.gene.txt", featureGenes, sep="\t", quote=F, row.names=F, col.names=F)




