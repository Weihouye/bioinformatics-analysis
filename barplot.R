######Video source: https://ke.biowolf.cn
######??????ѧ??: https://www.biowolf.cn/
######΢?Ź??ںţ?biowolf_cn
######???????䣺biowolf@foxmail.com
######????΢??: 18520221056

#install.packages("corrplot")


#???ð?
library(corrplot)
inputFile="CIBERSORT-Results.txt"     #????ϸ???????Ľ????ļ?
#setwd("C:\\biowolf\\Diagnostic\\19.barplot")     #???ù???Ŀ¼

#??ȡ????ϸ???????ļ?
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)

#????Ʒ????
con=grepl("_con", rownames(rt), ignore.case=T)
treat=grepl("_treat", rownames(rt), ignore.case=T)
conData=rt[con,]
treatData=rt[treat,]
conNum=nrow(conData)
treatNum=nrow(treatData)
data=t(rbind(conData,treatData))

#??????״ͼ
pdf(file="barplot2.pdf", width=16, height=7.5)
col=rainbow(nrow(data), s=0.7, v=0.7)
par(las=1,mar=c(8,5,4,16),mgp=c(3,0.1,0),cex.axis=1.5)
a1=barplot(data,col=col,yaxt="n",ylab="Relative Percent",xaxt="n",cex.lab=1.8)
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
par(srt=0,xpd=T)
rect(xleft = a1[1], ybottom = -0.01, xright = a1[conNum], ytop = -0.06,col="green")
text(a1[conNum]/2,-0.035,"Con",cex=1.8)
rect(xleft = a1[conNum], ybottom = -0.01, xright =a1[length(a1)] , ytop = -0.06,col="red")
text((a1[length(a1)]+a1[conNum])/2,-0.035,"Treat",cex=1.8)
ytick2 = cumsum(data[,ncol(data)])
ytick1 = c(0,ytick2[-length(ytick2)])
legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data),col=col,pch=15,bty="n",cex=1.2)
dev.off()

#?????Է???
pdf("corHeatmap.pdf", width=12, height=12)
corrplot(corr=cor(rt),
         method = "color",
         order = "hclust",
         tl.col="black",
         addCoef.col = "black",
         number.cex = 0.8,
         col=colorRampPalette(c("blue", "white", "red"))(50),
         )
dev.off()


######Video source: https://ke.biowolf.cn
######??????ѧ??: https://www.biowolf.cn/
######΢?Ź??ںţ?biowolf_cn
######???????䣺biowolf@foxmail.com
######????΢??: 18520221056

