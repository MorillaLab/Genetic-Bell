
#load('/Users/morilla/Desktop/Wainrib/LAGA/Bichat/Jean-Pierre/Epistasis/Bell/data/MRCA/correlation_analysis_MRCA_Bell.RData')
load('/Users/morilla/Desktop/Wainrib/LAGA/Bichat/Jean-Pierre/Epistasis/Bell/Algorithm/model/model_Ian/model1/data/liming_expression_rawdata/MRCA/correlation_analysis_MRCA_Bell.RData')

rm(res)

y <-log2(m)

# EM clustering for gaussian model mixture
library(mclust)

# fixing number of columns!
NCOL <- ncol(m)

# pre-allocating peaks
Npeaks_EM <- mat.or.vec(NCOL,1)
system.time(
for(i in 1:NCOL){
    x.gmm<-Mclust(y[,i])
    Npeaks_EM[i] <- length(x.gmm$parameters$mean)
    cat(i,"\n")
    #x.gmm <- NULL
    
}
)

save(Npeaks_EM, file="Npeaks_MRCA_Bell_EM.RData")

interesting <- which(Npeaks_EM==2)

s1<-sample(interesting,8)
s <- append(s1,641)


# depicting some few examples maybe significants!
pdf("some_examples_EM_top_mix_equal_2_MRCA.pdf")

par(mfrow=c(3,3))

for(i in s){
 hist(y[,i], prob=T, breaks=40)
 lines(density(y[,i]), col="blue", lwd=2)
}

#hist(data[,25809], prob=T)
#lines(density(data[,25809]), col="blue", lwd=2)

#hist(data[,22004], prob=T)
#lines(density(data[,22004]), col="blue", lwd=2)

#hist(data[,25851], prob=T)
#lines(density(data[,25851]), col="blue", lwd=2)

#hist(data[,22397], prob=T)
#lines(density(data[,22397]), col="blue", lwd=2)

#hist(data[,11313], prob=T)
#lines(density(data[,11313]), col="blue", lwd=2)

#hist(data[,4087], prob=T)
#lines(density(data[,4087]),  col="blue", lwd=2)

#hist(data[,26546], prob=T)
#lines(density(data[26546,]), col="blue", lwd=2)

#hist(data[,21191], prob=T)
#lines(density(data[,21191]), col="blue", lwd=2)

dev.off()





pdf("some_examples_MRCA_EM_top_mix_above_2.pdf")

interesting <- which(Npeaks_EM>2)

s<-sample(interesting,9)


par(mfrow=c(3,3))

for(i in s){
    hist(y[,i], prob=T, breaks=40, main=colnames(y)[i])
    lines(density(y[,i]), col="blue", lwd=2)
}

dev.off()




