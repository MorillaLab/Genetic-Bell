# Loading files into the workspace
#filesnames <- list.files("/Users/morilla/Desktop/Wainrib/LAGA/Bichat/Jean-Pierre/Epistasis/Bell/Algorithm/model/model_Ian/model1/Oyy_test_persample", pattern="*BP.txt", full.names=T)
load('/Users/morilla/Desktop/Wainrib/LAGA/Bichat/Jean-Pierre/Epistasis/Bell/Algorithm/model/model_Ian/model1/Oyy_test_persample/data/yliming_sGGMS.RData')

mlist <- list(m.1, m.2, m.3, m.4)
mlist <- lapply(mlist, function(x){rownames(x)<-paste0("pt",1:395); x})
mlist.scaled <- lapply(mlist, scale) #

library(mclust)
library(cluster)
library(fpc)

# Covariance matrices
mlist.s.cov <- lapply(lapply(mlist.scaled, t), cov)

# Combine the two lists
Mlist <- append(mlist.scaled, mlist.s.cov)

# Two clustering models to be compaired!

# Model Based Clustering
#fit.1 <- lapply(mlist.scaled, Mclust)
#par(mfrow = c(1,3))
#for(i in 1:4) {
#    plot(fit.1[[i]]) # plot results
#    summary(fit.1[[i]]) # display the best model
#}



# K-Means Clustering with 2 clusters
fit.2 <- lapply(Mlist, function(x){kmeans(x,2)})

# Cluster Plot against 1st 2 principal components


# vary parameters for most readable graph

functional.names <- c("Autophagy", "DNA Repair", "Fatty acid oxidation", "Response to toxic substance")

pdf('Vary parameter.pdf')

par(mfrow = c(2,2))

for(i in 1:4){
    clusplot(Mlist[[i]], fit.2[[i]]$cluster, color=TRUE, shade=TRUE,
    labels=0, lines=1, main=functional.names[i])
}

dev.off()

pdf('Vary parameter Cov.pdf')

par(mfrow = c(2,2))

for(i in 1:4){
    clusplot(Mlist[[i+4]], fit.2[[i+4]]$cluster, color=TRUE, shade=TRUE,
    labels=0, lines=1, main=functional.names[i])
}

dev.off()


# Centroid Plot against 1st 2nd discriminant functions

# functional.names <- c("Autophagy", "DNA Repair", "Fatty acid oxidation", "Response to toxic substance")

pdf('Centroid plot against 1st 2nd discr functions.pdf')

for(j in c("dc", "bc", "vbc", "mvdc", "adc", "awc", "arc", "nc", "wnc", "anc")){
    par(mfrow = c(2,2))
    for(i in 1:4){
        plotcluster(Mlist[[i]], fit.2[[i]]$cluster, main=functional.names[i], clnum=2, method=j)
    }
}

dev.off()

pdf('Centroid plot against 1st 2nd discr functions C.pdf')

for(j in c("dc", "adc")){# remaining picks are not singular!
    par(mfrow = c(2,2))
    for(i in 1:4){
        plotcluster(Mlist[[i+4]], fit.2[[i+4]]$cluster, main=functional.names[i], clnum=2, method=j)
    }
}

dev.off()

d <- date()

setwd('/Users/morilla/Desktop/Wainrib/LAGA/Bichat/Jean-Pierre/Epistasis/Bell/Algorithm/model/model_Ian/model1/Oyy_test_persample/history/')

savehistory(paste0(d,".history"))

# Determine number of clusters (I just wanna 2 of them!)
#wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
#for (i in 2:15) wss[i] <- sum(kmeans(mydata,centers=i)$withinss)
#plot(1:15, wss, type="b", xlab="Number of Clusters",
#ylab="Within groups sum of squares")