
# set the pathway to that we want to
setwd('/Users/morilla/Desktop/Wainrib/LAGA/Bichat/Jean-Pierre/Epistasis/Bell/Algorithm/model/model_Ian/model1/Oyy_test_persample/figures')

# Loading files into the workspace
#filesnames <- list.files("/Users/morilla/Desktop/Wainrib/LAGA/Bichat/Jean-Pierre/Epistasis/Bell/Algorithm/model/model_Ian/model1/Oyy_test_persample", pattern="*BP.txt", full.names=T)
load('/Users/morilla/Desktop/Wainrib/LAGA/Bichat/Jean-Pierre/Epistasis/Bell/Algorithm/model/model_Ian/model1/Oyy_test_persample/data/yliming_sGGMS.RData')

load('/Users/morilla/Desktop/Wainrib/LAGA/Bichat/Jean-Pierre/Epistasis/Bell/Algorithm/model/model_Ian/model1/Oyy_test_persample/Cleaning/CleanData_MRCA.RData')

m.1 <- CleanData[ , which(is.element(colnames(CleanData),colnames(m.1)))]
m.2 <- CleanData[ , which(is.element(colnames(CleanData),colnames(m.2)))]
m.3 <- CleanData[ , which(is.element(colnames(CleanData),colnames(m.3)))]
m.4 <- CleanData[ , which(is.element(colnames(CleanData),colnames(m.4)))]


mlist <- list(m.1, m.2, m.3, m.4)
mlist <- lapply(mlist, function(x){rownames(x)<-paste0("pt",1:395); x})
#mlist.scaled <- lapply(mlist, scale)
mlist.scaled <- lapply(mlist, t) # It must be transposed!


# getting a list from each column

mlist.s.4pat <- lapply(mlist.scaled, function(x){as.list(data.frame(x))})

# Difference matrices (colMeans)
D <- list()
for(i in 1:4){
    D[[i]] <- lapply(1:395, function(x) as.matrix(mlist.s.4pat[[i]][[x]])-mean(mlist.s.4pat[[i]][[x]]))
}

#b1<-apply(t(D),2,list)
# Covariance matrices multiplication
D2 <- list()
for(j in 1:4){
    D2[[j]] <- lapply(1:395, function(x) as.matrix(D[[j]][[x]])%*%t(as.matrix(D[[j]][[x]])))
}
#as.matrix(D[1,])%*%t(as.matrix(D[1,]))

# and finally, lets multiply by N to work the
# Covariance matrices out!
N <- 1/394
C <- list()
for(k in 1:4){
    C[[k]] <- lapply(D2[[k]], "*", N)
}



#mlist.s.4pat.D <- lapply(1:dim(mlist.scaled[[1]])[1], function(x){as.matrix(mlist.scaled[[]))-mean(x)})


# classical MDS
# N rows (objects) x p columns (variables)
mlist.s.d <- list()
for(t in 1:4){
    mlist.s.d[[t]] <- lapply(C[[t]], dist)
}
# projection/embedding onto 1-Dimension space
fit <- list()
for(s in 1:4){
    fit[[s]] <- lapply(mlist.s.d[[s]], function(x) cmdscale(x,eig=TRUE, k=1))
}

# Grand-Matrix Construction
Mlist.s.d <- list()
for(h in 1:4){
    Mlist.s.d[[h]] <- lapply(1:395, function(x) fit[[h]][[x]]$points)
}

Matrix <- list()
for(l in 1:4){
    Matrix[[l]] <- matrix(unlist(Mlist.s.d[[l]]), ncol=395, byrow=F) # previously T
}

# Combine the two lists
#Mlist <- append(mlist.scaled, mlist.s.cov)

# Two clustering models to be compaired!

# Model Based Clustering
#fit.1 <- lapply(mlist.scaled, Mclust)
#par(mfrow = c(1,3))
#for(i in 1:4) {
#    plot(fit.1[[i]]) # plot results
#    summary(fit.1[[i]]) # display the best model
#}




###### clustering core zone #######


library(mclust)
library(cluster)
library(fpc)

Matrix <- lapply(Matrix, t)

# K-Means Clustering with K clusters
K <- 2
fit.2 <- lapply(Matrix, function(x){kmeans(x,K)})

save(list=ls(), file=paste0("Oyy 4 patient ", K, " cleaned clusters.RData"))

# Cluster Plot against 1st 2 principal components


# vary parameters for most readable graph

functional.names <- c("Autophagy", "DNA Repair", "Fatty acid oxidation", "Response to toxic substance")

pdf(paste0('Vary parameter 4 patient ', K, ' cleaned clusters.pdf'))

par(mfrow = c(2,2))

for(i in 1:4){
    clusplot(Matrix[[i]], fit.2[[i]]$cluster, color=TRUE, shade=TRUE,
    labels=0, lines=1, main=functional.names[i])
}

dev.off()


# Centroid Plot against 1st 2nd discriminant functions

# functional.names <- c("Autophagy", "DNA Repair", "Fatty acid oxidation", "Response to toxic substance")
# bc does not work!

pdf(paste0('Centroid plot against 1st 2nd discr functions 4patient ', K, ' cleaned clusters.pdf'))

for(j in c("dc", "mvdc", "adc", "nc", "wnc")){
    par(mfrow = c(2,2))
    for(i in 1:4){
        plotcluster(Matrix[[i]], fit.2[[i]]$cluster, main=functional.names[i], clnum=K, method=j)
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