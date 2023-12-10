
load("/Users/morilla/Desktop/Wainrib/LAGA/Bichat/Jean-Pierre/Epistasis/Bell/Algorithm/model/model_Ian/model1/Oyy_test_persample/results/Oyy 4 patient 2 cleaned clusters_1.RData")

cluster.list <- list(fit.2[[1]]$cluster,fit.2[[2]]$cluster,fit.2[[3]]$cluster, fit.2[[4]]$cluster)

# patient cluster indices for the list
pts.cluster.2 <- lapply(cluster.list, function(x) which(x==2))
pts.cluster.1 <- lapply(pts.cluster.2, function(x) setdiff(1:395, x))

# c("Autophagy", "DNA Repair", "Fatty acid oxidation", "Response to toxic substance")
k<- 1

ptA1 <- m.1[pts.cluster.1[[k]],]
ptA2 <- m.1[pts.cluster.2[[k]],] # aut <- 1, DNA repair <- 2, etc.

rownames(m.1)<-paste0("sam_",1:395)

ptA1<-as.data.frame(t(ptA1))
ptA2<-as.data.frame(t(ptA2))

write.table(ptA1, file="GEF1.txt", quote=F, sep="\t")
write.table(ptA2, file="GEF2.txt", quote=F, sep="\t")

##
SNP<-read.table("/Users/morilla/Desktop/Wainrib/LAGA/Bichat/Jean-Pierre/Epistasis/Bell/Algorithm/model/model_Ian/model1/SCGGM_code/data/xtrain_IBD.txt")

SNP<-as.data.frame(t(SNP))
colnames(SNP)<-paste0("sam_", 1:395)
rownames(SNP)<-paste0("snp_", 1:218718)

SNP.F1 <- SNP[,pts.cluster.1[[k]]]
SNP.F2 <- SNP[,pts.cluster.2[[k]]]


write.table(SNP.F1, file="SNPF1.txt", quote=F, sep="\t")
write.table(SNP.F2, file="SNPF2.txt", quote=F, sep="\t")


##
cvrt <- read.table("/Users/morilla/Desktop/Wainrib/LAGA/Bichat/Jean-Pierre/Epistasis/Bell/Algorithm/model/model_Ian/model1/data/liming_expression_rawdata/MRCA/E-MTAB-1425.sdrf.txt",skip=1)

df<-data.frame(gender=cvrt$V10, age=cvrt$V12)

df$gender=ifelse(df$gender=="female", 0,1)
df<-as.data.frame(t(df))

colnames(df)<-paste0("sam_",1:395)


cvrt.F1 <- df[,pts.cluster.1[[k]]]
cvrt.F2 <- df[,pts.cluster.2[[k]]]
write.table(cvrt.F1, file="cvrtF1.txt", quote=F, sep="\t")
write.table(cvrt.F2, file="cvrtF2.txt", quote=F, sep="\t")

##

load("/Users/morilla/Desktop/Wainrib/LAGA/Bichat/Jean-Pierre/Epistasis/Bell/Algorithm/model/model_Ian/model1/Oyy_test_persample/figures/graphs/graph_networks/data/Oyy_clean_clustering_nets.RData")

#TBC Co-expression
# Element-wise mean over list of matrices

C.2 <- lapply(1:4, function(i) C[[i]][pts.cluster.2[[i]]])
C.1 <- lapply(1:4, function(i) C[[i]][pts.cluster.1[[i]]])

C.2.glued <- lapply(1:4, function(i) do.call(cbind, C.2[[i]]))
C.1.glued <- lapply(1:4, function(i) do.call(cbind, C.1[[i]]))

C.2.ar <- lapply(1:4, function(i) array(C.2.glued[[i]], dim=c(dim(C.2[[i]][[1]]), length(C.2[[i]]))))
C.1.ar <- lapply(1:4, function(i) array(C.1.glued[[i]], dim=c(dim(C.1[[i]][[1]]), length(C.1[[i]]))))

temp_autoph.1 <- C.1.ar[[1]]
temp_autoph.2 <- C.2.ar[[1]]
ngenes <- dim(temp_autoph.1)[1]
npatients <- append(dim(temp_autoph.1)[3], dim(temp_autoph.2)[3])

autoph.1.pca <- lapply(1:npatients[1], function(i) prcomp(temp_autoph.1[,,i]))
autoph.1 <- lapply(1:npatients[1], function(i) autoph.1.pca[[i]]$rotation[,1])
autoph.1 <- t(matrix(unlist(autoph.1), ncol=ngenes, byrow=T))
colnames(autoph.1) <- paste0("sam_", pts.cluster.1[[1]])
rownames(autoph.1) <- colnames(m.1)

autoph.2.pca <- lapply(1:npatients[2], function(i) prcomp(temp_autoph.2[,,i]))
autoph.2 <- lapply(1:npatients[2], function(i) autoph.2.pca[[i]]$rotation[,1])
autoph.2 <- t(matrix(unlist(autoph.2), ncol=ngenes, byrow=T))
colnames(autoph.2) <- paste0("sam_", pts.cluster.2[[1]])
rownames(autoph.2) <- colnames(m.1)

base.dir <- '/Users/morilla/Desktop/Wainrib/LAGA/Bichat/Jean-Pierre/Epistasis/Bell/Algorithm/model/model_Ian/model1/Oyy_test_persample/data/eQTL';

write.table(autoph.1, file=paste0(base.dir, "/GCEA1.txt"), quote=F, row.names=T, sep="\t")
write.table(autoph.2, file=paste0(base.dir, "/GCEA2.txt"), quote=F, row.names=T, sep="\t")



