

load("/Users/morilla/Desktop/Wainrib/LAGA/Bichat/Jean-Pierre/Epistasis/Bell/Algorithm/model/model_Ian/model1/Oyy_test_persample/figures/graphs/graph_networks/data/Oyy_clean_clustering_nets.RData")


load("/Users/morilla/Desktop/Wainrib/LAGA/Bichat/Jean-Pierre/Epistasis/Bell/Algorithm/model/model_Ian/model1/Oyy_test_persample/results/Oyy 4 patient 2 cleaned clusters_1.RData")

cluster.list <- list(fit.2[[1]]$cluster,fit.2[[2]]$cluster,fit.2[[3]]$cluster, fit.2[[4]]$cluster)

pts.cluster.2 <- lapply(cluster.list, function(x) which(x==2))#
pts.cluster.1 <- lapply(pts.cluster.2, function(x) setdiff(1:395, x))#

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

