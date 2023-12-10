
load("/Users/morilla/Desktop/Wainrib/LAGA/Bichat/Jean-Pierre/Epistasis/Bell/Algorithm/model/model_Ian/model1/Oyy_test_persample/results/Oyy 4 patient 2 cleaned clusters_1.RData")


# make a list from every single cluster entry
cluster.list <- list(fit.2[[1]]$cluster,fit.2[[2]]$cluster,fit.2[[3]]$cluster, fit.2[[4]]$cluster)

# patient cluster indices for the list
pts.cluster.2 <- lapply(cluster.list, function(x) which(x==2))
pts.cluster.1 <- lapply(pts.cluster.2, function(x) setdiff(1:395, x))

# Element-wise mean over list of matrices

C.2 <- lapply(1:4, function(i) C[[i]][pts.cluster.2[[i]]])
C.1 <- lapply(1:4, function(i) C[[i]][pts.cluster.1[[i]]])

C.2.glued <- lapply(1:4, function(i) do.call(cbind, C.2[[i]]))
C.1.glued <- lapply(1:4, function(i) do.call(cbind, C.1[[i]]))

C.2.ar <- lapply(1:4, function(i) array(C.2.glued[[i]], dim=c(dim(C.2[[i]][[1]]), length(C.2[[i]]))))
C.1.ar <- lapply(1:4, function(i) array(C.1.glued[[i]], dim=c(dim(C.1[[i]][[1]]), length(C.1[[i]]))))

Mavg.2 <- lapply(1:4, function(i) apply(C.2.ar[[i]], c(1, 2), mean, na.rm = TRUE))
Mavg.1 <- lapply(1:4, function(i) apply(C.1.ar[[i]], c(1, 2), mean, na.rm = TRUE))


save(Mavg.1, Mavg.2, file="/Users/morilla/Desktop/Wainrib/LAGA/Bichat/Jean-Pierre/Epistasis/Bell/Algorithm/model/model_Ian/model1/Oyy_test_persample/results/Oyy 4 patient 2 cleaned clusters_Mavg1_2.RData")

# plot patients networks (2 and 1) with weight matrices Mavg.2 and Mavg.1
library(igraph)
#nodes.2 <- lapply(1:4, function(i) colnames(mlist[[i]]))

# construction of network objects

graph.adj.2 <- lapply(1:4, function(i) graph.adjacency(Mavg.2[[i]], weighted=T))
df.graph.adj.2 <- lapply(1:4, function(i) get.data.frame(graph.adj.2[[i]]))

graph.adj.1 <- lapply(1:4, function(i) graph.adjacency(Mavg.1[[i]], weighted=T))
df.graph.adj.1 <- lapply(1:4, function(i) get.data.frame(graph.adj.1[[i]]))

nets.2 <- lapply(1:4, function(i) graph.data.frame(df.graph.adj.2[[i]], 1:dim(C.2[[i]][[1]])[1], directed=F))
nets.1 <- lapply(1:4, function(i) graph.data.frame(df.graph.adj.1[[i]], 1:dim(C.1[[i]][[1]])[1], directed=F))

nets.2 <- lapply(nets.2, function(x) simplify(x, remove.multiple=T, remove.loops=T))
nets.1 <- lapply(nets.1, function(x) simplify(x, remove.multiple=T, remove.loops=T))

#Color scaling function
#c_scale <- colorRamp(c('red','yellow','cyan','blue'))

for(j in 1:4){
    E(nets.2[[j]])$width <- 1+abs(E(nets.2[[j]])$weight)*5000
    E(nets.1[[j]])$width <- 1+abs(E(nets.1[[j]])$weight)*5000
    #Applying the color scale to edge weights.
    #rgb method is to convert colors to a character vector.
    #E(nets.2[[j]])$color <- apply(c_scale(E(nets.2[[j]])$weight), 1, function(x) rgb(x[1]/255,x[2]/255,x[3]/255) )
    #E(nets.1[[j]])$color <- apply(c_scale(E(nets.1[[j]])$weight), 1, function(x) rgb(x[1]/255,x[2]/255,x[3]/255) )
    #E(nets.2[[j]])$color[E(nets.2[[j]])$weight<0] <- "blue"
    #E(nets.1[[j]])$color[E(nets.1[[j]])$weight<0] <- "blue"
    #E(nets.2[[j]])$color[E(nets.2[[j]])$weight>=0] <- "red"
    #E(nets.1[[j]])$color[E(nets.1[[j]])$weight>=0] <- "red"
}


layouts.2 <- lapply(nets.2, function(x) layout.graphopt(x))
layouts.1 <- lapply(nets.1, function(x) layout.graphopt(x))

functional.names <- c("Autophagy", "DNA Repair", "Fatty acid oxidation", "Response to toxic substance")

#
pdf(paste0(functional.names[1], " cleaned clusters.pdf"), width=11)

par(mfrow=c(1,2))

plot(nets.2[[1]], edge.color=ifelse(E(nets.2[[1]])$weight<=0, "blue","red"), edge.arrow.mode=0, edge.curved=.1, layout=layouts.2[[1]], main=paste0(functional.names[1], " 2"))
plot(nets.1[[1]], edge.color=ifelse(E(nets.1[[1]])$weight<=0, "blue","red"), edge.arrow.mode=0, edge.curved=.1, layout=layouts.1[[1]], main=paste0(functional.names[1], " 1"))

dev.off()
#
pdf(paste0(functional.names[2], " cleaned clusters.pdf"), width=14)

par(mfrow=c(1,2))

plot(nets.2[[2]], edge.color=ifelse(E(nets.2[[2]])$weight<=0, "blue","red"), edge.arrow.mode=0, edge.curved=.1, layout=layouts.2[[2]], main=paste0(functional.names[2], " 2"))
plot(nets.1[[2]], edge.color=ifelse(E(nets.1[[2]])$weight<=0, "blue","red"), edge.arrow.mode=0, edge.curved=.1, layout=layouts.1[[2]], main=paste0(functional.names[2], " 1"))

dev.off()

#
pdf(paste0(functional.names[3], " cleaned clusters.pdf"), width=11)

par(mfrow=c(1,2))

plot(nets.2[[3]], edge.color=ifelse(E(nets.2[[3]])$weight<=0, "blue","red"), edge.arrow.mode=0, edge.curved=.1, layout=layouts.2[[3]], main=paste0(functional.names[3], " 2"))
plot(nets.1[[3]], edge.color=ifelse(E(nets.1[[3]])$weight<=0, "blue","red"), edge.arrow.mode=0, edge.curved=.1, layout=layouts.1[[3]], main=paste0(functional.names[3], " 1"))

dev.off()

#
pdf(paste0(functional.names[4], " cleaned clusters.pdf"), width=11)

par(mfrow=c(1,2))

plot(nets.2[[4]], edge.color=ifelse(E(nets.2[[4]])$weight<=0, "blue","red"), edge.arrow.mode=0, edge.curved=.1, layout=layouts.2[[4]], main=paste0(functional.names[4], " 2"))
plot(nets.1[[4]], edge.color=ifelse(E(nets.1[[4]])$weight<=0, "blue","red"), edge.arrow.mode=0, edge.curved=.1, layout=layouts.1[[4]], main=paste0(functional.names[4], " 1"))

dev.off()
#
pdf("Histogram class 2-1 cleaned clusters.pdf", height=11)

par(mfrow=c(4,2))

for(i in 1:4){
    hist(E(nets.2[[i]])$weight, xlab="weight", main=paste0(functional.names[i]," 2"))
    hist(E(nets.1[[i]])$weight, xlab="weight", main=paste0(functional.names[i]," 1"))
}

dev.off()

#pdf("Histogram class 1 cleaned clusters.pdf")

#par(mfrow=c(2,2))

#for(i in 1:4){
#    hist(df.graph.adj.1[[i]]$weight, main=functional.names[i])
#}

#dev.off()

