
rm(list=ls())

base.dir <- "/Users/morilla/Desktop/Wainrib/LAGA/Bichat/Jean-Pierre/Epistasis/Bell/Algorithm/model/model_Ian/model1/"

library(igraph)

load(paste0(base.dir, "Oyy_test_persample/figures/graphs/graph_networks/data/Oyy_clean_clustering_nets.RData"))

g <- nets.2[[2]]

# Set Output for each frame
png(file=paste(base.dir, paste0("Oyy_test_persample/figures/graphs/movies/frames/", "output_%03d.png"), sep=""), width=800,height=450)

dt <- 1
total_time <- length(E(g))
E(g)$time <- 1:total_time

#generate first layout using weights
layout.old <- layout.fruchterman.reingold(g,params=list(weights=E(g)$weight))

#Time loop starts
for(time in seq(1, total_time,dt)){   #dt is the defined interval between successive plots
    
    gt <- delete_edges(g,which(E(g)$time > time)) #remove absent edges
    
    layout.new <- layout_with_fr(gt,coords=layout.old,niter=10,start.temp=0.05,grid="nogrid") #jitter layout
    
    plot(gt,layout=layout.new,vertex.label=V(g)$name,vertex.size=1+2*log(degree(gt)),
    vertex.frame.color=V(g)$color,edge.width=1+abs(E(gt)$weight)*5000,asp=9/16,margin=-0.15, edge.color=ifelse(E(gt)$weight<=0, "blue","red")) #make plot
    
    layout.old <- layout.new #keep layout for next image
}

dev.off()


system("/Users/morilla/Desktop/Wainrib/LAGA/Bichat/Jean-Pierre/Epistasis/Bell/Algorithm/model/model_Ian/model1/functions/ffmpeg -r 10 -i /Users/morilla/Desktop/Wainrib/LAGA/Bichat/Jean-Pierre/Epistasis/Bell/Algorithm/model/model_Ian/model1/Oyy_test_persample/figures/graphs/movies/frames/output_%03d.png /Users/morilla/Desktop/Wainrib/LAGA/Bichat/Jean-Pierre/Epistasis/Bell/Algorithm/model/model_Ian/model1/Oyy_test_persample/figures/graphs/movies/output2.mp4", wait=T)

system("rm -r /Users/morilla/Desktop/Wainrib/LAGA/Bichat/Jean-Pierre/Epistasis/Bell/Algorithm/model/model_Ian/model1/Oyy_test_persample/figures/graphs/movies/frames/*.png", wait=T)
