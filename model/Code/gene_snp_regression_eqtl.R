# source("Matrix_eQTL_R/Matrix_eQTL_engine.r");
library(MatrixEQTL)
library(Biobase)
library(devtools)

## Location of the package with the data files.
base.dir <- '/Users/morilla/Desktop/Wainrib/LAGA/Bichat/Jean-Pierre/Epistasis/Bell/Algorithm/model/model_Ian/model1/Oyy_test_persample';
# base.dir = '.';

## Settings

k <- 2
# Genotype file name
SNP_file_name = paste(base.dir, paste0("/data/eQTL/SNPA", k,".txt"), sep="");

# Gene expression file name
expression_file_name = paste(base.dir, paste0("/data/eQTL/GCEA", k, ".txt"), sep="");

# Covariates file name
# Set to character() for no covariates
response_file_name = paste(base.dir, paste0("/results/eQTL/Autophagyco", k, "_all_eqtls.txt"), sep=""); # authopagy, response, etc.

# modelLinear graphically

snp <- read.table(SNP_file_name, sep="\t", header=T)
gene <- read.table(expression_file_name, sep="\t", header=T)
response <- read.table(response_file_name, sep="\t", header=T)

# all_eqtls_data
#indicesg <- match(intersect(rownames(gene),response$gene),rownames(gene))
#indicessnp <- match(intersect(rownames(snp),response$snps),rownames(snp))
genes_id <- response$gene[which(response$gene %in% rownames(gene))]
snps_id <- response$snps[which(response$snps %in% rownames(snp))]

dim <- dim(response)[1]
lm <- list()
lm<-lapply(1:dim, function(i)lm(as.numeric(gene[rownames(gene)==genes_id[i],])~as.numeric(snp[rownames(snp)==snps_id[i],])))


# plot fitting (co~coexpression) you must write /figures/eQTL/X etc. down each time!

pdf(paste(base.dir, paste0("/figures/eQTL/Autophagyco", k, "_all_eqtls_genotype.pdf"), sep=""), height=11)
tropical <- c("darkorange", "dodgerblue", "hotpink", "limegreen", "yellow")
palette(tropical)
par(pch=19)
d <- ifelse(sqrt(dim)%%1==0, sqrt(dim(response)[1]), ceiling(sqrt(dim)))
par(mfrow=c(d,d), oma = c( 0, 0, 2, 0 ))

for(i in 1:dim){
    plot(as.numeric(gene[rownames(gene)==genes_id[i],]) ~ jitter(as.numeric(snp[rownames(snp)==snps_id[i],])),
    col=(as.numeric(snp[rownames(snp)==snps_id[i],])+1), xaxt="n", xlab="Genotype", ylab="Co-Expression", main=paste0(as.character(response$snps[i]),"\n", as.character(response$gene[i])))
axis(1, at=c(0:2), labels=c("AA", "Aa", "aa"))
lines(lm[[i]]$fitted ~ as.numeric(snp[rownames(snp)==snps_id[i],]), type="b", pch=15, col="darkgrey")
title( paste0("Autophagyco", k,"_all_eqtls_genotype") , outer = TRUE )
}
dev.off()

##
