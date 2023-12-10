# loading onto workspace the selected 4 GO modules related to IBD
filesnames <- list.files("/Users/morilla/Desktop/Wainrib/LAGA/Bichat/Jean-Pierre/Epistasis/Bell/Algorithm/model/model_Ian/model1/genes_per_modules", pattern="*BP.txt", full.names=T)

ldf<-lapply(filesnames,read.delim,header=F)

ldf.gnames <- lapply(ldf, "[", "V2")
ldf.gnames <- unique(do.call(rbind,ldf.gnames))
colnames(ldf.gnames) <- "Composite.DB.Entry"


###
# Phenotypic layer's generation

# MRCA phenotype_probe correspondance (adf.clean)
load("/Users/morilla/Desktop/Wainrib/LAGA/Bichat/Jean-Pierre/Epistasis/Bell/Algorithm/model/model_Ian/model1/data/liming_expression_rawdata/MRCA/adf.MRCA.RData")

#ldf.gnames <- lapply(ldf.gnames, function(x){
#    colnames(x) <- "Composite.DB.Entry"
#    return(x)
#})

###
# Composite.DB.Entry - affy probe correspondance (89 individuals / some few probes released from duplicated DB.Entries)
gprob <- merge(ldf.gnames, adf.clean, by="Composite.DB.Entry")
rownames(gprob) <- gprob$Composite.Name
gprob <- gprob[-2]

# loading expression matrix (m, 395 x 54675)
load("/Users/morilla/Desktop/Wainrib/LAGA/Bichat/Jean-Pierre/Epistasis/Bell/Algorithm/model/model_Ian/model1/data/liming_expression_rawdata/MRCA/correlation_analysis_MRCA_Bell.RData")

# re-arrenging labels in m
colnames(m) <- adf.clean$Composite.Name
# phenotypic layer for the chosen 4 sub-groups
MAT <- m[,rownames(gprob)]

###

# Genotypic network construction
# loading of genotype (SNPs) data into workspace
library(trio)

MRCA_SNPs <- read.pedfile("/Users/morilla/Desktop/Wainrib/LAGA/Bichat/Jean-Pierre/Epistasis/Bell/Algorithm/model/model_Ian/model1/data/liming/ILMN100k_plink.ped", p2g = F)
# MRCA family information
MRCA_ped <- read.pedfile("/Users/morilla/Desktop/Wainrib/LAGA/Bichat/Jean-Pierre/Epistasis/Bell/Algorithm/model/model_Ian/model1/data/liming/MRCA_expression_DeIdentified_395.ped", p2g = F)

# Both Phenotype and Genotype layers match the headers
ped.infor <- data.frame(famid=MRCA_ped$famid, pid=MRCA_ped$pid, fatid= MRCA_ped$fatid, motid= MRCA_ped$motid, sex=MRCA_ped$sex, affected= -9)
rm(MRCA_ped)

###
# anchorage

# e-QTL links
eqtl.infor <- read.csv("/Users/morilla/Desktop/Wainrib/LAGA/Bichat/Jean-Pierre/Epistasis/Bell/Algorithm/model/model_Ian/model1/data/liming_expression_rawdata/database_eQTL/MRCA_5percentFDR.csv", header=T)
temp.eqtl.infor <- eqtl.infor[,c(1,4)]

gprob$ProbeID<-rownames(gprob)
anchors <- merge(gprob,temp.eqtl.infor,by="ProbeID") # rs_ids are in there as well! so not forced to make use of the snps.gz file

# unique rs ids from the merged data
anc_SNP <- data.frame(unique(anchors$SNP))
colnames(anc_SNP) <- "Name"
I100k_SNP_samples <- read.delim("/Users/morilla/Desktop/Wainrib/LAGA/Bichat/Jean-Pierre/Epistasis/Bell/Algorithm/model/model_Ian/model1/data/liming/ILMN100k_DeIdentified_ExpressionSamples.dat", header=F)
I100k_SNP_samples <- data.frame(I100k_SNP_samples$V2[-1])
colnames(I100k_SNP_samples) <- "Name"
anchors_ix <- match(anc_SNP$Name, I100k_SNP_samples$Name)
anchors_ix <- anchors_ix[!is.na(anchors_ix)]

###
# Matrix of SNPs plus eqtl-based anchors!
MATG <- MRCA_SNPs[,7:dim(MRCA_SNPs)[2]]
col_ix <- colnames(MATG[, anchors_ix])
col_ix.2 <- grep(".2$",col_ix)
col_ix.1 <- grep(".1$",col_ix)
ncol.2 <- anchors_ix[grep(".2$",col_ix)]-1
ncol.1 <- anchors_ix[grep(".1$",col_ix)]+1
ncol <- append(ncol.2,ncol.1)
total.col <- append(ncol,anchors_ix)

# Genotypic layer for the chosen 4 sub-groups
MATG <- MATG[, sort(total.col)]

###
# saving data in RData format

save (list=ls(), file="/Users/morilla/Desktop/Wainrib/LAGA/Bichat/Jean-Pierre/Epistasis/Bell/Algorithm/model/model_Ian/model1/data/liming/Phenot_Genoty_eqtl_graphs.RData")



