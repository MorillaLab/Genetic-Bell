
# Customize as needed for file locations

# Modify data.dir to indicate the location of the GWAStutorial files
# Intermediate data files will also be stored in this same location unless you set out.dir
data.dir <- '/Users/morilla/Desktop/Wainrib/LAGA/Bichat/Jean-Pierre/Epistasis/Bell/Algorithm/model/model_Ian/model1/data/liming/ILMN100k'

out.dir <- data.dir

gwas.fn <- lapply(c(bed='bed',bim='bim',fam='fam',gds='gds'), function(n) sprintf("%s/myPlinkBinaryData.%s", data.dir, n))

#Read in clinical file
clinical <- read.table("/Users/morilla/Desktop/Wainrib/LAGA/Bichat/Jean-Pierre/Epistasis/Bell/Algorithm/model/model_Ian/model1/data/liming_expression_rawdata/MRCA/E-MTAB-1425.sdrf.txt",skip=1)
clinical <- data.frame(FamID=clinical$V6,CAD=clinical$V11, sex=clinical$V10, age=clinical$V12)

i <- which(clinical$sex=="female")
clinical$sex<-1
clinical$sex[i]<-2
# j <- which(clinical$CAD=="asthma")
# clinical$CAD<-0
# clinical$CAD[j]<-1
# networks 1 or 2
load('/Users/morilla/Desktop/Wainrib/LAGA/Bichat/Jean-Pierre/Epistasis/Bell/Algorithm/model/model_Ian/model1/Oyy_test_persample/results/Autophagy_cluster_phenotype.RData')
clinical$CAD <- phenotype

# Output files
gwaa.fname <- sprintf("%s/GWASout.txt", out.dir)
gwaa.unadj.fname <- sprintf("%s/GWASoutUnadj.txt", out.dir)

working.data.fname <- function(num) { sprintf("%s/working.%s.Rdata", out.dir, num) }

# Run this once interactively to download and install BioConductor packages and other packages.

#source("http://bioconductor.org/biocLite.R")
#biocLite("snpStats")
#biocLite("SNPRelate")
#biocLite("rtracklayer")
#biocLite("biomaRt")
#install.packages(c('plyr', 'GenABEL', 'LDheatmap','doParallel', 'ggplot2', 'coin', 'igraph', 'devtools', 'downloader'))
lp <- c('plyr', 'GenABEL', 'LDheatmap','doParallel', 'ggplot2', 'coin', 'igraph', 'devtools', 'downloader', 'postgwas')
lapply(lp, require, character.only = TRUE)

#library(devtools)
#install_url("http://cran.r-project.org/src/contrib/Archive/postgwas/postgwas_1.11.tar.gz")

# Data preprocessing

library(snpStats)

# Read in PLINK files
geno <- read.plink(gwas.fn$bed, gwas.fn$bim, gwas.fn$fam, na.strings = ("-9"))

# Obtain the SnpMatrix object (genotypes) table from geno list
# Note: Phenotypes and covariates will be read from the clinical data file, below
genotype <- geno$genotype
rownames(genotype)<- geno$fam$pedigree
# Subset genotype for subject data
genotype <- genotype[clinical$FamID, ]

#Obtain the SNP information from geno list
genoBim <- geno$map
colnames(genoBim) <- c("chr", "SNP", "gen.dist", "position", "A1", "A2")


# Write genotype, genoBim, clinical for future use
save(genotype, genoBim, clinical, file = working.data.fname(1))

# SNP level filtering

# Create SNP summary statistics (MAF, call rate, etc.)
snpsum.col <- col.summary(genotype)

# Setting thresholds
call <- 0.95
minor <- 0.01

# Filter on MAF and call rate (Using these summary statistics, we keep the subset of SNPs that meet our criteria for minimum call rate and minor allele frequency.)
use <- with(snpsum.col, (!is.na(MAF) & MAF > minor) & Call.rate >= call)
use[is.na(use)] <- FALSE                # Remove NA's as well

cat(ncol(genotype)-sum(use),"SNPs will be removed due to low MAF or call rate.\n")


# Subset genotype and SNP summary data for SNPs that pass call rate and MAF criteria
genotype <- genotype[,use]
snpsum.col <- snpsum.col[use,]

# Write subsetted genotype data and derived results for future use
save(genotype, snpsum.col, genoBim, clinical, file=working.data.fname(2))

# Sample level filtering
# Sample level filtering

# source("globals.R")

# load data created in previous snippets
# load(working.data.fname(2))

library(snpStats)

# Basic sample filtering

library(SNPRelate)                      # LD pruning, relatedness, PCA
library(plyr)

# Create sample statistics (Call rate, Heterozygosity)
snpsum.row <- row.summary(genotype)

# Add the F stat (inbreeding coefficient) to snpsum.row
MAF <- snpsum.col$MAF
callmatrix <- !is.na(genotype)
hetExp <- callmatrix %*% (2*MAF*(1-MAF))
hetObs <- with(snpsum.row, Heterozygosity*(ncol(genotype))*Call.rate)
snpsum.row$hetF <- 1-(hetObs/hetExp)

# Setting thresholds
sampcall <- 0.95    # Sample call rate cut-off
hetcutoff <- 0.1    # Inbreeding coefficient cut-off

sampleuse <- with(snpsum.row, !is.na(Call.rate) & Call.rate > sampcall & abs(hetF) <= hetcutoff)
sampleuse[is.na(sampleuse)] <- FALSE    # remove NA's as well
cat(nrow(genotype)-sum(sampleuse), "subjects will be removed due to low sample call rate or inbreeding coefficient.\n")

# Subset genotype and clinical data for subjects who pass call rate and heterozygosity crtieria
genotype <- genotype[sampleuse,]
clinical <- clinical[ rownames(genotype), ]
clinical$FamID <- gsub("\\.[0-9]+","",rownames(clinical))



# IBD analysis

# Checking for Relatedness

ld.thresh <- 0.2    # LD cut-off
kin.thresh <- 0.1   # Kinship cut-off

# Create gds file, required for SNPRelate functions
snpgdsBED2GDS(gwas.fn$bed, gwas.fn$fam, gwas.fn$bim, gwas.fn$gds)

genofile <- openfn.gds(gwas.fn$gds, readonly = FALSE)

# Automatically added "-1" sample suffixes are removed
gds.ids <- read.gdsn(index.gdsn(genofile, "sample.id"))
#gds.ids <- sub("\\-[0-9]+", "", gds.ids)
#add.gdsn(genofile, "sample.id", gds.ids, replace = TRUE)
# add.gdsn(genofile, "snp.id", gds.ids, replace = TRUE)

#Prune SNPs for IBD analysis
set.seed(1000)
#geno.sample.ids <- rownames(genotype)
geno.sample.ids <- gds.ids[sampleuse]
snpSUB <- snpgdsLDpruning(genofile, ld.threshold = ld.thresh,
sample.id = geno.sample.ids, # Only analyze the filtered samples
snp.id = colnames(genotype)) # Only analyze the filtered SNPs

snpset.ibd <- unlist(snpSUB, use.names=FALSE)
cat(length(snpset.ibd),"will be used in IBD analysis\n")

# Find IBD coefficients using Method of Moments procedure.  Include pairwise kinship.
ibd <- snpgdsIBDMoM(genofile, kinship=TRUE,
sample.id = NULL,
snp.id = snpset.ibd,
num.thread = 1)

ibdcoeff <- snpgdsIBDSelection(ibd)     # Pairwise sample comparison
head(ibdcoeff)




# iteratively remove samples with high kinship starting with the sample with the most pairings
related.samples <- NULL
while ( nrow(ibdcoeff) > 0 ) {
    
    # count the number of occurrences of each and take the top one
    sample.counts <- arrange(count(c(ibdcoeff$ID1, ibdcoeff$ID2)), -freq)
    rm.sample <- sample.counts[1, 'x']
    cat("Removing sample", as.character(rm.sample), 'too closely related to', sample.counts[1, 'freq'],'other samples.\n')
    
    # remove from ibdcoeff and add to list
    ibdcoeff <- ibdcoeff[ibdcoeff$ID1 != rm.sample & ibdcoeff$ID2 != rm.sample,]
    related.samples <- c(as.character(rm.sample), related.samples)
}

# filter genotype and clinical to include only unrelated samples
genotype <- genotype[ !(rownames(genotype) %in% related.samples), ]
clinical <- clinical[ !(clinical$FamID %in% related.samples), ]

geno.sample.ids <- gds.ids[sampleuse]

cat(length(related.samples), "similar samples removed due to correlation coefficient >=", kin.thresh,"\n")


# Checking for ancestry
# To better understand ancestry, we plot the first two principal components of the genotype data. We are reasonably confident that our samples are homogeneous, coming from european ancestry. Therefore, given that there are no clear outliers, we fail to remove any samples.

# Find PCA matrix
pca <- snpgdsPCA(genofile, sample.id = geno.sample.ids,  snp.id = snpset.ibd, num.thread=1)

# Create data frame of first two principal comonents
pctab <- data.frame(sample.id = pca$sample.id,
PC1 = pca$eigenvect[,1],    # the first eigenvector
PC2 = pca$eigenvect[,2],    # the second eigenvector
stringsAsFactors = FALSE)

# Plot the first two principal comonents
pdf(paste0(out.dir,'/AncestryPlot_2_PCA.pdf'))
plot(pctab$PC2, pctab$PC1, xlab="Principal Component 2", ylab="Principal Component 1", main = "Ancestry Plot")
dev.off()
# Close GDS file
closefn.gds(genofile)

# Overwrite old genotype with new filtered version
save(genotype, genoBim, clinical, ibdcoeff, sampleuse, file=working.data.fname(3))

# Hardy-Weinberg SNP filtering on CAD controls

hardy <- 10^-6      # HWE cut-off

CADcontrols <- clinical[ clinical$CAD==1, 'FamID' ]
snpsum.colCont <- col.summary( genotype[CADcontrols,] )
HWEuse <- with(snpsum.colCont, !is.na(z.HWE) & ( abs(z.HWE) < abs( qnorm(hardy/2) ) ) )
rm(snpsum.colCont)

HWEuse[is.na(HWEuse)] <- FALSE          # Remove NA's as well
cat(ncol(genotype)-sum(HWEuse),"SNPs will be removed due to high HWE.\n")  # 1296 SNPs removed

# Subset genotype and SNP summary data for SNPs that pass HWE criteria
genotype <- genotype[,HWEuse]

# Overwrite old genotype with new filtered version
save(genotype, genoBim, clinical, file=working.data.fname(4))


### Data generation
genofile <- openfn.gds(gwas.fn$gds, readonly = FALSE)
ld.thresh <- 0.2

set.seed(1000)
geno.sample.ids <- gds.ids[sampleuse]
snpSUB <- snpgdsLDpruning(genofile, ld.threshold = ld.thresh,
sample.id = geno.sample.ids, # Only analyze the filtered samples
snp.id = colnames(genotype)) # Only analyze the filtered SNPs

snpset.pca <- unlist(snpSUB, use.names=FALSE)
cat(length(snpset.pca),"\n")

pca <- snpgdsPCA(genofile, sample.id = NULL,  snp.id = snpset.pca, num.thread=1)

pcs <- data.frame(FamID = pca$sample.id, pca$eigenvect[,1 : 10],
stringsAsFactors = FALSE)
colnames(pcs)[2:11]<-paste("pc", 1:10, sep = "")

print(head(pcs))

closefn.gds(genofile)

# Store pcs for future reference with the rest of the derived data
save(genotype, genoBim, clinical, pcs, file=working.data.fname(5))

# Parallel model fitting

# Genome-wide Association Analysis
# Parallel implementation of linear model fitting on each SNP

GWAA <- function(genodata=genotypes,  phenodata=phenotypes, family = gaussian, filename=NULL,
append=FALSE, workers=getOption("mc.cores",2L), flip=TRUE,
select.snps=NULL, hosts=NULL, nSplits=10)
{
    if (!require(doParallel)) { stop("Missing doParallel package") }
    
    #Check that a filename was specified
    if(is.null(filename)) stop("Must specify a filename for output.")
    
    #Check that the genotype data is of class 'SnpMatrix'
    if( class(genodata)!="SnpMatrix") stop("Genotype data must of class 'SnpMatrix'.")
    
    #Check that there is a variable named 'phenotype' in phenodata table
    if( !"phenotype" %in% colnames(phenodata))  stop("Phenotype data must have column named 'phenotype'")
    
    #Check that there is a variable named 'id' in phenodata table
    if( !"id" %in% colnames(phenodata)) stop("Phenotype data must have column named 'id'.")
    
    #If a vector of SNPs is given, subset genotype data for these SNPs
    if(!is.null(select.snps)) genodata<-genodata[,which(colnames(genodata)%in%select.snps)]
    
    #Check that there are still SNPs in 'SnpMatrix' object
    if(ncol(genodata)==0) stop("There are no SNPs in the 'SnpMatrix' object.")
    
    #Print the number of SNPs to be checked
    cat(paste(ncol(genodata), " SNPs included in analysis.\n"))
    
    #If append=FALSE than we will overwrite file with column names
    if(!isTRUE(append)) {
        columns<-c("SNP", "Estimate", "Std.Error", "t-value", "p-value")
        write.table(t(columns), filename, row.names=FALSE, col.names=FALSE, quote=FALSE)
    }
    
    # Check sample counts
    if (nrow(phenodata) != nrow(genodata)) {
        warning("Number of samples mismatch.  Using subset found in phenodata.")
    }
    
    # Order genodata rows to be the same as phenodata
    genodata <- genodata[phenodata$id,]
    
    cat(nrow(genodata), "samples included in analysis.\n")
    
    # Change which allele is counted (major or minor)
    flip.matrix<-function(x) {
        zero2 <- which(x==0)
        two0 <- which(x==2)
        x[zero2] <- 2
        x[two0] <- 0
        return(x)
    }
    
    nSNPs <- ncol(genodata)
    genosplit <- ceiling(nSNPs/nSplits) # number of SNPs in each subset
    
    snp.start <- seq(1, nSNPs, genosplit) # index of first SNP in group
    snp.stop <- pmin(snp.start+genosplit-1, nSNPs) # index of last SNP in group
    
    if (is.null(hosts)) {
        # On Unix this will use fork and mclapply.  On Windows it
        # will create multiple processes on localhost.
        cl <- makeCluster(workers)
    } else {
        # The listed hosts must be accessible by the current user using
        # password-less ssh with R installed on all hosts, all
        # packages installed, and "rscript" is in the default PATH.
        # See docs for makeCluster() for more information.
        cl <- makeCluster(hosts, "PSOCK")
    }
    show(cl)                            # report number of workers and type of parallel implementation
    registerDoParallel(cl)
    
    foreach (part=1:nSplits) %do% {
        # Returns a standar matrix of the alleles encoded as 0, 1 or 2
        genoNum <- as(genodata[,snp.start[part]:snp.stop[part]], "numeric")
        
        # Flip the numeric values of genotypes to count minor allele
        if (isTRUE(flip)) genoNum <- flip.matrix(genoNum)
        
        # For each SNP, concatenate the genotype column to the
        # phenodata and fit a generalized linear model
        rsVec <- colnames(genoNum)
        res <- foreach(snp.name=rsVec, .combine='rbind') %dopar% {
            a <- summary(glm(phenotype~ . - id, family=family, data=cbind(phenodata, snp=genoNum[,snp.name])))
            a$coefficients['snp',]
        }
        
        # write results so far to a file
        write.table(cbind(rsVec,res), filename, append=TRUE, quote=FALSE, col.names=FALSE, row.names=FALSE)
        
        cat(sprintf("GWAS SNPs %s-%s (%s%% finished)\n", snp.start[part], snp.stop[part], 100*part/nSplits))
    }
    
    stopCluster(cl)
    
    return(print("Done."))
}

# Phenotype data preparation
library(GenABEL)
# Data frame of phenotype features that is the concatenation of clinical features and the first ten principal components.

#source("GWAA.R")
pcs <- pcs[sampleuse,]
# Merge clincal data and principal components to create phenotype table
#pcs$FamID <- gsub("[0-9]-","",pcs$FamID)
#pcs$FamID <- clinical$FamID
clinical$FamID <- pcs$FamID
phenoSub <- merge(clinical,pcs)      # data.frame => [ FamID CAD sex age hdl pc1 pc2 ... pc10 ]

# We will do a rank-based inverse normal transformation of age
phenoSub$phenotype <- rntransform(phenoSub$age, family="gaussian")

# Show that the assumptions of normality met after transformation
pdf(paste0(out.dir,'/Normality_assumption_after_transformation.pdf'))
par(mfrow=c(1,2))
hist(phenoSub$age, main="Histogram of age", xlab="age")
hist(phenoSub$phenotype, main="Histogram of Tranformed age", xlab="Transformed age")
dev.off()

# Remove unnecessary columns from table
phenoSub$CAD <- NULL

# Rename columns to match names necessary for GWAS() function
phenoSub <- rename(phenoSub, replace=c(FamID="id"))

# Include only subjects with age data
phenoSub<-phenoSub[!is.na(phenoSub$phenotype),]
# 1309 subjects included with phenotype data

print(head(phenoSub))

rownames(genotype) <- phenoSub$id
# Run GWAS analysis
# Note: This function writes a file, but does not produce an R object
start <- Sys.time()
GWAA(genodata=genotype, phenodata=phenoSub, filename="gwaa.fname")
end <- Sys.time()
print(end-start)


# Model fitting of non-typed SNPs
# Carry out association testing for imputed SNPs using snp.rhs.tests()
# rownames(phenoSub) <- phenoSub$id

# imp <- snp.rhs.tests(phenotype ~ sex + age + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10,
# family = "Gaussian", data = phenoSub, snp.data = target, rules = rules)

# Obtain p values for imputed SNPs by calling methods on the returned GlmTests object.
#results <- data.frame(SNP = imp@snp.names, p.value = p.value(imp), stringsAsFactors = FALSE)
#results <- results[!is.na(results$p.value),]

#Write a file containing the results
#write.csv(results, impute.out.fname, row.names=FALSE)

# Merge imputation testing results with support to obtain coordinates
#imputeOut<-merge(results, support[, c("SNP", "position")])
#imputeOut$chr <- 16

#imputeOut$type <- "imputed"

# Find the -log_10 of the p-values
#imputeOut$Neg_logP <- -log10(imputeOut$p.value)

# Order by p-value
#imputeOut <- arrange(imputeOut, p.value)
#print(head(imputeOut))

# Add phenosub to saved results
# save(genotype, genoBim, clinical, pcs, imputed, target, rules, phenoSub, support, file=working.data.fname(7))


# Read in GWAS output that was produced by GWAA function
GWASout <- read.table(paste0(out.dir,'/gwaa.fname'), header=TRUE, colClasses=c("character", rep("numeric",4)))

# Find the -log_10 of the p-values
GWASout$Neg_logP <- -log10(GWASout$p.value)

# Merge output with genoBim by SNP name to add position and chromosome number
GWASout <- merge(GWASout, genoBim[,c("SNP", "chr", "position")])
rm(genoBim)

# Order SNPs by significance
GWASout <- arrange(GWASout, -Neg_logP)
print(head(GWASout))
GWASout$type <- "typed"

save(genotype, clinical, pcs, phenoSub, GWASout, file=working.data.fname(6))


# Manhattan plots are used to visual GWA significant results by chromosome location
# Receives a data.frame of SNPs with Neg_logP, chr, position, and type.
# Plots Manhattan plot with significant SNPs highlighted.
GWAS_Manhattan <- function(GWAS, col.snps=c("black","gray"),
col.detected=c("black"), col.imputed=c("blue"), col.text="black",
title="GWAS Tutorial Manhattan Plot", display.text=TRUE,
bonferroni.alpha=0.05, bonferroni.adjustment=1000000,
Lstringent.adjustment=10000) {
    
    bonferroni.thresh <- -log10(bonferroni.alpha / bonferroni.adjustment)
    Lstringent.thresh <- -log10(bonferroni.alpha / Lstringent.adjustment)
    xscale <- 1000000
    
    manhat <- GWAS[!grepl("[A-z]",GWAS$chr),]
    
    #sort the data by chromosome and then location
    manhat.ord <- manhat[order(as.numeric(manhat$chr),manhat$position),]
    manhat.ord <- manhat.ord[!is.na(manhat.ord$position),]
    
    ##Finding the maximum position for each chromosome
    max.pos <- sapply(1:21, function(i) { max(manhat.ord$position[manhat.ord$chr==i],0) })
    max.pos2 <- c(0, cumsum(max.pos))
    
    #Add spacing between chromosomes
    max.pos2 <- max.pos2 + c(0:21) * xscale * 10
    
    #defining the positions of each snp in the plot
    manhat.ord$pos <- manhat.ord$position + max.pos2[as.numeric(manhat.ord$chr)]
    
    # alternate coloring of chromosomes
    manhat.ord$col <- col.snps[1 + as.numeric(manhat.ord$chr) %% 2]
    
    # draw the chromosome label roughly in the middle of each chromosome band
    text.pos <- sapply(c(1:22), function(i) { mean(manhat.ord$pos[manhat.ord$chr==i]) })
    
    # Plot the data
    plot(manhat.ord$pos[manhat.ord$type=="typed"]/xscale, manhat.ord$Neg_logP[manhat.ord$type=="typed"],
    pch=20, cex=.3, col= manhat.ord$col[manhat.ord$type=="typed"], xlab=NA,
    ylab="Negative Log P-value", axes=F, ylim=c(0,max(manhat$Neg_logP)+1))
    #Add x-label so that it is close to axis
    mtext(side = 1, "Chromosome", line = 1.25)
    
    points(manhat.ord$pos[manhat.ord$type=="imputed"]/xscale, manhat.ord$Neg_logP[manhat.ord$type=="imputed"],
    pch=20, cex=.4, col = col.imputed)
    
    points(manhat.ord$pos[manhat.ord$type=="typed"]/xscale, manhat.ord$Neg_logP[manhat.ord$type=="typed"],
    pch=20, cex=.3, col = manhat.ord$col[manhat.ord$type=="typed"])
    
    axis(2)
    abline(h=0)
    
    SigNifSNPs <- as.character(GWAS[GWAS$Neg_logP > Lstringent.thresh & GWAS$type=="typed", "SNP"])
    
    #Add legend
    legend("topright",c("Bonferroni corrected threshold (p = 5E-8)", "Candidate threshold (p = 5E-6)"),
    border="black", col=c("gray60", "gray60"), pch=c(0, 0), lwd=c(1,1),
    lty=c(1,2), pt.cex=c(0,0), bty="o", cex=0.6)
    
    #Add chromosome number
    text(text.pos/xscale, -.3, seq(1,22,by=1), xpd=TRUE, cex=.8)
    
    #Add bonferroni line
    abline(h=bonferroni.thresh, untf = FALSE, col = "gray60")
    
    #Add "less stringent" line
    abline(h=Lstringent.thresh, untf = FALSE, col = "gray60", lty = 2 )
    
    #Plotting detected genes
    #Were any genes detected?
    if (length(SigNifSNPs)>0){
        
        sig.snps <- manhat.ord[,'SNP'] %in% SigNifSNPs
        
        points(manhat.ord$pos[sig.snps]/xscale,
        manhat.ord$Neg_logP[sig.snps],
        pch=20,col=col.detected, bg=col.detected,cex=0.5)
        
        text(manhat.ord$pos[sig.snps]/xscale,
        manhat.ord$Neg_logP[sig.snps],
        as.character(manhat.ord[sig.snps,1]), col=col.text, offset=1, adj=-.1, cex=.5)
    }
}

pdf(paste0(out.dir,'/GWAS_Manhattan.pdf'))
GWAS_Manhattan(GWASout)
dev.off()

# Rerun the GWAS using unadjusted model
phenoSub2 <- phenoSub[,c("id","phenotype")] # remove all extra factors, leave only phenotype

GWAA(genodata=genotype, phenodata=phenoSub2, filename="gwaa.unadj.fname")

GWASoutUnadj <- read.table(paste0(out.dir,'/gwaa.unadj.fname'), header=TRUE, colClasses=c("character", rep("numeric",4)))

# Create QQ plots for adjusted and unadjusted model outputs
pdf(paste0(out.dir,'/GWAS_QQplots.pdf'))
par(mfrow=c(1,2))
lambdaAdj <- estlambda(GWASout$t.value^2,plot=TRUE,method="median")
title('Adjusted lambda')
lambdaUnadj <- estlambda(GWASoutUnadj$t.value^2,plot=TRUE,method="median")
title('Unadjusted lambda')
dev.off()

cat(sprintf("Unadjusted lambda: %s\nAdjusted lambda: %s\n", lambdaUnadj$estimate, lambdaAdj$estimate))

# Calculate standardized lambda
lambdaAdj_1000<-1+(lambdaAdj$estimate-1)/nrow(phenoSub)*1000
lambdaUnadj_1000<-1+(lambdaUnadj$estimate-1)/nrow(phenoSub)*1000
cat(sprintf("Standardized unadjusted lambda: %s\nStandardized adjusted lambda: %s\n", lambdaUnadj_1000, lambdaAdj_1000))


#library(LDheatmap)
#library(rtracklayer)

# Add "rs247617" to CETP
# CETP <- rbind.fill(GWASout[GWASout$SNP == "rs247617",], CETP)

# Combine genotypes and imputed genotypes for CETP region
#subgen <- cbind(genotype[,colnames(genotype) %in% CETP$SNP], impCETPgeno)     # CETP subsets from typed and imputed SNPs

# Subset SNPs for only certain genotypes
#certain <- apply(as(subgen, 'numeric'), 2, function(x) { all(x %in% c(0,1,2,NA)) })
#subgen <- subgen[,certain]

# Subset and order CETP SNPs by position
#CETP <- CETP[CETP$SNP %in% colnames(subgen),]
#CETP <- arrange(CETP, position)
#subgen <- subgen[, order(match(colnames(subgen),CETP$SNP)) ]

# Create LDheatmap
#ld <- ld(subgen, subgen, stats="R.squared") # Find LD map of CETP SNPs

#ll <- LDheatmap(ld, CETP$position, flip=TRUE, name="myLDgrob", title=NULL)

# Add genes, recombination
#llplusgenes <- LDheatmap.addGenes(ll, chr = "chr16", genome = "hg19", genesLocation = 0.01)

# Add plot of -log(p)
#library(ggplot2)

#plot.new()
#llQplot2<-LDheatmap.addGrob(llplusgenes, rectGrob(gp = gpar(col = "white")),height = .34)
#pushViewport(viewport(x = 0.483, y= 0.76, width = .91 ,height = .4))

#grid.draw(ggplotGrob({
#   qplot(position, Neg_logP, data = CETP, xlab="", ylab = "Negative Log P-value", xlim = range(CETP$position),
#    asp = 1/10, color = factor(type), colour=c("#000000", "#D55E00")) +
#    theme(axis.text.x = element_blank(),
#    axis.title.y = element_text(size = rel(0.75)), legend.position = "none",
#    panel.background = element_blank(),
#    axis.line = element_line(colour = "black")) +
#    scale_color_manual(values = c("red", "black"))
#}))

# Create regional association plot
# Create data.frame of most significant SNP only
library(postgwas)

# Change column names necessary to run regionalplot function
GWAScomb <- rename(GWASout, c(p.value="P", chr="CHR", position="BP"))

# first candidate threshold
snps <- GWAScomb[GWAScomb$P < 1e-06,]


# Edit biomartConfigs so regionalplot function
# pulls from human genome build 37/hg19

# myconfig <- biomartConfigs$hsapiens
# myconfig$hsapiens$gene$host <- "grch37.ensembl.org"
# myconfig$hsapiens$gene$mart <- "ENSEMBL_MART_ENSEMBL"
# myconfig$hsapiens$snp$host <- "grch37.ensembl.org"
# myconfig$hsapiens$snp$mart <- "ENSEMBL_MART_SNP"

load(paste0(out.dir,'/results/myconfig.Hsapiens_biomart.RData'))


# Run regionalplot using HAPMAP data (pop = CEU)
pdf(paste0(out.dir,'/GWAS_typed_significant_Regional_plots.pdf'))
regionalplot(snps[1:2,], GWAScomb, biomart.config = myconfig, window.size = 400000, draw.snpname = data.frame(
snps = snps$SNP[1:2],
text = snps$SNP[1:2],
angle = c(20, 160),
length = c(1, 1),
cex = c(0.8)
),
ld.options = list(
gts.source = 2,
max.snps.per.window = 2000,
rsquare.min = 0.8,
show.rsquare.text = FALSE
),
out.format = list(file = NULL, panels.per.page = 2))
dev.off()

# Parallel coordinates plot for the top principal components
library(MASS)
pdf(paste0(out.dir,'/parallel_coords_12_PCA_plots.pdf'))
#datpop <- factor(pop_code)[match(pca$sample.id, sample.id)]
parcoord(pca$eigenvect[,1:16], col= 1 + (0:180)%/%50)
dev.off()

# To calculate the SNP correlations between eigenvactors and SNP genotypes:

# Get chromosome index
genofile <- openfn.gds(gwas.fn$gds, readonly = FALSE)
chr <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
CORR <- snpgdsPCACorr(pca, genofile, eig.which=1:4)
closefn.gds(genofile)

png(paste0(out.dir,'/SNP_CORR_eigen_vs_SNP_genotypes.png'))
savepar <- par(mfrow=c(3,1), mai=c(0.3, 0.55, 0.1, 0.25))
for (i in 1:3)
{
    plot(abs(CORR$snpcorr[i,]), ylim=c(0,1), xlab="", ylab=paste("PC", i),
    col=chr, pch="+")
}
par(savepar)
dev.off()

# To perform cluster analysis on the n×nn×n matrix of genome-wide IBS pairwise distances, and determine the groups by a permutation score:

set.seed(100)
genofile <- openfn.gds(gwas.fn$gds, readonly = FALSE)
state <- snpgdsIBS(genofile, num.thread=2)

library(zoo)
state$ibs <- na.approx(state$ibs)
state$ibs[which(is.na(state$ibs))] <- 0
ibs.hc <- snpgdsHCluster(state)

# Determine groups of individuals automatically
rv <- snpgdsCutTree(ibs.hc)
pdf(paste0(out.dir,'/cluster_analysis_GWAS_distances.pdf'))
plot(rv$dendrogram, leaflab="none", main="HapMap Phase II")
dev.off()
closefn.gds(genofile)

table(rv$samp.group)

save(state, CORR, chr, rv, file=working.data.fname(7))


# GWAS to Network
gwasgenes.prox <- snp2gene.prox(snps, biomart.config = myconfig)
# gwasgenes.prox$P <- GWASout$p.value[as.numeric(rownames(gwasgenes.prox))]

# we can now also use a GO network if desired
# species.db <- "org.Hs.eg.db"

network <- getInteractions.GO(
gwasgenes.prox$geneid,
GOpackagename = "org.Hs.eg.db",
toFile = NULL
)

net <- gwas2network(
gwasgenes.prox,
network = network,
vertexcolor.GO.overrep = "org.Hs.eg.db",
max.communities = -1,
biomart.config = myconfig
)

# plot the resulting network
gwas2network.plot(net, device = options("device")$device)

# when file.verbosity >= 2, we can restore graph data from text files to do further study
# v <- read.table("gwas2networkGraphVertices.csv", header = TRUE)
# e <- read.table("gwas2networkGraphEdges.csv", header = TRUE)
# g <- graph.data.frame(e, vertices = v, directed = FALSE)

