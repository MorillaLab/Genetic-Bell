# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
#
# Be sure to use an up to date version of R and Matrix eQTL.

# source("Matrix_eQTL_R/Matrix_eQTL_engine.r");
library(MatrixEQTL)
library(Biobase)
library(devtools)

## Location of the package with the data files.
base.dir <- '/Users/morilla/Desktop/Wainrib/LAGA/Bichat/Jean-Pierre/Epistasis/Bell/Algorithm/model/model_Ian/model1/Oyy_test_persample';
# base.dir = '.';

## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel <- modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

k <- 2 # 1 or 2
# Genotype file name
SNP_file_name <- paste(base.dir, paste0("/data/eQTL/SNPA", k,".txt"), sep="");

# Gene expression file name
expression_file_name <- paste(base.dir, paste0("/data/eQTL/GEA", k, ".txt"), sep="");

# Gene co-expression file name
coexpression_file_name <- paste(base.dir, paste0("/data/eQTL/GCEA", k, ".txt"), sep="");

# Covariates file name
# Set to character() for no covariates
covariates_file_name <- paste(base.dir, paste0("/data/eQTL/cvrtA", k, ".txt"), sep="");

# Output file name
output_file_name = tempfile();

# Only associations significant at this level will be saved
pvOutputThreshold = 1e-5; #1e-2, 1e-5 (Autoph.), 1e-6 (DNA rep.), 1e-7 (response),

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");


## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load gene co-expression data

cgene = SlicedData$new();
cgene$fileDelimiter = "\t";      # the TAB character
cgene$fileOmitCharacters = "NA"; # denote missing values;
cgene$fileSkipRows = 1;          # one row of column labels
cgene$fileSkipColumns = 1;       # one column of row labels
cgene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
cgene$LoadFile(coexpression_file_name);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
    cvrt$LoadFile(covariates_file_name);
}

## Run the analysis

me = Matrix_eQTL_engine(
snps = snps,
gene = cgene,
cvrt = cvrt,
output_file_name = output_file_name,
pvOutputThreshold = pvOutputThreshold,
useModel = useModel,
errorCovariance = errorCovariance,
verbose = TRUE,
pvalue.hist = TRUE,
min.pv.by.genesnp = TRUE,#FALSE default
noFDRsaveMemory = FALSE);

unlink(output_file_name);

## Results:
# c("Autophagy", "DNA Repair", "Fatty acid oxidation", "Response to toxic substance")

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected eQTLs:', '\n');
show(me$all$eqtls)

# autophagy, response, etc.
write.table(me$all$eqtls, file=paste(base.dir, paste0("/results/eQTL/Autophagyco", k, "_all_eqtls.txt"), sep=""), quote=F, row.names=F, sep="\t")
save(me, file=paste(base.dir, paste0("/results/eQTL/Autophagyco", k, "_all_eqtls.RData"), sep=""))

## Plot the histogram of all p-values

pdf(paste(base.dir, paste0("/figures/eQTL/Autophagyco", k, "_all_eqtls_pvalues.pdf"), sep=""))
plot(me)
dev.off()

