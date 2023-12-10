
## select() interface:
## Objects in this package can be accessed using the select() interface
## from the AnnotationDbi package. See ?select for details.


library(hgu133plus2.db) # shift microarray ids
## Bimap interface:
x <- hgu133plus2ACCNUM # 2 entrezid, emsembl, etc.
# Get the probe identifiers that are mapped to an ACCNUM
mapped_probes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_probes])
if(length(xx) > 0) {
    # Get the ACCNUM for the first five probes
    xx[1:5]
    # Get the first one
    xx[[1]]
}


# genes idendification per group
genes.table <- read.table("~/Desktop/Wainrib/LAGA/Bichat/Jean-Pierre/Epistasis/Bell/Algorithm/model/model_Ian/model1/Oyy_test_persample/figures/graphs/graph_networks/data/Fatty acid oxidation_gene_names.txt")

t <- as.character(as.matrix(genes.table))

LS.df <- as.data.frame(do.call(rbind,xx[t]))
colnames(LS.df) <- 'ACCNUM'

write.table(LS.df, file="GROUP NAME_probe_gene_correspondence.txt", quote=F, col.names=T, row.names=T, sep="\t")
