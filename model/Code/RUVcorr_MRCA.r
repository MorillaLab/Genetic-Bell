## ----style-knitr, eval=TRUE, echo=FALSE, results="asis"---------------------------------
BiocStyle::latex()

## ---------------------------------------------------------------------------------------
library(RUVcorr)

## ----sim1-------------------------------------------------------------------------------
#set.seed(400)
#Yind <- simulateGEdata(n=3000, m=1000, k=10, size.alpha=2,
#                       corr.strength=5, g=NULL, Sigma.eps=0.1,
#                       nc=2000, ne=1000, intercept=TRUE, check=TRUE)
#print(Yind)

## ----sim2-------------------------------------------------------------------------------
#set.seed(400)
#Ydep <- simulateGEdata(n=3000, m=1000, k=10, size.alpha=2,
#                       corr.strength=5, g=2, Sigma.eps=0.1,
#                       nc=2000, ne=1000, intercept=TRUE, check=TRUE)
#print(Ydep)

## ----design, fig.show='hide', fig.width=10, fig.height=6, results='hide'----------------
#library(bladderbatch)
#data(bladderdata)
#expr.meta <- pData(bladderEset)
#plotDesign(expr.meta, c("cancer", "outcome", "batch"),
#           c("Diagnosis", "Outcome", "Batch"), orderby="batch")

## ---------------------------------------------------------------------------------------
load('/Users/morilla/Desktop/Wainrib/LAGA/Bichat/Jean-Pierre/Epistasis/Bell/Algorithm/model/model_Ian/model1/Oyy_test_persample/data/yliming_sGGMS.RData')

prb_names <- colnames(m.1)
prb_names <- append(prb_names, colnames(m.2))
prb_names <- append(prb_names, colnames(m.3))
prb_names <- append(prb_names, colnames(m.4))
prb_names <- unique(prb_names)

load("/Users/morilla/Desktop/Wainrib/LAGA/Bichat/Jean-Pierre/Epistasis/Bell/Algorithm/model/model_Ian/model1/data/liming_expression_rawdata/MRCA/MAT.RData")

expr <- t(n)
# expr <- expr[,1:20000]

library(hgu133plus2.db)
x <- hgu133plus2SYMBOL
colnames(expr) <- ls(x)
xx <- as.list(x[colnames(expr)])


## ---------------------------------------------------------------------------------------
# na_genes <- c("SCN1A", "SCN3A", "SCN4A", "SCN5A", "SCN7A", "SCN8A", "SCN11A",
#            "SCN1B", "SCN2B", "SCN3B", "SCN4B")

## ---------------------------------------------------------------------------------------
# na_affy <- names(which(unlist(lapply(xx, function(x) is.element(x, na_genes)[1]))))
na_affy <- prb_names
na_index <- which(is.element(colnames(expr),na_affy))
nc_index <- empNegativeControls(expr, exclude=na_index, nc=3000)

## ----genePlot, fig.show='hide', fig.width=10, fig.height=6, results='hide'--------------
genePlot(expr, index=nc_index, 
         legend="Negative Control Genes", title="IQR-Mean Plot")

## ----warning=FALSE----------------------------------------------------------------------
library(snowfall)
k <- c(1,2,3,4)
nu <- c(0,500,1000,5000)
k.nu.matrix <- cbind(rep(k, each=4), rep(nu, 4))
k.nu.matrix <- as.list(as.data.frame(t(k.nu.matrix)))

sfInit(parallel=TRUE, cpus=4)
sfLibrary(RUVcorr)
sfExport("expr", "k.nu.matrix", "nc_index")
expr_AllRUV <- sfLapply(k.nu.matrix, function(x)
  RUVNaiveRidge(expr, center=TRUE, nc_index, x[2], x[1]))
sfStop()

## ----histogramPlotNA, fig.show='hide', fig.width=10, fig.height=10, results='hide', warning=FALSE----
cor_AllRUV_na <- lapply(expr_AllRUV, function(x) cor(x[,na_index]))
cor_Raw_na <- cor(expr[,na_index])
                      
par(mfrow=c(2,2))
lapply(1:4, function(i) histogramPlot(cor_AllRUV_na[seq(0,15,4)+i], cor_Raw_na,
                                      title=paste("nu=", nu[i]), 
                                      legend=c(paste("k=", k), "Raw")))

## ----histogramPlotBG, fig.show='hide', fig.width=10, fig.height=10, results='hide', warning=FALSE----
bg_index <- background(expr, nBG=100, exclude=na_index, nc_index=nc_index)

cor_AllRUV_bg <- lapply(expr_AllRUV, function(x) cor(x[,bg_index]))
cor_Raw_bg <- cor(expr[,bg_index])
                      
par(mfrow=c(2,2))
lapply(1:4, function(i) histogramPlot(cor_AllRUV_bg[seq(0,15,4)+i], cor_Raw_bg,
                                      title=paste("nu=", nu[i]), 
                                      legend=c(paste("k=", k), "Raw")))

## ----RLEplot1, fig.show='hide', fig.width=10, fig.height=10, results='hide', warning=FALSE----
par(mfrow=c(2,2))
lapply(1:4, function(i) RLEPlot(expr, expr_AllRUV[[12+i]],
                                name=c("Raw", "RUV"), title=paste("nu=", nu[i]),
                                method="IQR.boxplots"))

## ----RLEplot2, fig.show='hide', fig.width=10, fig.height=5, results='hide', warning=FALSE----
par(mfrow=c(1,1))
RLEPlot(expr, expr_AllRUV[[16]], name=c("Raw", "RUV"),
        title="Batches", method="IQR.points", anno=NULL,
        Factor="batch", numeric=TRUE)

## ---------------------------------------------------------------------------------------
CleanData <- expr_AllRUV[[16]]

## ---------------------------------------------------------------------------------------
#cand_genes <- c("CACNA1A", "CACNA1B", "SNAP25", "STX1A")
#cand_affy <- names(which(unlist(lapply(xx, function(x) is.element(x, cand_genes)[1]))))
#cand_index <- which(is.element(colnames(CleanData),cand_affy))

## ---------------------------------------------------------------------------------------
#Prop <- calculateThreshold(CleanData, exclude=c(nc_index, cand_index),
#                      index.ref=na_index, set.size=length(cand_index),
#                      Weights=NULL)
# threshold <- predict(Prop$loess.estimate, 0.3)
# threshold

## ----EFDR, fig.show='hide', fig.width=10, fig.height=6, results='hide', warning=FALSE----
# plotThreshold(Prop)

## ---------------------------------------------------------------------------------------
#prior<-prioritise(CleanData, na_index, cand_index, Weight=NULL, threshold=threshold)
#print(prior)
#xx[which(is.element(names(xx), prior[,1]))]

## ----sessionInfo, results = 'asis', eval = TRUE, echo = TRUE----------------------------
toLatex(sessionInfo())

