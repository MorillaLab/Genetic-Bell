
snps <- read.table("/Users/morilla/Desktop/Wainrib/LAGA/Bichat/Jean-Pierre/Epistasis/Bell/Algorithm/model/model_Ian/model1/Oyy_test_persample/data/eQTL/SNPA1.txt")


nsamples <- dim(snps)[1]

library('parallel')
mc <- getOption("mc.cores", 2)

ctable <- mclapply(1:nsamples, function(i)table(factor(snps[i,], levels=c("0","1","2"), mc)))
dims.ctable <- lapply(1:nsamples, function(i)dim(ctable[[i]]))

temp.out <- which(dims.ctable!=3)

if (length(temp.out)==0){
 df <- as.data.frame(ctable, stringsAsFactors=F)
 even.ind <- seq(2,dim(df)[2],2)
 df_evi <- df[,even.ind]
    RS <- rowSums(df_evi,na.rm=T)} else

{df <- as.data.frame(ctable[-1*temp.out], stringsAsFactors=F)
 even.ind <- seq(2,dim(df)[2],2)
 df_evi <- df[,even.ind]
 RS.temp <- rowSums(df_evi,na.rm=T)
 RS.temp[1] <- RS.temp[1]+length(temp.out)*dim(snps)[1]
    RS <- RS.temp}



snps.counts.1 <- data.frame(Genotype=c("AA", "Aa", "aa"), Frequency=RS, stringsAsFactors=F)

save(snps.counts.1, file="snps_counts_Autophagy1.RData")




axis(1, at=c(0:2), labels=c("AA", "Aa", "aa"))
