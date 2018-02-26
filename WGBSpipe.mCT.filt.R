#run in R-3.3.1
#a few lines of code from methylCtools bcall2beta by Hovestadt et al. 2014 were retained 
#CpG position handling and coverage and Beta calculations by Katarzyna Sikora
#set working directory
wdir<-commandArgs(trailingOnly=TRUE)[1]
#system(paste0('mkdir -p ',wdir)) #for debugging
setwd(wdir)
message(sprintf("working directory is %s",getwd()))

#methylation table
mpath<-commandArgs(trailingOnly=TRUE)[2]
mshort<-basename(mpath)
message(paste0("processing ",mpath))

## filtering thresholds
a.snp <- 0.25     # threshold for SNP filtereing (allelic frequency of illegal base)
a.cov <- 10        # min coverage on both strands (KS)

#read in methylation data; use data.table
require(data.table)

dm<-fread(mpath,header=FALSE,quote="",sep="\t")
colnames(dm)<-c("CHROM", "POS", "STRAND", "context", "snp", "M", "U")
dm

##filtering out 'orphan' CpGs
dm.posV<-vector("numeric",nrow(dm))
for(i in 1:nrow(dm)){
    ifelse(dm$STRAND[i]=="+",dm.posV[i]<-dm$POS[i],dm.posV[i]<-(dm$POS[i]-1))
}

#table(ave(dm.posV,factor(dm.posV),FUN=length))
dm.filtV2<-ave(dm.posV,factor(dm.posV),FUN=length)
dm<-dm[dm.filtV2!=1]

nori<-nrow(dm)/2
message(paste0(nori," CpGs were extracted"))

##summarize counts and clean up table
dm2<-data.table(matrix(nrow=table(dm$STRAND)["+"],ncol=9))
colnames(dm2)<-c("CHROM","START","END","Beta","M","U","Cov","snp","ms")
dm2$CHROM<-dm$CHROM[dm$STRAND=="+"]
dm2$START<-dm$POS[dm$STRAND=="+"]
dm2$END<-dm$POS[dm$STRAND=="+"]+1
dm2

##check if diff always 1
#dm2.diffV<-dm2$END-dm2$START
#table(dm2.diffV)

dm2$ms<-paste(dm2$CHROM,dm2$START,sep="_")
dm$ms<-"NA"
dm$ms[dm$STRAND=="+"]<-paste(dm$CHROM[dm$STRAND=="+"],dm$POS[dm$STRAND=="+"],sep="_")
dm$ms[dm$STRAND=="-"]<-paste(dm$CHROM[dm$STRAND=="-"],(dm$POS[dm$STRAND=="-"]-1),sep="_")

require(dplyr)
dmAgg<-summarize(group_by(dm,ms),Msum=sum(M),Usum=sum(U),SNPsum=max(snp))
dm2[,c("M","U","snp")]<-dmAgg[match(dm2$ms,dmAgg$ms),c("Msum","Usum","SNPsum")]


message("filtering SNPs ..", appendLF=FALSE)
m <- which(dm2$snp > a.snp)
dm2[m, "M"]<-NA
dm2[m, "U"]<-NA
message(" done")
dm2
nsnp<-length(m)
pct.snp<-nsnp/nori*100

message(sprintf("%i (%.2f%%) CpGs were filtered out due to SNP sites",nsnp,pct.snp))


rm(dm)
invisible(gc())
message(" done")

head(dm2)
dm2$Cov<-rowSums(dm2[,c("M","U")])
head(dm2)

## calculate beta values
message("calculating beta values ..", appendLF=FALSE)
dm2$Beta<- dm2$M/dm2$Cov
message(" done")

##filter for low coverage
message("filtering low coverage calls ..", appendLF=FALSE)
dm2$Beta[dm2$Cov < a.cov] <- NA
message(" done")
tot<-nori/2
filt.snp<-nsnp
filt.cov<-sum(is.na(dm2$Beta))-filt.snp
pct.cov<-filt.cov/tot*100

message(sprintf("%i (%.2f%%) CpGs were filtered out to coverage less than %i (sum on both strands)",filt.cov,pct.cov,a.cov))

dm2<-dm2[,!colnames(dm2) %in% "snp",with=FALSE]
## write output files
message("writing output files ..", appendLF=FALSE)

write.table(dm2,sep="\t", row.names=FALSE, quote=FALSE, file=paste0(wdir, "/", gsub(".CG.call",".CpG.filt2.bed",mshort, )))

message("done all")

