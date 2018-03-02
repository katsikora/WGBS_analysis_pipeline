#run in R-3.3.1
#set working directory
wdir<-commandArgs(trailingOnly=TRUE)[1]
#system(paste0('mkdir -p ',wdir)) #for debugging
setwd(wdir)
message(sprintf("working directory is %s",getwd()))

options(stringsAsFactors=FALSE,na.rm=TRUE)

bedF<-commandArgs(trailingOnly=TRUE)[2]
message(sprintf("processing %s",bedF))
bedshort<-gsub(".bed","",basename(bedF))

bedtab<-read.table(bedF,header=FALSE,sep="\t",as.is=TRUE,quote="")
colnames(bedtab)<-c("CHROM","START","END","qvalue","MeanDiff","NumCpGs","pMWU","p2DKS","meanA","meanB")
bedtab$Name<-paste(bedtab$CHROM,bedtab$START,bedtab$END,sep="_")

#load single CpG data, count NAs per row/CpG
sCpGF<-commandArgs(trailingOnly=TRUE)[3]
message(sprintf("loading %s",sCpGF))
require(data.table)
load(sCpGF)

noNA<-apply(limdat.LG,1,function(X)sum(is.na(X)))
NAf<-ifelse(noNA>1,1,0)
limdat.LG$NAf<-NAf

require(GenomicRanges)

bedGR<-GRanges(seqnames=bedtab$CHROM,ranges=IRanges(start=bedtab$START,end=bedtab$END,names=paste(bedtab$CHROM,bedtab$START,bedtab$END,sep="_")))
###NOT EDITED
limdatGR<-GRanges(seqnames=gsub("_.+","",limdat.LG$ms),ranges=IRanges(start=as.numeric(gsub(".+_","",limdat.LG$ms)),width=2,names=limdat.LG$ms),strand="*")
mtch <- as.data.frame(findOverlaps(query=limdatGR, subject=bedGR))

limdat.LG.inCGI<-limdat.LG[mtch$queryHits,]
limdat.LG.inCGI$IntID<-names(ranges(bedGR))[mtch$subjectHits]

NA.inCGI<-with(limdat.LG.inCGI,ave(NAf,factor(IntID),FUN=sum))
limdat.LG.inCGI$NA.inCGI<-NA.inCGI

to.m<-limdat.LG.inCGI[,c("IntID", "NA.inCGI"),with=FALSE]

CGI.map<-unique(to.m)
bedtab$N.CG.NA<-CGI.map$NA.inCGI[match(bedtab$Name,CGI.map$IntID)]

N.CG.tot<-with(limdat.LG.inCGI,ave(NAf,IntID,FUN=length))
bedtab$N.CG.tot<-N.CG.tot[match(bedtab$Name,limdat.LG.inCGI$IntID)]

bedtab$CGI.NAf<-ifelse(bedtab$N.CG.NA>(0.8*bedtab$N.CG.tot),NA,1)
bedtab.CC<-bedtab[complete.cases(bedtab),]

CGI.limdat<-as.data.frame(apply(limdat.LG.inCGI[,2:(ncol(limdat.LG.inCGI)-3)],2,function(X){ave(X,factor(limdat.LG.inCGI$IntID),FUN=function(X)mean(X,na.rm=TRUE))}),stringsAsFactors=FALSE)

CGI.limdat$IntID<-limdat.LG.inCGI$IntID
CGI.limdat<-unique(CGI.limdat)
rownames(CGI.limdat)<-CGI.limdat$IntID
CGI.limdat<-CGI.limdat[,-grep("IntID",colnames(CGI.limdat))]

CGI.limdat.CC<-CGI.limdat[bedtab.CC$Name,]



####for differential interval methylation
### limma + ebayes + BH p value adjustment

require("limma")
require("car")
require("FactoMineR")
require("reshape2")
require("ggplot2")
require("dplyr")

CGI.limdat.CC.logit<-logit(CGI.limdat.CC,percents=FALSE,adjust=0.025) 
x1<-PCA(CGI.limdat.CC,graph=FALSE)

pdf(paste0(bedshort,".CGI.limdat.CC.PCA.pdf"),paper="a4",bg="white") 
plot.PCA(x1,choix="var")
dev.off()

#calculate row means
spath<-commandArgs(trailingOnly=TRUE)[4]
sampleInfo<-read.table(spath,header=TRUE,sep="\t",as.is=TRUE)
#calculate and save row means
CGI.limdat.CC$IntID<-rownames(CGI.limdat.CC)
CGI.limdat.CC.L<-melt(CGI.limdat.CC,id.vars="IntID",value.name="Beta",variable.name="SampleID")
CGI.limdat.CC.L$Group<-sampleInfo$Group[match(CGI.limdat.CC.L$SampleID,sampleInfo$PlottingID)]
CGI.limdat.CC.Means<-data.table(summarize(group_by(CGI.limdat.CC.L,IntID,Group),Beta.Mean=mean(Beta)))

##density plots
ggplot(data=CGI.limdat.CC.Means,aes(x=Beta.Mean))+geom_density(aes(group=Group,colour=Group,fill=Group),alpha=0.3)+ggtitle("Differentially methylated regions")+
    theme(text = element_text(size=16),axis.text = element_text(size=12),axis.title = element_text(size=14))+xlab("Mean methylation ratio")
ggsave(paste0(bedshort,".Beta.MeanXgroup.metilene.dens.png"))

##violin plots
ggplot(data=CGI.limdat.CC.Means)+geom_violin(aes(x=Group,y=Beta.Mean,fill=Group))+geom_boxplot(aes(x=Group,y=Beta.Mean),width=0.1)+ggtitle("Differentially methylated regions")+
    theme(text = element_text(size=16),axis.text = element_text(size=12),axis.title = element_text(size=14))+xlab("Mean methylation ratio")
ggsave(paste0(bedshort,".Beta.MeanXgroup.metilene.violin.png"))

#differential methylation
design<-as.data.frame(matrix(ncol=2,nrow=(ncol(CGI.limdat.CC.logit))),stringsAsFactors=FALSE)
colnames(design)<-c("Intercept","Group")
rownames(design)<-colnames(CGI.limdat.CC.logit)
design$Group<-as.numeric(factor(sampleInfo$Group[match(colnames(CGI.limdat.CC.logit),sampleInfo$PlottingID)]))
design$Intercept<-1
design<-as.matrix(design)

fit<-lmFit(CGI.limdat.CC.logit,design)
fit.eB<-eBayes(fit)
tT.FDR5<-topTable(fit.eB,2,p.value=0.05,number=Inf)[,c("logFC","t","adj.P.Val","B")]
write.table(tT.FDR5,file=paste0(bedshort,".CGI.limdat.CC.tT.FDR5.txt"),sep="\t",quote=FALSE)

nrow(tT.FDR5)
nrow(CGI.limdat.CC.logit)
nrow(tT.FDR5)/nrow(CGI.limdat.CC.logit)

#CGI.limdat.CC.Diff<-summarize(group_by(CGI.limdat.CC.Means,IntID),Diff=(Beta.Mean[1]-Beta.Mean[2]))
#tT.FDR5.an<-as.data.frame(merge(x=tT.FDR5,y=CGI.limdat.CC.Diff,by.x="row.names",by.y="IntID",all.x=TRUE,sort=FALSE))

#annotate metilene output with this information
CGI.bed.intT<-as.data.frame(merge(x=bedtab,y=tT.FDR5,by.x="Name",by.y="row.names",sort=FALSE,all.x=TRUE))
CGI.bed.intT<-CGI.bed.intT[,!colnames(CGI.bed.intT) %in% "Name"]
write.table(CGI.bed.intT,file=paste0(bedshort,".metilene.limma.bed"),sep="\t",quote=FALSE)
save(CGI.bed.intT,file=paste0(bedshort,".metilene.limma.RData"))

