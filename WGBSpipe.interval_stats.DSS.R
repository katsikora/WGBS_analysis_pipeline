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
cnpool<-c("CHROM","START","END","STRAND","Name")
colnames(bedtab)<-cnpool[1:ncol(bedtab)]
if(!unique(grepl("STRAND",colnames(bedtab)))){bedtab$STRAND<-"*"}
if(!unique(grepl("Name",colnames(bedtab)))){bedtab$Name<-paste(bedtab$CHROM,bedtab$START,bedtab$END,sep="_")}

sCpGF<-commandArgs(trailingOnly=TRUE)[3] ### don't really neeed single CpG data
message(sprintf("loading %s",sCpGF))
load(sCpGF)

require(GenomicRanges)

bedGR<-GRanges(seqnames=bedtab$CHROM,ranges=IRanges(start=bedtab$START,end=bedtab$END,names=bedtab$Name),strand=bedtab$STRAND)


###################

require(DSS) 
require(data.table)
require(dplyr)

mpath<-commandArgs(trailingOnly=TRUE)[4]
mdir<-dir(mpath,pattern=".CpG.filt2.bed",full.names=TRUE)
limdat.LG.CC<-complete.cases(limdat.LG)

##intersect via Genomic Ranges, get aggregate values per interval

cC<-c("character","integer","NULL","NULL","integer","NULL","integer","NULL")

mlist<-vector("list",length(mdir))
names(mlist)<-mshort

for(i in seq_along(mlist)){
tabi<-fread(mdir[i],header=TRUE,sep="\t",colClasses=cC)
colnames(tabi)<-c("chr","pos","X","N")
tabiGR<-GRanges(seqnames=tabi$chr,ranges=IRanges(start=tabi$pos,width=2),strand="*")
mtch <- as.data.frame(findOverlaps(query=tabiGR, subject=bedGR))
tabi.inCGI<-tabi[mtch$queryHits,]
tabi.inCGI$IntID<-names(ranges(bedGR))[mtch$subjectHits]
CGI.dat<-summarize(group_by(tabi.inCGI,chr,pos),X=integer(mean(X)),N=integer(mean(N)))
mlist[[i]]<-as.data.frame(CGI.dat,stringsAsFactors=FALSE)[,match(c("chr","pos","N","X"),colnames(CGI.dat))]
}

if("Merge" %in% colnames(sampleInfo)){snames<-sampleInfo$PlottingID[match(mshort,sampleInfo$Merge)]}else{snames<-sampleInfo$PlottingID[match(mshort,sampleInfo$SampleID)]}

BS.dat.int<-makeBSseqData(mlist,snames) 
BS.dat.int

save(BS.dat.int,file=paste0(bedshort,".BS.dat.int.RData"))

BS.test<-DMLtest(BS.dat.int,group1=sampleInfo$PlottingID[sampleInfo$Group %in% unique(sampleInfo$Group)[1]],group2=sampleInfo$PlottingID[sampleInfo$Group %in% unique(sampleInfo$Group)[2]])
BS.res<-callDML(BS.test,p.threshold=0.05)

BS.res$ID<-paste(BS.res$chr,BS.res$pos,sep="_")

save(BS.test,BS.res,file=paste0(bedshort,"DSS.res.RData"))


