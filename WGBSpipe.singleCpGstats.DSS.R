#run in R-3.3.1
#set working directory
wdir<-commandArgs(trailingOnly=TRUE)[1]
#system(paste0('mkdir -p ',wdir)) #for debugging
setwd(wdir)
message(sprintf("working directory is %s",getwd()))

options(stringsAsFactors=FALSE,na.rm=TRUE)

spath<-commandArgs(trailingOnly=TRUE)[2]
sampleInfo<-read.table(spath,header=TRUE,sep="\t",as.is=TRUE)


require(DSS) 
require(data.table)

mpath<-commandArgs(trailingOnly=TRUE)[3]
mdir<-dir(mpath,pattern=".CpG.filt2.bed",full.names=TRUE)
mshort<-gsub(".CpG.filt2.bed","",basename(mdir))

#reformat the tables to DSS format
##define colClasses per extraction software
cC<-c("character","integer","NULL","NULL","integer","NULL","integer","NULL")

mlist<-vector("list",length(mdir))
names(mlist)<-mshort

for(i in seq_along(mlist)){
tabi<-fread(mdir[i],header=TRUE,sep="\t",colClasses=cC)
colnames(tabi)<-c("chr","pos","X","N")
mlist[[i]]<-as.data.frame(tabi,stringsAsFactors=FALSE)[,match(c("chr","pos","N","X"),colnames(tabi))]
}

if("Merge" %in% colnames(sampleInfo)){snames<-sampleInfo$PlottingID[match(mshort,sampleInfo$Merge)]}else{snames<-sampleInfo$PlottingID[match(mshort,sampleInfo$SampleID)]}

BS.dat<-makeBSseqData(mlist,snames)
save(BS.dat,file="singleCpG.RData")


message("testing without smoothing")
BS.test<-DMLtest(BS.dat,group1=sampleInfo$PlottingID[sampleInfo$Group %in% unique(sampleInfo$Group)[1]],group2=sampleInfo$PlottingID[sampleInfo$Group %in% unique(sampleInfo$Group)[2]])
BS.res<-callDML(BS.test,p.threshold=0.05)
BS.res$ms<-paste(BS.res$chr,BS.res$pos,sep="_")
save(BS.test,BS.res,file="DSS.res.RData")

dmrs<-callDMR(BS.test, p.threshold=0.05)
write.table(dmrs,file="DSS.DMR.tab.pval0.05.txt",row.names=FALSE,sep="\t",quote=FALSE)
