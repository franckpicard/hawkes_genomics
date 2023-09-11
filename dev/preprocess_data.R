rm(list=ls())

source("hawkesfun.R")

chrlist      = c("chrI","chrII","chrIII","chrIV","chrV","chrX")
path_data    = ""
path_results = ""
enh  = read.table(paste0(path_data,"enhancer.bed")
prom = read.table(paste0(path_data,"coding_prom.bed")
sin3 = read.table(paste0(path_data,"sin3peaks.bed")

enh[enh[,1]=="chrI",1]   = "chr1"
enh[enh[,1]=="chrII",1]  = "chr2"
enh[enh[,1]=="chrIII",1] = "chr3"
enh[enh[,1]=="chrIV",1]  = "chr4"
enh[enh[,1]=="chrV",1]   = "chr5"

sin3[sin3[,1]=="chrI",1]   = "chr1"
sin3[sin3[,1]=="chrII",1]  = "chr2"
sin3[sin3[,1]=="chrIII",1] = "chr3"
sin3[sin3[,1]=="chrIV",1]  = "chr4"
sin3[sin3[,1]=="chrV",1]   = "chr5"

prom[prom[,1]=="chrI",1]   = "chr1"
prom[prom[,1]=="chrII",1]  = "chr2"
prom[prom[,1]=="chrIII",1] = "chr3"
prom[prom[,1]=="chrIV",1]  = "chr4"
prom[prom[,1]=="chrV",1]   = "chr5"

maxlag    = 1e5

dd        = list(enh,prom,sin3)
names(dd) = c("enh","prom","sin3")
regions   = get_regions_from_bed(dd,maxlag,check.valid=FALSE)

z       = enh
z.bedr  = df2bedr(z)
z.bedr  = bedr.sort.region(z.bedr, check.valid=FALSE)  
write.table(convert2bed(z.bedr),file=paste0(path_results,"enh_preprocessed.bed"),quote=F,col.names=F,row.names=F,sep="\t")

z       = prom
z.bedr  = df2bedr(z)
z.bedr  = bedr.sort.region(z.bedr, check.valid=FALSE)  
write.table(convert2bed(z.bedr),file=paste0(path_results,"prom_preprocessed.bed"),quote=F,col.names=F,row.names=F,sep="\t")

z       = sin3
z.bedr  = df2bedr(z)
z.bedr  = bedr.sort.region(z.bedr, check.valid=FALSE)  
write.table(convert2bed(z.bedr),file=paste0(path_results,"sin3_preprocessed.bed"),quote=F,col.names=F,row.names=F,sep="\t")


