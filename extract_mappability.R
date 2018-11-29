library("GenomicRanges")
library("rtracklayer")
library("Rsamtools")

args <- commandArgs(trailingOnly = T)

fname_umap_bw <- args[1]
fname_gtf <- args[2]
fname_isunmappable <- args[3]
thresh <- as.numeric(args[4])

mygtf <- read.table(fname_gtf, sep="\t")
gene_num <- dim(mygtf)[1]

minmappabilityIntoFixedBins <- function(cov, binsize){
  sapply(seq(1,length(cov), by = binsize), function(x){
    min(cov[x:min((x+binsize-1),length(cov))], na.rm=T)
  })
}

for(g in 1:gene_num){
	if((g%%100)==0){
		print(paste("calculating ", g, "-th gene",sep=""))
	}

	mychr <- as.character(mygtf[g,1])
	mystart <- mygtf[g,4]
	myend <- mygtf[g,5]

	select <- GRanges(seqnames = mychr, ranges = IRanges(start = mystart, end = myend))

	binsize <- 100
	nbin <- floor((myend - mystart)/binsize)+1

	#if((myend-mystart) > 500000){
	#	binsize <- floor((myend-mystart+1)/5000) + 1
	#	nbin <- floor((myend - mystart)/binsize) + 1
	#}

	mappability <- rep(0, nbin)
	cov <- numeric(width(select))  

 	#select@seqinfo@seqnames %in% seqnames(seqinfo(BigWigFile(fname_umap_bw)))
	    
	cov <- import(fname_umap_bw, selection = select, as = 'NumericList')[[1]]
	mappability <- minmappabilityIntoFixedBins(cov, binsize)
	mappability[mappability>thresh]<- 1
	mappability[mappability<=thresh]<- 0

	mappability <- 1-mappability

	write.table(t(mappability), fname_isunmappable, sep="\t", col.names=F, row.names=F, append=T)
}
