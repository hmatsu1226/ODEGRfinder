library("GenomicRanges")
library("rtracklayer")
library("Rsamtools")

args <- commandArgs(trailingOnly = T)

fname_gtf <- args[1]
fname_bw_path <- args[2]
out_dir <- args[3]
init_gene_idx <- as.numeric(args[4])
end_gene_idx <- as.numeric(args[5])

gene_num <- end_gene_idx - init_gene_idx + 1

mygtf <- read.table(fname_gtf, sep="\t")

coverageIntoFixedBins <- function(cov, binsize){
  sapply(seq(1,length(cov), by = binsize), function(x){
    mean(cov[x:min((x+binsize-1),length(cov))], na.rm=T)
  })
}

for(g in 1:gene_num){
	if((g%%100)==0){
		print(paste("calculating ", g, "-th gene",sep=""))
	}

	mychr <- as.character(mygtf[init_gene_idx+g-1,1])
	mystart <- as.numeric(mygtf[init_gene_idx+g-1,4])
	myend <- as.numeric(mygtf[init_gene_idx+g-1,5])

	select <- GRanges(seqnames = mychr, ranges = IRanges(start = mystart, end = myend))

	path_bw_files <- as.character(read.table(fname_bw_path)[,1])
	binsize <- 100
	nbin <- floor((myend - mystart)/binsize)+1

	#if((myend-mystart) > 500000){
	#	binsize <- floor((myend-mystart+1)/5000) + 1
	#	nbin <- floor((myend - mystart)/binsize) + 1
	#}

	mat <- matrix(0, ncol = nbin, nrow = length(path_bw_files))
	cov <- numeric(width(select))
	  
	for(i in seq_along(path_bw_files)){
		path_bw <- path_bw_files[i]
	    
	    cov <- import(path_bw, selection = select, as = 'NumericList')[[1]]
	    mat[i,] <- coverageIntoFixedBins(cov, binsize)
	}

	write.table(mat, paste(out_dir,"/data_",init_gene_idx+g-1,".txt",sep=""), sep="\t", col.names=F, row.names=F)
}
