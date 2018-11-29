args <- commandArgs(trailingOnly = T)

fname_cell_label <- args[1]
fname_isoverlap <- args[2]
fname_isunmappable <- args[3]
data_dir <- args[4]
fname_ttest_result <- args[5]
gene_num <- as.numeric(args[6])

cell_label <- read.table(fname_cell_label)[,2]
cell_num <- length(cell_label)

tstatistics <- rep(0, gene_num)
log10pval <- rep(0, gene_num)

f_isoverlap <- file(fname_isoverlap, "r")
f_unmappable <- file(fname_isunmappable, "r")

for(g in 1:gene_num){
	if(g%%100==0){
		print(paste("calculating ", g, "-th gene",sep=""))
	}

	a <- readLines(con=f_isoverlap,1)
	isoverlap <- as.numeric(unlist(strsplit(a, "\t")))

	b <- readLines(con=f_unmappable,1)
	isunmappable <- as.numeric(unlist(strsplit(b, "\t")))

	isremove <- isoverlap | isunmappable

	if(inherits(try(data <- as.matrix(read.table(paste(data_dir,"/data_",g,".txt",sep=""))), silent=TRUE), "try-error")){
		next
	}

	data <- log10(data+1)

	if(dim(data)[2] <= 2){
		next
	}
	else if(length(which(apply(data, MARGIN=1, sum)!=0)) < 2){
		next
	}
	else if(sum(isremove==0) < 2){
		next
	}

	data <- data[,isremove==0]

	tmp_cell_idx <- c(1:cell_num)
	if(sum(apply(data, MARGIN=1, sum)!=0) < 2){
		next
	}
	#remove cells which do not countain mapped read count
	else if(sum(apply(data, MARGIN=1, sum)==0) != 0){
		tmp_cell_idx <- tmp_cell_idx[-which(apply(data, MARGIN=1, sum)==0)]
		data <- data[-which(apply(data, MARGIN=1, sum)==0),]
	}

	#remove bins which do not contain mapped read count
	if(sum(apply(data, MARGIN=2, sum)==0) != 0){
		data <- data[,-which(apply(data, MARGIN=2, sum)==0)]
	}

	if(is.matrix(data) == FALSE){
		next
	}

	meandata <- rep(0, cell_num)
	meandata[tmp_cell_idx] <- apply(data, MARGIN=1, mean)
	
	ttest_result <- t.test(meandata[cell_label==1], meandata[cell_label==2])

	if(is.na(ttest_result$statistic)){
		next
	}

	tstatistics[g] <- ttest_result$statistic
	log10pval[g] <- -log10(ttest_result$p.value)
}

write.table(cbind(tstatistics, log10pval), fname_ttest_result, sep="\t", col.names=F, row.names=F)
