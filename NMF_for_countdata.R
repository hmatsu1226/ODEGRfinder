library(NMF)

args <- commandArgs(trailingOnly = T)

fname_isoverlap <- args[1]
fname_isunmappable <- args[2]
data_dir <- args[3]
out_dir <- args[4]

init_gene_idx <- as.numeric(args[5])
end_gene_idx <- as.numeric(args[6])
gene_num <- end_gene_idx - init_gene_idx + 1

cell_num <- as.numeric(args[7])
K <- as.numeric(args[8])
myseed <- as.numeric(args[9])
nmfmethod <- args[10]

NMF_coef <- matrix(rep(0, cell_num*gene_num*K), nrow=cell_num, ncol=gene_num*K)

f_isoverlap <- file(fname_isoverlap, "r")
f_unmappable <- file(fname_isunmappable, "r")

if(init_gene_idx != 1){
	for(i in 1:(init_gene_idx-1)){
		readLines(con=f_isoverlap,1)
		readLines(con=f_unmappable,1)
	}
}

for(g in 1:gene_num){
	if(g%%100==0){
		print(paste("calculating ", g, "-th gene",sep=""))
	}

	a <- readLines(con=f_isoverlap,1)
	isoverlap <- as.numeric(unlist(strsplit(a, "\t")))

	b <- readLines(con=f_unmappable,1)
	isunmappable <- as.numeric(unlist(strsplit(b, "\t")))

	isremove <- isoverlap | isunmappable

	if(inherits(try(data <- as.matrix(read.table(paste(data_dir,"/data_",init_gene_idx+g-1,".txt",sep=""))), silent=TRUE), "try-error")){
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
	if(dim(data)[1] < 10){
		next
	}

	if(inherits(try(res <- nmf(t(data), rank=K, seed=myseed, method=nmfmethod), silent=TRUE), "try-error")){
		next
	}

	for(k in 1:K){
		NMF_coef[tmp_cell_idx, K*g-(K-k)] <- signif(coef(res)[k,], digits=3)
	}
}

write.table(NMF_coef, paste(out_dir,"/NMF_",K,"_coef_",init_gene_idx,"_",end_gene_idx,".txt",sep=""), sep="\t", col.names=F, row.names=F)
