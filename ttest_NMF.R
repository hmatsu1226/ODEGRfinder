args <- commandArgs(trailingOnly = T)

fname_cell_label <- args[1]
fname_nmf_coef <- args[2]
fname_ttest_result <- args[3]

gene_num <- as.numeric(args[4])
K <- as.numeric(args[5])

cell_label <- read.table(fname_cell_label)[,2]
cell_num <- length(cell_label)

NMF_coef <- as.matrix(read.table(fname_nmf_coef))

tstatistics <- matrix(rep(0, gene_num*K), nrow=gene_num, ncol=K)
log10pval <- matrix(rep(0, gene_num*K), nrow=gene_num, ncol=K)

for(g in 1:gene_num){
	for(k in 1:K){
		ttest_result <- t.test(NMF_coef[cell_label==1, g*K-(K-k)], NMF_coef[cell_label==2, g*K-(K-k)])

		if(is.na(ttest_result$statistic)){
			next
		}

		tstatistics[g,k] <- ttest_result$statistic
		log10pval[g,k] <- -log10(ttest_result$p.value)
	}
}

write.table(cbind(tstatistics, log10pval), fname_ttest_result, sep="\t", col.names=F, row.names=F)
