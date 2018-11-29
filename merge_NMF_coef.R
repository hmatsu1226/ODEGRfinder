args <- commandArgs(trailingOnly = T)

data_dir <- args[1]
cell_num <- as.numeric(args[2])
gene_num <- as.numeric(args[3])
K <- as.numeric(args[4])
init_gene_idx_list <- args[5]
end_gene_idx_list <- args[6]

file_num <- length(strsplit(init_gene_idx_list, ",")[[1]])

NMF_coef <- matrix(rep(0, cell_num*gene_num*K), nrow=cell_num, ncol=gene_num*K)
for(i in 1:file_num){
	init_gene_idx <- as.numeric(strsplit(init_gene_idx_list, ",")[[1]][i])
	end_gene_idx <- as.numeric(strsplit(end_gene_idx_list, ",")[[1]][i])

	tmp <- as.matrix(read.table(paste(data_dir,"/NMF_",K,"_coef_",init_gene_idx,"_",end_gene_idx,".txt",sep="")))

	NMF_coef[,(K*init_gene_idx-(K-1)):(K*end_gene_idx)] <- tmp
}
write.table(NMF_coef, paste(data_dir,"/NMF_",K,"_coef_all.txt",sep=""), sep="\t", col.names=F, row.names=F)
