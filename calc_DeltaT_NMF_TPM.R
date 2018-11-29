args <- commandArgs(trailingOnly = T)

nmf_dir <- args[1]
tpm_dir <- args[2]
fname_out <- args[3]
gene_num <- as.numeric(args[4])
thresh <- as.numeric(args[5])

Kset <- c(2,5,10)

myscore <- rep(-Inf, gene_num)

ttest_result_tpm <- as.matrix(read.table(paste(tpm_dir,"/ttest_result_TPM.txt",sep="")))

for(i in 1:length(Kset)){
	K <- Kset[i]
	ttest_result_nmf <- as.matrix(read.table(paste(nmf_dir,"/ttest_result_NMF_",K,".txt",sep="")))
	
	for(g in 1:gene_num){
		if(is.na(ttest_result_nmf[g,1])==TRUE){
			next
		}

		TPM_T_max <- max(0, as.numeric(ttest_result_tpm[g,1]))
		TPM_T_min <- min(0, as.numeric(ttest_result_tpm[g,2]))

		for(k in 1:K){
			if(ttest_result_nmf[g,k] >= 0){
				if(TPM_T_max > thresh){
					myscore[g] <- max(myscore[g], min(0,ttest_result_nmf[g,k]-TPM_T_max))
				}
				else{
					myscore[g] <- max(myscore[g], ttest_result_nmf[g,k]-TPM_T_max)
				}
			}
			else if(ttest_result_nmf[g,k] < 0){
				if((-TPM_T_min) > thresh){
					myscore[g] <- max(myscore[g], min(0,-(ttest_result_nmf[g,k]-TPM_T_min)))
				}
				else{
					myscore[g] <- max(myscore[g], -(ttest_result_nmf[g,k]-TPM_T_min))
				}
			}
		}
	}
}

myscore_rank <- rank(-myscore)

write.table(cbind(myscore, myscore_rank), fname_out, col.names=F, row.names=F, sep="\t")
