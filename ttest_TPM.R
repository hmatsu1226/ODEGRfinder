library("stringr")

args <- commandArgs(trailingOnly = T)

fname_tpm <- args[1]
fname_tpm_transcriptid <- args[2]
fname_cell_label <- args[3]
fname_gtf <- args[4]
fname_transcriptid_to_geneid <- args[5]

fname_ttest_result <- args[6]

TPM <- as.matrix(read.table(fname_tpm, header=T, row.names=1))
TPM_transcriptid <- read.table(fname_tpm_transcriptid)[,1]
cell_label <- read.table(fname_cell_label)[,2]
mygtf <- read.table(fname_gtf, sep="\t")
transcriptid_to_geneid = read.table(fname_transcriptid_to_geneid)

rownames(TPM) <- TPM_transcriptid
gene_num <- dim(mygtf)[1]

TPM_tstatistic_pos <- rep(0, gene_num)
TPM_tstatistic_neg <- rep(0, gene_num)
TPM_log10pval_pos <- rep(0, gene_num)
TPM_log10pval_neg <- rep(0, gene_num)
TPM_transcriptid_pos <- rep("-", gene_num)
TPM_transcriptid_neg <- rep("-", gene_num)

for(g in 1:gene_num){
	if(g%%100==0){
		print(paste("calculating ", g, "-th gene",sep=""))
	}

	gene_id = str_sub(strsplit(as.character(mygtf[g,9]), " ")[[1]][2], end=-2)
	transcript_idx_list = which(transcriptid_to_geneid[,2] == gene_id)

	for(t in 1:length(transcript_idx_list)){
		transcript_id <- as.character(transcriptid_to_geneid[transcript_idx_list[t],1])

		if(length(which(rownames(TPM)==transcript_id)) == 0){
			next
		}
		
		ttest_result <- t.test(log10(TPM[transcript_id, cell_label==1]+1), log10(TPM[transcript_id, cell_label==2]+1),var.equal=F,paired=F)

		if(is.nan(ttest_result$statistic)){
			next
		}

		if(ttest_result$statistic > 0){
			if(TPM_tstatistic_pos[g] < ttest_result$statistic){
				TPM_tstatistic_pos[g] <- ttest_result$statistic
				TPM_log10pval_pos[g] <- -log10(ttest_result$p.value)
				TPM_transcriptid_pos[g] <- transcript_id
			}
		}
		else{
			if(TPM_tstatistic_neg[g] > ttest_result$statistic){
				TPM_tstatistic_neg[g] <- ttest_result$statistic
				TPM_log10pval_neg[g] <- -log10(ttest_result$p.value)
				TPM_transcriptid_neg[g] <- transcript_id
			}
		}
	}
}

write.table(cbind(round(TPM_tstatistic_pos, digits=3), round(TPM_tstatistic_neg, digits=3), round(TPM_log10pval_pos, digits=3), round(TPM_log10pval_neg, digits=3), TPM_transcriptid_pos, TPM_transcriptid_neg), fname_ttest_result, sep="\t", quote=F, row.names=F, col.names=F)
