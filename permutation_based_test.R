args <- commandArgs(trailingOnly = T)

fname_score <- args[1]
fname_score_shuffle <- args[2]
fname_out <- args[3]

myscore <- read.table(fname_score)[,1]
myscore_shuffle <- read.table(fname_score_shuffle)[,1]

pval <- rep(0, length(myscore))
for(g in 1:length(myscore)){
	pval[g] <- -log10(pnorm(myscore[g], mean=mean(myscore_shuffle), sd=sd(myscore_shuffle), lower.tail=FALSE))
}

write.table(cbind(myscore, pval), fname_out, col.names=F, row.names=F, sep="\t")
