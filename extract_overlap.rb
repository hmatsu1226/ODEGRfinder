#!/usr/bin/ruby

fin_all_gtf = open(ARGV[0],"r")
fin_my_gtf = open(ARGV[1],"r")
fout = open(ARGV[2],"w")

my_gene_num = ARGV[3].to_i

all_gene_num = 0
all_gene_chr = Array.new()
all_gene_left = Array.new()
all_gene_right = Array.new()

while line=fin_all_gtf.gets
	tmp = line.split("\t")

	if tmp[2] != "gene"
		next
	end

	all_gene_chr.push(tmp[0])
	all_gene_left.push(tmp[3].to_i)
	all_gene_right.push(tmp[4].to_i)

	all_gene_num = all_gene_num + 1
end


for g in 0..my_gene_num-1
	line=fin_my_gtf.gets
	tmp = line.split("\t")

	left = tmp[3].to_i
	right = tmp[4].to_i

	binsize = 100
	nbin = ((right - left)/binsize).floor+1
	
	#if (right-left) > 500000
	#	binsize = ((right-left+1)/5000).floor + 1
	#	nbin = ((right - left)/binsize).floor + 1
	#end

	overlap_check = Array.new(nbin, 0)

	for i in 0..all_gene_num-1
		if all_gene_chr[i] != tmp[0]
			next
		end

		if all_gene_left[i] > right
			break
		end

		#same gene
		if all_gene_left[i] == left && all_gene_right[i] == right
			next
		end

		#fully cover
		if all_gene_right[i] >= right && all_gene_left[i] <= left
			for j in 0..nbin-1
				overlap_check[j] = 1
			end
			break
		#right overlapped
		elsif all_gene_right[i] >= left && all_gene_left[i] <= left
			for j in left..all_gene_right[i]
				binid = ((j - left)/binsize).floor
				overlap_check[binid] = 1
			end
		#left overlapped
		elsif all_gene_right[i] >= right && all_gene_left[i] <= right
			for j in all_gene_left[i]..right
				binid = ((j - left)/binsize).floor
				overlap_check[binid] = 1
			end
		#inside gene
		elsif all_gene_right[i] <= right && all_gene_left[i] >= left
			for j in all_gene_left[i]..all_gene_right[i]
				binid = ((j - left)/binsize).floor
				overlap_check[binid] = 1
			end
		end

	end

	out = ""
	for i in 0..nbin-1
		out += "#{overlap_check[i]}\t"
	end
	fout.puts(out.chop)
end
