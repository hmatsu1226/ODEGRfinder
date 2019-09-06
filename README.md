# ODEGRfinder

An NMF-based approach to discover overlooked differentially expressed gene regions from single-cell RNA-seq data

## Reference

Under submission

## Requirements

ODEGRfinder is mainly written with R and uses "NMF" library. We also use "GenomicRanges", "rtracklayer", "Rsamtools", and "stringr" library in data pre-processing. In the manuscript, the following package versions were used.
```
NMF: 0.20.5
GenomicRanges: 1.30.3
rtracklayer: 1.38.3
Rsamtools: 1.30.0
stringr: 1.2.0
```

## Download

```
git clone https://github.com/hmatsu1226/ODEGRfinder
cd ODEGRfinder
```
Or download from "Download ZIP" button and unzip it.

## Dataset in the manuscript
The mES-PrE dataset and hNSC-NC dataset are avaiable at https://doi.org/10.6084/m9.figshare.7410509.v1 and https://doi.org/10.6084/m9.figshare.7410512, respectively.

## A small dataset to demo the code
Running ODEGRfinder with a small dataset (20 genes and 185 cells).
The computational time is about 5 minutes with MacBook Pro (2.5GHz Intel Core i7 and 16GB 2133MHz LPDDR3). 
The detailed explanation of each step are described at [NMF](#nmf), [t-test](#ttest), and [Calculate Delta(Tnmf-Ttpm)](#deltaT).


#### NMF for read count matrix
```
Rscript NMF_for_countdata.R demo_data/isoverlap.txt demo_data/isunmappable.txt demo_count_data demo_out 1 20 185 2 123456 lee
Rscript NMF_for_countdata.R demo_data/isoverlap.txt demo_data/isunmappable.txt demo_count_data demo_out 1 20 185 5 123456 lee
Rscript NMF_for_countdata.R demo_data/isoverlap.txt demo_data/isunmappable.txt demo_count_data demo_out 1 20 185 10 123456 lee
```

#### t-test for NMF results
```
Rscript ttest_NMF.R demo_data/cell_label.txt demo_out/NMF_2_coef_1_20.txt demo_out/ttest_result_NMF_2.txt 20 2
Rscript ttest_NMF.R demo_data/cell_label.txt demo_out/NMF_5_coef_1_20.txt demo_out/ttest_result_NMF_5.txt 20 5
Rscript ttest_NMF.R demo_data/cell_label.txt demo_out/NMF_10_coef_1_20.txt demo_out/ttest_result_NMF_10.txt 20 10
```

#### t-test for TPM matrix
```
Rscript ttest_TPM.R demo_data/TPM_ES_PrE.txt demo_data/TPM_transcriptid.txt demo_data/cell_label.txt demo_data/mygtf_gene.txt demo_data/transcriptid_to_geneid.txt demo_out/ttest_result_TPM.txt
```

#### Calculate Delta(Tnmf-Ttpm)
```
Rscript calc_DeltaT_NMF_TPM.R demo_out demo_out demo_out/DeltaT_NMF_TPM.txt 20 10
```

# Pre-processing
## Make read count matrix
Make read count matrix for each gene region from bigwig files of scRNA-seq.

#### Usage
```
Rscript extract_count_data_from_bw.R <Input_file1> <Input_file2> <Output_dir> <idx1> <idx2>
```

* Input_file1 : annotation file of target gene regions
* Input_file2 : paths to bigwig files
* Output_dir : output directory
* idx1 : index 1
* idx2 : index 2

##### Example
```
Rscript extract_count_data_from_bw.R ES_PrE/data/mygtf_gene.txt ES_PrE/fbw_ES_PrE.txt ES_PrE/count_data 1 1000
```

#### Format of Input_file1
Input_file1 is the gtf file of length G, where G is the number of target gene regions.
Each row represents a region of a gene.

#### Example of Input_file1
```
chr1	HAVANA	gene	5070018	5162529	.	+	.	gene_id "ENSMUSG00000033793.12"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "Atp6v1h"; level 2; havana_gene "OTTMUSG00000050145.9";
chr1	HAVANA	gene	6206197	6276648	.	+	.	gene_id "ENSMUSG00000025907.14"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "Rb1cc1"; level 2; havana_gene "OTTMUSG00000033467.12";
chr1	HAVANA	gene	7088920	7173628	.	+	.	gene_id "ENSMUSG00000051285.17"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "Pcmtd1"; level 2; havana_gene "OTTMUSG00000043373.5";
...
```

#### Format of Input_file2
Input_file2 is the list of the paths for bigwig file of each cell.
Because of the large size of bigwig files, the example datasets in the figshare do not contain this file and bigwig files.

#### Example of Input_file2
```
ES_PrE/data/bw/RamDA_00h_A04.bw
ES_PrE/data/bw/RamDA_00h_A05.bw
ES_PrE/data/bw/RamDA_00h_A06.bw
...
```

#### Output_dir
The read count matrix for each gene region are saved in the Output_dir directory.
The Output_dir/data_n.txt represents the C x L count matrix for the n-th gene region in Input_file1, where C is the number of cells and L is the number of bins for n-th gene.

#### idx1 and idx2
The read count matrix from idx1-th row to idx2-th row in "Input_file1" are calculated.
The count matrix can be calculated in parallel by dividing the start and end index.


## Make "isoverlap.txt" file
Make "isoverlap.txt" file to check whether a bin in a gene region is overlapped with different genes.

#### Usage
```
ruby extract_overlap.rb <Input_file1> <Input_file2> <Output_file> <G>
```

* Input_file1 : GENCODE annotation file
* Input_file2 : annotation file of target gene regions
* Ouput_file : "isoverlap.txt"
* G : the number of target gene regions

#### Example
```
ruby extract_overlap.rb mm10/gencode.vM9.annotation.gtf ES_PrE/data/mygtf_gene.txt ES_PrE/data/isoverlap.txt 4965
```

#### Format of Input_file1
Input_file1 is the GENCODE gene annotation file and is downloaded from https://www.gencodegenes.org.

#### Format of Output_file ("isoverlap.txt")
Each row of Output_file contains tab-separated values (0 or 1) of length L, where L is the number of bins for a gene.
The value 1 indicates that the bin is overlapped with different genes.


## Make "isunmappable.txt" file
Make "isunmappable.txt" file to check whether a mappability of a bin in a gene region is low.

#### Usage
```
rm <Output_file>
Rscript extract_mappability.R <Input_file1> <Input_file2> <Output_file> <a>
```

* Input_file1 : Mappability file
* Input_file2 : Gene regions annotation file
* Output_file : "isunmappable.txt"
* a : threshold to determine isunmappable

#### Example
```
rm ES_PrE/data/isunmappable.txt
Rscript extract_mappability.R mm10/k24.umap.bw ES_PrE/data/mygtf_gene.txt ES_PrE/data/isunmappable.txt 0.5
```

#### Format of Input_file1
The bigwig file of mappability of 24-bp and is downloaded from https://bismap.hoffmanlab.org.

#### Format of Output_file ("isunmappable.txt")
Each row of Output_file contains tab-separated values (0 or 1) of length L, where L is the number of bins for a gene.
The value 1 indicates that the minimum of the mappability of the bin is less than "a".



# <a name="nmf"></a> NMF
## NMF for read count matrix in parallel
Run NMF for read count matrix of each gene regions.

#### Usage
```
Rscript NMF_for_countdata.R <Input_file1> <Input_file2> <Data_dir> <Output_dir> <idx1> <idx2> <C> <K> <seed> <method>
```

* Input_file1 : isoverlap.txt
* Input_file2 : isunmappable.txt
* Data_dir : directory of read count data
* Output_dir : output directory
* idx1 : start index
* idx2 : end index
* C : the number of cells
* K : the factorization rank of NMF
* seed : numerical seed
* method : NMF algorithm, such as "lee", "snmf/l", and "snmf/r" (see https://cran.r-project.org/web/packages/NMF/index.html)

#### Example
```
Rscript NMF_for_countdata.R ES_PrE/data/isoverlap.txt ES_PrE/data/isunmappable.txt ES_PrE/count_data ES_PrE/out1 1 2000 185 5 123456 lee
```

#### Data_dir and idx1,2
The script computes NMF for Data_dir/data_idx1.txt, Data_dir/data_(idx1+1).txt, ..., Data_dir/data_idx2.txt

#### Output_dir
The script output Output_dir/NMF_K_coef_idx1_idx2.txt, which is the C x ((idx2-idx1+1)\*K) matrix.
The coefficient matrix of i-th gene is saved column from (i-idx1)\*K+1 to (i-idx1)\*K+K of NMF_K_coef_idx1_idx2.txt.

## Combine NMF results
Combine Output_dir/NMF_K_coef_a_b.txt, Output_dir/NMF_K_coef_c_d.txt, .....

#### Usage
```
Rscript merge_NMF_coef.R <Input_dir> <C> <G> <K> <idxs1> <idxs2>
```

* Input_dir : the directory of Output_dir in NMF_for_countdata.R
* C : the number of cells
* G : the number of genes
* K : the factorization rank of NMF
* Idxs1 : list of index of idx1 in NMF_for_countdata.R (separated with ",")
* Idxs2 : list of index of idx2 in NMF_for_countdata.R (separated with ",")

#### Example
```
Rscript merge_NMF_coef.R ES_PrE/out1 185 4965 5 1,2001 2000,4965
```

#### Output
The script output Input_dir/NMF_K_coef_all.txt, which is the C x (G\*K) matrix.


# <a name="ttest"></a> t-test
## t-test for NMF results
#### Usage
```
Rscript ttest_NMF.R <Input_file1> <Input_file2> <Output_file> <G> <K>
```

* Input_file1 : cell labels
* Input_file2 : coefficients of NMF
* Output_file : result of t-test
* G : the number of genes
* K : the factorization rank of NMF

#### Example
```
Rscript ttest_NMF.R ES_PrE/data/cell_label.txt ES_PrE/out1/NMF_5_coef_all.txt ES_PrE/out1/ttest_result_NMF_5.txt 4965 5
```

#### Format of Input_file1
The Input_file1 is the C x 2 Cell label table.
The first column corresponds to the sample names, and the second column corresponds to the cluster labels (1 or 2).

#### Example of Input_file1
```
PDIS9061.001	1
PDIS9061.002	1
PDIS9061.003	2
```

#### Format of Input_file2
Input_file2 is the results of NMF and correspond to Input_dir/NMF_K_coef_all.txt.

#### Format of Output_file
Output_file is the G x (K\*2) matrix.
Each row represents the result of t-test for coeffient of NMF.
The 1st to Kth column represents the t-statisics for each coeffient elements, and the (K+1)th to 2\*K th column represent the -log10(p-value) of corresponding t-statistics.


## t-test for TPM matrix
##### Usage
```
Rscript ttest_TPM.R <Input_file1> <Input_file2> <Input_file3> <Input_file4> <Input_file5> <Output_file>
```

* Input_file1 : TPM matrix including header and rowname
* Input_file2 : list of transcript id
* Input_file3 : list of cell labels
* Input_file4 : gene regions annotation file
* Input_file5 : correspondence table between transcript id and gene id
* Output_file : result of t-test

##### Example
```
Rscript ttest_TPM.R ES_PrE/data/TPM_ES_PrE.txt ES_PrE/data/TPM_transcriptid.txt ES_PrE/data/cell_label.txt ES_PrE/data/mygtf_gene.txt ES_PrE/data/transcriptid_to_geneid.txt ES_PrE/out/ttest_result_TPM.txt
```

#### Format of Input_file1
The Input_file1 is the T x C TPM matrix including header and row name, where T is the number of transcripts and C is the number of cells.
The first row corresponds to the sample name, and the first column corresponds to the transcript id.

#### Format of Input_file2
The Input_file2 is the list of transcript id of length T.

#### Format of Input_file5
The Input_file5 is the correspondence table between transcript id and gene id.
The first column corresponds to the transcript id, and the second column corresponds to the gene id.
Each row represents the correspondence between a transcript id and a gene id.

#### Format of Output_file
The Output_file represents the result of t-test for each gene region.
The first and second column represents T+ and T-, respectively.
The third and fourth column represents the -log10(p-value) of T+ and T-, respectively.
The fifth and sixth column represents the transcript ids correspond to T+ and T-, respectively.


## t-test for Mean
#### Usage
```
Rscript ttest_Mean.R <Input_file1> <Input_file2> <Input_file3> <Input_dir> <Output_file> <G>
```

* Input_file1 : list of cell labels
* Input_file2 : isoverlap.txt
* Input_file3 : isunmappable.txt
* Input_dir : directory of read count data
* Output_file : result of t-test
* G : the number of genes

#### Example
```
Rscript ttest_Mean.R ES_PrE/data/cell_label.txt ES_PrE/data/isoverlap.txt ES_PrE/data/isunmappable.txt ES_PrE/count_data ES_PrE/out/ttest_result_Mean.txt 4965
```

#### Format of Output_file
The Output_file represents the result of t-test of the mean values of mapped count for each gene region.
The first column represents T and the second column represents -log10(p-value), respectively.


# <a name="deltaT"></a> Calculate Delta(Tnmf-Ttpm)

#### Usage
```
Rscript calc_DeltaT_NMF_TPM.R <Input_dir1> <Input_dir2> <Output_file> <G> <a>
```

* Input_dir1 : directory for t-test results of NMF
* Input_dir2 : directory for t-test result of TPM
* Output_file : output file
* G : the number of genes
* a : threshold to ignore DE when Ttpm is sufficiently large

#### Example
```
Rscript calc_DeltaT_NMF_TPM.R ES_PrE/out1 ES_PrE/out ES_PrE/out1/DeltaT_NMF_TPM.txt 4965 10
```

#### Format of Output_file
The Output_file represents Delta(Tnmf-Ttpm) values for each gene region.
The first column represents Delta(Tnmf-Ttpm) values and the second column represents the rank of the score in descending order.
Each row represents the result corresponding to the row of gene regions annotation file. 


# <a name="pval"></a> permutation-based test for Delta(Tnmf-Ttpm)
## t-test for NMF results with shuffled cell labels.
#### Usage
```
Rscript ttest_NMF_label_shuffling.R <Input_file1> <Input_file2> <Output_file> <G> <K>
```

* Input_file1 : cell labels
* Input_file2 : coefficients of NMF
* Output_file : result of t-test
* G : the number of genes
* K : the factorization rank of NMF

#### Example
```
Rscript ttest_NMF_label_shuffling.R ES_PrE/data/cell_label.txt ES_PrE/out1/NMF_2_coef_all.txt ES_PrE/out_shuffle1/ttest_result_NMF_2.txt 4965 2
Rscript ttest_NMF_label_shuffling.R ES_PrE/data/cell_label.txt ES_PrE/out1/NMF_5_coef_all.txt ES_PrE/out_shuffle1/ttest_result_NMF_5.txt 4965 5
Rscript ttest_NMF_label_shuffling.R ES_PrE/data/cell_label.txt ES_PrE/out1/NMF_10_coef_all.txt ES_PrE/out_shuffle1/ttest_result_NMF_10.txt 4965 10
```

## Calculate Delta(Tnmf-Ttpm) for shuffled data
#### Example
```
Rscript calc_DeltaT_NMF_TPM.R ES_PrE/out_shuffle1 ES_PrE/out ES_PrE/out_shuffle1/DeltaT_NMF_TPM.txt 4965 10
```

## Permutation-based test
#### Usage
```
Rscript permutation_based_test.R <Input_file1> <Input_file2> <Output_file>
```

* Input_file1 : Delta(Tnmf-Ttpm) for original data
* Input_file2 : Delta(Tnmf-Ttpm) for shuffled data
* Output_file : result of permutation-based test

#### Example
```
Rscript permutation_based_test.R ES_PrE/out1/DeltaT_NMF_TPM.txt ES_PrE/out_shuffle1/DeltaT_NMF_TPM.txt ES_PrE/out_shuffle1/result_permutation_test.txt
```

#### Format of Output_file
The Output_file represents Delta(Tnmf-Ttpm) values and -log10(p-value) for each gene region.
The first column represents Delta(Tnmf-Ttpm) values and the second column represents -log10(p-value).
Each row represents the result corresponding to the row of gene regions annotation file. 
