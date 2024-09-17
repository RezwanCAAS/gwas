
#after filtering the variant calling data under different parameters
MAF,missing data, depth etc.
###################################
#convert vcf file into plink format
#
###################################
plink --vcf accessions.vcf --make-bed --out accessions

#this generates three files
accessions.bed
accessions.bim
accessions.fam

#####################
#perform PCA to reveal the diversity among accessions
#
######################
plink --bfile accessions --pca 10 --out pca_result

#This will generate a file pca_result.eigenvec

#######################
#make a kinship matrix calculations
#
#######################
gemma -bfile accessions -gk 1 -o kinship_matrix

#####################
#run GWAS with GEMMA
#
#####################

gemma -bfile accessions -k kinship_matrix -c pca_result.eigenvec -lmm 4 -p phenotype_trait.csv -o gwas_results.txt

#visualize the GWAS results

library(qqman)
gwasResults <- read.table("gwas_results.txt", header=TRUE)
manhattan(gwasResults, chr="CHR", bp="pos", snp="SNP", p="significance", 
    col=c("blue", "red"),suggestiveline=-log10(1e-5), genomewideline=-log10(5e-8))
