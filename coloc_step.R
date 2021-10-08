library(data.table)
library(tidyverse)
library(coloc)
library(sys)
args = commandArgs(trailingOnly=TRUE)
sample=args[1]
cases=as.numeric(args[2])
controls=as.numeric(args[3])
sum_stats<-fread(paste0("/shared/shared/cdwillis/cancer_GWAS/data/sum_stats/formatted/37to38translate/",sample,"/formatted/",sample,"_hg38_positions_all_chr.txt"),stringsAsFactors = FALSE,header=T,fill=TRUE)
cis_qtl<-fread("/shared/shared/cdwillis/hiv_coloc/hg38_final_cis_qtl_output_for_co_loc.txt",stringsAsFactors=FALSE,header=T,fill=TRUE)
maf<-fread("/shared/shared/cdwillis/cancer_GWAS/coloc_analysis/co_loc_analysis_snp_MAF.txt",stringsAsFactors = FALSE,header=T)
gene_exp<-fread("/shared/shared/cdwillis/cancer_GWAS/data/HIV_acq_matrix_eqtl_input_count_normalized_gene_counts_filtered_dge_de_genes.txt",stringsAsFactors = FALSE,header=T)
gene_list<-unique(cis_qtl$gene)
source("/shared/shared/cdwillis/cancer_GWAS/coloc_analysis/running_coloc.R")
run_coloc<-lapply(gene_list,running_coloc.abf,
                  qtl_data=cis_qtl,sum_stat_data=sum_stats,maf_data=maf, gene_exp_data=gene_exp, 
                  maf_position_col="GRCh38_position", maf_chromosome_col="Chromosome", maf_val_col="MAF",
                  qtl_beta_col="beta", qtl_se_col="SE", qtl_allele_1_col="Allele_1",qtl_allele_2_col="Allele_2",
                  qtl_position_col="GRCH38_position", qtl_chromosome_col="Chromosome", qtl_gene_col="gene",
                  sum_stat_beta_col="Beta", sum_stat_se_col="SE",
                  sum_stat_position_col="GRCh38_position", sum_stat_chromosome_col="chromosome",
                  sum_stat_allele_1_col="effect_allele",sum_stat_allele_2_col="other_allele",
                  sum_stat_phenotype_type="cc",sum_stat_phenotype_prop=(cases)/(cases+controls),sum_stat_sample_size=(cases+controls))

coloc_results<-do.call("rbind",run_coloc)
coloc_results_ordered<-coloc_results%>%arrange(desc(PP.H4.abf))
head(coloc_results_ordered)
nrow(coloc_results)
write.table(coloc_results_ordered,file=paste0("/shared/shared/cdwillis/cancer_GWAS/coloc_analysis/results/results_colocalization_",sample,"_summary_stats_with_qtls_with_strand_flips.txt"),
            col.names=T,row.names=T,sep="\t",quote=F)


