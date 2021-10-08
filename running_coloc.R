running_coloc.abf<-function(selected_gene, qtl_data,sum_stat_data,maf_data, gene_exp_data, 
                            maf_position_col, maf_chromosome_col, maf_val_col,
                            qtl_beta_col, qtl_se_col,
                            qtl_position_col, qtl_chromosome_col, qtl_gene_col, qtl_allele_1_col, qtl_allele_2_col,
                            sum_stat_beta_col, sum_stat_se_col,
                            sum_stat_position_col, sum_stat_chromosome_col, sum_stat_allele_1_col, sum_stat_allele_2_col,
                            sum_stat_phenotype_type,sum_stat_phenotype_sd=NULL,sum_stat_phenotype_prop,sum_stat_sample_size_col=NULL,
                            sum_stat_sample_size){
  
  #The coloc.abf function only uses type="cc" or type="quant". Other entries to sum_stat_phenotype_type will produce an error.
  
  if(!(sum_stat_phenotype_type %in% c("cc","quant"))){ stop("Error: sum_stat_phenotype_type should only be 'cc' or 'quant'") }
  
  
  
  #Removing the rows that did not match successfully to a build 38 position. Merging may or may not take care of that already.
  
  sum_stat_data<-sum_stat_data[which(!is.na(sum_stat_data[[sum_stat_position_col]])),]
  
  qtl_data<-qtl_data[which(!is.na(qtl_data[[qtl_position_col]])),]
  
  
  
  #We filter the qtl_data to just one gene at a time because for the qtl_data, each gene can be thought of as a trait.
  #It seems like if we include all the genes at once, 
  #we would be trying to find evidence of co-localization between the phenotype of interest and all the genes at once.
  #That seems much less likely and ultimately not what we want. 
  #However, you can input multiple genes too, to filter the data down to more than one gene if needed.
  
  #Not sure why this has to be written like this and not using the x$col_name syntax because that doesn't work.
  
  qtl_filter<-qtl_data[qtl_data[[qtl_gene_col]] %in% selected_gene,]%>%rename(qtl_beta=all_of(qtl_beta_col),qtl_se=all_of(qtl_se_col),
                                                                              qtl_A1=all_of(qtl_allele_1_col),qtl_A2=all_of(qtl_allele_2_col))
  
  #We also have to take into account the possiblity of strand flips. Bluntly, we will find the complement of A1 and A2 in the sum stats dataset,
  #and we will do the same filtering as below once the data is merged (matching qtl alleles with the strand-flipped sum stat alleles)
  
  sum_stat_data<-sum_stat_data[sum_stat_data[[sum_stat_position_col]] %in% qtl_filter[[qtl_position_col]] & 
                                 sum_stat_data[[sum_stat_chromosome_col]] %in% qtl_filter[[qtl_chromosome_col]],]%>%
    rename(sum_stat_A1=all_of(sum_stat_allele_1_col),sum_stat_A2=all_of(sum_stat_allele_2_col))%>%
    mutate(strand_flip_sum_stat_A1=case_when(sum_stat_A1=="G" ~ "C",
                                             sum_stat_A1=="C" ~ "G",
                                             sum_stat_A1=="A" ~ "T",
                                             sum_stat_A1=="T" ~ "A"),
           strand_flip_sum_stat_A2=case_when(sum_stat_A2=="G" ~ "C",
                                             sum_stat_A2=="C" ~ "G",
                                             sum_stat_A2=="A" ~ "T",
                                             sum_stat_A2=="T" ~ "A"))
  
  maf_data<-maf_data%>%rename(MAF=all_of(maf_val_col))
  
  if (is.null(sum_stat_sample_size_col)) {
    
    #print("Single sample size for summary stats data")
    
    if(is.numeric(sum_stat_sample_size)){
      
      sum_stat_data<-sum_stat_data%>%rename(sum_stats_beta=all_of(sum_stat_beta_col),sum_stats_se=all_of(sum_stat_se_col))
      
      merged_data<-inner_join(qtl_filter,sum_stat_data, 
                              by=setNames(c(sum_stat_position_col,sum_stat_chromosome_col),c(qtl_position_col,qtl_chromosome_col)),
                              suffix=c(".qtl",".sum_stats"))
      
      #Make sure we only select SNPs to include where the SNPs from the qtl and the sum stats data have the same reference 
      #and alternate alleles. Also, match using the strand flipped summary stats alleles as well.
      
      
      merged_data<-merged_data[(merged_data[["qtl_A1"]]==merged_data[["sum_stat_A1"]] & merged_data[["qtl_A2"]]==merged_data[["sum_stat_A2"]]) 
                               | (merged_data[["qtl_A1"]]==merged_data[["sum_stat_A2"]] & merged_data[["qtl_A2"]]==merged_data[["sum_stat_A1"]])
                               | (merged_data[["qtl_A1"]]==merged_data[["strand_flip_sum_stat_A1"]] & merged_data[["qtl_A2"]]==merged_data[["strand_flip_sum_stat_A2"]])
                               | (merged_data[["qtl_A1"]]==merged_data[["strand_flip_sum_stat_A2"]] & merged_data[["qtl_A2"]]==merged_data[["strand_flip_sum_stat_A1"]]),]
      
      merged_data$sum_stats_n<-sum_stat_sample_size
      
      
    } else {
      
      stop("Error: You are specifying a single sample size for the summary stats data. Input a numeric value for 'sum_stat_sample_size'.")
      
    }
    
    
    
  } else if (is.null(sum_stat_sample_size_col)==0) {
    
    #print("Sample size column for summary stats data")
    
    sum_stat_data<-sum_stat_data%>%rename(sum_stats_beta=all_of(sum_stat_beta_col),sum_stats_se=all_of(sum_stat_se_col),
                                          sum_stats_n=all_of(sum_stat_sample_size_col))
    
    merged_data<-inner_join(qtl_filter,sum_stat_data, 
                            by=setNames(c(sum_stat_position_col,sum_stat_chromosome_col),c(qtl_position_col,qtl_chromosome_col)),
                            suffix=c(".qtl",".sum_stats"))
    
    #Make sure we only select SNPs to include where the SNPs from the qtl and the sum stats data have the same reference 
    #and alternate alleles. Also, match using the strand flipped summary stats alleles as well.
    
    merged_data<-merged_data[(merged_data[["qtl_A1"]]==merged_data[["sum_stat_A1"]] & merged_data[["qtl_A2"]]==merged_data[["sum_stat_A2"]]) 
                             | (merged_data[["qtl_A1"]]==merged_data[["sum_stat_A2"]] & merged_data[["qtl_A2"]]==merged_data[["sum_stat_A1"]])
                             | (merged_data[["qtl_A1"]]==merged_data[["strand_flip_sum_stat_A1"]] & merged_data[["qtl_A2"]]==merged_data[["strand_flip_sum_stat_A2"]])
                             | (merged_data[["qtl_A1"]]==merged_data[["strand_flip_sum_stat_A2"]] & merged_data[["qtl_A2"]]==merged_data[["strand_flip_sum_stat_A1"]]),]
    
  } 
  
  
  if(nrow(merged_data)==0){ 
    
    print(paste0("Warning: Inner-joining summary stat data and cis-qtl data for gene ", selected_gene, " resulted in 0 rows."))
    
    results<-as.data.frame(cbind(0,NA,NA,NA,NA,NA))
    
  } else if(nrow(merged_data)>0) {
    
    
    complete_merged_data<-inner_join(merged_data,maf_data,
                                     by=setNames(c(maf_position_col,maf_chromosome_col),c(qtl_position_col,qtl_chromosome_col)),
                                     suffix=c(".snp_data",".maf"))
    
    length_MAF<-length(complete_merged_data[["MAF"]][!is.na(complete_merged_data[["MAF"]])])
    
    length_sum_stats_n<-length(complete_merged_data[["sum_stats_n"]][!is.na(complete_merged_data[["sum_stats_n"]])])
    
    if (length_MAF!=length_sum_stats_n){
      
      print(paste0("Warning: For ", selected_gene, " the number of snp MAF values and summary stats sample size values are not equal."))
      print(paste0("The number of snp MAF values is ", length_MAF))
      print(paste0("The number of summary stats sample size values is ", length_sum_stats_n))
      
    }
    
    
    #filtered_maf<-maf_data[maf_data[[maf_position_col]] %in% merged_data[[qtl_position_col]]]
    
    # filtered_maf<-filtered_maf[order(match(filtered_maf[[maf_position_col]],merged_data[[qtl_position_col]],)),]
    
    # filtered_maf2<-filtered_maf
    
    gene_matrix<-gene_exp_data[,-1]
    
    rownames(gene_matrix)<-unlist(gene_exp_data[,1])
    
    filtered_gene_matrix<-gene_matrix[rownames(gene_matrix)==selected_gene]
    
    sd_gene_exp<-sd(filtered_gene_matrix)
    
    if (sum_stat_phenotype_type=="cc"){
      
      print(selected_gene)
      
      my.res <- coloc.abf(dataset1=list(beta=complete_merged_data$qtl_beta, varbeta=(complete_merged_data$qtl_se)**2, 
                                        N=ncol(filtered_gene_matrix),sdY=sd_gene_exp,type="quant"),
                          dataset2=list(beta=complete_merged_data$sum_stats_beta, 
                                        N=complete_merged_data$sum_stats_n,type=sum_stat_phenotype_type,s=sum_stat_phenotype_prop,pvalues=complete_merged_data$p_value),
                          MAF=complete_merged_data$MAF)
      
    } else if (sum_stat_phenotype_type=="quant"){
      
      if (is.null(sum_stat_phenotype_sd)){
        
        print(selected_gene)
        
        
        my.res <- coloc.abf(dataset1=list(beta=complete_merged_data$qtl_beta, varbeta=(complete_merged_data$qtl_se)**2, 
                                          N=ncol(filtered_gene_matrix),sdY=sd_gene_exp,type="quant"),
                            dataset2=list(beta=complete_merged_data$sum_stats_beta, varbeta=(complete_merged_data$sum_stats_se)**2, 
                                          N=complete_merged_data$sum_stats_n,type=sum_stat_phenotype_type),
                            MAF=complete_merged_data$MAF)
        
      } else if ( (is.null(sum_stat_phenotype_sd)==0) & (sum_stat_phenotype_sd>=0) ) {
        
        print(selected_gene)
        
        
        my.res <- coloc.abf(dataset1=list(beta=complete_merged_data$qtl_beta, varbeta=(complete_merged_data$qtl_se)**2, 
                                          N=ncol(filtered_gene_matrix),sdY=sd_gene_exp,type="quant"),
                            dataset2=list(beta=complete_merged_data$sum_stats_beta, varbeta=(complete_merged_data$sum_stats_se)**2, 
                                          N=complete_merged_data$sum_stats_n,type=sum_stat_phenotype_type,sdY=sum_stat_phenotype_sd),
                            MAF=complete_merged_data$MAF)
        
      } else {
        
        stop("Error: Check your entry for the summary statistics SD. Is it a number greater than or equal to 0?")
      }
      
      
    }
    
    results<-as.data.frame(cbind(my.res[[1]][1],my.res[[1]][2],my.res[[1]][3],my.res[[1]][4],my.res[[1]][5],my.res[[1]][6]))
    
    
    
    
  }
  
  rownames(results)<-selected_gene
  colnames(results)<-c("N_SNPS", "PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf")
  #test<-list(maf1=filtered_maf,maf2=filtered_maf2,maf_col=filtered_maf$MAF,n=sum_stats_n)
  return(results)
  
  
}
