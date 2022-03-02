

library(dplyr)
library(tidyverse)

newflower_fpkm_table = read.table("fpkm_all_leaf.txt",header = TRUE,sep = "\t")

newflower_C_vs_D = read.table("leaf/C_vs_D_gene_exp.diff",header = TRUE,sep = "\t")
newflower_C_vs_H = read.table("leaf/C_vs_H_gene_exp.diff",header = TRUE,sep = "\t")
newflower_C_vs_HD = read.table("leaf/C_vs_HD_gene_exp.diff",header = TRUE,sep = "\t")

combined_table <- full_join(newflower_fpkm_table,newflower_C_vs_D,by = "gene_id")
combined_table <- full_join(combined_table,newflower_C_vs_H,by = "gene_id")
combined_table <- full_join(combined_table,newflower_C_vs_HD,by = "gene_id")


#write.table(annotated_counts,"combined_new_leaf_results.txt",sep = "\t",row.names = TRUE,col.names = TRUE)
#comb_table = read.csv("combined_new_leaf_results.txt",header = TRUE,sep = "\t")

anno_table = read.csv("Gmax_275_Wm82.a2.v1.annotation_info.txt",header = TRUE,sep = "\t")

#colnames(anno_table)
#reduced_anno_table<-select(anno_table, c(Gene_ID,Best.hit.arabi.name,arabi.symbol,arabi.defline))
#colnames(reduced_anno_table)

#distinct_anno_table<-distinct(reduced_anno_table)
distinct_anno_table_geneID<-distinct(anno_table,gene_id, .keep_all= TRUE)

length(unique(anno_table$gene_id))

annotated_counts_exp <- left_join(combined_table,distinct_anno_table_geneID,by = "gene_id")

TAF_all_table = read.csv("Gma_TF_list/Gmax_TF_list_mod.txt",header = TRUE,sep = "\t")
annotated_counts_exp_2 <- left_join(annotated_counts_exp,TAF_all_table,by = "gene_id")
write.table(annotated_counts_exp_2,"annotated_TF_new_leaf_results_fpkm.txt",sep = "\t",row.names = TRUE,col.names = TRUE)



# write.table(annotated_counts_exp,"annotated_leaf_results_fpkm.txt",sep = "\t",row.names = TRUE,col.names = TRUE)
# 
# TF_leaf_HD_uniq = read.table("TF_leaf_HD_uniq.txt",header = TRUE,sep = "\t")
# TF_annotated_exp_leaf <- left_join(TF_leaf_HD_uniq,annotated_counts_exp,by = "gene_id")
# write.table(TF_annotated_exp_leaf,"TF_annotated_leaf.txt",sep = "\t",row.names = TRUE,col.names = TRUE)
