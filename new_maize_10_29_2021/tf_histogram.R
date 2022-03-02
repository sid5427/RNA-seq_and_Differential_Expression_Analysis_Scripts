# ####
# TF - histogram
# 
# pattern from 1 fold to 2 fold and distribution
# 
# (WWvsD, WWvsMS, DvsMS) X 
# 5 sections X 
# genotypes
# 
# 3 X 5 X 2
# 
# 30 plots
# ####

##B73_D_S1_vs_B73_WW_S1_ALL

library(ggplot2)
library(ggfortify)
library(dplyr)

tf_list = read.table("TF_V3_V4.txt",header=T,sep = "\t")
master_df_fc<-data.frame(FC=c(2.0,1.9,1.8,1.7,1.6,1.5,1.4,1.3,1.2,1.1))

dist_col_names<-c("FC")

file_names<-list.files(path = "./EdgeR_ALL")

for (names in file_names){
  print (names)
  name_of_obs<-substr(names,1,nchar(names)-8)
  dist_col_names<-append(dist_col_names,name_of_obs)
  # TF_genes<-inner_join(tf_list,
  #                      read.table(paste("EdgeR_ALL/",names,sep = ""),
  #                                 header=T,row.names=1,sep = ","),
  #                      by = c("V4_id" = "genes"))
  TF_genes<-read.table(paste("EdgeR_ALL/",names,sep = ""),
             header=T,row.names=1,sep = ",")
  num_list<- c()
  
  #2 FC
  num_list<-append(num_list,(nrow(TF_genes[(TF_genes$FDR < 0.05) & 
                                             ((TF_genes$logFC > 1)|
                                                (TF_genes$logFC < -1)),])))
  
  #1.9 FC
  num_list<-append(num_list,(nrow(TF_genes[(TF_genes$FDR < 0.05) & 
                                             ((TF_genes$logFC > 0.92)|
                                                (TF_genes$logFC < -0.92)),])))
  #1.8 FC
  num_list<-append(num_list,(nrow(TF_genes[(TF_genes$FDR < 0.05) & 
                                             ((TF_genes$logFC > 0.85)|
                                                (TF_genes$logFC < -0.85)),])))
  #1.7 FC
  num_list<-append(num_list,(nrow(TF_genes[(TF_genes$FDR < 0.05) & 
                                             ((TF_genes$logFC > 0.76)|
                                                (TF_genes$logFC < -0.76)),])))
  
  #1.6 FC
  num_list<-append(num_list,(nrow(TF_genes[(TF_genes$FDR < 0.05) & 
                                             ((TF_genes$logFC > 0.67)|
                                                (TF_genes$logFC < -0.67)),])))
  
  #1.5 FC
  num_list<-append(num_list,(nrow(TF_genes[(TF_genes$FDR < 0.05) & 
                                             ((TF_genes$logFC > 0.58)|
                                                (TF_genes$logFC < -0.58)),])))
  
  #1.4 FC
  num_list<-append(num_list,(nrow(TF_genes[(TF_genes$FDR < 0.05) & 
                                             ((TF_genes$logFC > 0.48)|
                                                (TF_genes$logFC < -0.48)),])))
  
  #1.3 FC
  num_list<-append(num_list,(nrow(TF_genes[(TF_genes$FDR < 0.05) & 
                                             ((TF_genes$logFC > 0.38)|
                                                (TF_genes$logFC < -0.38)),])))
  
  #1.2 FC
  num_list<-append(num_list,(nrow(TF_genes[(TF_genes$FDR < 0.05) & 
                                             ((TF_genes$logFC > 0.26)|
                                                (TF_genes$logFC < -0.26)),])))
  
  #1.1 FC
  num_list<-append(num_list,(nrow(TF_genes[(TF_genes$FDR < 0.05) & 
                                             ((TF_genes$logFC > 0.13)|
                                                (TF_genes$logFC < -0.13)),])))
  
  master_df_fc<-cbind(master_df_fc,num_list)
  
}

colnames(master_df_fc)<-dist_col_names

write.csv(master_df_fc,file = "master_DE_FC_01_21_2022.csv",col.names = T)

