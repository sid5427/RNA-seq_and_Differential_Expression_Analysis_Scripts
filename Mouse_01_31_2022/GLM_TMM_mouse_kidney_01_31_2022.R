library(edgeR)
library(limma)
library(ggplot2)
#packageVersion('ggplot2')
library(ggfortify)
#install.packages("pheatmap")
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(magrittr)
library(tibble)
#packageVersion('dplyr') #old was 0.8.3
#library(dendextend)

#setwd("C:/Users/Sid/OneDrive/rna-seq documents/NSF-drought/FR697_EXTENDED_ASSEMBLY/SuperDuper/final_rep_set/rna-seq_files/DE_calc")

####load workspace image####

load('All_DE_02_02_2022.RData')

######

##load data matrix/ table
#ReadCounts = read.table("combined_rna_seq_counts.txt",header=T,row.names=1)
ReadCounts = read.table("rna_seq_counts_MK_01_30_2022_main.txt",header=T,row.names=1)


##remove WD.PTKO.2, WD.PTKO.6, CD.PTKO.3, CD.PTWT.3, CD.PTWT.4, 

#full_sample_groups_list
sample_groups_list = read.table("sample_groups_MK.txt",header=T)

######groups####

sample_groups_list$groups

#batch <- sample_groups_list$batch no batches in this analysis
groups <- sample_groups_list$groups
replicate<-sample_groups_list$replicate
groups

colnames(ReadCounts)
head(ReadCounts)
tail(ReadCounts)


######Generate normalization factors from TMM method####

#Put the counts and other information into a DGEList object:
##original
#d <- DGEList(counts=ReadCounts, group=groups)
##modified for GLM
row.names(ReadCounts)
d <- DGEList(counts=ReadCounts, group=groups, genes = row.names(ReadCounts)) ##maybe not use this?
dim(d)
table(groups)
#table(batch)

###extra filtering - this can be changed
keep <- rowSums(cpm(d)>1) >= 3
#keep <- rowSums(cpm(d)>1) >= 6
dim(d[keep,])
d_filt <- d[keep,]
dim(d_filt)
#boxplot(d_filt$counts, col="gray", las=3)
#number filtered
dim(d)[1]-dim(d_filt)[1] ##14484
#%filtered out
(dim(d)[1]-dim(d_filt)[1])/dim(d)[1] ##0.3370175
#library size update ----- IMP STEP
d_filt$samples$lib.size <- colSums(d_filt$counts)

##Get normalized data from TMM via CPM method
#TMM <- calcNormFactors(d, method="TMM")
dim(d_filt)
TMM <- calcNormFactors(d_filt, method="TMM") ##if using extra filtering
normalized_data <- cpm(TMM)
d_filt <- cpm(TMM)
d_filt <- DGEList(counts=d_filt, group=groups, genes = row.names(d_filt))

head(normalized_data)
dim(normalized_data)
dim(ReadCounts)

colnames(normalized_data)
colnames(ReadCounts)

zero_rows<-subset(ReadCounts, !(row.names(ReadCounts) %in% row.names(normalized_data)))

dim(zero_rows)
norm_zero_combined_data<-rbind(normalized_data,zero_rows)
colnames(norm_zero_combined_data)

####write to CSV files - comment out after use to prevent overwriting.
# normalized_data_tb <- tibble::rownames_to_column(as.data.frame(normalized_data), "genes")
# normalized_data_tb<-left_join(normalized_data_tb,mouse_mart_export,by="genes")
# write.csv(normalized_data_tb,file = "norm_filtered_MK_only_01_31_2022.csv")

norm_zero_combined_data_tb <- tibble::rownames_to_column(as.data.frame(norm_zero_combined_data), "genes")
norm_zero_combined_data_tb<-left_join(norm_zero_combined_data_tb,mouse_mart_export,by="genes")
#write.csv(norm_zero_combined_data_tb,file = "norm_zero_ALL_MK_01_31_2022.csv")

####PLOTS ####
normalized_data_t <- t(normalized_data) #### ALL

normalized_data_t[1:10,1:10]
dim(normalized_data_t)

####dendogram plot####
d <- dist(normalized_data_t)
hc <- hclust(d)
dhc <- as.dendrogram(hc)
#specific_leaf <- dhc[[1]][[1]][[1]]
#attributes(specific_leaf)
plot(hc)

normalized_data_t <- data.frame(normalized_data_t)
fpkm_table_data_groups <- normalized_data_t
dim(fpkm_table_data_groups)

FULL_table_data_groups<-as.data.frame(t(ReadCounts))
dim(FULL_table_data_groups)
groups

fpkm_table_data_groups['groups'] <- groups
FULL_table_data_groups['groups'] <- groups
dim(fpkm_table_data_groups)
fpkm_table_data_groups['sample_names'] <- sample_groups_list$Sample
FULL_table_data_groups['sample_names'] <- sample_groups_list$Sample
dim(fpkm_table_data_groups)

fpkm_table_data_groups['replicate'] <- replicate

##raw counts PCA
dim(FULL_table_data_groups)
autoplot(prcomp(FULL_table_data_groups[1:25706]), 
         data = FULL_table_data_groups, colour = 'groups' , 
         label = TRUE, label.size = 3,shape = FALSE)

##normalized counts PCA
fpkm_table_data_groups$groups
autoplot(prcomp(fpkm_table_data_groups[1:12358]), data = fpkm_table_data_groups, colour = 'groups' , label = TRUE, label.size = 2,shape = FALSE)
#autoplot(prcomp(fpkm_table_data_groups[1:12895]), data = fpkm_table_data_groups, colour = 'sample_names' , label = TRUE, label.size = 3)

####other pca####

##remember to update the row(column) no.s for the PCAs

##filter down to groups which contain 18-27m
fpkm_table_data_groups_18_27m<- fpkm_table_data_groups %>% filter(grepl('18_27m', groups))
autoplot(prcomp(fpkm_table_data_groups_18_27m[1:12530]), 
         data = fpkm_table_data_groups_18_27m, 
         colour = 'groups' , 
         label = TRUE, 
         label.size = 3,shape = FALSE)

fpkm_table_data_groups_CD<- fpkm_table_data_groups %>% filter(grepl('CD', groups))
autoplot(prcomp(fpkm_table_data_groups_CD[1:12895]), 
         data = fpkm_table_data_groups_CD, 
         colour = 'groups' , 
         label = TRUE, 
         label.size = 3,shape = FALSE,shape.size=5)

fpkm_table_data_groups_CDWT<- fpkm_table_data_groups %>% filter(grepl('WT', groups)) %>% filter(grepl('CD', groups))
autoplot(prcomp(fpkm_table_data_groups_CDWT[1:12530]), 
         data = fpkm_table_data_groups_CDWT, 
         colour = 'groups' , 
         label = TRUE, 
         label.size = 3,shape = FALSE,shape.size=5)

fpkm_table_data_groups_WD<- fpkm_table_data_groups %>% filter(grepl('WD', groups))
autoplot(prcomp(fpkm_table_data_groups_WD[1:12895]), 
         data = fpkm_table_data_groups_WD, 
         colour = 'groups' , 
         label = TRUE, 
         label.size = 3,shape = FALSE,shape.size=5)

fpkm_table_data_groups_WDWT<- fpkm_table_data_groups %>% filter(grepl('WT', groups)) %>% filter(grepl('WD', groups))
autoplot(prcomp(fpkm_table_data_groups_WDWT[1:12530]), 
         data = fpkm_table_data_groups_WDWT, 
         colour = 'groups' , 
         label = TRUE, 
         label.size = 3,shape = FALSE,shape.size=5)


fpkm_table_data_groups_KO_only<- fpkm_table_data_groups %>% filter(grepl('KO', groups))
autoplot(prcomp(fpkm_table_data_groups_KO_only[1:12895]), 
         data = fpkm_table_data_groups_KO_only, 
         colour = 'groups' , 
         label = TRUE, 
         label.size = 3,shape = FALSE,shape.size=5)

fpkm_table_data_groups_WTall<- fpkm_table_data_groups %>% filter(grepl('WT', groups))
autoplot(prcomp(fpkm_table_data_groups_WTall[1:12530]), 
         data = fpkm_table_data_groups_WTall, 
         colour = 'groups' , 
         label = TRUE, 
         label.size = 3,shape = FALSE,shape.size=5)

fpkm_table_data_groups_WT_18_27m<- fpkm_table_data_groups %>% filter(grepl('WT_18-27m', groups))
autoplot(prcomp(fpkm_table_data_groups_WT_18_27m[1:12895]), 
         data = fpkm_table_data_groups_WT_18_27m, 
         colour = 'groups' , 
         label = TRUE, 
         label.size = 3,shape = FALSE,shape.size=5)

fpkm_table_data_groups_6m<- fpkm_table_data_groups %>% filter(grepl('WT_6m', groups))
autoplot(prcomp(fpkm_table_data_groups_6m[1:12895]), 
         data = fpkm_table_data_groups_6m, 
         colour = 'groups' , 
         label = TRUE, 
         label.size = 3,shape = FALSE,shape.size=5)

fpkm_table_data_groups_16m<- fpkm_table_data_groups %>% filter(grepl('WT_16m', groups))
autoplot(prcomp(fpkm_table_data_groups_16m[1:12895]), 
         data = fpkm_table_data_groups_16m, 
         colour = 'groups' , 
         label = TRUE, 
         label.size = 3,shape = FALSE,shape.size=5)

fpkm_table_data_groups_KO<- fpkm_table_data_groups %>% filter(grepl('KO', groups))
autoplot(prcomp(fpkm_table_data_groups_KO[1:12895]), 
         data = fpkm_table_data_groups_KO, 
         colour = 'groups' , 
         label = TRUE, 
         label.size = 3,shape = FALSE,shape.size=5)
####3d PCA####
#install.packages("pca3d")
library( pca3d )

pca.data <- prcomp(fpkm_table_data_groups[1:12895])
scores = as.data.frame(pca.data$x)
plot(pca.data, type="lines")
summary(pca.data)

#plot3d(pca.data$scores[,1:3])


pca3d(pca.data,group=sample_groups_list$groups)
all_3d_plot<-pca3d(pca.data,group=sample_groups_list$groups)

all_3d_plot$colors
all_3d_plot[colors] <- lapply(all_3d_plot[colors], colhex2col)

list_of_cols = list()
for (cols in all_3d_plot$colors){
  print(colhex2col(cols))
  }


pca3d(pca.data,group=sample_groups_list$rep_n)
pca3d(pca.data,group=sample_groups_list$treatment,col = sample_groups_list$Color_grp,legend = "topright")
pca3d(pca.data,group=sample_groups_list$treatment,col = sample_groups_list$section_color)
?pca3d

sample_groups_list
#####

# dend <- normalized_data_t %>%  scale %>% 
#   dist %>% hclust %>% as.dendrogram
# dend %>% plot
# 
# dend <- dend %>% set("labels_colors", as.numeric(sample_groups_list$Tip_region), order_value = TRUE) %>%
#   set("labels_cex", 1)
# 
# #par(mar = c(4,1,0,8))
# plot(dend)

####adding group info to dataset -- optional - needed for colored plots atm####
# normalized_data_dt <- data.frame(normalized_data_t)
# normalized_data_groups <- normalized_data_dt
# dim(normalized_data_groups)
# #add groups
# dim(groups)
# normalized_data_groups['groups'] <- groups
# dim(normalized_data_groups)
# normalized_data_groups['years'] <- years
# dim(normalized_data_groups)
# normalized_data_groups['region'] <- region
# dim(normalized_data_groups)
# normalized_data_groups['genotype'] <- genotype
# dim(normalized_data_groups)
# norm_data_BACKUP <- normalized_data_groups
# ##LOAD BACKUP DATA
# #normalized_data_groups <- norm_data_BACKUP
# 
# ##add the row names a starting column for coloring
# normalized_data_groups <- cbind(sample = rownames(normalized_data_groups), normalized_data_groups)
# #normalized_data_groups[,1] <- rownames(normalized_data_groups)
# dim(normalized_data_groups)
# normalized_data_groups[,28078]
# normalized_data_groups[,28079]
# normalized_data_groups[,28080]




#### DGE computation####

summary(TMM$samples$norm.factors)
summary(TMM$samples$lib.size)
TMM$samples

#table(years) ##replace years with batch
table(batch) ##replace years with batch
table(groups) 

##design matrix
#design=model.matrix(~0+groups+years) #<--main ##replace years with batch
#design=model.matrix(~0+groups+batch) #<--main ##replace years with batch
design=model.matrix(~0+groups) #<--main ##replace years with batch
#design=model.matrix(~0+years) #<--main
#design=model.matrix(~0+years+groups) #<--main
design

colnames(d_filt)
rownames(design)=colnames(d_filt)
#rownames(design)=colnames(ReadCounts)
#rownames(design)=og_row_names
#rownames(design)
#write.csv(design,file = "design_groups_only.csv")


#EdgeR DE
dge_final <- estimateGLMCommonDisp(d_filt, design, verbose=TRUE) 
dge_final <- estimateGLMTrendedDisp(dge_final, design,verbose=TRUE)
####disp####
# Disp = 0 , BCV = 3e-04 
# Disp = 0 , BCV = 6e-04 
# Disp = 0 , BCV = 3e-04 
# Disp = 0 , BCV = 3e-04 
# Disp = 0 , BCV = 3e-04 
# Disp = 0 , BCV = 5e-04 
# Disp = 0 , BCV = 3e-04 
# Disp = 0 , BCV = 3e-04 
# Disp = 0 , BCV = 1e-04 
# Disp = 0 , BCV = 4e-04 
# Disp = 0 , BCV = 3e-04 
# Disp = 0 , BCV = 3e-04 
# Disp = 0 , BCV = 4e-04 
# Disp = 0 , BCV = 4e-04 
# Disp = 0 , BCV = 5e-04 
# Disp = 0 , BCV = 0 
# Disp = 0 , BCV = 4e-04 
# Disp = 0 , BCV = 4e-04 
# Disp = 0 , BCV = 3e-04 
# Disp = 0 , BCV = 4e-04 
# Disp = 0 , BCV = 6e-04 
# Disp = 0 , BCV = 3e-04 
# Disp = 0 , BCV = 1e-04 
# Disp = 0.01949 , BCV = 0.1396 
# Disp = 0.00838 , BCV = 0.0915 
# Disp = 0.00707 , BCV = 0.0841 
# Disp = 0.01258 , BCV = 0.1121 
# Disp = 0.01691 , BCV = 0.1301 
# Disp = 0.02747 , BCV = 0.1657 
# Disp = 0.02946 , BCV = 0.1716 
# Disp = 0.03332 , BCV = 0.1825 
# Disp = 0.02967 , BCV = 0.1722 
# Disp = 0.02621 , BCV = 0.1619 
# Disp = 0.03789 , BCV = 0.1947 
# Disp = 0.03609 , BCV = 0.19 
# Disp = 0.04476 , BCV = 0.2116 
# Disp = 0.04812 , BCV = 0.2194 
# Disp = 0.05124 , BCV = 0.2264 
# Disp = 0.03529 , BCV = 0.1879 
# Disp = 0.04039 , BCV = 0.201 
# Disp = 0.04072 , BCV = 0.2018 
# Disp = 0.03477 , BCV = 0.1865 
# Disp = 0.03094 , BCV = 0.1759
######
dge_final <- estimateGLMTagwiseDisp(dge_final, design)

fit <- glmFit(dge_final, design)
#fit$samples
fit

table(groups)
list_groups <- as.data.frame(table(groups))

# groupsB73_Ssa groupsB73_SSb groupsB73_SSc 
#groupsB73_Wwa groupsB73_Wwb groupsB73_Wwc 
#groupsFR697_Ssa groupsFR697_Ssb groupsFR697_SSc 
#groupsFR697_Wwa groupsFR697_Wwb groupsFR697_Wwc 
#groupsLab_FR697_Msa groupsLab_FR697_Msb groupsLab_FR697_Msc 
#groupsLab_FR697_Ssa groupsLab_FR697_Ssb groupsLab_FR697_Ssc 
#groupsLab_FR697_Wwa groupsLab_FR697_Wwb groupsLab_FR697_Wwc

#save.image(file='All_DE_02_02_2022.RData')

####CONTRASTS####

#mouse_mart_export = read.csv("MOUSE_IDS_mart_export.txt",header=T,sep = "\t")
#head(mouse_mart_export)


#### NEW COMBINATION - 4th feb 2022####
# a2. WD_PTWT_18_27m vs CD_PTWT_18_27m

lrt.WD_PTWT_18_27m_vs_CD_PTWT_18_27m= glmLRT(fit,contrast=makeContrasts(groupsWD_PTWT_18_27m-groupsCD_PTWT_18_27m, levels=design))
WD_PTWT_18_27m_vs_CD_PTWT_18_27m<- as.data.frame(topTags(lrt.WD_PTWT_18_27m_vs_CD_PTWT_18_27m, n=Inf))
#add info
WD_PTWT_18_27m_vs_CD_PTWT_18_27m<-left_join(WD_PTWT_18_27m_vs_CD_PTWT_18_27m,mouse_mart_export,by="genes")
#pvalue 0.05 & FC = 2.0
dim(WD_PTWT_18_27m_vs_CD_PTWT_18_27m[(WD_PTWT_18_27m_vs_CD_PTWT_18_27m$PValue < 0.05) & 
                                       ((WD_PTWT_18_27m_vs_CD_PTWT_18_27m$logFC > 1)|
                                          (WD_PTWT_18_27m_vs_CD_PTWT_18_27m$logFC < -1)),])
write.csv(WD_PTWT_18_27m_vs_CD_PTWT_18_27m[(WD_PTWT_18_27m_vs_CD_PTWT_18_27m$PValue < 0.05) & 
                                             ((WD_PTWT_18_27m_vs_CD_PTWT_18_27m$logFC > 1)|
                                                (WD_PTWT_18_27m_vs_CD_PTWT_18_27m$logFC < -1)),],
          file= "Results/EdgeR_PVal0.5_2fold/WD_PTWT_18_27m_vs_CD_PTWT_18_27m_2fold.csv",row.names=F)
#pvalue 0.05 & FC = 1.5
dim(WD_PTWT_18_27m_vs_CD_PTWT_18_27m[(WD_PTWT_18_27m_vs_CD_PTWT_18_27m$PValue < 0.05) & 
                                       ((WD_PTWT_18_27m_vs_CD_PTWT_18_27m$logFC > 0.58)|
                                          (WD_PTWT_18_27m_vs_CD_PTWT_18_27m$logFC < -0.58)),])
write.csv(WD_PTWT_18_27m_vs_CD_PTWT_18_27m[(WD_PTWT_18_27m_vs_CD_PTWT_18_27m$PValue < 0.05) & 
                                             ((WD_PTWT_18_27m_vs_CD_PTWT_18_27m$logFC > 0.58)|
                                                (WD_PTWT_18_27m_vs_CD_PTWT_18_27m$logFC < -0.58)),],
          file= "Results/EdgeR_PVal0.5_2fold/WD_PTWT_18_27m_vs_CD_PTWT_18_27m_1.5fold.csv",row.names=F)


# a4. WD_PTWT_18_27m vs WD_PTKO_18_27m

lrt.WD_PTWT_18_27m_vs_WD_PTKO_18_27m= glmLRT(fit,contrast=makeContrasts(groupsWD_PTWT_18_27m-groupsWD_PTKO_18_27m, levels=design))
WD_PTWT_18_27m_vs_WD_PTKO_18_27m<- as.data.frame(topTags(lrt.WD_PTWT_18_27m_vs_WD_PTKO_18_27m, n=Inf))
#add info
WD_PTWT_18_27m_vs_WD_PTKO_18_27m<-left_join(WD_PTWT_18_27m_vs_WD_PTKO_18_27m,mouse_mart_export,by="genes")
#pvalue 0.05 & FC = 2.0
dim(WD_PTWT_18_27m_vs_WD_PTKO_18_27m[(WD_PTWT_18_27m_vs_WD_PTKO_18_27m$PValue < 0.05) & 
                                       ((WD_PTWT_18_27m_vs_WD_PTKO_18_27m$logFC > 1)|
                                          (WD_PTWT_18_27m_vs_WD_PTKO_18_27m$logFC < -1)),])
write.csv(WD_PTWT_18_27m_vs_WD_PTKO_18_27m[(WD_PTWT_18_27m_vs_WD_PTKO_18_27m$PValue < 0.05) & 
                                             ((WD_PTWT_18_27m_vs_WD_PTKO_18_27m$logFC > 1)|
                                                (WD_PTWT_18_27m_vs_WD_PTKO_18_27m$logFC < -1)),],
          file= "Results/EdgeR_PVal0.5_2fold/WD_PTWT_18_27m_vs_WD_PTKO_18_27m_2fold.csv",row.names=F)
#pvalue 0.05 & FC = 1.5
dim(WD_PTWT_18_27m_vs_WD_PTKO_18_27m[(WD_PTWT_18_27m_vs_WD_PTKO_18_27m$PValue < 0.05) & 
                                       ((WD_PTWT_18_27m_vs_WD_PTKO_18_27m$logFC > 0.58)|
                                          (WD_PTWT_18_27m_vs_WD_PTKO_18_27m$logFC < -0.58)),])
write.csv(WD_PTWT_18_27m_vs_WD_PTKO_18_27m[(WD_PTWT_18_27m_vs_WD_PTKO_18_27m$PValue < 0.05) & 
                                             ((WD_PTWT_18_27m_vs_WD_PTKO_18_27m$logFC > 0.58)|
                                                (WD_PTWT_18_27m_vs_WD_PTKO_18_27m$logFC < -0.58)),],
          file= "Results/EdgeR_PVal0.5_2fold/WD_PTWT_18_27m_vs_WD_PTKO_18_27m_2fold.csv",row.names=F)

# a5. WD_PTKO_18_27m vs CD_PTKO_18_27m <--- added later not in original comparisons

lrt.WD_PTKO_18_27m_vs_CD_PTKO_18_27m= glmLRT(fit,contrast=makeContrasts(groupsWD_PTKO_18_27m-groupsCD_PTKO_18_27m, levels=design))
WD_PTKO_18_27m_vs_CD_PTKO_18_27m<- as.data.frame(topTags(lrt.WD_PTKO_18_27m_vs_CD_PTKO_18_27m, n=Inf))
#add info
WD_PTKO_18_27m_vs_CD_PTKO_18_27m<-left_join(WD_PTKO_18_27m_vs_CD_PTKO_18_27m,mouse_mart_export,by="genes")
#pvalue 0.05 & FC = 2.0
dim(WD_PTKO_18_27m_vs_CD_PTKO_18_27m[(WD_PTKO_18_27m_vs_CD_PTKO_18_27m$PValue < 0.05) & 
                                       ((WD_PTKO_18_27m_vs_CD_PTKO_18_27m$logFC > 1)|
                                          (WD_PTKO_18_27m_vs_CD_PTKO_18_27m$logFC < -1)),])
write.csv(WD_PTKO_18_27m_vs_CD_PTKO_18_27m[(WD_PTKO_18_27m_vs_CD_PTKO_18_27m$PValue < 0.05) & 
                                             ((WD_PTKO_18_27m_vs_CD_PTKO_18_27m$logFC > 1)|
                                                (WD_PTKO_18_27m_vs_CD_PTKO_18_27m$logFC < -1)),],
          file= "Results/EdgeR_PVal0.5_2fold/WD_PTKO_18_27m_vs_CD_PTKO_18_27m_2fold.csv",row.names=F)
#pvalue 0.05 & FC = 1.5
dim(WD_PTKO_18_27m_vs_CD_PTKO_18_27m[(WD_PTKO_18_27m_vs_CD_PTKO_18_27m$PValue < 0.05) & 
                                       ((WD_PTKO_18_27m_vs_CD_PTKO_18_27m$logFC > 0.58)|
                                          (WD_PTKO_18_27m_vs_CD_PTKO_18_27m$logFC < -0.58)),])
write.csv(WD_PTKO_18_27m_vs_CD_PTKO_18_27m[(WD_PTKO_18_27m_vs_CD_PTKO_18_27m$PValue < 0.05) & 
                                             ((WD_PTKO_18_27m_vs_CD_PTKO_18_27m$logFC > 0.58)|
                                                (WD_PTKO_18_27m_vs_CD_PTKO_18_27m$logFC < -0.58)),],
          file= "Results/EdgeR_PVal0.5_2fold/WD_PTKO_18_27m_vs_CD_PTKO_18_27m_1.5fold.csv",row.names=F)


# c2. WD_PTWT_18_27m vs WD_PTWT_6m

lrt.WD_PTWT_18_27m_vs_WD_PTWT_6m= glmLRT(fit,contrast=makeContrasts(groupsWD_PTWT_18_27m-groupsWD_PTWT_6m, levels=design))
WD_PTWT_18_27m_vs_WD_PTWT_6m<- as.data.frame(topTags(lrt.WD_PTWT_18_27m_vs_WD_PTWT_6m, n=Inf))
#add info
WD_PTWT_18_27m_vs_WD_PTWT_6m<-left_join(WD_PTWT_18_27m_vs_WD_PTWT_6m,mouse_mart_export,by="genes")
#pvalue 0.05 & FC = 2.0
dim(WD_PTWT_18_27m_vs_WD_PTWT_6m[(WD_PTWT_18_27m_vs_WD_PTWT_6m$PValue < 0.05) & 
                                   ((WD_PTWT_18_27m_vs_WD_PTWT_6m$logFC > 1)|
                                      (WD_PTWT_18_27m_vs_WD_PTWT_6m$logFC < -1)),])
write.csv(WD_PTWT_18_27m_vs_WD_PTWT_6m[(WD_PTWT_18_27m_vs_WD_PTWT_6m$PValue < 0.05) & 
                                         ((WD_PTWT_18_27m_vs_WD_PTWT_6m$logFC > 1)|
                                            (WD_PTWT_18_27m_vs_WD_PTWT_6m$logFC < -1)),],
          file= "Results/EdgeR_PVal0.5_2fold/WD_PTWT_18_27m_vs_WD_PTWT_6m_2fold.csv",row.names=F)
#pvalue 0.05 & FC = 1.5
dim(WD_PTWT_18_27m_vs_WD_PTWT_6m[(WD_PTWT_18_27m_vs_WD_PTWT_6m$PValue < 0.05) & 
                                   ((WD_PTWT_18_27m_vs_WD_PTWT_6m$logFC > 0.58)|
                                      (WD_PTWT_18_27m_vs_WD_PTWT_6m$logFC < -0.58)),])
write.csv(WD_PTWT_18_27m_vs_WD_PTWT_6m[(WD_PTWT_18_27m_vs_WD_PTWT_6m$PValue < 0.05) & 
                                         ((WD_PTWT_18_27m_vs_WD_PTWT_6m$logFC > 0.58)|
                                            (WD_PTWT_18_27m_vs_WD_PTWT_6m$logFC < -0.58)),],
          file= "Results/EdgeR_PVal0.5_2fold/WD_PTWT_18_27m_vs_WD_PTWT_6m_1.5fold.csv",row.names=F)







########## A ############

# a1 CD_PTKO_18_27m vs CD_PTWT_18_27m

lrt.CD_PTKO_18_27m_vs_CD_PTWT_18_27m= glmLRT(fit,contrast=makeContrasts(groupsCD_PTKO_18_27m-groupsCD_PTWT_18_27m, levels=design))
CD_PTKO_18_27m_vs_CD_PTWT_18_27m<- as.data.frame(topTags(lrt.CD_PTKO_18_27m_vs_CD_PTWT_18_27m, n=Inf))
#add info
CD_PTKO_18_27m_vs_CD_PTWT_18_27m<-left_join(CD_PTKO_18_27m_vs_CD_PTWT_18_27m,mouse_mart_export,by="genes")
#
dim(CD_PTKO_18_27m_vs_CD_PTWT_18_27m[CD_PTKO_18_27m_vs_CD_PTWT_18_27m$PValue < 0.05,])
write.csv(CD_PTKO_18_27m_vs_CD_PTWT_18_27m[CD_PTKO_18_27m_vs_CD_PTWT_18_27m$PValue < 0.05,], 
          file= "Results/EdgeR_Pval/CD_PTKO_18_27m_vs_CD_PTWT_18_27m_Pval.csv",row.names=F)
dim(CD_PTKO_18_27m_vs_CD_PTWT_18_27m[CD_PTKO_18_27m_vs_CD_PTWT_18_27m$FDR < 0.5,])
write.csv(CD_PTKO_18_27m_vs_CD_PTWT_18_27m[CD_PTKO_18_27m_vs_CD_PTWT_18_27m$FDR < 0.05,], 
          file= "Results/EdgeR_FDR/CD_PTKO_18_27m_vs_CD_PTWT_18_27m_FDR.csv",row.names=F)
dim(CD_PTKO_18_27m_vs_CD_PTWT_18_27m[(CD_PTKO_18_27m_vs_CD_PTWT_18_27m$FDR < 0.05) & 
                               ((CD_PTKO_18_27m_vs_CD_PTWT_18_27m$logFC > 1)|
                                  (CD_PTKO_18_27m_vs_CD_PTWT_18_27m$logFC < -1)),])

dim(CD_PTKO_18_27m_vs_CD_PTWT_18_27m[(CD_PTKO_18_27m_vs_CD_PTWT_18_27m$PValue < 0.05) & 
                                       ((CD_PTKO_18_27m_vs_CD_PTWT_18_27m$logFC > 1)|
                                          (CD_PTKO_18_27m_vs_CD_PTWT_18_27m$logFC < -1)),])


write.csv(CD_PTKO_18_27m_vs_CD_PTWT_18_27m[(CD_PTKO_18_27m_vs_CD_PTWT_18_27m$FDR < 0.05) & 
                                     ((CD_PTKO_18_27m_vs_CD_PTWT_18_27m$logFC > 1)|
                                        (CD_PTKO_18_27m_vs_CD_PTWT_18_27m$logFC < -1)),],
          file= "Results/EdgeR_2fold/CD_PTKO_18_27m_vs_CD_PTWT_18_27m_2fold.csv",row.names=F)
write.csv(CD_PTKO_18_27m_vs_CD_PTWT_18_27m, file= "Results/EdgeR_ALL/CD_PTKO_18_27m_vs_CD_PTWT_18_27m_all.csv",row.names=F)
##fdr <0.1 and 2 fold change
dim(CD_PTKO_18_27m_vs_CD_PTWT_18_27m[(CD_PTKO_18_27m_vs_CD_PTWT_18_27m$FDR < 0.5) & 
                                       ((CD_PTKO_18_27m_vs_CD_PTWT_18_27m$logFC > 1)|
                                          (CD_PTKO_18_27m_vs_CD_PTWT_18_27m$logFC < -1)),])
write.csv(CD_PTKO_18_27m_vs_CD_PTWT_18_27m[(CD_PTKO_18_27m_vs_CD_PTWT_18_27m$FDR < 0.5) & 
                                             ((CD_PTKO_18_27m_vs_CD_PTWT_18_27m$logFC > 1)|
                                                (CD_PTKO_18_27m_vs_CD_PTWT_18_27m$logFC < -1)),],
          file= "Results/EdgeR_0.1/CD_PTKO_18_27m_vs_CD_PTWT_18_27m_2fold.csv",row.names=F)
##fdr <0.05 and 1.5 fold change
dim(CD_PTKO_18_27m_vs_CD_PTWT_18_27m[(CD_PTKO_18_27m_vs_CD_PTWT_18_27m$FDR < 0.05) & 
                                       ((CD_PTKO_18_27m_vs_CD_PTWT_18_27m$logFC > 0.26)|
                                          (CD_PTKO_18_27m_vs_CD_PTWT_18_27m$logFC < -0.26)),])
write.csv(CD_PTKO_18_27m_vs_CD_PTWT_18_27m[(CD_PTKO_18_27m_vs_CD_PTWT_18_27m$FDR < 0.05) & 
                                             ((CD_PTKO_18_27m_vs_CD_PTWT_18_27m$logFC > 0.58)|
                                                (CD_PTKO_18_27m_vs_CD_PTWT_18_27m$logFC < -0.58)),],
          file= "Results/EdgeR_FC_1.5/CD_PTKO_18_27m_vs_CD_PTWT_18_27m_FC_1.5.csv",row.names=F)
#fdr < 0.1 and 1.5 fold change
dim(CD_PTKO_18_27m_vs_CD_PTWT_18_27m[(CD_PTKO_18_27m_vs_CD_PTWT_18_27m$FDR < 0.1) & 
                                       ((CD_PTKO_18_27m_vs_CD_PTWT_18_27m$logFC > 0.58)|
                                          (CD_PTKO_18_27m_vs_CD_PTWT_18_27m$logFC < -0.58)),])
write.csv(CD_PTKO_18_27m_vs_CD_PTWT_18_27m[(CD_PTKO_18_27m_vs_CD_PTWT_18_27m$FDR < 0.1) & 
                                             ((CD_PTKO_18_27m_vs_CD_PTWT_18_27m$logFC > 0.58)|
                                                (CD_PTKO_18_27m_vs_CD_PTWT_18_27m$logFC < -0.58)),],
          file= "Results/EdgeR_FDR0.1_FC_1.5/CD_PTKO_18_27m_vs_CD_PTWT_18_27m_FDR0.1_FC_1.5.csv",row.names=F)
####
##

lrt.CD_PTWT_18_27m_vs_CD_PTKO_18_27m= glmLRT(fit,contrast=makeContrasts(groupsCD_PTWT_18_27m-groupsCD_PTKO_18_27m, levels=design))
CD_PTWT_18_27m_vs_CD_PTKO_18_27m<- as.data.frame(topTags(lrt.CD_PTWT_18_27m_vs_CD_PTKO_18_27m, n=Inf))
dim(CD_PTWT_18_27m_vs_CD_PTKO_18_27m[CD_PTWT_18_27m_vs_CD_PTKO_18_27m$PValue < 0.05,])
write.csv(CD_PTWT_18_27m_vs_CD_PTKO_18_27m[CD_PTWT_18_27m_vs_CD_PTKO_18_27m$PValue < 0.05,], 
          file= "Results_2/EdgeR_Pval/CD_PTWT_18_27m_vs_CD_PTKO_18_27m_Pval.csv",row.names=F)
dim(CD_PTWT_18_27m_vs_CD_PTKO_18_27m[CD_PTWT_18_27m_vs_CD_PTKO_18_27m$FDR < 0.05,])
write.csv(CD_PTWT_18_27m_vs_CD_PTKO_18_27m[CD_PTWT_18_27m_vs_CD_PTKO_18_27m$FDR < 0.05,], 
          file= "Results_2/EdgeR_FDR/CD_PTWT_18_27m_vs_CD_PTKO_18_27m_FDR.csv",row.names=F)
dim(CD_PTWT_18_27m_vs_CD_PTKO_18_27m[(CD_PTWT_18_27m_vs_CD_PTKO_18_27m$FDR < 0.05) & 
                                       ((CD_PTWT_18_27m_vs_CD_PTKO_18_27m$logFC > 1)|
                                          (CD_PTWT_18_27m_vs_CD_PTKO_18_27m$logFC < -1)),])
write.csv(CD_PTWT_18_27m_vs_CD_PTKO_18_27m[(CD_PTWT_18_27m_vs_CD_PTKO_18_27m$FDR < 0.05) & 
                                             ((CD_PTWT_18_27m_vs_CD_PTKO_18_27m$logFC > 1)|
                                                (CD_PTWT_18_27m_vs_CD_PTKO_18_27m$logFC < -1)),],
          file= "Results_2/EdgeR_2fold/CD_PTWT_18_27m_vs_CD_PTKO_18_27m_2fold.csv",row.names=F)
write.csv(CD_PTWT_18_27m_vs_CD_PTKO_18_27m, file= "Results_2/EdgeR_ALL/CD_PTWT_18_27m_vs_CD_PTKO_18_27m_all.csv",row.names=F)


# a2. WD_PTWT_18_27m vs CD_PTWT_18_27m

lrt.WD_PTWT_18_27m_vs_CD_PTWT_18_27m= glmLRT(fit,contrast=makeContrasts(groupsWD_PTWT_18_27m-groupsCD_PTWT_18_27m, levels=design))
WD_PTWT_18_27m_vs_CD_PTWT_18_27m<- as.data.frame(topTags(lrt.WD_PTWT_18_27m_vs_CD_PTWT_18_27m, n=Inf))
#add info
WD_PTWT_18_27m_vs_CD_PTWT_18_27m<-left_join(WD_PTWT_18_27m_vs_CD_PTWT_18_27m,mouse_mart_export,by="genes")
#
dim(WD_PTWT_18_27m_vs_CD_PTWT_18_27m[WD_PTWT_18_27m_vs_CD_PTWT_18_27m$PValue < 0.05,])
write.csv(WD_PTWT_18_27m_vs_CD_PTWT_18_27m[WD_PTWT_18_27m_vs_CD_PTWT_18_27m$PValue < 0.05,], 
          file= "Results/EdgeR_Pval/WD_PTWT_18_27m_vs_CD_PTWT_18_27m_Pval.csv",row.names=F)
dim(WD_PTWT_18_27m_vs_CD_PTWT_18_27m[WD_PTWT_18_27m_vs_CD_PTWT_18_27m$FDR < 0.05,])
write.csv(WD_PTWT_18_27m_vs_CD_PTWT_18_27m[WD_PTWT_18_27m_vs_CD_PTWT_18_27m$FDR < 0.05,], 
          file= "Results/EdgeR_FDR/WD_PTWT_18_27m_vs_CD_PTWT_18_27m_FDR.csv",row.names=F)
dim(WD_PTWT_18_27m_vs_CD_PTWT_18_27m[(WD_PTWT_18_27m_vs_CD_PTWT_18_27m$FDR < 0.05) & 
                                       ((WD_PTWT_18_27m_vs_CD_PTWT_18_27m$logFC > 1)|
                                          (WD_PTWT_18_27m_vs_CD_PTWT_18_27m$logFC < -1)),])
write.csv(WD_PTWT_18_27m_vs_CD_PTWT_18_27m[(WD_PTWT_18_27m_vs_CD_PTWT_18_27m$FDR < 0.05) & 
                                             ((WD_PTWT_18_27m_vs_CD_PTWT_18_27m$logFC > 1)|
                                                (WD_PTWT_18_27m_vs_CD_PTWT_18_27m$logFC < -1)),],
          file= "Results/EdgeR_2fold/WD_PTWT_18_27m_vs_CD_PTWT_18_27m_2fold.csv",row.names=F)
write.csv(WD_PTWT_18_27m_vs_CD_PTWT_18_27m, file= "Results/EdgeR_ALL/WD_PTWT_18_27m_vs_CD_PTWT_18_27m_all.csv",row.names=F)
dim(WD_PTWT_18_27m_vs_CD_PTWT_18_27m[(WD_PTWT_18_27m_vs_CD_PTWT_18_27m$FDR < 0.1) & 
                                       ((WD_PTWT_18_27m_vs_CD_PTWT_18_27m$logFC > 1)|
                                          (WD_PTWT_18_27m_vs_CD_PTWT_18_27m$logFC < -1)),])
write.csv(WD_PTWT_18_27m_vs_CD_PTWT_18_27m[(WD_PTWT_18_27m_vs_CD_PTWT_18_27m$FDR < 0.1) & 
                                             ((WD_PTWT_18_27m_vs_CD_PTWT_18_27m$logFC > 1)|
                                                (WD_PTWT_18_27m_vs_CD_PTWT_18_27m$logFC < -1)),],
          file= "Results/EdgeR_0.1/WD_PTWT_18_27m_vs_CD_PTWT_18_27m_2fold.csv",row.names=F)

##fdr <0.05 and 1.5 fold change
dim(WD_PTWT_18_27m_vs_CD_PTWT_18_27m[(WD_PTWT_18_27m_vs_CD_PTWT_18_27m$FDR < 0.05) & 
                                       ((WD_PTWT_18_27m_vs_CD_PTWT_18_27m$logFC > 0.58)|
                                          (WD_PTWT_18_27m_vs_CD_PTWT_18_27m$logFC < -0.58)),])
write.csv(WD_PTWT_18_27m_vs_CD_PTWT_18_27m[(WD_PTWT_18_27m_vs_CD_PTWT_18_27m$FDR < 0.05) & 
                                             ((WD_PTWT_18_27m_vs_CD_PTWT_18_27m$logFC > 0.58)|
                                                (WD_PTWT_18_27m_vs_CD_PTWT_18_27m$logFC < -0.58)),],
          file= "Results/EdgeR_FC_1.5/WD_PTWT_18_27m_vs_CD_PTWT_18_27m_FC_1.5.csv",row.names=F)
#fdr < 0.1 and 1.5 fold change
dim(WD_PTWT_18_27m_vs_CD_PTWT_18_27m[(WD_PTWT_18_27m_vs_CD_PTWT_18_27m$FDR < 0.1) & 
                                       ((WD_PTWT_18_27m_vs_CD_PTWT_18_27m$logFC > 0.58)|
                                          (WD_PTWT_18_27m_vs_CD_PTWT_18_27m$logFC < -0.58)),])
write.csv(WD_PTWT_18_27m_vs_CD_PTWT_18_27m[(WD_PTWT_18_27m_vs_CD_PTWT_18_27m$FDR < 0.1) & 
                                             ((WD_PTWT_18_27m_vs_CD_PTWT_18_27m$logFC > 0.58)|
                                                (WD_PTWT_18_27m_vs_CD_PTWT_18_27m$logFC < -0.58)),],
          file= "Results/EdgeR_FDR0.1_FC_1.5/WD_PTWT_18_27m_vs_CD_PTWT_18_27m_FDR0.1_FC_1.5.csv",row.names=F)
####




lrt.CD_PTWT_18_27m_vs_WD_PTWT_18_27m= glmLRT(fit,contrast=makeContrasts(groupsCD_PTWT_18_27m-groupsWD_PTWT_18_27m, levels=design))
CD_PTWT_18_27m_vs_WD_PTWT_18_27m<- as.data.frame(topTags(lrt.CD_PTWT_18_27m_vs_WD_PTWT_18_27m, n=Inf))
dim(CD_PTWT_18_27m_vs_WD_PTWT_18_27m[CD_PTWT_18_27m_vs_WD_PTWT_18_27m$PValue < 0.05,])
write.csv(CD_PTWT_18_27m_vs_WD_PTWT_18_27m[CD_PTWT_18_27m_vs_WD_PTWT_18_27m$PValue < 0.05,], 
          file= "Results_2/EdgeR_Pval/CD_PTWT_18_27m_vs_WD_PTWT_18_27m_Pval.csv",row.names=F)
dim(CD_PTWT_18_27m_vs_WD_PTWT_18_27m[CD_PTWT_18_27m_vs_WD_PTWT_18_27m$FDR < 0.05,])
write.csv(CD_PTWT_18_27m_vs_WD_PTWT_18_27m[CD_PTWT_18_27m_vs_WD_PTWT_18_27m$FDR < 0.05,], 
          file= "Results_2/EdgeR_FDR/CD_PTWT_18_27m_vs_WD_PTWT_18_27m_FDR.csv",row.names=F)
dim(CD_PTWT_18_27m_vs_WD_PTWT_18_27m[(CD_PTWT_18_27m_vs_WD_PTWT_18_27m$FDR < 0.05) & 
                                       ((CD_PTWT_18_27m_vs_WD_PTWT_18_27m$logFC > 1)|
                                          (CD_PTWT_18_27m_vs_WD_PTWT_18_27m$logFC < -1)),])
write.csv(CD_PTWT_18_27m_vs_WD_PTWT_18_27m[(CD_PTWT_18_27m_vs_WD_PTWT_18_27m$FDR < 0.05) & 
                                             ((CD_PTWT_18_27m_vs_WD_PTWT_18_27m$logFC > 1)|
                                                (CD_PTWT_18_27m_vs_WD_PTWT_18_27m$logFC < -1)),],
          file= "Results_2/EdgeR_2fold/CD_PTWT_18_27m_vs_WD_PTWT_18_27m_2fold.csv",row.names=F)
write.csv(CD_PTWT_18_27m_vs_WD_PTWT_18_27m, file= "Results_2/EdgeR_ALL/CD_PTWT_18_27m_vs_WD_PTWT_18_27m_all.csv",row.names=F)


# a3. WD_PTKO_18_27m vs CD_PTWT_18_27m

lrt.WD_PTKO_18_27m_vs_CD_PTWT_18_27m= glmLRT(fit,contrast=makeContrasts(groupsWD_PTKO_18_27m-groupsCD_PTWT_18_27m, levels=design))
WD_PTKO_18_27m_vs_CD_PTWT_18_27m<- as.data.frame(topTags(lrt.WD_PTKO_18_27m_vs_CD_PTWT_18_27m, n=Inf))
#add info
WD_PTKO_18_27m_vs_CD_PTWT_18_27m<-left_join(WD_PTKO_18_27m_vs_CD_PTWT_18_27m,mouse_mart_export,by="genes")
#
dim(WD_PTKO_18_27m_vs_CD_PTWT_18_27m[WD_PTKO_18_27m_vs_CD_PTWT_18_27m$PValue < 0.05,])
write.csv(WD_PTKO_18_27m_vs_CD_PTWT_18_27m[WD_PTKO_18_27m_vs_CD_PTWT_18_27m$PValue < 0.05,], 
          file= "Results/EdgeR_Pval/WD_PTKO_18_27m_vs_CD_PTWT_18_27m_Pval.csv",row.names=F)
dim(WD_PTKO_18_27m_vs_CD_PTWT_18_27m[WD_PTKO_18_27m_vs_CD_PTWT_18_27m$FDR < 0.05,])
write.csv(WD_PTKO_18_27m_vs_CD_PTWT_18_27m[WD_PTKO_18_27m_vs_CD_PTWT_18_27m$FDR < 0.05,], 
          file= "Results/EdgeR_FDR/WD_PTKO_18_27m_vs_CD_PTWT_18_27m_FDR.csv",row.names=F)
dim(WD_PTKO_18_27m_vs_CD_PTWT_18_27m[(WD_PTKO_18_27m_vs_CD_PTWT_18_27m$FDR < 0.05) & 
                                       ((WD_PTKO_18_27m_vs_CD_PTWT_18_27m$logFC > 1)|
                                          (WD_PTKO_18_27m_vs_CD_PTWT_18_27m$logFC < -1)),])
write.csv(WD_PTKO_18_27m_vs_CD_PTWT_18_27m[(WD_PTKO_18_27m_vs_CD_PTWT_18_27m$FDR < 0.05) & 
                                             ((WD_PTKO_18_27m_vs_CD_PTWT_18_27m$logFC > 1)|
                                                (WD_PTKO_18_27m_vs_CD_PTWT_18_27m$logFC < -1)),],
          file= "Results/EdgeR_2fold/WD_PTKO_18_27m_vs_CD_PTWT_18_27m_2fold.csv",row.names=F)
write.csv(WD_PTKO_18_27m_vs_CD_PTWT_18_27m, file= "Results/EdgeR_ALL/WD_PTKO_18_27m_vs_CD_PTWT_18_27m_all.csv",row.names=F)
dim(WD_PTKO_18_27m_vs_CD_PTWT_18_27m[(WD_PTKO_18_27m_vs_CD_PTWT_18_27m$FDR < 0.1) & 
                                       ((WD_PTKO_18_27m_vs_CD_PTWT_18_27m$logFC > 1)|
                                          (WD_PTKO_18_27m_vs_CD_PTWT_18_27m$logFC < -1)),])
write.csv(WD_PTKO_18_27m_vs_CD_PTWT_18_27m[(WD_PTKO_18_27m_vs_CD_PTWT_18_27m$FDR < 0.1) & 
                                             ((WD_PTKO_18_27m_vs_CD_PTWT_18_27m$logFC > 1)|
                                                (WD_PTKO_18_27m_vs_CD_PTWT_18_27m$logFC < -1)),],
          file= "Results/EdgeR_0.1/WD_PTKO_18_27m_vs_CD_PTWT_18_27m_2fold.csv",row.names=F)
##fdr <0.05 and 1.5 fold change
dim(WD_PTWT_18_27m_vs_CD_PTWT_18_27m[(WD_PTWT_18_27m_vs_CD_PTWT_18_27m$FDR < 0.05) & 
                                       ((WD_PTWT_18_27m_vs_CD_PTWT_18_27m$logFC > 0.58)|
                                          (WD_PTWT_18_27m_vs_CD_PTWT_18_27m$logFC < -0.58)),])
write.csv(WD_PTWT_18_27m_vs_CD_PTWT_18_27m[(WD_PTWT_18_27m_vs_CD_PTWT_18_27m$FDR < 0.05) & 
                                             ((WD_PTWT_18_27m_vs_CD_PTWT_18_27m$logFC > 0.58)|
                                                (WD_PTWT_18_27m_vs_CD_PTWT_18_27m$logFC < -0.58)),],
          file= "Results/EdgeR_FC_1.5/WD_PTWT_18_27m_vs_CD_PTWT_18_27m_FC_1.5.csv",row.names=F)
#fdr < 0.1 and 1.5 fold change
dim(WD_PTWT_18_27m_vs_CD_PTWT_18_27m[(WD_PTWT_18_27m_vs_CD_PTWT_18_27m$FDR < 0.1) & 
                                       ((WD_PTWT_18_27m_vs_CD_PTWT_18_27m$logFC > 0.58)|
                                          (WD_PTWT_18_27m_vs_CD_PTWT_18_27m$logFC < -0.58)),])
write.csv(WD_PTWT_18_27m_vs_CD_PTWT_18_27m[(WD_PTWT_18_27m_vs_CD_PTWT_18_27m$FDR < 0.1) & 
                                             ((WD_PTWT_18_27m_vs_CD_PTWT_18_27m$logFC > 0.58)|
                                                (WD_PTWT_18_27m_vs_CD_PTWT_18_27m$logFC < -0.58)),],
          file= "Results/EdgeR_FDR0.1_FC_1.5/WD_PTWT_18_27m_vs_CD_PTWT_18_27m_FDR0.1_FC_1.5.csv",row.names=F)
####





lrt.CD_PTWT_18_27m_vs_WD_PTKO_18_27m= glmLRT(fit,contrast=makeContrasts(groupsCD_PTWT_18_27m-groupsWD_PTKO_18_27m, levels=design))
CD_PTWT_18_27m_vs_WD_PTKO_18_27m<- as.data.frame(topTags(lrt.CD_PTWT_18_27m_vs_WD_PTKO_18_27m, n=Inf))
dim(CD_PTWT_18_27m_vs_WD_PTKO_18_27m[CD_PTWT_18_27m_vs_WD_PTKO_18_27m$PValue < 0.05,])
write.csv(CD_PTWT_18_27m_vs_WD_PTKO_18_27m[CD_PTWT_18_27m_vs_WD_PTKO_18_27m$PValue < 0.05,], 
          file= "Results_2/EdgeR_Pval/CD_PTWT_18_27m_vs_WD_PTKO_18_27m_Pval.csv",row.names=F)
dim(CD_PTWT_18_27m_vs_WD_PTKO_18_27m[CD_PTWT_18_27m_vs_WD_PTKO_18_27m$FDR < 0.05,])
write.csv(CD_PTWT_18_27m_vs_WD_PTKO_18_27m[CD_PTWT_18_27m_vs_WD_PTKO_18_27m$FDR < 0.05,], 
          file= "Results_2/EdgeR_FDR/CD_PTWT_18_27m_vs_WD_PTKO_18_27m_FDR.csv",row.names=F)
dim(CD_PTWT_18_27m_vs_WD_PTKO_18_27m[(CD_PTWT_18_27m_vs_WD_PTKO_18_27m$FDR < 0.05) & 
                                       ((CD_PTWT_18_27m_vs_WD_PTKO_18_27m$logFC > 1)|
                                          (CD_PTWT_18_27m_vs_WD_PTKO_18_27m$logFC < -1)),])
write.csv(CD_PTWT_18_27m_vs_WD_PTKO_18_27m[(CD_PTWT_18_27m_vs_WD_PTKO_18_27m$FDR < 0.05) & 
                                             ((CD_PTWT_18_27m_vs_WD_PTKO_18_27m$logFC > 1)|
                                                (CD_PTWT_18_27m_vs_WD_PTKO_18_27m$logFC < -1)),],
          file= "Results_2/EdgeR_2fold/CD_PTWT_18_27m_vs_WD_PTKO_18_27m_2fold.csv",row.names=F)
write.csv(CD_PTWT_18_27m_vs_WD_PTKO_18_27m, file= "Results_2/EdgeR_ALL/CD_PTWT_18_27m_vs_WD_PTKO_18_27m_all.csv",row.names=F)

# a4. WD_PTWT_18_27m vs WD_PTKO_18_27m

lrt.WD_PTWT_18_27m_vs_WD_PTKO_18_27m= glmLRT(fit,contrast=makeContrasts(groupsWD_PTWT_18_27m-groupsWD_PTKO_18_27m, levels=design))
WD_PTWT_18_27m_vs_WD_PTKO_18_27m<- as.data.frame(topTags(lrt.WD_PTWT_18_27m_vs_WD_PTKO_18_27m, n=Inf))
#add info
WD_PTWT_18_27m_vs_WD_PTKO_18_27m<-left_join(WD_PTWT_18_27m_vs_WD_PTKO_18_27m,mouse_mart_export,by="genes")
#
dim(WD_PTWT_18_27m_vs_WD_PTKO_18_27m[WD_PTWT_18_27m_vs_WD_PTKO_18_27m$PValue < 0.05,])
write.csv(WD_PTWT_18_27m_vs_WD_PTKO_18_27m[WD_PTWT_18_27m_vs_WD_PTKO_18_27m$PValue < 0.05,], 
          file= "Results/EdgeR_Pval/WD_PTWT_18_27m_vs_WD_PTKO_18_27m_Pval.csv",row.names=F)
dim(WD_PTWT_18_27m_vs_WD_PTKO_18_27m[WD_PTWT_18_27m_vs_WD_PTKO_18_27m$FDR < 0.05,])
write.csv(WD_PTWT_18_27m_vs_WD_PTKO_18_27m[WD_PTWT_18_27m_vs_WD_PTKO_18_27m$FDR < 0.05,], 
          file= "Results/EdgeR_FDR/WD_PTWT_18_27m_vs_WD_PTKO_18_27m_FDR.csv",row.names=F)
dim(WD_PTWT_18_27m_vs_WD_PTKO_18_27m[(WD_PTWT_18_27m_vs_WD_PTKO_18_27m$FDR < 0.05) & 
                                       ((WD_PTWT_18_27m_vs_WD_PTKO_18_27m$logFC > 1)|
                                          (WD_PTWT_18_27m_vs_WD_PTKO_18_27m$logFC < -1)),])
write.csv(WD_PTWT_18_27m_vs_WD_PTKO_18_27m[(WD_PTWT_18_27m_vs_WD_PTKO_18_27m$FDR < 0.05) & 
                                             ((WD_PTWT_18_27m_vs_WD_PTKO_18_27m$logFC > 1)|
                                                (WD_PTWT_18_27m_vs_WD_PTKO_18_27m$logFC < -1)),],
          file= "Results/EdgeR_2fold/WD_PTWT_18_27m_vs_WD_PTKO_18_27m_2fold.csv",row.names=F)
write.csv(WD_PTWT_18_27m_vs_WD_PTKO_18_27m, file= "Results/EdgeR_ALL/WD_PTWT_18_27m_vs_WD_PTKO_18_27m_all.csv",row.names=F)
dim(WD_PTWT_18_27m_vs_WD_PTKO_18_27m[(WD_PTWT_18_27m_vs_WD_PTKO_18_27m$FDR < 0.1) & 
                                       ((WD_PTWT_18_27m_vs_WD_PTKO_18_27m$logFC > 1)|
                                          (WD_PTWT_18_27m_vs_WD_PTKO_18_27m$logFC < -1)),])
write.csv(WD_PTWT_18_27m_vs_WD_PTKO_18_27m[(WD_PTWT_18_27m_vs_WD_PTKO_18_27m$FDR < 0.1) & 
                                             ((WD_PTWT_18_27m_vs_WD_PTKO_18_27m$logFC > 1)|
                                                (WD_PTWT_18_27m_vs_WD_PTKO_18_27m$logFC < -1)),],
          file= "Results/EdgeR_0.1/WD_PTWT_18_27m_vs_WD_PTKO_18_27m_2fold.csv",row.names=F)
##no vice-versa here


############ B ##############

# b1 CD_PTWT_16m vs CD_PTWT_6m

lrt.CD_PTWT_16m_vs_CD_PTWT_6m= glmLRT(fit,contrast=makeContrasts(groupsCD_PTWT_16m-groupsCD_PTWT_6m, levels=design))
CD_PTWT_16m_vs_CD_PTWT_6m<- as.data.frame(topTags(lrt.CD_PTWT_16m_vs_CD_PTWT_6m, n=Inf))
#add info
CD_PTWT_16m_vs_CD_PTWT_6m<-left_join(CD_PTWT_16m_vs_CD_PTWT_6m,mouse_mart_export,by="genes")
#
dim(CD_PTWT_16m_vs_CD_PTWT_6m[CD_PTWT_16m_vs_CD_PTWT_6m$PValue < 0.05,])
write.csv(CD_PTWT_16m_vs_CD_PTWT_6m[CD_PTWT_16m_vs_CD_PTWT_6m$PValue < 0.05,], 
          file= "Results/EdgeR_Pval/CD_PTWT_16m_vs_CD_PTWT_6m_Pval.csv",row.names=F)
dim(CD_PTWT_16m_vs_CD_PTWT_6m[CD_PTWT_16m_vs_CD_PTWT_6m$FDR < 0.05,])
write.csv(CD_PTWT_16m_vs_CD_PTWT_6m[CD_PTWT_16m_vs_CD_PTWT_6m$FDR < 0.05,], 
          file= "Results/EdgeR_FDR/CD_PTWT_16m_vs_CD_PTWT_6m_FDR.csv",row.names=F)
dim(CD_PTWT_16m_vs_CD_PTWT_6m[(CD_PTWT_16m_vs_CD_PTWT_6m$FDR < 0.05) & 
                                       ((CD_PTWT_16m_vs_CD_PTWT_6m$logFC > 1)|
                                          (CD_PTWT_16m_vs_CD_PTWT_6m$logFC < -1)),])
write.csv(CD_PTWT_16m_vs_CD_PTWT_6m[(CD_PTWT_16m_vs_CD_PTWT_6m$FDR < 0.05) & 
                                             ((CD_PTWT_16m_vs_CD_PTWT_6m$logFC > 1)|
                                                (CD_PTWT_16m_vs_CD_PTWT_6m$logFC < -1)),],
          file= "Results/EdgeR_2fold/CD_PTWT_16m_vs_CD_PTWT_6m_2fold.csv",row.names=F)
write.csv(CD_PTWT_16m_vs_CD_PTWT_6m, file= "Results/EdgeR_ALL/CD_PTWT_16m_vs_CD_PTWT_6m_all.csv",row.names=F)
dim(CD_PTWT_16m_vs_CD_PTWT_6m[(CD_PTWT_16m_vs_CD_PTWT_6m$FDR < 0.1) & 
                                ((CD_PTWT_16m_vs_CD_PTWT_6m$logFC > 1)|
                                   (CD_PTWT_16m_vs_CD_PTWT_6m$logFC < -1)),])
write.csv(CD_PTWT_16m_vs_CD_PTWT_6m[(CD_PTWT_16m_vs_CD_PTWT_6m$FDR < 0.1) & 
                                      ((CD_PTWT_16m_vs_CD_PTWT_6m$logFC > 1)|
                                         (CD_PTWT_16m_vs_CD_PTWT_6m$logFC < -1)),],
          file= "Results/EdgeR_0.1/CD_PTWT_16m_vs_CD_PTWT_6m_2fold.csv",row.names=F)


lrt.CD_PTWT_6m_vs_CD_PTWT_16m= glmLRT(fit,contrast=makeContrasts(groupsCD_PTWT_6m-groupsCD_PTWT_16m, levels=design))
CD_PTWT_6m_vs_CD_PTWT_16m<- as.data.frame(topTags(lrt.CD_PTWT_6m_vs_CD_PTWT_16m, n=Inf))
dim(CD_PTWT_6m_vs_CD_PTWT_16m[CD_PTWT_6m_vs_CD_PTWT_16m$PValue < 0.05,])
write.csv(CD_PTWT_6m_vs_CD_PTWT_16m[CD_PTWT_6m_vs_CD_PTWT_16m$PValue < 0.05,], 
          file= "Results_2/EdgeR_Pval/CD_PTWT_6m_vs_CD_PTWT_16m_Pval.csv",row.names=F)
dim(CD_PTWT_6m_vs_CD_PTWT_16m[CD_PTWT_6m_vs_CD_PTWT_16m$FDR < 0.05,])
write.csv(CD_PTWT_6m_vs_CD_PTWT_16m[CD_PTWT_6m_vs_CD_PTWT_16m$FDR < 0.05,], 
          file= "Results_2/EdgeR_FDR/CD_PTWT_6m_vs_CD_PTWT_16m_FDR.csv",row.names=F)
dim(CD_PTWT_6m_vs_CD_PTWT_16m[(CD_PTWT_6m_vs_CD_PTWT_16m$FDR < 0.05) & 
                                       ((CD_PTWT_6m_vs_CD_PTWT_16m$logFC > 1)|
                                          (CD_PTWT_6m_vs_CD_PTWT_16m$logFC < -1)),])
write.csv(CD_PTWT_6m_vs_CD_PTWT_16m[(CD_PTWT_6m_vs_CD_PTWT_16m$FDR < 0.05) & 
                                             ((CD_PTWT_6m_vs_CD_PTWT_16m$logFC > 1)|
                                                (CD_PTWT_6m_vs_CD_PTWT_16m$logFC < -1)),],
          file= "Results_2/EdgeR_2fold/CD_PTWT_6m_vs_CD_PTWT_16m_2fold.csv",row.names=F)
write.csv(CD_PTWT_6m_vs_CD_PTWT_16m, file= "Results_2/EdgeR_ALL/CD_PTWT_6m_vs_CD_PTWT_16m_all.csv",row.names=F)


# b2. CD_PTWT_18_27m vs CD_PTWT_6m

lrt.CD_PTWT_18_27m_vs_CD_PTWT_6m= glmLRT(fit,contrast=makeContrasts(groupsCD_PTWT_18_27m-groupsCD_PTWT_6m, levels=design))
CD_PTWT_18_27m_vs_CD_PTWT_6m<- as.data.frame(topTags(lrt.CD_PTWT_18_27m_vs_CD_PTWT_6m, n=Inf))
#add info
CD_PTWT_18_27m_vs_CD_PTWT_6m<-left_join(CD_PTWT_18_27m_vs_CD_PTWT_6m,mouse_mart_export,by="genes")
#
dim(CD_PTWT_18_27m_vs_CD_PTWT_6m[CD_PTWT_18_27m_vs_CD_PTWT_6m$PValue < 0.05,])
write.csv(CD_PTWT_18_27m_vs_CD_PTWT_6m[CD_PTWT_18_27m_vs_CD_PTWT_6m$PValue < 0.05,], 
          file= "Results/EdgeR_Pval/CD_PTWT_18_27m_vs_CD_PTWT_6m_Pval.csv",row.names=F)
dim(CD_PTWT_18_27m_vs_CD_PTWT_6m[CD_PTWT_18_27m_vs_CD_PTWT_6m$FDR < 0.05,])
write.csv(CD_PTWT_18_27m_vs_CD_PTWT_6m[CD_PTWT_18_27m_vs_CD_PTWT_6m$FDR < 0.05,], 
          file= "Results/EdgeR_FDR/CD_PTWT_18_27m_vs_CD_PTWT_6m_FDR.csv",row.names=F)
dim(CD_PTWT_18_27m_vs_CD_PTWT_6m[(CD_PTWT_18_27m_vs_CD_PTWT_6m$FDR < 0.05) & 
                                       ((CD_PTWT_18_27m_vs_CD_PTWT_6m$logFC > 1)|
                                          (CD_PTWT_18_27m_vs_CD_PTWT_6m$logFC < -1)),])
write.csv(CD_PTWT_18_27m_vs_CD_PTWT_6m[(CD_PTWT_18_27m_vs_CD_PTWT_6m$FDR < 0.05) & 
                                             ((CD_PTWT_18_27m_vs_CD_PTWT_6m$logFC > 1)|
                                                (CD_PTWT_18_27m_vs_CD_PTWT_6m$logFC < -1)),],
          file= "Results/EdgeR_2fold/CD_PTWT_18_27m_vs_CD_PTWT_6m_2fold.csv",row.names=F)
write.csv(CD_PTWT_18_27m_vs_CD_PTWT_6m, file= "Results/EdgeR_ALL/CD_PTWT_18_27m_vs_CD_PTWT_6m_all.csv",row.names=F)
dim(CD_PTWT_18_27m_vs_CD_PTWT_6m[(CD_PTWT_18_27m_vs_CD_PTWT_6m$FDR < 0.1) & 
                                   ((CD_PTWT_18_27m_vs_CD_PTWT_6m$logFC > 1)|
                                      (CD_PTWT_18_27m_vs_CD_PTWT_6m$logFC < -1)),])
write.csv(CD_PTWT_18_27m_vs_CD_PTWT_6m[(CD_PTWT_18_27m_vs_CD_PTWT_6m$FDR < 0.1) & 
                                         ((CD_PTWT_18_27m_vs_CD_PTWT_6m$logFC > 1)|
                                            (CD_PTWT_18_27m_vs_CD_PTWT_6m$logFC < -1)),],
          file= "Results/EdgeR_0.1/CD_PTWT_18_27m_vs_CD_PTWT_6m_2fold.csv",row.names=F)



lrt.CD_PTWT_6m_vs_CD_PTWT_18_27m= glmLRT(fit,contrast=makeContrasts(groupsCD_PTWT_6m-groupsCD_PTWT_18_27m, levels=design))
CD_PTWT_6m_vs_CD_PTWT_18_27m<- as.data.frame(topTags(lrt.CD_PTWT_6m_vs_CD_PTWT_18_27m, n=Inf))
dim(CD_PTWT_6m_vs_CD_PTWT_18_27m[CD_PTWT_6m_vs_CD_PTWT_18_27m$PValue < 0.05,])
write.csv(CD_PTWT_6m_vs_CD_PTWT_18_27m[CD_PTWT_6m_vs_CD_PTWT_18_27m$PValue < 0.05,], 
          file= "Results_2/EdgeR_Pval/CD_PTWT_6m_vs_CD_PTWT_18_27m_Pval.csv",row.names=F)
dim(CD_PTWT_6m_vs_CD_PTWT_18_27m[CD_PTWT_6m_vs_CD_PTWT_18_27m$FDR < 0.05,])
write.csv(CD_PTWT_6m_vs_CD_PTWT_18_27m[CD_PTWT_6m_vs_CD_PTWT_18_27m$FDR < 0.05,], 
          file= "Results_2/EdgeR_FDR/CD_PTWT_6m_vs_CD_PTWT_18_27m_FDR.csv",row.names=F)
dim(CD_PTWT_6m_vs_CD_PTWT_18_27m[(CD_PTWT_6m_vs_CD_PTWT_18_27m$FDR < 0.05) & 
                                       ((CD_PTWT_6m_vs_CD_PTWT_18_27m$logFC > 1)|
                                          (CD_PTWT_6m_vs_CD_PTWT_18_27m$logFC < -1)),])
write.csv(CD_PTWT_6m_vs_CD_PTWT_18_27m[(CD_PTWT_6m_vs_CD_PTWT_18_27m$FDR < 0.05) & 
                                             ((CD_PTWT_6m_vs_CD_PTWT_18_27m$logFC > 1)|
                                                (CD_PTWT_6m_vs_CD_PTWT_18_27m$logFC < -1)),],
          file= "Results_2/EdgeR_2fold/CD_PTWT_6m_vs_CD_PTWT_18_27m_2fold.csv",row.names=F)
write.csv(CD_PTWT_6m_vs_CD_PTWT_18_27m, file= "Results_2/EdgeR_ALL/CD_PTWT_6m_vs_CD_PTWT_18_27m_all.csv",row.names=F)

############ C ############

# C1 WD_PTWT_16m vs WD_PTWT_6m

lrt.WD_PTWT_16m_vs_WD_PTWT_6m= glmLRT(fit,contrast=makeContrasts(groupsWD_PTWT_16m-groupsWD_PTWT_6m, levels=design))
WD_PTWT_16m_vs_WD_PTWT_6m<- as.data.frame(topTags(lrt.WD_PTWT_16m_vs_WD_PTWT_6m, n=Inf))
#add info
WD_PTWT_16m_vs_WD_PTWT_6m<-left_join(WD_PTWT_16m_vs_WD_PTWT_6m,mouse_mart_export,by="genes")
#
dim(WD_PTWT_16m_vs_WD_PTWT_6m[WD_PTWT_16m_vs_WD_PTWT_6m$PValue < 0.05,])
write.csv(WD_PTWT_16m_vs_WD_PTWT_6m[WD_PTWT_16m_vs_WD_PTWT_6m$PValue < 0.05,], 
          file= "Results/EdgeR_Pval/WD_PTWT_16m_vs_WD_PTWT_6m_Pval.csv",row.names=F)
dim(WD_PTWT_16m_vs_WD_PTWT_6m[WD_PTWT_16m_vs_WD_PTWT_6m$FDR < 0.05,])
write.csv(WD_PTWT_16m_vs_WD_PTWT_6m[WD_PTWT_16m_vs_WD_PTWT_6m$FDR < 0.05,], 
          file= "Results/EdgeR_FDR/WD_PTWT_16m_vs_WD_PTWT_6m_FDR.csv",row.names=F)
dim(WD_PTWT_16m_vs_WD_PTWT_6m[(WD_PTWT_16m_vs_WD_PTWT_6m$FDR < 0.05) & 
                                ((WD_PTWT_16m_vs_WD_PTWT_6m$logFC > 1)|
                                   (WD_PTWT_16m_vs_WD_PTWT_6m$logFC < -1)),])
write.csv(WD_PTWT_16m_vs_WD_PTWT_6m[(WD_PTWT_16m_vs_WD_PTWT_6m$FDR < 0.05) & 
                                      ((WD_PTWT_16m_vs_WD_PTWT_6m$logFC > 1)|
                                         (WD_PTWT_16m_vs_WD_PTWT_6m$logFC < -1)),],
          file= "Results/EdgeR_2fold/WD_PTWT_16m_vs_WD_PTWT_6m_2fold.csv",row.names=F)
write.csv(WD_PTWT_16m_vs_WD_PTWT_6m, file= "Results/EdgeR_ALL/WD_PTWT_16m_vs_WD_PTWT_6m_all.csv",row.names=F)
dim(WD_PTWT_16m_vs_WD_PTWT_6m[(WD_PTWT_16m_vs_WD_PTWT_6m$FDR < 0.1) & 
                                ((WD_PTWT_16m_vs_WD_PTWT_6m$logFC > 1)|
                                   (WD_PTWT_16m_vs_WD_PTWT_6m$logFC < -1)),])
write.csv(WD_PTWT_16m_vs_WD_PTWT_6m[(WD_PTWT_16m_vs_WD_PTWT_6m$FDR < 0.1) & 
                                      ((WD_PTWT_16m_vs_WD_PTWT_6m$logFC > 1)|
                                         (WD_PTWT_16m_vs_WD_PTWT_6m$logFC < -1)),],
          file= "Results/EdgeR_0.1/WD_PTWT_16m_vs_WD_PTWT_6m_2fold.csv",row.names=F)


lrt.WD_PTWT_6m_vs_WD_PTWT_16m= glmLRT(fit,contrast=makeContrasts(groupsWD_PTWT_6m-groupsWD_PTWT_16m, levels=design))
WD_PTWT_6m_vs_WD_PTWT_16m<- as.data.frame(topTags(lrt.WD_PTWT_6m_vs_WD_PTWT_16m, n=Inf))
dim(WD_PTWT_6m_vs_WD_PTWT_16m[WD_PTWT_6m_vs_WD_PTWT_16m$PValue < 0.05,])
write.csv(WD_PTWT_6m_vs_WD_PTWT_16m[WD_PTWT_6m_vs_WD_PTWT_16m$PValue < 0.05,], 
          file= "Results_2/EdgeR_Pval/WD_PTWT_6m_vs_WD_PTWT_16m_Pval.csv",row.names=F)
dim(WD_PTWT_6m_vs_WD_PTWT_16m[WD_PTWT_6m_vs_WD_PTWT_16m$FDR < 0.05,])
write.csv(WD_PTWT_6m_vs_WD_PTWT_16m[WD_PTWT_6m_vs_WD_PTWT_16m$FDR < 0.05,], 
          file= "Results_2/EdgeR_FDR/WD_PTWT_6m_vs_WD_PTWT_16m_FDR.csv",row.names=F)
dim(WD_PTWT_6m_vs_WD_PTWT_16m[(WD_PTWT_6m_vs_WD_PTWT_16m$FDR < 0.05) & 
                                ((WD_PTWT_6m_vs_WD_PTWT_16m$logFC > 1)|
                                   (WD_PTWT_6m_vs_WD_PTWT_16m$logFC < -1)),])
write.csv(WD_PTWT_6m_vs_WD_PTWT_16m[(WD_PTWT_6m_vs_WD_PTWT_16m$FDR < 0.05) & 
                                      ((WD_PTWT_6m_vs_WD_PTWT_16m$logFC > 1)|
                                         (WD_PTWT_6m_vs_WD_PTWT_16m$logFC < -1)),],
          file= "Results_2/EdgeR_2fold/WD_PTWT_6m_vs_WD_PTWT_16m_2fold.csv",row.names=F)
write.csv(WD_PTWT_6m_vs_WD_PTWT_16m, file= "Results_2/EdgeR_ALL/WD_PTWT_6m_vs_WD_PTWT_16m_all.csv",row.names=F)


# c2. WD_PTWT_18_27m vs WD_PTWT_6m

lrt.WD_PTWT_18_27m_vs_WD_PTWT_6m= glmLRT(fit,contrast=makeContrasts(groupsWD_PTWT_18_27m-groupsWD_PTWT_6m, levels=design))
WD_PTWT_18_27m_vs_WD_PTWT_6m<- as.data.frame(topTags(lrt.WD_PTWT_18_27m_vs_WD_PTWT_6m, n=Inf))
#add info
WD_PTWT_18_27m_vs_WD_PTWT_6m<-left_join(WD_PTWT_18_27m_vs_WD_PTWT_6m,mouse_mart_export,by="genes")
#
dim(WD_PTWT_18_27m_vs_WD_PTWT_6m[WD_PTWT_18_27m_vs_WD_PTWT_6m$PValue < 0.05,])
write.csv(WD_PTWT_18_27m_vs_WD_PTWT_6m[WD_PTWT_18_27m_vs_WD_PTWT_6m$PValue < 0.05,], 
          file= "Results/EdgeR_Pval/WD_PTWT_18_27m_vs_WD_PTWT_6m_Pval.csv",row.names=F)
dim(WD_PTWT_18_27m_vs_WD_PTWT_6m[WD_PTWT_18_27m_vs_WD_PTWT_6m$FDR < 0.05,])
write.csv(WD_PTWT_18_27m_vs_WD_PTWT_6m[WD_PTWT_18_27m_vs_WD_PTWT_6m$FDR < 0.05,], 
          file= "Results/EdgeR_FDR/WD_PTWT_18_27m_vs_WD_PTWT_6m_FDR.csv",row.names=F)
dim(WD_PTWT_18_27m_vs_WD_PTWT_6m[(WD_PTWT_18_27m_vs_WD_PTWT_6m$FDR < 0.05) & 
                                   ((WD_PTWT_18_27m_vs_WD_PTWT_6m$logFC > 1)|
                                      (WD_PTWT_18_27m_vs_WD_PTWT_6m$logFC < -1)),])
write.csv(WD_PTWT_18_27m_vs_WD_PTWT_6m[(WD_PTWT_18_27m_vs_WD_PTWT_6m$FDR < 0.05) & 
                                         ((WD_PTWT_18_27m_vs_WD_PTWT_6m$logFC > 1)|
                                            (WD_PTWT_18_27m_vs_WD_PTWT_6m$logFC < -1)),],
          file= "Results/EdgeR_2fold/WD_PTWT_18_27m_vs_WD_PTWT_6m_2fold.csv",row.names=F)
write.csv(WD_PTWT_18_27m_vs_WD_PTWT_6m, file= "Results/EdgeR_ALL/WD_PTWT_18_27m_vs_WD_PTWT_6m_all.csv",row.names=F)
dim(WD_PTWT_18_27m_vs_WD_PTWT_6m[(WD_PTWT_18_27m_vs_WD_PTWT_6m$FDR < 0.1) & 
                                   ((WD_PTWT_18_27m_vs_WD_PTWT_6m$logFC > 1)|
                                      (WD_PTWT_18_27m_vs_WD_PTWT_6m$logFC < -1)),])
write.csv(WD_PTWT_18_27m_vs_WD_PTWT_6m[(WD_PTWT_18_27m_vs_WD_PTWT_6m$FDR < 0.1) & 
                                         ((WD_PTWT_18_27m_vs_WD_PTWT_6m$logFC > 1)|
                                            (WD_PTWT_18_27m_vs_WD_PTWT_6m$logFC < -1)),],
          file= "Results/EdgeR_0.1/WD_PTWT_18_27m_vs_WD_PTWT_6m_2fold.csv",row.names=F)

##
lrt.WD_PTWT_6m_vs_WD_PTWT_18_27m= glmLRT(fit,contrast=makeContrasts(groupsWD_PTWT_6m-groupsWD_PTWT_18_27m, levels=design))
WD_PTWT_6m_vs_WD_PTWT_18_27m<- as.data.frame(topTags(lrt.WD_PTWT_6m_vs_WD_PTWT_18_27m, n=Inf))
dim(WD_PTWT_6m_vs_WD_PTWT_18_27m[WD_PTWT_6m_vs_WD_PTWT_18_27m$PValue < 0.05,])
write.csv(WD_PTWT_6m_vs_WD_PTWT_18_27m[WD_PTWT_6m_vs_WD_PTWT_18_27m$PValue < 0.05,], 
          file= "Results/EdgeR_Pval/WD_PTWT_6m_vs_WD_PTWT_18_27m_Pval.csv",row.names=F)
dim(WD_PTWT_6m_vs_WD_PTWT_18_27m[WD_PTWT_6m_vs_WD_PTWT_18_27m$FDR < 0.05,])
write.csv(WD_PTWT_6m_vs_WD_PTWT_18_27m[WD_PTWT_6m_vs_WD_PTWT_18_27m$FDR < 0.05,], 
          file= "Results/EdgeR_FDR/WD_PTWT_6m_vs_WD_PTWT_18_27m_FDR.csv",row.names=F)
dim(WD_PTWT_6m_vs_WD_PTWT_18_27m[(WD_PTWT_6m_vs_WD_PTWT_18_27m$FDR < 0.05) & 
                                   ((WD_PTWT_6m_vs_WD_PTWT_18_27m$logFC > 1)|
                                      (WD_PTWT_6m_vs_WD_PTWT_18_27m$logFC < -1)),])
write.csv(WD_PTWT_6m_vs_WD_PTWT_18_27m[(WD_PTWT_6m_vs_WD_PTWT_18_27m$FDR < 0.05) & 
                                         ((WD_PTWT_6m_vs_WD_PTWT_18_27m$logFC > 1)|
                                            (WD_PTWT_6m_vs_WD_PTWT_18_27m$logFC < -1)),],
          file= "Results/EdgeR_2fold/WD_PTWT_6m_vs_WD_PTWT_18_27m_2fold.csv",row.names=F)
write.csv(WD_PTWT_6m_vs_WD_PTWT_18_27m, file= "Results/EdgeR_ALL/WD_PTWT_6m_vs_WD_PTWT_18_27m_all.csv",row.names=F)

################# D ######################
# d1 CD_PTWT_6m vs WD_PTWT_6m

lrt.CD_PTWT_6m_vs_WD_PTWT_6m= glmLRT(fit,contrast=makeContrasts(groupsCD_PTWT_6m-groupsWD_PTWT_6m, levels=design))
CD_PTWT_6m_vs_WD_PTWT_6m<- as.data.frame(topTags(lrt.CD_PTWT_6m_vs_WD_PTWT_6m, n=Inf))
#add info
CD_PTWT_6m_vs_WD_PTWT_6m<-left_join(CD_PTWT_6m_vs_WD_PTWT_6m,mouse_mart_export,by="genes")
#
dim(CD_PTWT_6m_vs_WD_PTWT_6m[CD_PTWT_6m_vs_WD_PTWT_6m$PValue < 0.05,])
write.csv(CD_PTWT_6m_vs_WD_PTWT_6m[CD_PTWT_6m_vs_WD_PTWT_6m$PValue < 0.05,], 
          file= "Results/EdgeR_Pval/CD_PTWT_6m_vs_WD_PTWT_6m_Pval.csv",row.names=F)
dim(CD_PTWT_6m_vs_WD_PTWT_6m[CD_PTWT_6m_vs_WD_PTWT_6m$FDR < 0.05,])
write.csv(CD_PTWT_6m_vs_WD_PTWT_6m[CD_PTWT_6m_vs_WD_PTWT_6m$FDR < 0.05,], 
          file= "Results/EdgeR_FDR/CD_PTWT_6m_vs_WD_PTWT_6m_FDR.csv",row.names=F)
dim(CD_PTWT_6m_vs_WD_PTWT_6m[(CD_PTWT_6m_vs_WD_PTWT_6m$FDR < 0.05) & 
                                ((CD_PTWT_6m_vs_WD_PTWT_6m$logFC > 1)|
                                   (CD_PTWT_6m_vs_WD_PTWT_6m$logFC < -1)),])
write.csv(CD_PTWT_6m_vs_WD_PTWT_6m[(CD_PTWT_6m_vs_WD_PTWT_6m$FDR < 0.05) & 
                                      ((CD_PTWT_6m_vs_WD_PTWT_6m$logFC > 1)|
                                         (CD_PTWT_6m_vs_WD_PTWT_6m$logFC < -1)),],
          file= "Results/EdgeR_2fold/CD_PTWT_6m_vs_WD_PTWT_6m_2fold.csv",row.names=F)
write.csv(CD_PTWT_6m_vs_WD_PTWT_6m, file= "Results/EdgeR_ALL/CD_PTWT_6m_vs_WD_PTWT_6m_all.csv",row.names=F)
dim(CD_PTWT_6m_vs_WD_PTWT_6m[(CD_PTWT_6m_vs_WD_PTWT_6m$FDR < 0.1) & 
                               ((CD_PTWT_6m_vs_WD_PTWT_6m$logFC > 1)|
                                  (CD_PTWT_6m_vs_WD_PTWT_6m$logFC < -1)),])
write.csv(CD_PTWT_6m_vs_WD_PTWT_6m[(CD_PTWT_6m_vs_WD_PTWT_6m$FDR < 0.1) & 
                                     ((CD_PTWT_6m_vs_WD_PTWT_6m$logFC > 1)|
                                        (CD_PTWT_6m_vs_WD_PTWT_6m$logFC < -1)),],
          file= "Results/EdgeR_0.1/CD_PTWT_6m_vs_WD_PTWT_6m_2fold.csv",row.names=F)



##
lrt.WD_PTWT_6m_vs_CD_PTWT_6m= glmLRT(fit,contrast=makeContrasts(groupsWD_PTWT_6m-groupsCD_PTWT_6m, levels=design))
WD_PTWT_6m_vs_CD_PTWT_6m<- as.data.frame(topTags(lrt.WD_PTWT_6m_vs_CD_PTWT_6m, n=Inf))
dim(WD_PTWT_6m_vs_CD_PTWT_6m[WD_PTWT_6m_vs_CD_PTWT_6m$PValue < 0.05,])
write.csv(WD_PTWT_6m_vs_CD_PTWT_6m[WD_PTWT_6m_vs_CD_PTWT_6m$PValue < 0.05,], 
          file= "Results_2/EdgeR_Pval/WD_PTWT_6m_vs_CD_PTWT_6m_Pval.csv",row.names=F)
dim(WD_PTWT_6m_vs_CD_PTWT_6m[WD_PTWT_6m_vs_CD_PTWT_6m$FDR < 0.05,])
write.csv(WD_PTWT_6m_vs_CD_PTWT_6m[WD_PTWT_6m_vs_CD_PTWT_6m$FDR < 0.05,], 
          file= "Results_2/EdgeR_FDR/WD_PTWT_6m_vs_CD_PTWT_6m_FDR.csv",row.names=F)
dim(WD_PTWT_6m_vs_CD_PTWT_6m[(WD_PTWT_6m_vs_CD_PTWT_6m$FDR < 0.05) & 
                                ((WD_PTWT_6m_vs_CD_PTWT_6m$logFC > 1)|
                                   (WD_PTWT_6m_vs_CD_PTWT_6m$logFC < -1)),])
write.csv(WD_PTWT_6m_vs_CD_PTWT_6m[(WD_PTWT_6m_vs_CD_PTWT_6m$FDR < 0.05) & 
                                      ((WD_PTWT_6m_vs_CD_PTWT_6m$logFC > 1)|
                                         (WD_PTWT_6m_vs_CD_PTWT_6m$logFC < -1)),],
          file= "Results_2/EdgeR_2fold/WD_PTWT_6m_vs_CD_PTWT_6m_2fold.csv",row.names=F)
write.csv(WD_PTWT_6m_vs_CD_PTWT_6m, file= "Results_2/EdgeR_ALL/WD_PTWT_6m_vs_CD_PTWT_6m_all.csv",row.names=F)


# d2. CD_PTWT_16m vs WD_PTWT_16m

lrt.CD_PTWT_16m_vs_WD_PTWT_16m= glmLRT(fit,contrast=makeContrasts(groupsWD_PTWT_18_27m-groupsWD_PTWT_16m, levels=design))
CD_PTWT_16m_vs_WD_PTWT_16m<- as.data.frame(topTags(lrt.CD_PTWT_16m_vs_WD_PTWT_16m, n=Inf))
#add info
CD_PTWT_16m_vs_WD_PTWT_16m<-left_join(CD_PTWT_16m_vs_WD_PTWT_16m,mouse_mart_export,by="genes")
#
dim(CD_PTWT_16m_vs_WD_PTWT_16m[CD_PTWT_16m_vs_WD_PTWT_16m$PValue < 0.05,])
write.csv(CD_PTWT_16m_vs_WD_PTWT_16m[CD_PTWT_16m_vs_WD_PTWT_16m$PValue < 0.05,], 
          file= "Results/EdgeR_Pval/CD_PTWT_16m_vs_WD_PTWT_16m_Pval.csv",row.names=F)
dim(CD_PTWT_16m_vs_WD_PTWT_16m[CD_PTWT_16m_vs_WD_PTWT_16m$FDR < 0.05,])
write.csv(CD_PTWT_16m_vs_WD_PTWT_16m[CD_PTWT_16m_vs_WD_PTWT_16m$FDR < 0.05,], 
          file= "Results/EdgeR_FDR/CD_PTWT_16m_vs_WD_PTWT_16m_FDR.csv",row.names=F)
dim(CD_PTWT_16m_vs_WD_PTWT_16m[(CD_PTWT_16m_vs_WD_PTWT_16m$FDR < 0.05) & 
                                   ((CD_PTWT_16m_vs_WD_PTWT_16m$logFC > 1)|
                                      (CD_PTWT_16m_vs_WD_PTWT_16m$logFC < -1)),])
write.csv(CD_PTWT_16m_vs_WD_PTWT_16m[(CD_PTWT_16m_vs_WD_PTWT_16m$FDR < 0.05) & 
                                         ((CD_PTWT_16m_vs_WD_PTWT_16m$logFC > 1)|
                                            (CD_PTWT_16m_vs_WD_PTWT_16m$logFC < -1)),],
          file= "Results/EdgeR_2fold/CD_PTWT_16m_vs_WD_PTWT_16m_2fold.csv",row.names=F)
write.csv(CD_PTWT_16m_vs_WD_PTWT_16m, file= "Results/EdgeR_ALL/CD_PTWT_16m_vs_WD_PTWT_16m_all.csv",row.names=F)
dim(CD_PTWT_16m_vs_WD_PTWT_16m[(CD_PTWT_16m_vs_WD_PTWT_16m$FDR < 0.1) & 
                                 ((CD_PTWT_16m_vs_WD_PTWT_16m$logFC > 1)|
                                    (CD_PTWT_16m_vs_WD_PTWT_16m$logFC < -1)),])
write.csv(CD_PTWT_16m_vs_WD_PTWT_16m[(CD_PTWT_16m_vs_WD_PTWT_16m$FDR < 0.1) & 
                                       ((CD_PTWT_16m_vs_WD_PTWT_16m$logFC > 1)|
                                          (CD_PTWT_16m_vs_WD_PTWT_16m$logFC < -1)),],
          file= "Results/EdgeR_0.1/CD_PTWT_16m_vs_WD_PTWT_16m_0.1FDR.csv",row.names=F)

##
lrt.WD_PTWT_16m_vs_CD_PTWT_16m= glmLRT(fit,contrast=makeContrasts(groupsWD_PTWT_16m-groupsCD_PTWT_16m, levels=design))
WD_PTWT_16m_vs_CD_PTWT_16m<- as.data.frame(topTags(lrt.WD_PTWT_16m_vs_CD_PTWT_16m, n=Inf))
dim(WD_PTWT_16m_vs_CD_PTWT_16m[WD_PTWT_16m_vs_CD_PTWT_16m$PValue < 0.05,])
write.csv(WD_PTWT_16m_vs_CD_PTWT_16m[WD_PTWT_16m_vs_CD_PTWT_16m$PValue < 0.05,], 
          file= "Results_2/EdgeR_Pval/WD_PTWT_16m_vs_CD_PTWT_16m_Pval.csv",row.names=F)
dim(WD_PTWT_16m_vs_CD_PTWT_16m[WD_PTWT_16m_vs_CD_PTWT_16m$FDR < 0.05,])
write.csv(WD_PTWT_16m_vs_CD_PTWT_16m[WD_PTWT_16m_vs_CD_PTWT_16m$FDR < 0.05,], 
          file= "Results_2/EdgeR_FDR/WD_PTWT_16m_vs_CD_PTWT_16m_FDR.csv",row.names=F)
dim(WD_PTWT_16m_vs_CD_PTWT_16m[(WD_PTWT_16m_vs_CD_PTWT_16m$FDR < 0.05) & 
                                   ((WD_PTWT_16m_vs_CD_PTWT_16m$logFC > 1)|
                                      (WD_PTWT_16m_vs_CD_PTWT_16m$logFC < -1)),])
write.csv(WD_PTWT_16m_vs_CD_PTWT_16m[(WD_PTWT_16m_vs_CD_PTWT_16m$FDR < 0.05) & 
                                         ((WD_PTWT_16m_vs_CD_PTWT_16m$logFC > 1)|
                                            (WD_PTWT_16m_vs_CD_PTWT_16m$logFC < -1)),],
          file= "Results_2/EdgeR_2fold/WD_PTWT_16m_vs_CD_PTWT_16m_2fold.csv",row.names=F)
write.csv(WD_PTWT_16m_vs_CD_PTWT_16m, file= "Results_2/EdgeR_ALL/WD_PTWT_16m_vs_CD_PTWT_16m_all.csv",row.names=F)


# d3. CD_PTWT_18_27m vs WD_PTWT_18_27m

lrt.CD_PTWT_18_27m_vs_WD_PTWT_18_27m= glmLRT(fit,contrast=makeContrasts(groupsCD_PTWT_18_27m-groupsWD_PTWT_18_27m, levels=design))
CD_PTWT_18_27m_vs_WD_PTWT_18_27m<- as.data.frame(topTags(lrt.CD_PTWT_18_27m_vs_WD_PTWT_18_27m, n=Inf))
#add info
CD_PTWT_18_27m_vs_WD_PTWT_18_27m<-left_join(CD_PTWT_18_27m_vs_WD_PTWT_18_27m,mouse_mart_export,by="genes")
#
dim(CD_PTWT_18_27m_vs_WD_PTWT_18_27m[CD_PTWT_18_27m_vs_WD_PTWT_18_27m$PValue < 0.05,])
write.csv(CD_PTWT_18_27m_vs_WD_PTWT_18_27m[CD_PTWT_18_27m_vs_WD_PTWT_18_27m$PValue < 0.05,], 
          file= "Results/EdgeR_Pval/CD_PTWT_18_27m_vs_WD_PTWT_18_27m_Pval.csv",row.names=F)
dim(CD_PTWT_18_27m_vs_WD_PTWT_18_27m[CD_PTWT_18_27m_vs_WD_PTWT_18_27m$FDR < 0.05,])
write.csv(CD_PTWT_18_27m_vs_WD_PTWT_18_27m[CD_PTWT_18_27m_vs_WD_PTWT_18_27m$FDR < 0.05,], 
          file= "Results/EdgeR_FDR/CD_PTWT_18_27m_vs_WD_PTWT_18_27m_FDR.csv",row.names=F)
dim(CD_PTWT_18_27m_vs_WD_PTWT_18_27m[(CD_PTWT_18_27m_vs_WD_PTWT_18_27m$FDR < 0.05) & 
                                 ((CD_PTWT_18_27m_vs_WD_PTWT_18_27m$logFC > 1)|
                                    (CD_PTWT_18_27m_vs_WD_PTWT_18_27m$logFC < -1)),])
write.csv(CD_PTWT_18_27m_vs_WD_PTWT_18_27m[(CD_PTWT_18_27m_vs_WD_PTWT_18_27m$FDR < 0.05) & 
                                       ((CD_PTWT_18_27m_vs_WD_PTWT_18_27m$logFC > 1)|
                                          (CD_PTWT_18_27m_vs_WD_PTWT_18_27m$logFC < -1)),],
          file= "Results/EdgeR_2fold/CD_PTWT_18_27m_vs_WD_PTWT_18_27m_2fold.csv",row.names=F)
write.csv(CD_PTWT_18_27m_vs_WD_PTWT_18_27m, file= "Results/EdgeR_ALL/CD_PTWT_18_27m_vs_WD_PTWT_18_27m_all.csv",row.names=F)
dim(CD_PTWT_18_27m_vs_WD_PTWT_18_27m[(CD_PTWT_18_27m_vs_WD_PTWT_18_27m$FDR < 0.1) & 
                                       ((CD_PTWT_18_27m_vs_WD_PTWT_18_27m$logFC > 1)|
                                          (CD_PTWT_18_27m_vs_WD_PTWT_18_27m$logFC < -1)),])
write.csv(CD_PTWT_18_27m_vs_WD_PTWT_18_27m[(CD_PTWT_18_27m_vs_WD_PTWT_18_27m$FDR < 0.1) & 
                                             ((CD_PTWT_18_27m_vs_WD_PTWT_18_27m$logFC > 1)|
                                                (CD_PTWT_18_27m_vs_WD_PTWT_18_27m$logFC < -1)),],
          file= "Results/EdgeR_0.1/CD_PTWT_18_27m_vs_WD_PTWT_18_27m_0.1FDR.csv",row.names=F)


##
lrt.WD_PTWT_18_27m_vs_CD_PTWT_18_27m= glmLRT(fit,contrast=makeContrasts(groupsWD_PTWT_18_27m-groupsCD_PTWT_18_27m, levels=design))
WD_PTWT_18_27m_vs_CD_PTWT_18_27m<- as.data.frame(topTags(lrt.WD_PTWT_18_27m_vs_CD_PTWT_18_27m, n=Inf))
dim(WD_PTWT_18_27m_vs_CD_PTWT_18_27m[WD_PTWT_18_27m_vs_CD_PTWT_18_27m$PValue < 0.05,])
write.csv(WD_PTWT_18_27m_vs_CD_PTWT_18_27m[WD_PTWT_18_27m_vs_CD_PTWT_18_27m$PValue < 0.05,], 
          file= "Results_2/EdgeR_Pval/WD_PTWT_18_27m_vs_CD_PTWT_18_27m_Pval.csv",row.names=F)
dim(WD_PTWT_18_27m_vs_CD_PTWT_18_27m[WD_PTWT_18_27m_vs_CD_PTWT_18_27m$FDR < 0.05,])
write.csv(WD_PTWT_18_27m_vs_CD_PTWT_18_27m[WD_PTWT_18_27m_vs_CD_PTWT_18_27m$FDR < 0.05,], 
          file= "Results_2/EdgeR_FDR/WD_PTWT_18_27m_vs_CD_PTWT_18_27m_FDR.csv",row.names=F)
dim(WD_PTWT_18_27m_vs_CD_PTWT_18_27m[(WD_PTWT_18_27m_vs_CD_PTWT_18_27m$FDR < 0.05) & 
                                 ((WD_PTWT_18_27m_vs_CD_PTWT_18_27m$logFC > 1)|
                                    (WD_PTWT_18_27m_vs_CD_PTWT_18_27m$logFC < -1)),])
write.csv(WD_PTWT_18_27m_vs_CD_PTWT_18_27m[(WD_PTWT_18_27m_vs_CD_PTWT_18_27m$FDR < 0.05) & 
                                       ((WD_PTWT_18_27m_vs_CD_PTWT_18_27m$logFC > 1)|
                                          (WD_PTWT_18_27m_vs_CD_PTWT_18_27m$logFC < -1)),],
          file= "Results_2/EdgeR_2fold/WD_PTWT_18_27m_vs_CD_PTWT_18_27m_2fold.csv",row.names=F)
write.csv(WD_PTWT_18_27m_vs_CD_PTWT_18_27m, file= "Results_2/EdgeR_ALL/WD_PTWT_18_27m_vs_CD_PTWT_18_27m_all.csv",row.names=F)


########END of main pipeline########

##
#extract groupwise columns
dim(norm_zero_combined_data_tb)
colnames(norm_zero_combined_data_tb)

KO_columns <- norm_zero_combined_data_tb %>% select(genes, 
                                                    CD.PTKO.1, CD.PTKO.2, CD.PTKO.3, CD.PTKO.4, CD.PTKO.5,
                                                    WD.PTKO.1, WD.PTKO.3, WD.PTKO.4, WD.PTKO.5, WD.PTKO.7)
                                                    #Gene_name, Gene_description)

kegg_gene_list<-read.csv("Results/Kegg_gene_list.txt",sep="\t")
filtered_KO_columns<-filter(KO_columns, genes %in% kegg_gene_list$Ensembl_IDS)

#write.csv(KO_columns,file="Results/KO_columns.csv",row.names = FALSE)
#write.csv(filtered_KO_columns,file="Results/filtered_KO_columns.csv",row.names = FALSE)


###### clusterprofile####
library(clusterProfiler)

dim(WD_PTWT_18_27m_vs_WD_PTKO_18_27m[(WD_PTWT_18_27m_vs_WD_PTKO_18_27m$PValue < 0.05) & 
                                       ((WD_PTWT_18_27m_vs_WD_PTKO_18_27m$logFC > 0.58)|
                                          (WD_PTWT_18_27m_vs_WD_PTKO_18_27m$logFC < -0.58)),])

sigGenes <- WD_PTWT_18_27m_vs_WD_PTKO_18_27m[(WD_PTWT_18_27m_vs_WD_PTKO_18_27m$PValue < 0.05) & 
                                               ((WD_PTWT_18_27m_vs_WD_PTKO_18_27m$logFC > 0.58)|
                                                  (WD_PTWT_18_27m_vs_WD_PTKO_18_27m$logFC < -0.58)),]
## switch to entrez id? 
kk <- enrichKEGG(gene = sigGenes, organism = 'mmu') ## wil not work - no entrez ids

ensembl_entrez_mouse_IDs<-read.csv("ensembl_entrez_mouse_IDs.txt",sep="\t")

sigGenes_entrez<-left_join(sigGenes,ensembl_entrez_mouse_IDs,by="genes") #<- adds entrez id to table with significant genes

sigGenes_entrez_onlyids <- na.exclude(sigGenes_entrez$Entrezgene_ID)
head(sigGenes_entrez)
kk <- enrichKEGG(gene = sigGenes_entrez_onlyids, organism = "mmu")

WD_PTWT_18_27m_vs_WD_PTKO_18_27m_KEGG_PTS<-head(kk, n=50)
View(WD_PTWT_18_27m_vs_WD_PTKO_18_27m_KEGG_PTS)

browseKEGG(kk, 'mmu05020')

##show and save pathway

library(pathview)

logFC <- sigGenes_entrez$logFC #<- uses sigGenes_entrez table
names(logFC) <- sigGenes_entrez$Entrezgene_ID
pathview(gene.data = logFC, 
         pathway.id = "mmu05020", 
         species = "mmu", 
         limit = list(gene=5, cpd=1))

?pathview()

#######
#extract logfc for plot
colnames(WD_PTWT_18_27m_vs_WD_PTKO_18_27m)

KO_genes_Metabolic_pathways<-c("ENSMUSG00000064345","ENSMUSG00000024843","ENSMUSG00000034639","ENSMUSG00000064367","ENSMUSG00000064351","ENSMUSG00000064370","ENSMUSG00000064357","ENSMUSG00000006345","ENSMUSG00000064341","ENSMUSG00000021957","ENSMUSG00000028307","ENSMUSG00000020456","ENSMUSG00000024644","ENSMUSG00000021577","ENSMUSG00000003477","ENSMUSG00000018770","ENSMUSG00000029455","ENSMUSG00000025479","ENSMUSG00000031818","ENSMUSG00000069805","ENSMUSG00000032350","ENSMUSG00000032527","ENSMUSG00000027227","ENSMUSG00000023456","ENSMUSG00000025428","ENSMUSG00000026473","ENSMUSG00000064354","ENSMUSG00000035637","ENSMUSG00000029063","ENSMUSG00000042096","ENSMUSG00000034707","ENSMUSG00000030541","ENSMUSG00000030246","ENSMUSG00000026170","ENSMUSG00000068874","ENSMUSG00000023120","ENSMUSG00000031231","ENSMUSG00000074207","ENSMUSG00000022994","ENSMUSG00000027406","ENSMUSG00000026617","ENSMUSG00000028961","ENSMUSG00000030268","ENSMUSG00000037260","ENSMUSG00000020774","ENSMUSG00000061838","ENSMUSG00000003865","ENSMUSG00000032373","ENSMUSG00000002010","ENSMUSG00000032047","ENSMUSG00000027665","ENSMUSG00000004880","ENSMUSG00000000088","ENSMUSG00000030884","ENSMUSG00000032437","ENSMUSG00000038296","ENSMUSG00000024899","ENSMUSG00000059534","ENSMUSG00000028124","ENSMUSG00000027332","ENSMUSG00000027870","ENSMUSG00000062908","ENSMUSG00000028032","ENSMUSG00000024726","ENSMUSG00000025176","ENSMUSG00000029632","ENSMUSG00000036427","ENSMUSG00000028405","ENSMUSG00000066097","ENSMUSG00000064363","ENSMUSG00000039648","ENSMUSG00000020988","ENSMUSG00000027984","ENSMUSG00000025041","ENSMUSG00000089678","ENSMUSG00000007038","ENSMUSG00000001467","ENSMUSG00000025465","ENSMUSG00000075232","ENSMUSG00000035936","ENSMUSG00000029680","ENSMUSG00000076441","ENSMUSG00000011179","ENSMUSG00000022683","ENSMUSG00000038690","ENSMUSG00000023827","ENSMUSG00000021646","ENSMUSG00000039105","ENSMUSG00000031059","ENSMUSG00000052520","ENSMUSG00000000171","ENSMUSG00000020051")
KO_genes_Oxidative_phosphorylation<-c("ENSMUSG00000064345","ENSMUSG00000064367","ENSMUSG00000064351","ENSMUSG00000064370","ENSMUSG00000064357","ENSMUSG00000064341","ENSMUSG00000021577","ENSMUSG00000018770","ENSMUSG00000031818","ENSMUSG00000025428","ENSMUSG00000064354","ENSMUSG00000031231","ENSMUSG00000000088","ENSMUSG00000030884","ENSMUSG00000059534","ENSMUSG00000029632","ENSMUSG00000064363","ENSMUSG00000038690","ENSMUSG00000039105","ENSMUSG00000031059","ENSMUSG00000000171")
KO_genes_Carbon_metabolism<-c("ENSMUSG00000021957","ENSMUSG00000028307","ENSMUSG00000020456","ENSMUSG00000021577","ENSMUSG00000069805","ENSMUSG00000032527","ENSMUSG00000023456","ENSMUSG00000030541","ENSMUSG00000027406","ENSMUSG00000028961","ENSMUSG00000061838","ENSMUSG00000002010","ENSMUSG00000032047","ENSMUSG00000027870","ENSMUSG00000036427","ENSMUSG00000028405","ENSMUSG00000025465","ENSMUSG00000000171")
KO_genes_Thermogenesis<-c("ENSMUSG00000064345","ENSMUSG00000064367","ENSMUSG00000064351","ENSMUSG00000064370","ENSMUSG00000064357","ENSMUSG00000064341","ENSMUSG00000021577","ENSMUSG00000018770","ENSMUSG00000031818","ENSMUSG00000025428","ENSMUSG00000064354","ENSMUSG00000031231","ENSMUSG00000022994","ENSMUSG00000000088","ENSMUSG00000030884","ENSMUSG00000059534","ENSMUSG00000029632","ENSMUSG00000064363","ENSMUSG00000062825","ENSMUSG00000038690","ENSMUSG00000031059","ENSMUSG00000000171","ENSMUSG00000005469")

KO_genes_Citrate_cycle<-c("ENSMUSG00000061787","ENSMUSG00000032518","ENSMUSG00000062647","ENSMUSG00000032399","ENSMUSG00000059291","ENSMUSG00000028861","ENSMUSG00000008682","ENSMUSG00000029614","ENSMUSG00000047215","ENSMUSG00000090862","ENSMUSG00000025290","ENSMUSG00000047675","ENSMUSG00000067274","ENSMUSG00000007892","ENSMUSG00000012848","ENSMUSG00000063457","ENSMUSG00000046330")
KO_genes_Citrate_cycle<-c("ENSMUSG00000029455","ENSMUSG00000032527","ENSMUSG00000030268","ENSMUSG00000032047","ENSMUSG00000027332","ENSMUSG00000062908","ENSMUSG00000027984","ENSMUSG00000089678","ENSMUSG00000025465","ENSMUSG00000021646")


filtered_KO_columns<-filter(WD_PTWT_18_27m_vs_WD_PTKO_18_27m, genes %in% KO_genes_Citrate_cycle)

write.table(filtered_KO_columns,
            "filtered_KO_columns/Citrate_cycle_FC.txt",
            sep = "\t",
            row.names = FALSE)

colnames(filtered_KO_columns)

# df <- data.frame(gene = filtered_KO_columns$Gene_name, 
#                  FoldChange = filtered_KO_columns$logFC, 
#                  pvalue = filtered_KO_columns$PValue)
# ggplot(data=df, aes(x=gene, y=FoldChange, fill = pvalue)) +
#   geom_bar(stat="identity") + coord_flip()

df <- data.frame(gene = filtered_KO_columns$Gene_name, 
                 FoldChange = filtered_KO_columns$logFC)

plot_1<-ggplot(data=df, aes(x=gene, y=FoldChange, fill = FoldChange)) +
  geom_bar(stat="identity") +
  coord_flip() + 
  ggtitle("KEGG: Citrate cycle") +
  labs(y = "Fold Change", x = "Mouse Gene Names", fill = "Fold Change") +
  # scale_fill_gradientn(colours = c("green", "black", "red"),
  #                       values = scales::rescale(c(-2, -0.5, 0, 0.5, 2)))
  scale_fill_gradientn(colours = c("black", "red"),
                       values = scales::rescale(c(0, 0.5, 2)))


plot_1

# ggsave(
#   "X_1",
#   plot = plot_1,
#   device = jpeg,
#   path = "filtered_KO_columns/",
#   scale = 1,
#   width = 3,
#   height = 4,
#   units = c("in"),
#   dpi = 300,
#   limitsize = TRUE,
#   bg = NULL,
# )

?ggsave()

# plot_1+scale_fill_gradient(low="green", high="red")
# 
# plot_1+scale_fill_gradientn(colours = c("green", "black", "red"),
#                      values = scales::rescale(c(-0.5, -0.05, 0, 0.05, 0.5)))
# 
# plot_1+scale_fill_gradientn(colours = c("green", "black", "red"),
#                             values = scales::rescale(c(-2, -0.25, 0, 0.25, 2)))
# 
# plot_1 + ggtitle("Metabolic_pathways")


# ####convert DGE to table####
#                                                   #-1     0     1   total
# ##Field_B73_SSvsWW (cuffdiff_1)
# table(decideTestsDGE(lrt.Field_B73_SSa_vs_WWa))  #695 25499  3323   4018
# table(decideTestsDGE(lrt.Field_B73_SSb_vs_WWb))  #4194 16660  8663  12857
# table(decideTestsDGE(lrt.Field_B73_SSc_vs_WWc))  #2812 20511  6194  9006
# 
# ##Field_FR697_SSvsWW (cuffiff_4)
# table(decideTestsDGE(lrt.Field_FR697_SSa_vs_WWa)) #971 27082  2064  3035 
# table(decideTestsDGE(lrt.Field_FR697_SSb_vs_WWb)) #3197 20475  5845 9042
# table(decideTestsDGE(lrt.Field_FR697_SSc_vs_WWc)) #2383 20993 6141  8524
# 
# ##Lab_FR697_SSvsWW (cuffiff_8)
# table(decideTestsDGE(lrt.Lab_FR697_SSa_vs_WWa))   #99 29299   119 298
# table(decideTestsDGE(lrt.Lab_FR697_SSb_vs_WWb))   #931  28073 513 1444
# table(decideTestsDGE(lrt.Lab_FR697_SSc_vs_WWc))   #559  28450 508 1067
# 
# ####Write Full table to file_MOD####
# 
# ####volcano plot####
# Lab_FR697_SSc_vs_WWc <- as.data.frame(topTags(lrt.Lab_FR697_SSc_vs_WWc, n=Inf))
# library(EnhancedVolcano)
# 
# colnames(Lab_FR697_SSc_vs_WWc)
# rownames(Lab_FR697_SSc_vs_WWc)
# 
# ##og##dif <- data.frame(fc =fold_changes,pv =pvalues)
# dif <- data.frame(genes = Lab_FR697_SSc_vs_WWc$genes,
#                   logfc = Lab_FR697_SSc_vs_WWc$logFC,
#                   fdr = Lab_FR697_SSc_vs_WWc$FDR)
# dif$threshold <- ifelse(dif$logfc > 1 & dif$fdr < 0.05, "upreg",
#                         ifelse(dif$logfc < -1 & dif$fdr < 0.05,"downreg","NotSign" ))
# head(dif)
# 
# vol_plot_2<-ggplot(data=dif, aes(x=logfc, y=-log10(fdr))) +
#   geom_point( size=1 ,aes(color=as.factor(threshold))) +
#   theme(legend.position = "none") +
#   xlim(c(-10, 10)) + ylim(c(0, 15)) +
#   xlab("log2 fold change") + ylab("-log10 p-value")  + theme_bw()+
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),)
# 
# vol_plot_2
# vol_plot<-ggplot(Lab_FR697_SSc_vs_WWc, aes(x = logFC, y=-log10(FDR), 
#                                            col=(FDR < 0.05 & ((logFC > 1)|(logFC < -1))))) + geom_point()
# vol_plot
# 
# #ggplot_build(vol_plot)
# #ggplot_build(vol_plot)$data[[1]]$x
# #ggplot_build(vol_plot)$data[[1]]$y
# ggplot_build(vol_plot_2)$data[[1]]$x
# ggplot_build(vol_plot_2)$data[[1]]$y
# 
# 
# volc_plot_XY<-cbind(vol_plot[["data"]][["genes"]],ggplot_build(vol_plot)$data[[1]]$x,ggplot_build(vol_plot)$data[[1]]$y,dif$threshold)
# dim(volc_plot_XY)
# head(volc_plot_XY)
# 
# 
# volcano_fxn <- function(inputFile){
#   #add classificationn for significant genes
#   # dif <- data.frame(genes = Lab_FR697_SSc_vs_WWc$genes,
#   #                   logfc = Lab_FR697_SSc_vs_WWc$logFC,
#   #                   fdr = -log10(Lab_FR697_SSc_vs_WWc$FDR))
#   
#   # dif$threshold <- ifelse(dif$logfc > 1 & dif$fdr < 0.05, "upreg",
#   #                         ifelse(dif$logfc < -1 & dif$fdr < 0.05,"downreg","NotSign" ))
#   
#   # vol_plot<-ggplot(Lab_FR697_SSc_vs_WWc, aes(x = logFC, y=-log10(FDR), 
#   #                                            col=(FDR < 0.05 & ((logFC > 1)|(logFC < -1))))) + geom_point()
#   
#   volc_plot_XY<-data.frame(inputFile$genes,
#                            volc.X = inputFile$logFC,
#                            volc.y = -log10(inputFile$FDR),
#                            threshold = ifelse(inputFile$logfc > 1 & inputFile$fdr < 0.05, "upreg",
#                                               ifelse(inputFile$logfc < -1 & inputFile$fdr < 0.05,"downreg","NotSign" ))
#                            )
#   return(volc_plot_XY)
#   
# }
# 
# ######
# 
# # EnhancedVolcano(Lab_FR697_SSc_vs_WWc,
# #                 lab = rownames(Lab_FR697_SSc_vs_WWc),
# #                 x = 'logFC',
# #                 y = 'FDR',
# #                 pCutoff = 0.05,
# #                 FCcutoff = 1,
# #                 # ylab = bquote("FDR"),
# #                 xlim = c(-15,15))
# 
# ##ggplot volcano plot
# #ggplot(Lab_FR697_SSc_vs_WWc, aes(x = logFC, y=-log10(FDR), col=FDR < 0.05)) + geom_point()
# 
# 
# # [(Lab_FR697_SSc_vs_WWc$FDR < 0.05) & 
# #     ((Lab_FR697_SSc_vs_WWc$logFC > 1)|
# #        (Lab_FR697_SSc_vs_WWc$logFC < -1)),]
# 
# 
# #-log10(0.05)
# 
# ######
# ######filter counts tible by fold2 dataframe gene list
# ##counts_Field_B73_SSa_vs_WWa (tibble)
# ##fold2_Field_B73_SSa_vs_WWa (data.frame)
# 
# #colnames(fold2_Field_B73_SSa_vs_WWa)
# 
# library(tidyverse)
#  
# tb_normalized_data <- as_tibble(normalized_data, rownames = "rowname")
# rownames_to_column(.data, var = "rowname")
# 
# ## SELECT COLUMNS SPECIFIC TO DE COMPARISON
# #Field_B73
# #1
# counts_Field_B73_SSa_vs_WWa = select(tb_normalized_data,rowname,
#                                      S1031_2017:S6031_2017,s1041_2018,s2041_2018,s3011_2018,s4031_2018,s5031_2018,  #ssa
#                                      S1011_2017:S6011_2017,s1021_2018,s2011_2018,s3041_2018,s4041_2018,s5021_2018)  #wwa
# #2
# counts_Field_B73_SSb_vs_WWb = select(tb_normalized_data,rowname,
#                                      S1032_2017:S6032_2017,s1042_2018,s2042_2018,s3012_2018,s4032_2018,s5032_2018,  #ssb
#                                      S1012_2017:S6012_2017,s1022_2018,s2012_2018,s3042_2018,s4042_2018,s5022_2018)  #wwb
# #3
# counts_Field_B73_SSc_vs_WWc = select(tb_normalized_data,rowname,
#                                      S1033_2017:S6033_2017,s1043_2018,s2043_2018,s3013_2018,s4033_2018,s5033_2018,  #ssb
#                                      S1013_2017:S6013_2017,s1023_2018,s2013_2018,s3043_2018,s4043_2018,s5023_2018)  #wwb
# 
# 
# #Field_FR697
# #4
# counts_Field_FR697_SSa_vs_WWa = select(tb_normalized_data,rowname,
#                                        S1021_2017:S6021_2017,s2021_2018,s3031_2018,s4021_2018,s5011_2018,  #ssa
#                                        S1041_2017:S6041_2017,s1011_2018,s2031_2018,s3021_2018,s4011_2018,s5041_2018)  #wwa
# #5
# counts_Field_FR697_SSb_vs_WWb = select(tb_normalized_data,rowname,
#                                        S1022_2017:S6022_2017,s2022_2018,s3032_2018,s4022_2018,s5012_2018,  #ssb
#                                        S2012_2017:S6042_2017,s1012_2018,s2032_2018,s3022_2018,s4012_2018,s5042_2018)  #wwb
# #6
# counts_Field_FR697_SSc_vs_WWc = select(tb_normalized_data,rowname,
#                                        S1023_2017:S6023_2017,s2023_2018,s3033_2018,s4023_2018,s5013_2018,  #ssb
#                                        S1043_2017:S6043_2017,s1013_2018,s2033_2018,s3023_2018,s4013_2018,s5043_2018)  #wwb
# 
# #Lab_FR697
# #7
# counts_Lab_FR697_SSa_vs_WWa = select(tb_normalized_data,rowname,
#                                      S8021_2017:S8051_2017,  #ssa
#                                      S9011_2017:S9051_2017)  #wwa
# #8
# counts_Lab_FR697_SSb_vs_WWb = select(tb_normalized_data,rowname,
#                                      S8012_2017:S8052_2017,  #ssb
#                                      S9012_2017:S9052_2017)  #wwb
# #9
# counts_Lab_FR697_SSc_vs_WWc = select(tb_normalized_data,rowname,
#                                      S8013_2017:S8053_2017,  #ssb
#                                      S9013_2017:S9053_2017)  #wwb
# ##initialize this fxn first
# counts_list_filter <- function(counts_file,de_file,file_name){
#   fold2_file <- as.data.frame(topTags(de_file, n=Inf))
#   filt_fold2_file <- fold2_file[which((fold2_file$FDR < 0.05) & 
#                                                              ((fold2_file$logFC > 1)|
#                                                                 (fold2_file$logFC < -1))),]
#   filter_fold2_counts_file <- filter(counts_file, rowname %in% filt_fold2_file$genes)
#   
#   write.table(filter_fold2_counts_file,
#               file = sprintf("C:/Users/Sid/OneDrive/rna-seq documents/NSF-drought/rnaseq_counts_combined/Field_lab/counts_gene_filtered/%s.txt",file_name),
#               sep = "\t",
#               row.names = FALSE
#   )
# }
# 
# #counts_file = counts_Field_B73_SSa_vs_WWa
# #de_file = lrt.Field_B73_SSa_vs_WWa
# #file name = filter_fold2_Field_B73_SSa_vs_WWa
# 
# counts_list_filter(counts_Field_B73_SSa_vs_WWa,lrt.Field_B73_SSa_vs_WWa,"filter_fold2_Field_B73_SSa_vs_WWa")
# counts_list_filter(counts_Field_B73_SSb_vs_WWb,lrt.Field_B73_SSb_vs_WWb,"filter_fold2_Field_B73_SSb_vs_WWb")
# counts_list_filter(counts_Field_B73_SSc_vs_WWc,lrt.Field_B73_SSc_vs_WWc,"filter_fold2_Field_B73_SSc_vs_WWc")
# 
# counts_list_filter(counts_Field_FR697_SSa_vs_WWa,lrt.Field_FR697_SSa_vs_WWa,"filter_fold2_Field_FR697_SSa_vs_WWa")
# counts_list_filter(counts_Field_FR697_SSb_vs_WWb,lrt.Field_FR697_SSb_vs_WWb,"filter_fold2_Field_FR697_SSb_vs_WWb")
# counts_list_filter(counts_Field_FR697_SSc_vs_WWc,lrt.Field_FR697_SSc_vs_WWc,"filter_fold2_Field_FR697_SSc_vs_WWc")
# 
# counts_list_filter(counts_Lab_FR697_SSa_vs_WWa,lrt.Lab_FR697_SSa_vs_WWa,"filter_fold2_Lab_FR697_SSa_vs_WWa")
# counts_list_filter(counts_Lab_FR697_SSb_vs_WWb,lrt.Lab_FR697_SSb_vs_WWb,"filter_fold2_Lab_FR697_SSb_vs_WWb")
# counts_list_filter(counts_Lab_FR697_SSc_vs_WWc,lrt.Lab_FR697_SSc_vs_WWc,"filter_fold2_Lab_FR697_SSc_vs_WWc")
# 
# 
# # fold2_Field_B73_SSa_vs_WWa <- Field_B73_SSa_vs_WWa[which((Field_B73_SSa_vs_WWa$FDR < 0.05) & 
# #                                                            ((Field_B73_SSa_vs_WWa$logFC > 1)|
# #                                                               (Field_B73_SSa_vs_WWa$logFC < -1))),]
# # filter_fold2_Field_B73_SSa_vs_WWa <- filter(counts_Field_B73_SSa_vs_WWa, rowname %in% fold2_Field_B73_SSa_vs_WWa$genes)
# # colnames(filter_fold2_Field_B73_SSa_vs_WWa)
# # 
# # write.table(filter_fold2_Field_B73_SSa_vs_WWa,
# #           "C:/Users/Sid/OneDrive/rna-seq documents/NSF-drought/rnaseq_counts_combined/Field_lab/counts_gene_filtered/filter_fold2_Field_B73_SSa_vs_WWa.txt",
# #           sep = "\t",
# #           row.names = FALSE
# #           )
# 
# ####end####