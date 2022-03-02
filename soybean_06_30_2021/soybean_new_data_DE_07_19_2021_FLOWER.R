####
#Author: Sidharth Sen Ph.D.

#date: 19 July 2021

#University of Missouri - Division of Plant Sciences and IPG

##data used: soybean leaf data from Ron Mittler's lab

####

library(edgeR)
library(limma)
library(ggplot2)
library(ggfortify)
library(pheatmap)
library(RColorBrewer)
library(dplyr)

##load.image("DE_soybean_07_12_2021.RData")


##load data matrix/ table
ReadCounts_og = read.table("soybean_leaf_raw_counts_7_19_2021.txt",header=T,row.names=1,check.names=FALSE)

colnames(ReadCounts_og)
dim(ReadCounts_og)

##define groups
groups = c(rep("CT",3),
           rep("WD",3),
           rep("HS",3),
           rep("CS",3))

##define sample names
sample_names = c("CT_1","CT_2","CT_3",
                 "WD_1","WD_2","WD_3",
                 "HS_1","HS_2","HS_3",
                 "CS_1","CS_2","CS_3")

##modified for GLM
row.names(ReadCounts_og)
d <- DGEList(counts=ReadCounts_og, group=groups, genes = row.names(ReadCounts_og))
dim(d)
table(groups)
#table(years)


#### section 2 extra filtering to remove genes with very low expression - ####
keep <- rowSums(cpm(d)>1) >= 3 ## this removes any gene who's mean count across all replicates is less than 3
#keep <- rowSums(cpm(d)>1) >= 6
dim(d[keep,])
d_filt <- d[keep,]
dim(d_filt)

#boxplot(d_filt$counts, col="gray", las=3)
#number filtered
dim(d)[1]-dim(d_filt)[1] ##42866 <-- verifies the number of genes dropped from the original counts matrix. i.e. no. of low expressed genes dropped.
#%filtered out
(dim(d)[1]-dim(d_filt)[1])/dim(d)[1] ##0.7668748 ## <-- %age of genes dropped.
colSums(d_filt$counts)
#library size update ----- IMP STEP
d_filt$samples$lib.size <- colSums(d_filt$counts)
### section 2 end ###

#### section 3: normalization of replicates ####
##Get normalized data from TMM via CPM method
#TMM <- calcNormFactors(d, method="TMM") #if no filtering is used.
TMM <- calcNormFactors(d_filt, method="TMM") ##if using extra filtering (i.e. if section 2 is active)
normalized_data <- cpm(TMM)

head(normalized_data)
dim(normalized_data)

##write normalized couunts to file -->
#write.csv(normalized_data,file = "norm_filtered_leaf_counts_07_19_2021.csv")

## section 3 end

####section 4: PLOTS ####
normalized_data_t <- t(normalized_data) #### <- transposes columns and rows, required for plots

normalized_data_t[,1:10]
dim(normalized_data_t) 

###dendogram plot###
d <- dist(normalized_data_t)
hc <- hclust(d)
dhc <- as.dendrogram(hc)
#specific_leaf <- dhc[[1]][[1]][[1]]
#attributes(specific_leaf)
plot(hc)

###adding group info to dataset -- optional - needed for colored plots atm###
normalized_data_dt <- data.frame(normalized_data_t)
normalized_data_groups <- normalized_data_dt
dim(normalized_data_groups)
###add groups (i.e. groups of annotations, not just the "group" of replicates)

normalized_data_groups['groups'] <- groups
dim(normalized_data_groups)
normalized_data_groups['sample_names'] <- sample_names
dim(normalized_data_groups)


##plot pca
autoplot(prcomp(normalized_data_t))

##plot pca with group colors
autoplot(prcomp(normalized_data_groups[1:13031]), data = normalized_data_groups, colour = 'groups' , label = TRUE, label.size = 3,shape = FALSE)
autoplot(prcomp(normalized_data_groups[1:13031]), data = normalized_data_groups, colour = 'sample_names' , label = TRUE, label.size = 3)

##remove sample 5
# no5_normalized_data_t<-as.data.frame(t(select(as.data.frame(normalized_data), -"5_S30")))
# no5_normalized_data_t['groups'] <- groups[-5]
# rownames(no5_normalized_data_t)
# dim(no5_normalized_data_t)
# autoplot(prcomp(no5_normalized_data_t[1:30484]), data = no5_normalized_data_t, colour = 'groups' , label = TRUE, label.size = 5)

##pca plot from raw data
ReadCounts_og_t<-as.data.frame(t(ReadCounts_og))
dim(ReadCounts_og_t)
ReadCounts_og_t['groups'] <- groups
dim(ReadCounts_og_t)
autoplot(prcomp(ReadCounts_og_t[1:55897]), data = ReadCounts_og_t, colour = 'groups' , label = TRUE, label.size = 3)

#save.image("DE_soybean_07_12_2021.RData")

## section 4 end

#### section 5: DGE computation####

summary(TMM$samples$norm.factors)
summary(TMM$samples$lib.size)
TMM$samples
table(groups)

##design matrix
design=model.matrix(~0+groups) #<--main
rownames(design)=colnames(d_filt)
design

#EdgeR DE
TMM
#dge_final <- estimateGLMCommonDisp(d_filt, design, verbose=TRUE) ##this takes some time
dge_final <- estimateGLMCommonDisp(TMM, design, verbose=TRUE) ##this takes some time
#Disp = 0.44291 , BCV = 0.6655

dge_final <- estimateGLMTrendedDisp(dge_final, design,verbose=TRUE)
####results####
# Disp = 1.5136 , BCV = 1.2303 
# Disp = 3.5862 , BCV = 1.8937 
# Disp = 3.43774 , BCV = 1.8541 
# Disp = 3.24789 , BCV = 1.8022 
# Disp = 3.48429 , BCV = 1.8666 
# Disp = 3.32772 , BCV = 1.8242 
# Disp = 2.73995 , BCV = 1.6553 
# Disp = 3.72376 , BCV = 1.9297 
# Disp = 2.76874 , BCV = 1.664 
# Disp = 3.13113 , BCV = 1.7695 
# Disp = 2.66822 , BCV = 1.6335 
# Disp = 2.52812 , BCV = 1.59 
# Disp = 2.03328 , BCV = 1.4259 
# Disp = 2.01439 , BCV = 1.4193 
# Disp = 2.07529 , BCV = 1.4406 
# Disp = 1.95699 , BCV = 1.3989 
# Disp = 1.68478 , BCV = 1.298 
# Disp = 1.63272 , BCV = 1.2778 
# Disp = 1.55174 , BCV = 1.2457 
# Disp = 1.3591 , BCV = 1.1658 
# Disp = 1.28364 , BCV = 1.133 
# Disp = 1.18997 , BCV = 1.0909 
# Disp = 1.08587 , BCV = 1.042 
# Disp = 0.9417 , BCV = 0.9704 
# Disp = 0.95999 , BCV = 0.9798 
# Disp = 0.83458 , BCV = 0.9136 
# Disp = 0.79147 , BCV = 0.8896 
# Disp = 0.69154 , BCV = 0.8316 
# Disp = 0.67933 , BCV = 0.8242 
# Disp = 0.62135 , BCV = 0.7883 
# Disp = 0.53878 , BCV = 0.734 
# Disp = 0.48841 , BCV = 0.6989 
# Disp = 0.49609 , BCV = 0.7043 
# Disp = 0.41017 , BCV = 0.6404 
# Disp = 0.36763 , BCV = 0.6063 
# Disp = 0.35042 , BCV = 0.592 
# Disp = 0.35754 , BCV = 0.5979 
# Disp = 0.31637 , BCV = 0.5625 
# Disp = 0.29067 , BCV = 0.5391 
# Disp = 0.27971 , BCV = 0.5289 
# Disp = 0.24817 , BCV = 0.4982 
# Disp = 0.22915 , BCV = 0.4787 
# Disp = 0.17015 , BCV = 0.4125 
# Disp = 0.15236 , BCV = 0.3903
######

dge_final <- estimateGLMTagwiseDisp(dge_final, design)## this takes time as well
fit <- glmFit(dge_final, design)
fit

####CONTRASTS####

#CT_vs_WD
lrt.leaf_CT_vs_WD= glmLRT(fit,contrast=makeContrasts(groupsCT-groupsWD, levels=design)) #a

#CT_vs_HS
lrt.leaf_CT_vs_HS= glmLRT(fit,contrast=makeContrasts(groupsCT-groupsHS, levels=design)) #b

#CT_vs_CS
lrt.leaf_CT_vs_CS= glmLRT(fit,contrast=makeContrasts(groupsCT-groupsCS, levels=design)) #c

#SAVE SESSION
#save.image(file="soybean_leaf_DE_07_21_2021.RSession")

######leaf_CT_vs_WD##########################################################################################
leaf_CT_vs_WD <- as.data.frame(topTags(lrt.leaf_CT_vs_WD, n=Inf))
dim(leaf_CT_vs_WD[leaf_CT_vs_WD$PValue < 0.05,]) #[1] 8908    6
write.csv(leaf_CT_vs_WD[leaf_CT_vs_WD$PValue < 0.05,], file= "EdgeR_Pval/leaf_CT_vs_WD_Pval.csv",row.names=T)
dim(leaf_CT_vs_WD[leaf_CT_vs_WD$FDR < 0.05,]) #[1] 6236    6
write.csv(leaf_CT_vs_WD[leaf_CT_vs_WD$FDR < 0.05,], file= "EdgeR_FDR/leaf_CT_vs_WD_FDR.csv",row.names=T)
dim(leaf_CT_vs_WD[(leaf_CT_vs_WD$FDR < 0.05) & 
                           ((leaf_CT_vs_WD$logFC > 1)|
                              (leaf_CT_vs_WD$logFC < -1)),]) #[1] 3963    6
write.csv(leaf_CT_vs_WD[(leaf_CT_vs_WD$FDR < 0.05) & 
                                 ((leaf_CT_vs_WD$logFC > 1)|
                                    (leaf_CT_vs_WD$logFC < -1)),], 
          file= "EdgeR_2fold/leaf_CT_vs_WD_2fold.csv",row.names=T)
write.csv(leaf_CT_vs_WD, file= "EdgeR_ALL/leaf_CT_vs_WD_all.csv",row.names=T)

######leaf_CT_vs_HS##########################################################################################
leaf_CT_vs_HS <- as.data.frame(topTags(lrt.leaf_CT_vs_HS, n=Inf))
dim(leaf_CT_vs_HS[leaf_CT_vs_HS$PValue < 0.05,]) #[1] 8908    6
write.csv(leaf_CT_vs_HS[leaf_CT_vs_HS$PValue < 0.05,], file= "EdgeR_Pval/leaf_CT_vs_HS_Pval.csv",row.names=T)
dim(leaf_CT_vs_HS[leaf_CT_vs_HS$FDR < 0.05,]) #[1] 6236    6
write.csv(leaf_CT_vs_HS[leaf_CT_vs_HS$FDR < 0.05,], file= "EdgeR_FDR/leaf_CT_vs_HS_FDR.csv",row.names=T)
dim(leaf_CT_vs_HS[(leaf_CT_vs_HS$FDR < 0.05) & 
                    ((leaf_CT_vs_HS$logFC > 1)|
                       (leaf_CT_vs_HS$logFC < -1)),]) #[1] 3963    6
write.csv(leaf_CT_vs_HS[(leaf_CT_vs_HS$FDR < 0.05) & 
                          ((leaf_CT_vs_HS$logFC > 1)|
                             (leaf_CT_vs_HS$logFC < -1)),], 
          file= "EdgeR_2fold/leaf_CT_vs_HS_2fold.csv",row.names=T)
write.csv(leaf_CT_vs_HS, file= "EdgeR_ALL/leaf_CT_vs_HS_all.csv",row.names=T)

######leaf_CT_vs_CS##########################################################################################
leaf_CT_vs_CS <- as.data.frame(topTags(lrt.leaf_CT_vs_CS, n=Inf))
dim(leaf_CT_vs_CS[leaf_CT_vs_CS$PValue < 0.05,]) #[1] 8908    6
write.csv(leaf_CT_vs_CS[leaf_CT_vs_CS$PValue < 0.05,], file= "EdgeR_Pval/leaf_CT_vs_CS_Pval.csv",row.names=T)
dim(leaf_CT_vs_CS[leaf_CT_vs_CS$FDR < 0.05,]) #[1] 6236    6
write.csv(leaf_CT_vs_CS[leaf_CT_vs_CS$FDR < 0.05,], file= "EdgeR_FDR/leaf_CT_vs_CS_FDR.csv",row.names=T)
dim(leaf_CT_vs_CS[(leaf_CT_vs_CS$FDR < 0.05) & 
                    ((leaf_CT_vs_CS$logFC > 1)|
                       (leaf_CT_vs_CS$logFC < -1)),]) #[1] 3963    6
write.csv(leaf_CT_vs_CS[(leaf_CT_vs_CS$FDR < 0.05) & 
                          ((leaf_CT_vs_CS$logFC > 1)|
                             (leaf_CT_vs_CS$logFC < -1)),], 
          file= "EdgeR_2fold/leaf_CT_vs_CS_2fold.csv",row.names=T)
write.csv(leaf_CT_vs_CS, file= "EdgeR_ALL/leaf_CT_vs_CS_all.csv",row.names=T)

## section 5 end

#### section 6: Figures####

##volcano plots##
