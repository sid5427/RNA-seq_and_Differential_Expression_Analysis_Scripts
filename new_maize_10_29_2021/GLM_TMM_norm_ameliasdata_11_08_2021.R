library(edgeR)
library(limma)
library(ggplot2)
#packageVersion('ggplot2')
library(ggfortify)
#install.packages("pheatmap")
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(dendextend)
library(rgl)
#install.packages("pca3d")
library( pca3d )

#setwd("C:/Users/Sid/OneDrive/post_doc_research/amelia_data/DE_analysis_11_07_2021")
#pwd

####load workspace image####

load('All_DE_11_09_2021.RData')

######

##load data matrix/ table
#ReadCounts = read.table("combined_rna_seq_counts.txt",header=T,row.names=1)
ReadCounts = read.table("rna_seq_counts_11_7_2021.txt",header=T,row.names=1)

#sample_groups_list = read.table("combined_sample_codes_11_7_2021.txt",header=T,row.names=1,sep = "\t")
sample_groups_list = read.csv("combined_sample_codes_11_7_2021.txt",header=T,row.names=1,sep = "\t")


#batch <- sample_groups_list$batch
groups <- sample_groups_list$Groups
section <- sample_groups_list$S_n
treatment <- sample_groups_list$Treatment
Genotype <- sample_groups_list$Genotype

write.table(sample_groups_list,"combined_sample_codes_11_30_2021.txt",sep = "\t")

groups
colnames(ReadCounts)
head(ReadCounts)
tail(ReadCounts)

######groups####
# groups = c("B73_Ssa","B73_Ssa","B73_Ssa","B73_Ssa","B73_Ssa","B73_Ssa",
#            "B73_SSb","B73_SSb","B73_SSb","B73_SSb","B73_SSb","B73_SSb",
#            "B73_SSc","B73_SSc","B73_SSc","B73_SSc","B73_SSc","B73_SSc",
#            "B73_Wwa","B73_Wwa","B73_Wwa","B73_Wwa","B73_Wwa","B73_Wwa",
#            "B73_Wwb","B73_Wwb","B73_Wwb","B73_Wwb","B73_Wwb","B73_Wwb",
#            "B73_Wwc","B73_Wwc","B73_Wwc","B73_Wwc","B73_Wwc",
#            "FR697_Ssa","FR697_Ssa","FR697_Ssa","FR697_Ssa","FR697_Ssa","FR697_Ssa",
#            "FR697_Ssb","FR697_Ssb","FR697_Ssb","FR697_Ssb","FR697_Ssb","FR697_Ssb",
#            "FR697_SSc","FR697_SSc","FR697_SSc","FR697_SSc","FR697_SSc","FR697_SSc",
#            "FR697_Wwa","FR697_Wwa","FR697_Wwa","FR697_Wwa","FR697_Wwa","FR697_Wwa",
#            "FR697_Wwb","FR697_Wwb","FR697_Wwb","FR697_Wwb","FR697_Wwb",
#            "FR697_Wwc","FR697_Wwc","FR697_Wwc","FR697_Wwc","FR697_Wwc",
#            "Lab_FR697_Msa","Lab_FR697_Msa","Lab_FR697_Msa","Lab_FR697_Msa","Lab_FR697_Msa",
#            "Lab_FR697_Msb","Lab_FR697_Msb","Lab_FR697_Msb","Lab_FR697_Msb","Lab_FR697_Msb",
#            "Lab_FR697_Msc","Lab_FR697_Msc","Lab_FR697_Msc","Lab_FR697_Msc","Lab_FR697_Msc",
#            "Lab_FR697_Ssa","Lab_FR697_Ssa","Lab_FR697_Ssa","Lab_FR697_Ssa","Lab_FR697_Ssa",
#            "Lab_FR697_Ssb","Lab_FR697_Ssb","Lab_FR697_Ssb","Lab_FR697_Ssb","Lab_FR697_Ssb",
#            "Lab_FR697_Ssc","Lab_FR697_Ssc","Lab_FR697_Ssc","Lab_FR697_Ssc","Lab_FR697_Ssc",
#            "Lab_FR697_Wwa","Lab_FR697_Wwa","Lab_FR697_Wwa","Lab_FR697_Wwa","Lab_FR697_Wwa",
#            "Lab_FR697_Wwb","Lab_FR697_Wwb","Lab_FR697_Wwb","Lab_FR697_Wwb","Lab_FR697_Wwb",
#            "Lab_FR697_Wwc","Lab_FR697_Wwc","Lab_FR697_Wwc","Lab_FR697_Wwc","Lab_FR697_Wwc",
#            "FR697_Wwa","FR697_Wwb","FR697_Wwc",
#            "B73_Wwa","B73_Wwb","B73_Wwc",
#            "B73_Ssa","B73_SSb","B73_SSc",
#            "B73_Wwa","B73_Wwb","B73_Wwc",
#            "FR697_Ssa","FR697_Ssb","FR697_SSc",
#            "FR697_Wwa","FR697_Wwb","FR697_Wwc",
#            "B73_Ssa","B73_SSb","B73_SSc",
#            "B73_Ssa","B73_SSb","B73_SSc",
#            "FR697_Wwa","FR697_Wwb","FR697_Wwc",
#            "FR697_Ssa","FR697_Ssb","FR697_SSc",
#            "B73_Wwa","B73_Wwb","B73_Wwc",
#            "FR697_Wwa","FR697_Wwb","FR697_Wwc",
#            "FR697_Ssa","FR697_Ssb","FR697_SSc",
#            "B73_Ssa","B73_SSb","B73_SSc",
#            "B73_Wwa","B73_Wwb","B73_Wwc",
#            "FR697_Ssa","FR697_Ssb","FR697_SSc",
#            "B73_Wwa","B73_Wwb","B73_Wwc",
#            "B73_Ssa","B73_SSb","B73_SSc",
#            "FR697_Wwa","FR697_Wwb","FR697_Wwc")
# 
#years = c(rep("2017" , 114) , rep("2018" , 57))
# region = c(rep("A",6),rep("B",6),rep("C",6),
#            rep("A",6),rep("B",6),rep("C",5),
#            rep("A",6),rep("B",6),rep("C",6),
#            rep("A",6),rep("B",5),rep("C",5),
#            rep("A",5),rep("B",5),rep("C",5),
#            rep("A",5),rep("B",5),rep("C",5),
#            rep("A",5),rep("B",5),rep("C",5),
#            "A","B","C","A","B","C","A","B","C","A","B","C","A","B","C",
#            "A","B","C","A","B","C","A","B","C","A","B","C","A","B","C",
#            "A","B","C","A","B","C","A","B","C","A","B","C","A","B","C",
#            "A","B","C","A","B","C","A","B","C","A","B","C")
# 
# genotype = c(rep("B73",35),rep("FR697",82),rep("B73",9),rep("FR697",6),
#              rep("B73",6),rep("FR697",6),rep("B73",3),rep("FR697",6),
#              rep("B73",6),rep("FR697",3),rep("B73",6),rep("FR697",3))
#            
# dim(region)
# dim(genotype)
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

###extra filtering
keep <- rowSums(cpm(d)>1) >= 3
#keep <- rowSums(cpm(d)>1) >= 6
dim(d[keep,])
d_filt <- d[keep,]
dim(d_filt)
#boxplot(d_filt$counts, col="gray", las=3)
#number filtered
dim(d)[1]-dim(d_filt)[1] ##18197
#%filtered out
(dim(d)[1]-dim(d_filt)[1])/dim(d)[1] ##0.3932616
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

#write.csv(normalized_data,file = "norm_filtered_counts_3_batches_10_07_21.csv")

#write.csv(norm_zero_combined_data,file = "norm_zero_combined_counts_11_08_2021.csv")

colnames(normalized_data)<- sample_groups_list$combined


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
normalized_data_groups['section'] <- section
dim(normalized_data_groups)
normalized_data_groups['treatment'] <- treatment
dim(normalized_data_groups)
normalized_data_groups['Genotype'] <- Genotype
dim(normalized_data_groups)

##plot pca
autoplot(prcomp(normalized_data_t))

##plot pca with group colors
#autoplot(prcomp(normalized_data_groups[1:27865]), data = normalized_data_groups, colour = 'groups' , label = TRUE, label.size = 3,shape = FALSE)
autoplot(prcomp(normalized_data_groups[1:27865]), data = normalized_data_groups, 
         colour = 'section', 
         label = FALSE,
         shape = 'treatment',
         size = 3,
         #shape.size = 5,
         label.size = 5)

autoplot(prcomp(normalized_data_groups[1:27865]), data = normalized_data_groups, 
         colour = 'Genotype', 
         label = FALSE,
         shape = 'treatment',
         size = 3,
         #shape.size = 5,
         label.size = 5)

####3d PCA

pca.data <- prcomp(normalized_data_groups[1:27865])
scores = as.data.frame(pca.data$x)
plot(pca.data, type="lines")
summary(pca.data)

#plot3d(pca.data$scores[,1:3])

pca3d(pca.data,group=sample_groups_list$Groups)
pca3d(pca.data,group=sample_groups_list$rep_n)
pca3d(pca.data,group=sample_groups_list$treatment,col = sample_groups_list$Color_grp,legend = "topright")
pca3d(pca.data,group=sample_groups_list$treatment,col = sample_groups_list$section_color)
?pca3d

sample_groups_list

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

d_filt

#EdgeR DE
dge_final <- estimateGLMCommonDisp(d_filt, design, verbose=TRUE) 
dge_final <- estimateGLMTrendedDisp(dge_final, design,verbose=TRUE)
####disp####
# Disp = 0 , BCV = 2e-04 
# Disp = 0 , BCV = 3e-04 
# Disp = 0 , BCV = 3e-04 
# Disp = 0 , BCV = 3e-04 
# Disp = 0 , BCV = 4e-04 
# Disp = 0 , BCV = 2e-04 
# Disp = 0 , BCV = 3e-04 
# Disp = 0 , BCV = 4e-04 
# Disp = 0 , BCV = 4e-04 
# Disp = 0 , BCV = 4e-04 
# Disp = 0 , BCV = 4e-04 
# Disp = 0 , BCV = 5e-04 
# Disp = 0 , BCV = 5e-04 
# Disp = 0 , BCV = 4e-04 
# Disp = 0 , BCV = 1e-04 
# Disp = 0 , BCV = 6e-04 
# Disp = 0.00162 , BCV = 0.0402 
# Disp = 0 , BCV = 6e-04 
# Disp = 0.00483 , BCV = 0.0695 
# Disp = 0.0091 , BCV = 0.0954 
# Disp = 0.00485 , BCV = 0.0696 
# Disp = 0.00403 , BCV = 0.0635 
# Disp = 0.00925 , BCV = 0.0962 
# Disp = 0.00746 , BCV = 0.0864 
# Disp = 0.00476 , BCV = 0.069 
# Disp = 0.00287 , BCV = 0.0536 
# Disp = 0.00399 , BCV = 0.0632 
# Disp = 0.00374 , BCV = 0.0611 
# Disp = 0.00537 , BCV = 0.0733 
# Disp = 0.00597 , BCV = 0.0773 
# Disp = 0.00419 , BCV = 0.0647 
# Disp = 0.00399 , BCV = 0.0631 
# Disp = 0.00799 , BCV = 0.0894 
# Disp = 0.01055 , BCV = 0.1027 
# Disp = 0.00595 , BCV = 0.0771 
# Disp = 0.00678 , BCV = 0.0823 
# Disp = 0.00492 , BCV = 0.0701 
# Disp = 0.01136 , BCV = 0.1066 
# Disp = 0.01205 , BCV = 0.1098 
# Disp = 0.00645 , BCV = 0.0803 
# Disp = 0.00788 , BCV = 0.0887 
# Disp = 0.0094 , BCV = 0.097 
# Disp = 0.00979 , BCV = 0.099 
# Disp = 0.01328 , BCV = 0.1153 
# Disp = 0.01445 , BCV = 0.1202 
# Disp = 0.01212 , BCV = 0.1101 
# Disp = 0.01255 , BCV = 0.112 
# Disp = 0.01406 , BCV = 0.1186 
# Disp = 0.01271 , BCV = 0.1128 
# Disp = 0.0178 , BCV = 0.1334 
# Disp = 0.01679 , BCV = 0.1296 
# Disp = 0.01837 , BCV = 0.1355 
# Disp = 0.02157 , BCV = 0.1469 
# Disp = 0.02054 , BCV = 0.1433 
# Disp = 0.02484 , BCV = 0.1576 
# Disp = 0.02627 , BCV = 0.1621 
# Disp = 0.03162 , BCV = 0.1778 
# Disp = 0.03462 , BCV = 0.1861 
# Disp = 0.04022 , BCV = 0.2005 
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

save.image(file='All_DE_11_09_2021.RData')

####CONTRASTS & DE Tables####

comaparisons = read.table("comparisons_to_run_11_08_2021.txt",header=T)

numbers_df <- data.frame(comparison=character(),Pval=integer(),FDR=integer(),twofold=integer())

# for (row in 1:nrow(comaparisons)){
#   lrt_value = toString(comaparisons[row,4])
#   grps_value = toString(comaparisons[row,5])
#   comp_value = toString(comaparisons[row,3])
# 
#   lrt_value= glmLRT(fit,contrast=makeContrasts(grps_value, levels=design))
#   comp_value<- as.data.frame(topTags(lrt_value, n=Inf))
#   #print(dim(comp_value[comp_value$PValue < 0.05,]))
#   print(nrow(comp_value[comp_value$PValue < 0.05,]))
#   # write.csv(comp_value[comp_value$PValue < 0.05,],
#   #           paste("EdgeR_Pval/",toString(comaparisons[row,3]),"_Pval.csv",sep = ''),row.names=T)
# 
#   #print(dim(comp_value[comp_value$FDR < 0.05,]))
#   # write.csv(comp_value[comp_value$FDR < 0.05,],
#   #           paste("EdgeR_FDR/",toString(comaparisons[row,3]),"_FDR.csv",sep = ''),row.names=T)
# 
#   # print(dim(comp_value[(comp_value$FDR < 0.05) & 
#   #                        ((comp_value$logFC > 1)|
#   #                           (comp_value$logFC < -1)),]))
#   # write.csv(comp_value[(comp_value$FDR < 0.05) & 
#   #                        ((comp_value$logFC > 1)|
#   #                           (comp_value$logFC < -1)),],
#   #           paste("EdgeR_2fold/",toString(comaparisons[row,3]),"_2fold.csv",sep = ''),row.names=T)
# 
#   add_row(comp_value, comparison=, ######################
#   #print(dim(comp_value))
#   # write.csv(comp_value,
#   #           paste("EdgeR_ALL/",toString(comaparisons[row,3]),"_ALL.csv",sep = ''),row.names=T)
# }


#lrt.B73_WW_A_vs_B73_WW_B= glmLRT(fit,contrast=makeContrasts(groupsB73_WW_A-groupsB73_WW_B, levels=design))
#B73_WW_A_vs_B73_WW_B<- as.data.frame(topTags(lrt.B73_WW_A_vs_B73_WW_B, n=Inf))
#dim(B73_WW_A_vs_B73_WW_B[B73_WW_A_vs_B73_WW_B$PValue < 0.05,])
# write.csv(B73_WW_A_vs_B73_WW_B[B73_WW_A_vs_B73_WW_B$PValue < 0.05,], 
#           file= "EdgeR_Pval/B73_WW_A_vs_B73_WW_B_Pval.csv",row.names=T)
# dim(B73_WW_A_vs_B73_WW_B[B73_WW_A_vs_B73_WW_B$FDR < 0.05,])
# write.csv(B73_WW_A_vs_B73_WW_B[B73_WW_A_vs_B73_WW_B$FDR < 0.05,], 
#           file= "EdgeR_FDR/B73_WW_A_vs_B73_WW_B_FDR.csv",row.names=T)

#dim(B73_WW_A_vs_B73_WW_B[(B73_WW_A_vs_B73_WW_B$FDR < 0.05) & ((B73_WW_A_vs_B73_WW_B$logFC > 1)|(B73_WW_A_vs_B73_WW_B$logFC < -1)),])

# write.csv(B73_WW_A_vs_B73_WW_B[(B73_WW_A_vs_B73_WW_B$FDR < 0.05) & 
#                                  ((B73_WW_A_vs_B73_WW_B$logFC > 1)|
#                                     (B73_WW_A_vs_B73_WW_B$logFC < -1)),],
#           file= "EdgeR_2fold/B73_WW_A_vs_B73_WW_B_2fold.csv",row.names=T)
# 
# write.csv(B73_WW_A_vs_B73_WW_B, file= "EdgeR_ALL/B73_WW_A_vs_B73_WW_B_all.csv",row.names=T)
























####volcano plot####
Lab_FR697_SSc_vs_WWc <- as.data.frame(topTags(lrt.Lab_FR697_SSc_vs_WWc, n=Inf))
library(EnhancedVolcano)

colnames(Lab_FR697_SSc_vs_WWc)
rownames(Lab_FR697_SSc_vs_WWc)

##og##dif <- data.frame(fc =fold_changes,pv =pvalues)
dif <- data.frame(genes = Lab_FR697_SSc_vs_WWc$genes,
                  logfc = Lab_FR697_SSc_vs_WWc$logFC,
                  fdr = Lab_FR697_SSc_vs_WWc$FDR)
dif$threshold <- ifelse(dif$logfc > 1 & dif$fdr < 0.05, "upreg",
                        ifelse(dif$logfc < -1 & dif$fdr < 0.05,"downreg","NotSign" ))
head(dif)

vol_plot_2<-ggplot(data=dif, aes(x=logfc, y=-log10(fdr))) +
  geom_point( size=1 ,aes(color=as.factor(threshold))) +
  theme(legend.position = "none") +
  xlim(c(-10, 10)) + ylim(c(0, 15)) +
  xlab("log2 fold change") + ylab("-log10 p-value")  + theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),)

vol_plot_2
vol_plot<-ggplot(Lab_FR697_SSc_vs_WWc, aes(x = logFC, y=-log10(FDR), 
                                           col=(FDR < 0.05 & ((logFC > 1)|(logFC < -1))))) + geom_point()
vol_plot

#ggplot_build(vol_plot)
#ggplot_build(vol_plot)$data[[1]]$x
#ggplot_build(vol_plot)$data[[1]]$y
ggplot_build(vol_plot_2)$data[[1]]$x
ggplot_build(vol_plot_2)$data[[1]]$y


volc_plot_XY<-cbind(vol_plot[["data"]][["genes"]],ggplot_build(vol_plot)$data[[1]]$x,ggplot_build(vol_plot)$data[[1]]$y,dif$threshold)
dim(volc_plot_XY)
head(volc_plot_XY)


volcano_fxn <- function(inputFile){
  #add classificationn for significant genes
  # dif <- data.frame(genes = Lab_FR697_SSc_vs_WWc$genes,
  #                   logfc = Lab_FR697_SSc_vs_WWc$logFC,
  #                   fdr = -log10(Lab_FR697_SSc_vs_WWc$FDR))
  
  # dif$threshold <- ifelse(dif$logfc > 1 & dif$fdr < 0.05, "upreg",
  #                         ifelse(dif$logfc < -1 & dif$fdr < 0.05,"downreg","NotSign" ))
  
  # vol_plot<-ggplot(Lab_FR697_SSc_vs_WWc, aes(x = logFC, y=-log10(FDR), 
  #                                            col=(FDR < 0.05 & ((logFC > 1)|(logFC < -1))))) + geom_point()
  
  volc_plot_XY<-data.frame(inputFile$genes,
                           volc.X = inputFile$logFC,
                           volc.y = -log10(inputFile$FDR),
                           threshold = ifelse(inputFile$logfc > 1 & inputFile$fdr < 0.05, "upreg",
                                              ifelse(inputFile$logfc < -1 & inputFile$fdr < 0.05,"downreg","NotSign" ))
                           )
  return(volc_plot_XY)
  
}

######

# EnhancedVolcano(Lab_FR697_SSc_vs_WWc,
#                 lab = rownames(Lab_FR697_SSc_vs_WWc),
#                 x = 'logFC',
#                 y = 'FDR',
#                 pCutoff = 0.05,
#                 FCcutoff = 1,
#                 # ylab = bquote("FDR"),
#                 xlim = c(-15,15))

##ggplot volcano plot
#ggplot(Lab_FR697_SSc_vs_WWc, aes(x = logFC, y=-log10(FDR), col=FDR < 0.05)) + geom_point()


# [(Lab_FR697_SSc_vs_WWc$FDR < 0.05) & 
#     ((Lab_FR697_SSc_vs_WWc$logFC > 1)|
#        (Lab_FR697_SSc_vs_WWc$logFC < -1)),]


#-log10(0.05)

######
######filter counts tible by fold2 dataframe gene list
##counts_Field_B73_SSa_vs_WWa (tibble)
##fold2_Field_B73_SSa_vs_WWa (data.frame)

#colnames(fold2_Field_B73_SSa_vs_WWa)

library(tidyverse)
 
tb_normalized_data <- as_tibble(normalized_data, rownames = "rowname")
rownames_to_column(.data, var = "rowname")

## SELECT COLUMNS SPECIFIC TO DE COMPARISON
#Field_B73
#1
counts_Field_B73_SSa_vs_WWa = select(tb_normalized_data,rowname,
                                     S1031_2017:S6031_2017,s1041_2018,s2041_2018,s3011_2018,s4031_2018,s5031_2018,  #ssa
                                     S1011_2017:S6011_2017,s1021_2018,s2011_2018,s3041_2018,s4041_2018,s5021_2018)  #wwa
#2
counts_Field_B73_SSb_vs_WWb = select(tb_normalized_data,rowname,
                                     S1032_2017:S6032_2017,s1042_2018,s2042_2018,s3012_2018,s4032_2018,s5032_2018,  #ssb
                                     S1012_2017:S6012_2017,s1022_2018,s2012_2018,s3042_2018,s4042_2018,s5022_2018)  #wwb
#3
counts_Field_B73_SSc_vs_WWc = select(tb_normalized_data,rowname,
                                     S1033_2017:S6033_2017,s1043_2018,s2043_2018,s3013_2018,s4033_2018,s5033_2018,  #ssb
                                     S1013_2017:S6013_2017,s1023_2018,s2013_2018,s3043_2018,s4043_2018,s5023_2018)  #wwb


#Field_FR697
#4
counts_Field_FR697_SSa_vs_WWa = select(tb_normalized_data,rowname,
                                       S1021_2017:S6021_2017,s2021_2018,s3031_2018,s4021_2018,s5011_2018,  #ssa
                                       S1041_2017:S6041_2017,s1011_2018,s2031_2018,s3021_2018,s4011_2018,s5041_2018)  #wwa
#5
counts_Field_FR697_SSb_vs_WWb = select(tb_normalized_data,rowname,
                                       S1022_2017:S6022_2017,s2022_2018,s3032_2018,s4022_2018,s5012_2018,  #ssb
                                       S2012_2017:S6042_2017,s1012_2018,s2032_2018,s3022_2018,s4012_2018,s5042_2018)  #wwb
#6
counts_Field_FR697_SSc_vs_WWc = select(tb_normalized_data,rowname,
                                       S1023_2017:S6023_2017,s2023_2018,s3033_2018,s4023_2018,s5013_2018,  #ssb
                                       S1043_2017:S6043_2017,s1013_2018,s2033_2018,s3023_2018,s4013_2018,s5043_2018)  #wwb

#Lab_FR697
#7
counts_Lab_FR697_SSa_vs_WWa = select(tb_normalized_data,rowname,
                                     S8021_2017:S8051_2017,  #ssa
                                     S9011_2017:S9051_2017)  #wwa
#8
counts_Lab_FR697_SSb_vs_WWb = select(tb_normalized_data,rowname,
                                     S8012_2017:S8052_2017,  #ssb
                                     S9012_2017:S9052_2017)  #wwb
#9
counts_Lab_FR697_SSc_vs_WWc = select(tb_normalized_data,rowname,
                                     S8013_2017:S8053_2017,  #ssb
                                     S9013_2017:S9053_2017)  #wwb
##initialize this fxn first
counts_list_filter <- function(counts_file,de_file,file_name){
  fold2_file <- as.data.frame(topTags(de_file, n=Inf))
  filt_fold2_file <- fold2_file[which((fold2_file$FDR < 0.05) & 
                                                             ((fold2_file$logFC > 1)|
                                                                (fold2_file$logFC < -1))),]
  filter_fold2_counts_file <- filter(counts_file, rowname %in% filt_fold2_file$genes)
  
  write.table(filter_fold2_counts_file,
              file = sprintf("C:/Users/Sid/OneDrive/rna-seq documents/NSF-drought/rnaseq_counts_combined/Field_lab/counts_gene_filtered/%s.txt",file_name),
              sep = "\t",
              row.names = FALSE
  )
}

#counts_file = counts_Field_B73_SSa_vs_WWa
#de_file = lrt.Field_B73_SSa_vs_WWa
#file name = filter_fold2_Field_B73_SSa_vs_WWa

counts_list_filter(counts_Field_B73_SSa_vs_WWa,lrt.Field_B73_SSa_vs_WWa,"filter_fold2_Field_B73_SSa_vs_WWa")
counts_list_filter(counts_Field_B73_SSb_vs_WWb,lrt.Field_B73_SSb_vs_WWb,"filter_fold2_Field_B73_SSb_vs_WWb")
counts_list_filter(counts_Field_B73_SSc_vs_WWc,lrt.Field_B73_SSc_vs_WWc,"filter_fold2_Field_B73_SSc_vs_WWc")

counts_list_filter(counts_Field_FR697_SSa_vs_WWa,lrt.Field_FR697_SSa_vs_WWa,"filter_fold2_Field_FR697_SSa_vs_WWa")
counts_list_filter(counts_Field_FR697_SSb_vs_WWb,lrt.Field_FR697_SSb_vs_WWb,"filter_fold2_Field_FR697_SSb_vs_WWb")
counts_list_filter(counts_Field_FR697_SSc_vs_WWc,lrt.Field_FR697_SSc_vs_WWc,"filter_fold2_Field_FR697_SSc_vs_WWc")

counts_list_filter(counts_Lab_FR697_SSa_vs_WWa,lrt.Lab_FR697_SSa_vs_WWa,"filter_fold2_Lab_FR697_SSa_vs_WWa")
counts_list_filter(counts_Lab_FR697_SSb_vs_WWb,lrt.Lab_FR697_SSb_vs_WWb,"filter_fold2_Lab_FR697_SSb_vs_WWb")
counts_list_filter(counts_Lab_FR697_SSc_vs_WWc,lrt.Lab_FR697_SSc_vs_WWc,"filter_fold2_Lab_FR697_SSc_vs_WWc")


# fold2_Field_B73_SSa_vs_WWa <- Field_B73_SSa_vs_WWa[which((Field_B73_SSa_vs_WWa$FDR < 0.05) & 
#                                                            ((Field_B73_SSa_vs_WWa$logFC > 1)|
#                                                               (Field_B73_SSa_vs_WWa$logFC < -1))),]
# filter_fold2_Field_B73_SSa_vs_WWa <- filter(counts_Field_B73_SSa_vs_WWa, rowname %in% fold2_Field_B73_SSa_vs_WWa$genes)
# colnames(filter_fold2_Field_B73_SSa_vs_WWa)
# 
# write.table(filter_fold2_Field_B73_SSa_vs_WWa,
#           "C:/Users/Sid/OneDrive/rna-seq documents/NSF-drought/rnaseq_counts_combined/Field_lab/counts_gene_filtered/filter_fold2_Field_B73_SSa_vs_WWa.txt",
#           sep = "\t",
#           row.names = FALSE
#           )




####end####