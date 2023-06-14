#!/nfs/sw/R/R-3.6.1/bin/Rscript
library(dplyr)
require(NanoStringQCPro)
require(ggplot2)
require(EnvStats)
library('pheatmap')

library(RUVSeq)
library(DESeq2)
library(limma)
library(matrixStats)

library(survival)
library(stringr)
library(survminer)
library(limma)

library(lsa)
library(tidyr)
library(biomaRt)
library(ggrepel)

library("xlsx")
library("readxl")
library(writexl)

 
library(ComplexHeatmap)
library(circlize)
library(rstatix)

library(ggpubr)


# setwd("/pathTo/analysis")
source('./nanostring_RUV_functions.r')
source("./nanostring_functions.r")


load("brainmet_081522.RData")
# save.image("brainmet_081522.RData")


####
primary_prostate <- read.table("primary_protate.txt", header=TRUE, sep='\t')
primary_prostate$Met_or_prostate <- primary_prostate$GeneralTissueSite
rownames(primary_prostate) <- primary_prostate$NanostringID
primary_prostate_codeset375 <- na.omit(primary_prostate[colnames(nanoData$codeset375)[c(-1,-2)] ,])
primary_prostate_codeset375_179 <- primary_prostate_codeset375

write_xlsx(primary_prostate_codeset375[, -6], "primary_prostate_codeset375.xlsx")

# 140/179, after removing the denovo small cell cases, and the cell line cases from Thompson)
primary_prostate_codeset375 <- read.table("primary_protate_140_updated_MS.txt", header=TRUE, sep='\t')
primary_prostate_codeset375$Met_or_prostate <- primary_prostate_codeset375$GeneralTissueSite
rownames(primary_prostate_codeset375) <- primary_prostate_codeset375$NanostringID

# anno <- read.table("annoData.txt", header=TRUE, sep='\t')
# anno_old <- anno
# # "20150918_NEPC34_NEPC34_11" / "20180925_NEPC76_NEPC76_05" is in codeset163 and excluded because less gene numbers.
# anno_codeset375 <- subset(anno, Nanostring_ID!="20150918_NEPC34_NEPC34_11")
anno <- read_xlsx("annoData_updated.xlsx")
anno <- as.data.frame(anno)
anno_codeset375 <- subset(anno, Nanostring_ID!="20180925_NEPC76_NEPC76_05")

table(smp_anno_all_matched_31anno$NanostringID == anno_codeset375$Nanostring_ID)

# TRUE 
  # 31 

anno_codeset375 <- cbind(anno_codeset375, smp_anno_all_matched_31anno[, c(2, 3)])

anno_codeset375_Brain_nonBrain <- subset(anno_codeset375, Brain_vs_nonBrain != "Prostate")

###########

# load("../RData_raw_andrea/nanoData_Raw_20220815.RData")
names <- c(colnames(nanoData$codeset89)[c(-1, -2)], colnames(nanoData$codeset163)[c(-1, -2)], colnames(nanoData$codeset375)[c(-1, -2)], colnames(nanoData$codeset262)[c(-1, -2)], colnames(nanoData$codeset784)[c(-1, -2)], colnames(nanoData$codeset431)[c(-1, -2)])

# old anno with "20150918_NEPC34_NEPC34_11", "20180925_NEPC76_NEPC76_05"
all(anno$Nanostring_ID %in% colnames(nanoData$codeset89))
# [1] FALSE
all(anno$Nanostring_ID %in% colnames(nanoData$codeset163))
# [1] FALSE
all(anno$Nanostring_ID %in% colnames(nanoData$codeset375))
# [1] FALSE
all(anno$Nanostring_ID %in% colnames(nanoData$codeset262))
# [1] FALSE
all(anno$Nanostring_ID %in% colnames(nanoData$codeset784))
# [1] FALSE
all(anno$Nanostring_ID %in% colnames(nanoData$codeset431))
# [1] FALSE
all(anno$Nanostring_ID %in% names)
# [1] FALSE


# 4 ids are with a space in name, delete the space because there is no space in the names of nanoData
# [1] "20181120_NEPC 83_NEPC83_01" "20181120_NEPC 83_NEPC83_10" "20181120_NEPC 83_NEPC83_11"
# [4] "20181120_NEPC 83_NEPC83_12"
# anno_199 <- read.table("anno_test.txt", header=TRUE, sep='\t')
anno_199 <- read_xlsx("anno_32wgs_52npec_140prostate_195nanoIds_081122.xlsx", sheet=1)
anno_199 <- head(anno_199, n=199)

# 51 in final, 1  20160623_NEPC49_NEPC49_02 is from unknown tissue, so dropped.
anno_52_nepc <- as.data.frame(read_xlsx("../from_researcher_or_DL/CRPC-NEPC cases for brain met analysis.xlsx", sheet=1))

anno_199 <- as.data.frame(anno_199)


Primary_Prostate_Untrt <- filter(anno_199, anno_199$Group == "Primary Prostate (Untreated)" & anno_199$`Met site or Prostate?`=="Prostate")
dim(Primary_Prostate_Untrt)
# [1] 139  28

met_trt <- filter(anno_199, anno_199$`Treated vs non-treated` == "Treated" & anno_199$`Metastatic or non-metastatic?`=="Metastatic")
dim(met_trt)

all_met <- filter(anno_199, anno_199$`Metastatic or non-metastatic?`=="Metastatic")
dim(all_met)

brainMet_trt <- filter(anno_199, anno_199$`Treated vs non-treated` == "Treated" & anno_199$Group=="Met/Brain")
dim(brainMet_trt)

nonbrainMet_trt <- filter(anno_199, anno_199$`Treated vs non-treated` == "Treated" & anno_199$Group=="Met/Non-Brain")
dim(nonbrainMet_trt)

allProstate <- filter(anno_199, anno_199$`Met site or Prostate?`=="Prostate")
dim(allProstate)

prostate_trt <- filter(anno_199, anno_199$Group == "Prostate (Treated)")
dim(prostate_trt)

dural <- filter(anno_199, anno_199$Brain_Met_Type=="Dural")
dim(dural)

parenchymal <- filter(anno_199, anno_199$Brain_Met_Type=="Parenchymal")
dim(parenchymal)

prostate_brainMet_nonbrainMet_trt <- filter(anno_199, anno_199$`Treated vs non-treated` == "Treated")
dim(prostate_brainMet_nonbrainMet_trt)

# Primary Prostate (Untreated): 139+1 failed. $Group
# Met [Met/Brain(treated) + Met/Non-brain(treated)]: 19+17. Metastatic or non-metastatic? and Treated vs non-treated
# Met/Brain(treated): 19. Group  and Treated vs non-treated
# Met/Non-brain(treated): 17. Group
# Prostate (untreated, treated or unknown treatment): 161. Met site or Prostate? 
# Met [Met/Brain(treated) + Met/Non-brain(treated)]: 36. 
# Met [Met/Brain(treated) + Met/Brain(untreated) + Met/Non-brain(treated)]: 37. Metastatic or non-metastatic?
# treated prostate: 16. Group == Prostate (Treated)
# Dural: 9. Brain_Met_Type==Dural 
# parenchimal: 10. Brain_Met_Type=="Parenchymal", including 1 case of untreated brain met
# treated prostate+Met/brain(treated)+Met/Non-brain(treated): 52. Treated vs non-treated=="Treated"

# 9 comparisons
# 1. Primary prostate (untreated) VS Met [Met/Brain(treated) + Met/Non-brain(treated)]
# 2. Primary prostate  (untreated)VS Met/Brain(treated)
# 3. Primary prostate  (untreated) VS Met/Non-brain(treated)
# 4. Met/Brain(treated) VS Met/Non-brain(treated)
# 5. Prostate (untreated, treated or unknown treatment) VS Met [Met/Brain(treated) + Met/Non-brain(treated)]
# 6. Prostate (untreated, treated or unknown treatment) VS Met [Met/Brain(treated) + Met/Brain(untreated) + Met/Non-brain(treated)]
# 7. Primary prostate VS treated prostate
# 8. Dural VS parenchimal
# 9. Primary prostate VS treated prostate+Met/brain(treated)+Met/Non-brain(treated): primary VS treated (brain, non-brain, and prostate)

# merge codesets because of difference in codeset size
nanoData_brainMet_codset375_262 <- mergeCodesets(nanoData[c(3,4)], nanoData$codeset262$Gene)

nanoData_brainMet_codset375_262_199ids <- nanoData_brainMet_codset375_262[, c(1,2, which(colnames(nanoData_brainMet_codset375_262) %in% anno_199$Nanostring_ID))]

# RUV_seq

# pData_match: 198 rows because one is failed in anno_199
pData_match <- anno_199[match(colnames(nanoData_brainMet_codset375_262)[c(-1,-2)], anno_199$Nanostring_ID)[!is.na(match(colnames(nanoData_brainMet_codset375_262)[c(-1,-2)], anno_199$Nanostring_ID))], ]
table(pData_match$Nanostring_ID == colnames(nanoData_brainMet_codset375_262_199ids)[c(-1,-2)])
# TRUE 
 # 198

raw = nanoData_brainMet_codset375_262_199ids[,-c(1:2)]
fData = nanoData_brainMet_codset375_262_199ids[,c(1:2)]
rownames(raw) = fData$Gene
cIdx <- fData$Gene[fData$Class == "Housekeeping"]
negative_controls <-  fData$Gene[fData$Class == "Negative"]
pData_match$HK_Gene_Miss = colSums(raw[cIdx,] == 0)
rownames(fData) = fData$Gene
rownames(pData_match) = colnames(raw)

# 261 genes: codeset375 and 262
write_xlsx(fData, "261genes_241_EndogenousGenes.xlsx")


#### CHECK IF HK ARE ASSOCIATED WITH PRIMARY PHENO
hk_raw = raw[cIdx,]
pval = vector(length = nrow(hk_raw))

require(MASS)

for (i in 1:nrow(hk_raw)){
  
  reg = glm.nb(as.numeric(hk_raw[i,]) ~ as.factor(pData_match$Group))
  pval[i] = coef(summary(reg))[2,4]
  
}

sum(pval <= .05)

#### samples with high numbers of endogenous and housekeeping genes below the limit of detection

neg_raw = raw[negative_controls,]
lod = colMeans(neg_raw) - apply(neg_raw,2,sd)
num_endogenous_blod = colSums(nanoData_brainMet_codset375_262_199ids[nanoData_brainMet_codset375_262_199ids$Class == 'Endogenous',-c(1:2)] < lod) 

# no samples with HK genes below LOD
num_hk_blod = colSums(nanoData_brainMet_codset375_262_199ids[nanoData_brainMet_codset375_262_199ids$Class == 'Housekeeping',-c(1:2)] < lod)


### grid of datasets
library(RUVSeq)
library(DESeq2)
library(limma)
library(matrixStats)

# norm.dat <- RUV_total(raw, pData_match[, c(1:4, 8:13)], fData, k=1)

norm.dat <- RUV_total(raw, pData_match[, c(1:4, 6:29)], fData, k=1)
# got error as below when adding "5"
# converting counts to integer mode
# Error in lmFit(x, cbind(design, X.batch), ...) : 
  # design must be a numeric matrix




# 1 st comparison
# Primary prostate (untreated) VS Met [Met/Brain(treated) + Met/Non-brain(treated)]
countData_1st = counts(norm.dat$set)[1:241, c(met_trt$Nanostring_ID, Primary_Prostate_Untrt$Nanostring_ID)]
colData_1st <- pData(norm.dat$set)[colnames(countData_1st), ]
table(rownames(colData_1st) == colnames(countData_1st))
# TRUE 
 # 175 
# not necessary to give trt1
# colData_1st$trt1 <- c(rep("brMet_nonBrMet_trt", 36), rep("primary_prostate_Untrt", 139))

colData_1st$Group2 <- colData_1st$`Treated vs non-treated`
 
# Primary_Prostate_Untrt_VS_metTrt_dds <- DESeqDataSetFromMatrix(countData = countData_1st, colData = colData_1st, design = ~ W_1 + "Treated vs non-treated")

Primary_Prostate_Untrt_VS_metTrt_dds <- DESeqDataSetFromMatrix(countData = countData_1st, colData = colData_1st, design = ~ W_1 + Group2)


Primary_Prostate_Untrt_VS_metTrt_dds <- DESeq(Primary_Prostate_Untrt_VS_metTrt_dds)

res_1st <- results(Primary_Prostate_Untrt_VS_metTrt_dds,contrast = c('Group2','Treated','Non-treated'))
res_1st.df <- as.data.frame(res_1st)

res_0.05_1st <- subset(res_1st, padj <0.05 & abs(log2FoldChange) >1) 
res_0.01_1st <- subset(res_1st, padj <0.01 & abs(log2FoldChange) >1) 

# 2 nd comparison
# Primary prostate (untreated)VS Met/Brain(treated)

countData_2nd = counts(norm.dat$set)[1:241, c(brainMet_trt$Nanostring_ID, Primary_Prostate_Untrt$Nanostring_ID)]
colData_2nd <- pData(norm.dat$set)[colnames(countData_2nd), ]
table(rownames(colData_2nd) == colnames(countData_2nd))
# TRUE 
 # 158
colData_2nd$Group2 <- colData_2nd$`Treated vs non-treated`

# 3 rd comparison
# Primary prostate  (untreated) VS Met/Non-brain(treated)
countData_3rd = counts(norm.dat$set)[1:241, c(nonbrainMet_trt$Nanostring_ID, Primary_Prostate_Untrt$Nanostring_ID)]
colData_3rd <- pData(norm.dat$set)[colnames(countData_3rd), ]
table(rownames(colData_3rd) == colnames(countData_3rd))
# TRUE 
 # 156 
colData_3rd$Group2 <- colData_3rd$`Treated vs non-treated`

# 4 th comparison
# Met/Brain(treated) VS Met/Non-brain(treated)
countData_4th = counts(norm.dat$set)[1:241, c(nonbrainMet_trt$Nanostring_ID, brainMet_trt$Nanostring_ID)]
colData_4th <- pData(norm.dat$set)[colnames(countData_4th), ]
table(rownames(colData_4th) == colnames(countData_4th))
# TRUE 
  # 36
colData_4th$Group2 <- colData_4th$Group

# 5 th comparison
# all Prostate (untreated, treated or unknown treatment) VS Met [Met/Brain(treated) + Met/Non-brain(treated)]
countData_5th = counts(norm.dat$set)[1:241, c(allProstate$Nanostring_ID, met_trt$Nanostring_ID)]
colData_5th <- pData(norm.dat$set)[colnames(countData_5th), ]
table(rownames(colData_5th) == colnames(countData_5th))
# TRUE 
  # 197
colData_5th$Group2 <- colData_5th$`Metastatic or non-metastatic?`

# 6 th comparison
# Prostate (untreated, treated or unknown treatment) VS Met [Met/Brain(treated) + Met/Brain(untreated) + Met/Non-brain(treated)]
countData_6th = counts(norm.dat$set)[1:241, c(allProstate$Nanostring_ID, all_met$Nanostring_ID)]
colData_6th <- pData(norm.dat$set)[colnames(countData_6th), ]
table(rownames(colData_6th) == colnames(countData_6th))
# TRUE 
  # 198
colData_6th$Group2 <- colData_6th$`Metastatic or non-metastatic?`

# 7 th comparison
# Primary prostate VS treated prostate
countData_7th = counts(norm.dat$set)[1:241, c(Primary_Prostate_Untrt$Nanostring_ID, prostate_trt$Nanostring_ID)]
colData_7th <- pData(norm.dat$set)[colnames(countData_7th), ]
table(rownames(colData_7th) == colnames(countData_7th))
# TRUE 
  # 155
colData_7th$Group2 <- colData_7th$`Treated vs non-treated`


# 8 th comparison
# Dural VS parenchimal
countData_8th = counts(norm.dat$set)[1:241, c(dural$Nanostring_ID, parenchymal$Nanostring_ID)]
colData_8th <- pData(norm.dat$set)[colnames(countData_8th), ]
table(rownames(colData_8th) == colnames(countData_8th))
# TRUE 
  # 155
colData_8th$Group2 <- colData_8th$Brain_Met_Type

# 9 th comparison
# Primary prostate VS treated [prostate+Met/brain(treated)+Met/Non-brain(treated)]: primary VS treated (brain, non-brain, and prostate)
countData_9th = counts(norm.dat$set)[1:241, c(Primary_Prostate_Untrt$Nanostring_ID, prostate_brainMet_nonbrainMet_trt$Nanostring_ID)]
colData_9th <- pData(norm.dat$set)[colnames(countData_9th), ]
table(rownames(colData_9th) == colnames(countData_9th))
# TRUE 
  # 191
colData_9th$Group2 <- colData_9th$`Treated vs non-treated`


countData_ls <- list(countData_1st, countData_2nd, countData_3rd, countData_4th, countData_5th, countData_6th, countData_7th, countData_8th, countData_9th)
colData_ls <- list(colData_1st, colData_2nd, colData_3rd, colData_4th, colData_5th, colData_6th, colData_7th, colData_8th, colData_9th)
names(countData_ls) <- c("comp1", "comp2", "comp3", "comp4", "comp5", "comp6", "comp7", "comp8", "comp9")
names(colData_ls) <- c("comp1", "comp2", "comp3", "comp4", "comp5", "comp6", "comp7", "comp8", "comp9")

lapply(countData_ls, dim)

dds_ls <- list()


for (i in 1:9){
	dds_ls[[i]] <- DESeqDataSetFromMatrix(countData = countData_ls[[i]], colData = colData_ls[[i]], design = ~ W_1 + Group2)
	dds_ls[[i]] <- DESeq(dds_ls[[i]])

}

names(dds_ls) <- c("comp1", "comp2", "comp3", "comp4", "comp5", "comp6", "comp7", "comp8", "comp9")

# drop the sample from NEPC
dds_ls[[8]] <- DESeqDataSetFromMatrix(countData = countData_ls[[8]][, -13], colData = colData_ls[[8]][-13,], design = ~ W_1 + Group2)
dds_ls[[8]] <- DESeq(dds_ls[[8]])

res_1st <- results(dds_ls[[1]],contrast = c('Group2','Treated','Non-treated'))
res_2nd <- results(dds_ls[[2]],contrast = c('Group2','Treated','Non-treated'))
res_3rd <- results(dds_ls[[3]],contrast = c('Group2','Treated','Non-treated'))
res_4th <- results(dds_ls[[4]],contrast = c('Group2','Met/Brain','Met/Non-Brain'))
res_5th <- results(dds_ls[[5]],contrast = c('Group2','Metastatic','Non-metastatic'))
res_6th <- results(dds_ls[[6]],contrast = c('Group2','Metastatic','Non-metastatic'))
res_7th <- results(dds_ls[[7]],contrast = c('Group2','Treated','Non-treated'))
res_8th <- results(dds_ls[[8]],contrast = c('Group2','Parenchymal','Dural'))
res_9th <- results(dds_ls[[9]],contrast = c('Group2','Treated','Non-treated'))

for (i in 1:9){
	print(table(colData_ls[[i]]$Group2))
}

res_ls <- list(res_1st, res_2nd, res_3rd, res_4th, res_5th, res_6th, res_7th, res_8th, res_9th) 
names(res_ls) <- c("comp1", "comp2", "comp3", "comp4", "comp5", "comp6", "comp7", "comp8", "comp9")

res_0.05_ls <- list()
res_0.01_ls <- list()
res_0.05_ls_1.5logFC <- list()


for (i in 1:9){
	res_0.05_ls[[i]] <- subset(res_ls[[i]], padj <0.05 & abs(log2FoldChange) >1) 
	res_0.01_ls[[i]] <- subset(res_ls[[i]], padj <0.01 & abs(log2FoldChange) >1) 
	res_0.05_ls_1.5logFC[[i]] <- subset(res_ls[[i]], padj <0.05 & abs(log2FoldChange) >=1.5) 
}

names(res_0.05_ls) <- c("comp1", "comp2", "comp3", "comp4", "comp5", "comp6", "comp7", "comp8", "comp9")
names(res_0.01_ls) <- c("comp1", "comp2", "comp3", "comp4", "comp5", "comp6", "comp7", "comp8", "comp9")
names(res_0.05_ls_1.5logFC) <- c("comp1", "comp2", "comp3", "comp4", "comp5", "comp6", "comp7", "comp8", "comp9")


# write_xlsx(cbind(" "=rownames(res_0.05_ls[[4]]), as.data.frame(res_0.05_ls[[4]])), "brainMets_vs_nonBrainMets_19vs17_comp4.xlsx")

write.xlsx(cbind(" "=rownames(res_0.05_ls[[4]]), as.data.frame(res_0.05_ls[[4]])), "brainMets_vs_nonBrainMets_19vs17_comp4.xlsx", row.names=FALSE, sheetName="DEG")
write.xlsx(cbind(" "=rownames(res_ls[[4]]), as.data.frame(res_ls[[4]])), "brainMets_vs_nonBrainMets_19vs17_comp4.xlsx", row.names=FALSE, sheetName="allGenes", append=TRUE)



# expression of some genes? As we are seeing some of these genes amplified in brain mets,
write.table(rownames(res_ls[[4]]), "genes_241.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
geneList <- c("PPP2RA", "BRCA2", "CHEK2", "RAD54L", "CHEK1", "ATM", "BRCA1", "RAD51B", "CDK12", "RAD51C", "FANCL", "PALB2", "BARD1", "BRIP1")

res_ls[[4]][geneList[which(geneList %in% rownames(res_ls[[4]]))] , ]
# log2 fold change (MLE): Group2 Met/Brain vs Met/Non-Brain 
# Wald test p-value: Group2 Met/Brain vs Met/Non-Brain 
# DataFrame with 4 rows and 6 columns
       # baseMean log2FoldChange     lfcSE        stat     pvalue      padj
      # <numeric>      <numeric> <numeric>   <numeric>  <numeric> <numeric>
# BRCA2   85.7212   -0.020364545  0.309791 -0.06573639 0.94758770 0.9790774
# ATM    119.4644   -0.000232294  0.204159 -0.00113781 0.99909216 0.9990922
# FANCL  263.0963   -0.661845958  0.238898 -2.77041281 0.00559853 0.0449748
# BRIP1  168.7405   -0.959159433  0.329923 -2.90722571 0.00364650 0.0349856

res_ls[[1]][geneList[which(geneList %in% rownames(res_ls[[1]]))] , ]
# log2 fold change (MLE): Group2 Treated vs Non-treated 
# Wald test p-value: Group2 Treated vs Non.treated 
# DataFrame with 4 rows and 6 columns
       # baseMean log2FoldChange     lfcSE      stat      pvalue        padj
      # <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
# BRCA2   51.7840       0.909605  0.192133  4.734243 2.19874e-06 7.06529e-06
# ATM    257.7853      -0.542925  0.114479 -4.742585 2.11009e-06 6.87204e-06
# FANCL  326.4946      -0.119174  0.119782 -0.994926 3.19772e-01 3.72295e-01
# BRIP1   84.8725       1.616352  0.230060  7.025777 2.12879e-12 2.57703e-11

res_ls[[2]][geneList[which(geneList %in% rownames(res_ls[[2]]))] , ]
# log2 fold change (MLE): Group2 Treated vs Non-treated 
# Wald test p-value: Group2 Treated vs Non.treated 
# DataFrame with 4 rows and 6 columns
       # baseMean log2FoldChange     lfcSE      stat      pvalue        padj
      # <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
# BRCA2   46.3407       1.065506  0.242622   4.39162 1.12509e-05 6.02546e-05
# ATM    272.7819      -0.626905  0.149484  -4.19378 2.74341e-05 1.32233e-04
# FANCL  321.7822      -0.444040  0.141020  -3.14876 1.63964e-03 4.34234e-03
# BRIP1   63.6109       1.138702  0.293018   3.88611 1.01862e-04 3.83574e-04


# RLE plot

countData_allsmp = counts(norm.dat$set)[1:241, ]
colData_allsmp <- pData(norm.dat$set)

dds_allsmp <- DESeqDataSetFromMatrix(countData = countData_allsmp, colData = colData_allsmp, design = ~ W_1)

dds_allsmp <- DESeq(dds_allsmp)
dds_allsmp_vsd <- assay(varianceStabilizingTransformation(dds_allsmp, blind=FALSE))

table(pData(norm.dat$set)$Nanostring_ID == colnames(dds_allsmp_vsd))
# TRUE 
 # 198 

table(pData(norm.dat$set)$Nanostring_ID == colnames(raw))
# TRUE 
 # 198 
table(rownames(pData(norm.dat$set))== colnames(raw))
# TRUE 
 # 198 
colData_20_35_140_ordered <- colData_20_35_140[order(colData_20_35_140$Brain_vs_nonBrain2),]

ordered_group <- pData(norm.dat$set)[order(pData(norm.dat$set)$Group),]
ordered_raw <- raw[, ordered_group$Nanostring_ID]
table(rownames(ordered_group)== colnames(ordered_raw))
# TRUE 
 # 198 
 
ordered_dds_allsmp_vsd <- dds_allsmp_vsd[, ordered_group$Nanostring_ID]
table(rownames(ordered_group)== colnames(ordered_dds_allsmp_vsd))
# TRUE 
 # 198 
   
raw_rle <- newSeqExpressionSet(as.matrix(ordered_raw), phenoData=AnnotatedDataFrame(data.frame(conditions= factor(ordered_group $Group), row.names =rownames(ordered_group))))

 
rle_col<-rep.int("#fe7e0c", length(ordered_group $Group))
rle_col[ordered_group$Group =="Met/Non-Brain"] <- "#2390fd"
rle_col[ordered_group$Group =="Primary Prostate (Untreated)"] <- "#8eec91"
rle_col[ordered_group$Group =="Prostate - Treatment Unknown"] <- "#8ca3c7"
rle_col[ordered_group$Group =="Prostate (Treated)"] <- "#f3e00b"

pdf("RLE_2.pdf", height = 16, width=24)
# par(mfrow = c(2,1), mar = c(4,4,4,9) , xpd=TRUE) # Add extra space to right of plot area; change clipping to figure
par(mfrow = c(2,1), mar = c(9,5,3,4), mgp=c(3, 1.5, 0))
plotRLE(raw_rle, outline=FALSE, col= rle_col, main="Raw", cex.axis=2.5, cex.main=2.6, names=NULL, ylim=c(-6, 6))
plotRLE(ordered_dds_allsmp_vsd, col= rle_col, outline=FALSE, main="RUVseq", cex.axis=2.6, cex.main=3, names=NULL, ylim=c(-1, 1))
legend("bottom", inset = c(0, -0.3), legend=c("Brain","Non-brain", "Prostate (Untreated)", "Prostate (Treatment Unknown)", "Prostate (Treated)"), col=c("#fe7e0c", "#2390fd", "#8eec91", "#8ca3c7", "#f3e00b"), pch=c(20,20), cex=1.7, box.lty=0, xpd = TRUE, horiz=TRUE, bg="#efefef")
dev.off()

# volcano plot 

for_volcano_dds_ls <- list()
up_sig_indices_ls <- list()
down_sig_indices_ls <- list()

for_volcano_dds_ls_1.5logFC <- list()
up_sig_indices_ls_1.5logFC <- list()
down_sig_indices_ls_1.5logFC <- list()

for (i in 1:9){
	for_volcano_dds_ls[[i]] <- data.frame('log2FoldChange'= res_ls[[i]]$log2FoldChange,
						  'padj'= res_ls[[i]]$padj,
						  'Trend' = rep('Not sig', length(res_ls[[i]]$log2FoldChange)))
	up_sig_indices_ls[[i]] <- intersect(which(for_volcano_dds_ls[[i]]$log2FoldChange >1), which(for_volcano_dds_ls[[i]]$padj <0.05))
	down_sig_indices_ls[[i]] <- intersect(which(for_volcano_dds_ls[[i]]$log2FoldChange < (-1)), which(for_volcano_dds_ls[[i]]$padj <0.05))
	for_volcano_dds_ls[[i]][up_sig_indices_ls[[i]], 'Trend'] <- 'Up'
	for_volcano_dds_ls[[i]][down_sig_indices_ls[[i]], 'Trend'] <- 'Down'

	for_volcano_dds_ls[[i]]$Trend <- as.factor(for_volcano_dds_ls[[i]]$Trend)
	for_volcano_dds_ls[[i]]$Lg_padj <- -log10(for_volcano_dds_ls[[i]]$padj)
	for_volcano_dds_ls[[i]]$gene_symbol <- rownames(res_ls[[i]])

	# differentially expressed  label
	for_volcano_dds_ls[[i]]$delabel <- NA
	for_volcano_dds_ls[[i]]$delabel[for_volcano_dds_ls[[i]]$Trend != "Not sig"] <- for_volcano_dds_ls[[i]]$gene_symbol[for_volcano_dds_ls[[i]]$Trend != "Not sig"]
	
	# log2FC 1.5 cutoff
	for_volcano_dds_ls_1.5logFC[[i]] <- data.frame('log2FoldChange'= res_ls[[i]]$log2FoldChange,
						  'padj'= res_ls[[i]]$padj,
						  'Trend' = rep('Not sig', length(res_ls[[i]]$log2FoldChange)))
	up_sig_indices_ls_1.5logFC[[i]] <- intersect(which(for_volcano_dds_ls_1.5logFC[[i]]$log2FoldChange >1.5), which(for_volcano_dds_ls_1.5logFC[[i]]$padj <0.05))
	down_sig_indices_ls_1.5logFC[[i]] <- intersect(which(for_volcano_dds_ls_1.5logFC[[i]]$log2FoldChange < (-1.5)), which(for_volcano_dds_ls_1.5logFC[[i]]$padj <0.05))
	for_volcano_dds_ls_1.5logFC[[i]][up_sig_indices_ls_1.5logFC[[i]], 'Trend'] <- 'Up'
	for_volcano_dds_ls_1.5logFC[[i]][down_sig_indices_ls_1.5logFC[[i]], 'Trend'] <- 'Down'

	for_volcano_dds_ls_1.5logFC[[i]]$Trend <- as.factor(for_volcano_dds_ls_1.5logFC[[i]]$Trend)
	for_volcano_dds_ls_1.5logFC[[i]]$Lg_padj <- -log10(for_volcano_dds_ls_1.5logFC[[i]]$padj)
	for_volcano_dds_ls_1.5logFC[[i]]$gene_symbol <- rownames(res_ls[[i]])

	# differentially expressed  label
	for_volcano_dds_ls_1.5logFC[[i]]$delabel <- NA
	for_volcano_dds_ls_1.5logFC[[i]]$delabel[for_volcano_dds_ls_1.5logFC[[i]]$Trend != "Not sig"] <-for_volcano_dds_ls_1.5logFC[[i]]$gene_symbol[for_volcano_dds_ls_1.5logFC[[i]]$Trend != "Not sig"]	
}

names(for_volcano_dds_ls) <- c("comp1", "comp2", "comp3", "comp4", "comp5", "comp6", "comp7", "comp8", "comp9")
names(up_sig_indices_ls) <- c("comp1", "comp2", "comp3", "comp4", "comp5", "comp6", "comp7", "comp8", "comp9")
names(down_sig_indices_ls) <- c("comp1", "comp2", "comp3", "comp4", "comp5", "comp6", "comp7", "comp8", "comp9")


names(for_volcano_dds_ls_1.5logFC) <- c("comp1", "comp2", "comp3", "comp4", "comp5", "comp6", "comp7", "comp8", "comp9")
names(up_sig_indices_ls_1.5logFC) <- c("comp1", "comp2", "comp3", "comp4", "comp5", "comp6", "comp7", "comp8", "comp9")
names(down_sig_indices_ls_1.5logFC) <- c("comp1", "comp2", "comp3", "comp4", "comp5", "comp6", "comp7", "comp8", "comp9")


library(ggrepel)
# show.legend=FALSE to remove a in legend

volcano_ls <- list()


for (i in c(4, 8)){
# for (i in 1:9){
for (i in 8){	
	volcano_ls[[i]] <- ggplot(data= for_volcano_dds_ls[[i]], aes(x=log2FoldChange, y= Lg_padj, col= Trend, label=delabel)) + 
    geom_point() + 
    theme_minimal() +
    theme(
    legend.title=element_blank(),
    axis.title.x = element_text( size = 14, hjust=0.5, vjust=-2),  # change the position of titels
    axis.title.y = element_text( size = 14, vjust=2),
    axis.text=element_text(size= 14),
    plot.margin = unit(c(1,1,0.6,1), "cm"),
    legend.text = element_text(size=12),
    # legend.key.size = unit(0.1, 'cm'), #change legend key size
    legend.key.height = unit(0.3, 'cm'), #change legend key height
    legend.key.width = unit(0.8, 'cm'), #change legend key width
    legend.position="bottom",
    # Background
    legend.background = element_rect(fill = "#f7f4f4", color="#f7f4f4")       # Border
    )+    # scale_color_manual(values=c('No'='black', 'Up'='red', 'Down'='blue','UpSec'='orange', 'DownSec'='green'))+
    scale_color_manual(values=c('Up'='red', 'Down'='blue', 'Not sig'='#8c8888'))+
    # geom_text()	+
        # labs(x="log2FoldChange", y='-log10padj')+
        labs(x=substitute(paste({log[2]}, "FC")), y=substitute(paste({-log[10]}, "FDR")))+
        geom_vline(xintercept = c(-1, 0, 1), lty=2, size=I(0.4), colour='grey11')+
	    # geom_hline(yintercept = c(-log(x=0.05, base=10), -log(x=0.01, base=10)), lty=2, size=I(0.4), colour='grey11')+
	    geom_hline(yintercept = c(-log(x=0.05, base=10)), lty=2, size=I(0.4), colour='grey11')+
	    geom_text_repel(box.padding = 0.6, max.overlaps = 40, show.legend=FALSE) +
	    # xlim(-6, 6)	+   
	    ylim(0, 5)	+
	    scale_x_continuous(breaks=c(-6, -4, -2, 0, 2, 4, 6))+
	    annotate("text", x = 5, y = 1.1, label = "FDR = 0.05", size = 2.5) + # Label to horizontal cut-off line. y=1 (4 and 8), otherwise  0.3
        annotate("text", x = 1.2, y = 4.5, label = "2FC", size = 4, srt = 90) + # Label to vertical cut-off line. y=7 for i in 4 and 8, otherwise 40
        annotate("text", x = -1.2, y = 4.5, label = "-2FC", size = 4, srt = 90) # Label to vertical cut-off line.   
	ggsave(volcano_ls[[i]], file=paste0("volplot_comp", i, ".pdf"), width = 18, height = 18, units = "cm")	
} 



# heatmap

# transform ids

trans_ids <- as.data.frame(read_xlsx("../from_researcher_or_DL/Transformed IDs_Fig5.xlsx", sheet=2))
trans_ids_match <- trans_ids[match(colData_ls[[4]]$NewSampleName, trans_ids$NewSampleName),]
table(trans_ids_match$NewSampleName == colData_ls[[4]]$NewSampleName)

# TRUE 
  # 36 


rld_ls <- list()
expr_vst_ls <- list()
choose_gene_ls <- list() 
choose_gene_ls_1.5logFC <- list()
expr_vst_ls_DEG <- list()
expr_vst_ls_DEG_1.5logFC <- list()
z_score_matrix_ls <- list()
choose_matrix_ls <- list()
choose_matrix_ls_1.5logFC <- list()
topAnno_ls <- list()
heatmap_ls <- list()

# DEG
for (i in 1:9){
	choose_gene_ls[[i]] <- filter(for_volcano_dds_ls[[i]], for_volcano_dds_ls[[i]]$Trend != "Not sig")$gene_symbol
	choose_gene_ls_1.5logFC[[i]] <- filter(for_volcano_dds_ls_1.5logFC[[i]], for_volcano_dds_ls_1.5logFC[[i]]$Trend != "Not sig")$gene_symbol

	rld_ls[[i]] <- varianceStabilizingTransformation(dds_ls[[i]], blind=FALSE)
	expr_vst_ls[[i]] <- assay(rld_ls[[i]])
	expr_vst_ls_DEG[[i]] <- expr_vst_ls[[i]][choose_gene_ls[[i]], ]
	expr_vst_ls_DEG_1.5logFC[[i]] <- expr_vst_ls[[i]][choose_gene_ls_1.5logFC[[i]], ]

	z_score_matrix_ls[[i]] <- t(scale(t(expr_vst_ls[[i]])))
	choose_matrix_ls[[i]] <- z_score_matrix_ls[[i]][choose_gene_ls[[i]], ]
	choose_matrix_ls[[i]][choose_matrix_ls[[i]] > 2] =2
	choose_matrix_ls[[i]][choose_matrix_ls[[i]] < -2] =-2
	
	choose_matrix_ls_1.5logFC[[i]] <- z_score_matrix_ls[[i]][choose_gene_ls_1.5logFC[[i]], ]
	choose_matrix_ls_1.5logFC[[i]][choose_matrix_ls_1.5logFC[[i]] > 2] =2
	choose_matrix_ls_1.5logFC[[i]][choose_matrix_ls_1.5logFC[[i]] < -2] =-2

}

for (i in 1:9){
	print(table(colData_ls[[i]]$Group2))
	
	}


# add sample names
table(colData_ls[[4]]$Nanostring_ID == colnames(choose_matrix_ls[[4]]))

# TRUE 
  # 36 
choose_matrix_ls_comp4 <- choose_matrix_ls[[4]]
colnames(choose_matrix_ls_comp4) <- colData_ls[[4]]$NewSampleName

choose_matrix_ls_comp4_2 <- choose_matrix_ls[[4]]
colnames(choose_matrix_ls_comp4_2) <- trans_ids_match$`New ID`

par(mfrow = c(2,1), mar = c(4,8,4,9) , xpd=TRUE) # Add extra space to right of plot area; change clipping to figure
  
pdf (paste0("heatmap_comp", 4, "_3clusters_smpNmNew.pdf"), height = 7, width=11)
	# Heatmap(choose_matrix_ls[[4]], name='z-score', show_column_names=FALSE, show_row_names=TRUE, top_annotation= topAnno4, column_split=3, heatmap_legend_param = list(legend_height = unit(5.36, "cm"), legend_width = unit(5, "cm")))
	# Heatmap(choose_matrix_ls_comp4, name='z-score', show_column_names=TRUE, column_names_side="bottom", column_names_rot=45, show_row_names=TRUE, top_annotation= topAnno4, bottom_annotation= topAnno4_bottom, column_split=3, heatmap_legend_param = list(legend_height = unit(6.8, "cm"), legend_width = unit(5, "cm")))
	Heatmap(choose_matrix_ls_comp4_2, name='z-score', show_column_names=TRUE, column_names_side="bottom", column_names_rot=65, show_row_names=TRUE, top_annotation= topAnno4, bottom_annotation= topAnno4_bottom_2, column_split=3, heatmap_legend_param = list(legend_height = unit(6.8, "cm"), legend_width = unit(5, "cm")))

dev.off()



##### boxplot for TMPRSS2

TMPRSS2=c("TMPRSS2", "TMPRSS2-ERG")
norm_counts_dds_ls <- list()
norm_counts_dds_ls_TMPRSS2 <- list()
vsd_ls <- list()
vsd_ls_TERG <- list()
norm_counts_dds_ls_TMPRSS2_df <- list()
norm_counts_dds_ls_TMPRSS2_df_mod <- list()
stat.test.TMPRSS2_ls <- list()

# for (i in 1:9){
# for (i in c(1, 2, 3, 7, 9)){
# for (i in c(4)){
# for (i in c(5, 6)){	
for (i in c(8)){
	norm_counts_dds_ls[[i]] <- counts(dds_ls[[i]], normalized=TRUE)
	norm_counts_dds_ls_TMPRSS2[[i]] <- norm_counts_dds_ls[[i]][TMPRSS2,]
	vsd_ls[[i]] <- assay(rld_ls[[i]])
	vsd_ls_TERG[[i]] <- vsd_ls[[i]][TMPRSS2, ]
	# print(table(colnames(vsd_ls_TERG[[i]]) == colData_ls[[i]]$Nanostring_ID))
	# print(res_ls[[i]][TMPRSS2,])
	norm_counts_dds_ls_TMPRSS2_df[[i]] <- gather(cbind(t(norm_counts_dds_ls_TMPRSS2[[i]]), colData_ls[[i]][-13,]), TMPRSS2, norm_counts, 1:2) # colData_ls[[i]][-13,] only for comp8
	norm_counts_dds_ls_TMPRSS2_df[[i]]$TMPRSS2 <- as.factor(norm_counts_dds_ls_TMPRSS2_df[[i]]$TMPRSS2)
	norm_counts_dds_ls_TMPRSS2_df[[i]]$LogNormCounts <- log2(norm_counts_dds_ls_TMPRSS2_df[[i]]$norm_counts)
	norm_counts_dds_ls_TMPRSS2_df_mod[[i]] <- norm_counts_dds_ls_TMPRSS2_df[[i]] %>% mutate(Gene = recode(TMPRSS2, "TMPRSS2" = "TMPRSS2", "TMPRSS2-ERG" = "TMPRSS2-ERG"))
	
	stat.test.TMPRSS2_ls[[i]] <- norm_counts_dds_ls_TMPRSS2_df_mod[[i]] %>%
   		group_by(Gene) %>%
   		# t_test(data =., norm_counts ~ Group2, ref.group = "Non-treated") %>%
   		# t_test(data =., norm_counts ~ Group2, ref.group = "Met/Non-Brain") %>%
   		# t_test(data =., norm_counts ~ Group2, ref.group = "Non-metastatic") %>%
   		t_test(data =., norm_counts ~ Group2, ref.group = "Dural") %>%
   		adjust_pvalue(method = "BH") %>%
   		add_significance("p.adj")
     
}

# replace t test padj value by deseq2 padj for TMPRSS2 and TMPRSS2D, don't replace the stats and others
stat.test.TMPRSS2_ls[[1]][1,10] <- 0.58
stat.test.TMPRSS2_ls[[1]][2,10] <- 0.33 

stat.test.TMPRSS2_ls[[2]][1,10] <- 0.85
stat.test.TMPRSS2_ls[[2]][2,10] <- 0.88 

stat.test.TMPRSS2_ls[[3]][1,10] <- 0.79
stat.test.TMPRSS2_ls[[3]][2,10] <- 0.34

stat.test.TMPRSS2_ls[[4]][1,10] <- 0.79
stat.test.TMPRSS2_ls[[4]][2,10] <- 0.32

stat.test.TMPRSS2_ls[[5]][1,10] <- 0.79
stat.test.TMPRSS2_ls[[5]][2,10] <- 0.61

stat.test.TMPRSS2_ls[[6]][1,10] <- 0.64
stat.test.TMPRSS2_ls[[6]][2,10] <- 0.36

stat.test.TMPRSS2_ls[[7]][1,10] <- 0.27
stat.test.TMPRSS2_ls[[7]][2,10] <- 0.22

stat.test.TMPRSS2_ls[[8]][1,10] <- 0.99
stat.test.TMPRSS2_ls[[8]][2,10] <- 0.90

stat.test.TMPRSS2_ls[[9]][1,10] <- 0.94
stat.test.TMPRSS2_ls[[9]][2,10] <- 0.50

table(colnames(vsd_ls_TERG[[1]]) == colData_ls[[1]]$Nanostring_ID)
# TRUE 
 # 175
table(colnames(norm_counts_dds_ls_TMPRSS2[[1]]) == colData_ls[[1]]$Nanostring_ID)
# TRUE 
 # 175

res_ls[[i]][TMPRSS2,]

# log2 fold change (MLE): Group2 Treated vs Non-treated 
# Wald test p-value: Group2 Treated vs Non.treated 
# DataFrame with 2 rows and 6 columns
             # baseMean log2FoldChange     lfcSE      stat    pvalue      padj
            # <numeric>      <numeric> <numeric> <numeric> <numeric> <numeric>
# TMPRSS2      3022.416       0.210160  0.349090  0.602022  0.547160  0.583476
# TMPRSS2-ERG   110.759       0.635579  0.590094  1.077080  0.281445  0.330869
# log2 fold change (MLE): Group2 Treated vs Non-treated 
# Wald test p-value: Group2 Treated vs Non.treated 
# DataFrame with 2 rows and 6 columns
             # baseMean log2FoldChange     lfcSE      stat    pvalue      padj
            # <numeric>      <numeric> <numeric> <numeric> <numeric> <numeric>
# TMPRSS2      2977.345       0.117235  0.420598  0.278734  0.780449  0.851078
# TMPRSS2-ERG   105.444       0.154606  0.770832  0.200570  0.841035  0.885107
# log2 fold change (MLE): Group2 Treated vs Non-treated 
# Wald test p-value: Group2 Treated vs Non.treated 
# DataFrame with 2 rows and 6 columns
             # baseMean log2FoldChange     lfcSE      stat    pvalue      padj
            # <numeric>      <numeric> <numeric> <numeric> <numeric> <numeric>
# TMPRSS2      2885.380       0.150221  0.450862  0.333186  0.738994  0.795078
# TMPRSS2-ERG   111.172       0.803589  0.744055  1.080013  0.280136  0.349808
# log2 fold change (MLE): Group2 Met/Brain vs Met/Non-Brain 
# Wald test p-value: Group2 Met/Brain vs Met/Non-Brain 
# DataFrame with 2 rows and 6 columns
             # baseMean log2FoldChange     lfcSE      stat    pvalue      padj
            # <numeric>      <numeric> <numeric> <numeric> <numeric> <numeric>
# TMPRSS2      3609.045       0.298604  0.625795   0.47716  0.633248  0.798908
# TMPRSS2-ERG   125.345       1.357058  0.923380   1.46966  0.141653  0.325127
# log2 fold change (MLE): Group2 Metastatic vs Non-metastatic 
# Wald test p-value: Group2 Metastatic vs Non-metastatic 
# DataFrame with 2 rows and 6 columns
             # baseMean log2FoldChange     lfcSE      stat    pvalue      padj
            # <numeric>      <numeric> <numeric> <numeric> <numeric> <numeric>
# TMPRSS2      2943.029       0.108866  0.348339  0.312529  0.754638  0.790730
# TMPRSS2-ERG   116.269       0.351296  0.559322  0.628075  0.529955  0.614034
# log2 fold change (MLE): Group2 Metastatic vs Non-metastatic 
# Wald test p-value: Group2 Metastatic vs Non-metastatic 
# DataFrame with 2 rows and 6 columns
             # baseMean log2FoldChange     lfcSE      stat    pvalue      padj
            # <numeric>      <numeric> <numeric> <numeric> <numeric> <numeric>
# TMPRSS2      2939.894       0.193762  0.340032  0.569835   0.56879  0.640553
# TMPRSS2-ERG   118.872       0.585763  0.546977  1.070910   0.28421  0.360498
# log2 fold change (MLE): Group2 Treated vs Non-treated 
# Wald test p-value: Group2 Treated vs Non.treated 
# DataFrame with 2 rows and 6 columns
             # baseMean log2FoldChange     lfcSE      stat    pvalue      padj
            # <numeric>      <numeric> <numeric> <numeric> <numeric> <numeric>
# TMPRSS2     2812.8324      -0.517931  0.434391  -1.19232  0.233137  0.276779
# TMPRSS2-ERG   97.1328      -0.939369  0.716820  -1.31047  0.190038  0.228996
# log2 fold change (MLE): Group2 Parenchymal vs Dural 
# Wald test p-value: Group2 Parenchymal vs Dural 
# DataFrame with 2 rows and 6 columns
             # baseMean log2FoldChange     lfcSE      stat    pvalue      padj
            # <numeric>      <numeric> <numeric> <numeric> <numeric> <numeric>
# TMPRSS2      3720.011      0.0296402  0.525555  0.056398  0.955025  0.994793
# TMPRSS2-ERG   113.454     -1.0309386  1.332131 -0.773902  0.438989  0.904761
# log2 fold change (MLE): Group2 Treated vs Non-treated 
# Wald test p-value: Group2 Treated vs Non.treated 
# DataFrame with 2 rows and 6 columns
             # baseMean log2FoldChange     lfcSE      stat    pvalue      padj
            # <numeric>      <numeric> <numeric> <numeric> <numeric> <numeric>
# TMPRSS2      2980.026      0.0194022  0.309392 0.0627109  0.949997  0.949997
# TMPRSS2-ERG   103.157      0.3515165  0.500453 0.7023966  0.482432  0.503316


#####
#box plot 
# switch the color for levels
norm_counts_dds_ls_TMPRSS2_df_mod[[1]]$Group2 <- factor(norm_counts_dds_ls_TMPRSS2_df_mod[[1]]$Group2, levels=c("Treated", "Non-treated"))
norm_counts_dds_ls_TMPRSS2_df_mod[[2]]$Group2 <- factor(norm_counts_dds_ls_TMPRSS2_df_mod[[2]]$Group2, levels=c("Treated", "Non-treated"))
norm_counts_dds_ls_TMPRSS2_df_mod[[3]]$Group2 <- factor(norm_counts_dds_ls_TMPRSS2_df_mod[[3]]$Group2, levels=c("Treated", "Non-treated"))

norm_counts_dds_ls_TMPRSS2_df_mod[[4]]$Group2 <- factor(norm_counts_dds_ls_TMPRSS2_df_mod[[4]]$Group2, levels=c("Met/Brain", "Met/Non-Brain"))

norm_counts_dds_ls_TMPRSS2_df_mod[[5]]$Group2 <- factor(norm_counts_dds_ls_TMPRSS2_df_mod[[5]]$Group2, levels=c("Metastatic", "Non-metastatic"))
norm_counts_dds_ls_TMPRSS2_df_mod[[6]]$Group2 <- factor(norm_counts_dds_ls_TMPRSS2_df_mod[[6]]$Group2, levels=c("Metastatic", "Non-metastatic"))

norm_counts_dds_ls_TMPRSS2_df_mod[[7]]$Group2 <- factor(norm_counts_dds_ls_TMPRSS2_df_mod[[7]]$Group2, levels=c("Treated", "Non-treated"))

norm_counts_dds_ls_TMPRSS2_df_mod[[8]]$Group2 <- factor(norm_counts_dds_ls_TMPRSS2_df_mod[[8]]$Group2, levels=c("Parenchymal", "Dural"))

norm_counts_dds_ls_TMPRSS2_df_mod[[9]]$Group2 <- factor(norm_counts_dds_ls_TMPRSS2_df_mod[[9]]$Group2, levels=c("Treated", "Non-treated"))



TMPRSS2_normCounts_boxplot_ls <- list()

# # First specify the y.position of each comparison, different y position for each facet
# stat.test.TMPRSS2 <- stat.test.TMPRSS2 %>%
 # mutate(y.position = c(7500, 45, 2))
   
# norm_counts
# Expanding the y-axis using the ggplot2 function expansion()

for (i in 1:9){

TMPRSS2_normCounts_boxplot_ls[[i]] <- ggboxplot(norm_counts_dds_ls_TMPRSS2_df_mod[[i]], x = "Group2", y = "LogNormCounts",
          color = "Group2", palette = "jco",
          add = "jitter",
          facet.by = "Gene", short.panel.labs = FALSE, scales="free_y")
          
TMPRSS2_normCounts_boxplot_ls[[i]] <- TMPRSS2_normCounts_boxplot_ls[[i]] +

			# stat_compare_means(label = "p.format", ref.group="unCor",  vjust=3, hjust=0)+ # two-sides test
		# stat_compare_means(label = "p.format", ref.group="Non_Brain", vjust=3, hjust=-0.2, size=7)+
		# change p position
		stat_pvalue_manual( 
    				stat.test.TMPRSS2_ls[[i]], 
    				y.position = 15,
    				# y.position = 13,    				
    				label.size = 6,
  		 		bracket.size = 0.5,
    				label = "p.adj = {p.adj}"
    			) +
	   scale_y_continuous(expand = expansion(mult = c(0.2, 0.2)))  +  			
		labs(x ="", y = "mRNA expression log2(Normalized counts)") +
		theme_light()+
		theme(
		strip.placement = "outside",                      # Place facet labels outside x axis labels.
	    strip.text=element_text(size=15, color="black"),	
        strip.background = element_rect(fill = "grey90"),  # Make facet label background white.
        # axis.text.x = element_text(margin=margin(5,5,0,5,"pt")),
        # axis.text.y = element_text(margin=margin(5,5,5,5,"pt")),
        # plot.margin = unit(c(1,1,1,1), "cm"),
        # axis.title = element_blank())+   # Remove x and y axis titles.
       # axis.title.x = element_text( size = 20),  # change the position of titles
       axis.title.y = element_text( size = 22),
       axis.text=element_text(size=20), # axis ticks word size
       plot.title = element_text(hjust = 0.5, vjust=2, size=30, face="bold"),
       legend.key.size = unit(1.6, 'cm'), #change legend key size
        # legend.key.height = unit(1, 'cm'), #change legend key height
        # legend.key.width = unit(1, 'cm'), #change legend key width
        # legend.title = element_text(size=20), #change legend title font size
        legend.title=element_blank(),
        legend.text = element_text(size=15),
        legend.position="right"
        # legend.direction="horizontal"
       )  
ggsave(TMPRSS2_normCounts_boxplot_ls[[i]], file=paste0("TM_boxplot_comp", i, "_18smps.pdf"), width = 35, height = 20, units = "cm")	
}	
# # # # # # # # # #########



#################

          
############################

# computing AR score
### 0.0 function to compute the score
computeScore = function( expr, referenceCol= 2 , startDataCol = 4) {
  cor.res = apply( expr[,-c(1:(startDataCol-1))], 2, cor.test, expr[, referenceCol] ) ## 
  cor.p = unlist(lapply(cor.res, "[", "p.value"))
  cor.r = unlist(lapply(cor.res, "[", "estimate"))
  names(cor.p) = gsub( ".p.value", "", names(cor.p) )
  names(cor.r) = gsub( ".estimate.cor", "", names(cor.r))
  return( list( score=cor.r, pval=cor.p))
}

df.expr.ls <- list()
work.xprl.ls <- list()
ar.score.res.ls <- list()
df.ls <- list()
df.expr.nepc.ls <- list()
work.xprl.nepc.ls <- list()
nepc.score.res.ls <- list()
nepc.score.res.df.ls <- list()
df.annotation.ls <- list()

## 0.3 reference datasets
lncap.ref = read.table("../../CALGB_Alliance/nanostring_Andrea/Resources/LnCaP.reference.AR_score.csv", sep="\t")
nepc.ref = read.table("../../CALGB_Alliance/nanostring_Andrea/Resources/Integrated_NEPC_Score_reference_vector.txt", sep="\t", skip = 1)

# should be normalized counts
          
# table(rownames(colData_31_match) == rownames(t(normalized_counts_dds_violin)))
# TRUE 
  # 31
# table(rownames(colData_25) == rownames(t(normalized_counts_dds_violin)))
# TRUE 
  # 25  
  
for (i in 1:9){
	df.expr.ls[[i]] <- norm_counts_dds_ls[[i]]
	df.expr.ls[[i]] <- cbind(head(fData, n=241), df.expr.ls[[i]])
	## must be TRUE (if not true, merge anyway, older codeset missed 1 gene, "GLRA2")
	# print(all( lncap.ref$V1 %in% df.expr.ls[[i]]$Gene & !is.null( df.expr.ls[[i]]$Gene )))
	work.xprl.ls[[i]] = merge( lncap.ref, df.expr.ls[[i]], by.x="V1", by.y="Gene")
	# print(dim(work.xprl.ls[[i]]))
	# ar scores
	ar.score.res.ls[[i]] = computeScore( work.xprl.ls[[i]] )
	
	if( all(names(ar.score.res.ls[[i]]$score)==names(ar.score.res.ls[[i]]$pval))) df.ls[[i]] = as.data.frame( cbind( ar.score.res.ls[[i]]$score, ar.score.res.ls[[i]]$pval)) else stop("not the same names")
	colnames(df.ls[[i]]) = c("AR.score", "AR.score.pvalue")
	# print(head(df.ls[[i]]))
	# nepc
	# print(all( nepc.ref$V1 %in% df.expr.ls[[i]]$Gene & !is.null( df.expr.ls[[i]]$Gene ) )) ## must be TRUE (but proceed and see what happens)
	
	df.expr.nepc.ls[[i]] <- norm_counts_dds_ls[[i]]
	df.expr.nepc.ls[[i]] <- cbind(head(fData, n=241), df.expr.nepc.ls[[i]])
	work.xprl.nepc.ls[[i]] = merge( nepc.ref, df.expr.nepc.ls[[i]], by.x="V1", by.y="Gene")
	print(dim(work.xprl.nepc.ls[[i]]))
	
	startData = 4 ## column index of starting data
	nepc.score.res.ls[[i]] = computeScore( work.xprl.nepc.ls[[i]], startDataCol = startData )
	nepc.score.res.df.ls[[i]] = as.data.frame( cbind(nepc.score.res.ls[[i]]$score, nepc.score.res.ls[[i]]$pval) )
	colnames(nepc.score.res.df.ls[[i]])=c("NEPC.score", "NEPC.score.pvalue") ## creating a dataframe from list

	if( nrow(df.ls[[i]]) != nrow(nepc.score.res.df.ls[[i]] ) ) {
  		df.ls[[i]]=merge( df.ls[[i]], nepc.score.res.df.ls[[i]], by="row.names", all.x=T)
		} else {
  			if( all( rownames(df.ls[[i]])==names(nepc.score.res.ls[[i]]$score) )) df.ls[[i]]$NEPC.score        = nepc.score.res.ls[[i]]$score else stop("score: not the same names")
  			if( all( rownames(df.ls[[i]])==names(nepc.score.res.ls[[i]]$pval ) )) df.ls[[i]]$NEPC.score.pvalue = nepc.score.res.ls[[i]]$pval  else stop( "pval: not the same names") 
	} 
	# print(head(df.ls[[i]]))
	df.ls[[i]]$Nanostring_ID=rownames(df.ls[[i]])
	#colnames(df)[1]="NanostringID"
	# print(head(df.ls[[i]]))

	df.annotation.ls[[i]] <- colData_ls[[i]]
	df.ls[[i]] = merge(df.ls[[i]], df.annotation.ls[[i]], by="Nanostring_ID", all.x=T)
	print(head(df.ls[[i]]))
	print(duplicated(colnames(df.ls[[i]])))

}  


# AR NEPC score boxplot
df_AR_NEPC_ls <- list()
df_AR_NEPC_mod_ls <- list()
AR_NEPC_boxplot_ls <- list()

# for (i in 1:9){
# for (i in c(1, 2, 3, 7, 9)){
for (i in c(4)){
# for (i in c(5, 6)){	
# for (i in c(8)){	
df_AR_NEPC_ls[[i]] <- gather(df.ls[[i]], Score, Scores, c(2, 4))
df_AR_NEPC_mod_ls[[i]] <- df_AR_NEPC_ls[[i]] %>% mutate(Score = recode(Score, "AR.score" = "AR", "NEPC.score" = "NEPC"))

# switch the color for levels
# df_AR_NEPC_mod_ls[[i]]$Group2 <- factor(df_AR_NEPC_mod_ls[[i]]$Group2, levels=c("Treated", "Non-treated"))
df_AR_NEPC_mod_ls[[i]]$Group2 <- factor(df_AR_NEPC_mod_ls[[i]]$Group2, levels=c("Met/Brain", "Met/Non-Brain"))
# df_AR_NEPC_mod_ls[[i]]$Group2 <- factor(df_AR_NEPC_mod_ls[[i]]$Group2, levels=c("Metastatic", "Non-metastatic"))
# df_AR_NEPC_mod_ls[[i]]$Group2 <- factor(df_AR_NEPC_mod_ls[[i]]$Group2, levels=c("Parenchymal", "Dural"))

AR_NEPC_boxplot_ls[[i]] <- ggboxplot(df_AR_NEPC_mod_ls[[i]], x = "Group2", y = "Scores",
          fill = "Group2", palette = c("#fe7e0c", "#2390fd"),
          add = "jitter", size=0.3, add.params = list(size = 0.7), 
          facet.by = "Score", short.panel.labs = FALSE, scales="free_y")
          
AR_NEPC_boxplot_ls[[i]] <- AR_NEPC_boxplot_ls[[i]] +

			# stat_compare_means(label = "p.format", ref.group="unCor",  vjust=3, hjust=0)+ # two-sides test
		# stat_compare_means(label = "p.format", ref.group="Non-treated", vjust=-1.8, hjust=1, size=6)+
		stat_compare_means(label = "p.format", ref.group="Met/Non-Brain", vjust=-1.8, hjust=-0.5, size=6)+
		# stat_compare_means(label = "p.format", ref.group="Non-metastatic", vjust=-1.8, hjust=1, size=6)+
		# stat_compare_means(label = "p.format", ref.group="Dural", vjust=-1.8, hjust=1, size=6)+
		# stat_compare_means(ref.group="Dural", vjust=-1.8, hjust=0, size=4)+
		# change p position
		# stat_pvalue_manual( 
    				# stat.test.TMPRSS2, 
    				# y.position = 14,
    				# # y.position = 13,    				
    				# label.size = 6,
  		 			# bracket.size = 0.5,
    				# label = "p.adj"
    			# ) +
	   scale_y_continuous(expand = expansion(mult = c(0.2, 0.2)))  + 
	   scale_x_discrete(breaks=c("Met/Brain", "Met/Non-Brain"), labels=c("Brain", "Non-brain"))+			
		labs(x ="", y = "Scores") +
		theme_light()+
		theme(
		strip.placement = "outside",                      # Place facet labels outside x axis labels.
	    strip.text=element_text(size=15, color="black"),	
        strip.background = element_rect(fill = "grey90"),  # Make facet label background white.
        # axis.text.x = element_text(margin=margin(5,5,0,5,"pt")),
        # axis.text.y = element_text(margin=margin(5,5,5,5,"pt")),
        # plot.margin = unit(c(1,1,1,1), "cm"),
        # axis.title = element_blank())+   # Remove x and y axis titles.
        # axis.title.x = element_text( size = 20),  # change the position of titles
        axis.title.y = element_text( size = 22),
        axis.text=element_text(size=20), # axis ticks word size
        plot.title = element_text(hjust = 0.5, vjust=2, size=30, face="bold"),
        legend.key.size = unit(0.8, 'cm'), #change legend key size
        # legend.key.height = unit(1, 'cm'), #change legend key height
        # legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=12, face="bold"), #change legend title font size
        # legend.title=element_blank(),
        legend.text = element_text(size=12),
        legend.position="right"
        # legend.direction="horizontal"
       )  

AR_NEPC_boxplot_ls[[i]] <- AR_NEPC_boxplot_ls[[i]]+ guides(fill=guide_legend(title="Meta_site"))
ggsave(AR_NEPC_boxplot_ls[[i]], file=paste0("AR_NEPC_boxplot_comp", i, "_2.pdf"), width = 30, height = 15, units = "cm")	
	
}

# comp8, i=8
compare_means(Scores ~ Group2, data = df_AR_NEPC_mod_ls[[i]][df_AR_NEPC_mod_ls[[i]]$Score=="AR",])
# A tibble: 1 × 8
  # .y.    group1      group2     p p.adj p.format p.signif method  
  # <chr>  <chr>       <chr>  <dbl> <dbl> <chr>    <chr>    <chr>   
# 1 Scores Parenchymal Dural  0.863  0.86 0.86     ns       Wilcoxon
compare_means(Scores ~ Group2, data = df_AR_NEPC_mod_ls[[i]][df_AR_NEPC_mod_ls[[i]]$Score=="NEPC",])
# A tibble: 1 × 8
  # .y.    group1      group2     p p.adj p.format p.signif method  
  # <chr>  <chr>       <chr>  <dbl> <dbl> <chr>    <chr>    <chr>   
# 1 Scores Parenchymal Dural      1     1 1        ns       Wilcoxon

# Correlate specific genes with IHC staining analysis

# all samples of IHC information are brain mets except 20180920_NEPC75_NEPC75_05


IHC <- c("NKX3-1", "AR", "FOLH1 (PSMA)", "CHGA", "SYP", "RB1", "ERG")

IHC_normCounts <- norm_counts_dds_ls[[4]][IHC,]
table(rownames(colData_ls[[4]]) == rownames(t(IHC_normCounts)))

# TRUE 
# 36 

IHC_normCounts_df <- gather(cbind(t(IHC_normCounts), colData_ls[[4]]), IHC, norm_counts, 1:7)

# table(rownames(colData_25) == IHC_normCounts_df[c(26:50),]$Nanostring_ID)
table(rownames(colData_ls[[4]]) == IHC_normCounts_df[,]$Nanostring_ID)

# TRUE 
# 252 
IHC_normCounts_df$IHC <- as.factor(IHC_normCounts_df$IHC)

#box plot 
IHC_normCounts_df$LogNormCounts <- log2(IHC_normCounts_df$norm_counts)

ls_counts_IHC_7 <- list()

colnames(IHC_normCounts_df)[19] <- "NKX3.1"

for (i in 1:7){
  ls_counts_IHC_7[[i]] <- IHC_normCounts_df[c((36*(i-1)+1):(36*(i-1)+ 36)), ][!is.na(IHC_normCounts_df[c((36*(i-1)+1):(36*(i-1)+ 36)), 18+i]), ]
}
names(ls_counts_IHC_7) <- IHC

ls_plot_IHC_7 <- list()

for (i in 1:7){
  ls_plot_IHC_7[[i]] <- ggboxplot(ls_counts_IHC_7_reLevel[[i]], x = colnames(ls_counts_IHC_7_reLevel[[i]])[18+i], y = "LogNormCounts",
                                color = colnames(ls_counts_IHC_7_reLevel[[i]])[18+i], palette = "jco",
                                add = "jitter",
                                add.params = list(size=7), 
                                size =1.5,
                                # facet.by = "clade", 
                                short.panel.labs = FALSE, 
                                scales="free_y")
  

  ls_plot_IHC_7[[i]] <- ls_plot_IHC_7[[i]] +
    
    # stat_compare_means(label = "p.format", ref.group="unCor",  vjust=3, hjust=0)+ # two-sides test
    # stat_compare_means(label = "p.format", ref.group="Non_Brain", vjust=3, hjust=-0.2, size=7)+
    # change p position
    # stat_pvalue_manual( 
    #   stat.test.TMPRSS2, 
    #   y.position = 14,
    #   # y.position = 13,    				
    #   label.size = 6,
    #   bracket.size = 0.5,
    #   label = "p.adj = {p.adj}"
  # ) +
  # scale_y_continuous(expand = expansion(mult = c(0.2, 0.2)))  +  			
  labs(x="", y = "mRNA expression log2(Normalized counts)", title=paste0("IHC expression plot (", IHC[i], ")")) +
    theme_light()+
    # theme_linedraw()+
    theme(
      # strip.placement = "outside",                      # Place facet labels outside x axis labels.
      # strip.text=element_text(size=15, color="black"),	
      # strip.background = element_rect(fill = "grey90"),  # Make facet label background white.
      # axis.text.x = element_text(margin=margin(5,5,0,5,"pt")),
      # axis.text.y = element_text(margin=margin(5,5,5,5,"pt")),
      plot.margin = unit(c(3,1,2,3), "cm"),
      # axis.title = element_blank())+   # Remove x and y axis titles.
      # axis.title.x = element_text( size = 20),  # change the position of titles
      axis.title.y = element_text( size = 40),
      axis.text=element_text(size=40), # axis ticks word size
      # axis.text.x = element_text(angle=45, vjust=1, hjust=1),
      plot.title = element_text(hjust = 0.5, vjust=2, size=45, face="bold"),
      legend.key.size = unit(3, 'cm'), #change legend key size
      # legend.key.height = unit(1, 'cm'), #change legend key height
      # legend.key.width = unit(1, 'cm'), #change legend key width
      legend.title = element_text(size=35), #change legend title font size
      # legend.title=element_blank(),
      legend.text = element_text(size=30),
      legend.position="right"
      # legend.direction="horizontal"
    )  
  ggsave(paste0("IHC_normCounts_boxplot_comp4_", IHC[i], ".pdf"), ls_plot_IHC_7[[i]], height=15, width=22, limitsize = FALSE)
  
}


# need to re-level, switch the color for level

IHC_reLevel <- c( "CHGA", "SYP", "ERG")

ls_counts_IHC_7_reLevel <- ls_counts_IHC_7

ls_counts_IHC_7_reLevel[[4]]$Chromogranin <- factor(ls_counts_IHC_7_reLevel[[4]]$Chromogranin, levels=c("Positive", "Negative", "Focal"))

ls_counts_IHC_7_reLevel[[5]]$Synaptophysin <- factor(ls_counts_IHC_7_reLevel[[5]]$Synaptophysin, levels=c("Positive", "Negative", "Focal"))

ls_counts_IHC_7_reLevel[[7]]$ERG <- factor(ls_counts_IHC_7_reLevel[[7]]$ERG, levels=c("Positive", "Negative"))



# write colData 
colData_ls_out <- list()

for (i in 1:9){
	colData_ls_out[[i]]<- pData_match[match(colData_ls[[i]]$Nanostring_ID, pData_match$Nanostring_ID), -29]
}

names(colData_ls_out) <- c("comp1", "comp2", "comp3", "comp4", "comp5", "comp6", "comp7", "comp8", "comp9")


for (i in 1:9){
	
	write.xlsx(colData_ls_out[[i]], "anno_8cate_9comps_199NanoIds_092722.xlsx", sheetName=names(colData_ls_out)[i], append = TRUE, showNA = FALSE, row.names = FALSE)
}
	
	




