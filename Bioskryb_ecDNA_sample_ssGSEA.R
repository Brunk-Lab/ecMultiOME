###This notebook
library(GSVA)
packageVersion("GSVA")
library(msigdbr)
library(dplyr)
library(tidyr)


msig_hallmark <- msigdbr(species = "Homo sapiens", category = "H")
msig_hallmark <- msig_hallmark[, c("gs_name", "gene_symbol")]
gene_sets_list <- split(msig_hallmark$gene_symbol, msig_hallmark$gs_name)

#ssGSEA for Harper's/2nd plate NCIH2170
seurat_2170_h <- readRDS('/home/yue1118/Bioskryb_data_analysis/NCIH2170_plate2_Bioskryb_RNA_logTPM_ecDNA_counts_rds_03252025.rds')
X_2170_h <- as.matrix(seurat_2170_h@assays$RNA@data)

gsvaPar_2170_h <- gsvaParam(X_2170_h, gene_sets_list)
es_2170_h <- gsva(gsvaPar_2170_h, verbose=FALSE)
es_2170_h <- as.data.frame(t(es_2170_h))
View(es_2170_h)


rownames(es_2170_h) == colnames(seurat_2170_h)
es_2170_h$Cell <- seurat_2170_h$Cell
es_2170_h$ecDNA_CN <- seurat_2170_h$chr8_copy_number

plot(es_2170_h$ecDNA_CN, es_2170_h$HALLMARK_MYC_TARGETS_V1)
cor.test(es_2170_h$ecDNA_CN, es_2170_h$HALLMARK_MYC_TARGETS_V1)
plot(es_2170_h$ecDNA_CN, es_2170_h$HALLMARK_MYC_TARGETS_V2)
plot(es_2170_h$ecDNA_CN, es_2170_h$HALLMARK_MTORC1_SIGNALING)
cor.test(es_2170_h$ecDNA_CN, es_2170_h$HALLMARK_MTORC1_SIGNALING)

write.csv(es_2170_h, "/home/yue1118/Bioskryb_data_analysis/NCIH2170_Bioskryb_RNA_Harper_2nd_plate_ssGSEA_score_08292025.csv")

#ssGSEA for NCIH2170
seurat_2170 <- readRDS('/home/yue1118/Bioskryb_data_analysis/NCIH2170_Bioskryb_RNA_logTPM_ecDNA_counts_rds_03082025.rds')
X_2170 <- as.matrix(seurat_2170@assays$RNA@data)

gsvaPar_2170 <- gsvaParam(X_2170, gene_sets_list)
es_2170 <- gsva(gsvaPar_2170, verbose=FALSE)
es_2170 <- as.data.frame(t(es_2170))
View(es_2170)


rownames(es_2170) == colnames(seurat_2170)
es_2170$Cell <- seurat_2170$Cell
es_2170$ecDNA_CN <- seurat_2170$chr8_copy_number

plot(es_2170$HALLMARK_MTORC1_SIGNALING, es_2170$ecDNA_CN)
plot(es_2170$HALLMARK_GLYCOLYSIS, es_2170$ecDNA_CN)
plot(es_2170$HALLMARK_MYC_TARGETS_V2, es_2170$ecDNA_CN)
plot(es_2170$HALLMARK_MYC_TARGETS_V1, es_2170$ecDNA_CN)
plot(es_2170$HALLMARK_G2M_CHECKPOINT, es_2170$ecDNA_CN)

cor_results_2170 <- as.data.frame(sapply(es_2170[, -ncol(es_2170)], function(x) cor(x, es_2170[[ncol(es_2170)]], use = "complete.obs")))
View(cor_results_2170)

###To combine NCIH2170 PVT1 splice junction information###
df_SJ_2170high <- read.delim("/datastore/lbcfs/labs/brunk_lab/private/Bioskryb_2170_RNA/2170High_RNA/STAR_alignment_for_single_cells/H2170High_multiple_high_PVT1_chr8_SJs.tsv")
df_SJ_2170med <- read.delim("/datastore/lbcfs/labs/brunk_lab/private/Bioskryb_2170_RNA/2170Med_RNA/STAR_alignment_for_single_cells/H2170Med_multiple_high_PVT1_chr8_SJs.tsv")
df_SJ_2170low <- read.delim("/datastore/lbcfs/labs/brunk_lab/private/Bioskryb_2170_RNA/2170Low_RNA/STAR_alignment_for_single_cells/H2170Low_multiple_high_PVT1_chr8_SJs.tsv")
df_SJ_2170all <- rbind(df_SJ_2170high, df_SJ_2170med, df_SJ_2170low)

df_SJ_2170all$splice_isoform <- paste(df_SJ_2170all$chr, df_SJ_2170all$start, df_SJ_2170all$end, sep = "_")
df_SJ_2170all <- df_SJ_2170all[, c("Sample", "numUnique", "splice_isoform")]
df_SJ_2170wide <- pivot_wider(df_SJ_2170all,  names_from = splice_isoform, values_from = numUnique)
df_SJ_2170wide$Cell <- sub(".*-R-(.*)-_S.*", "\\1", df_SJ_2170wide$Sample)

View(df_SJ_2170wide)

df_merge_2170 <- merge(es_2170, df_SJ_2170wide, by.x="Cell", by.y="Cell", all.x=TRUE)
View(df_merge_2170)
write.csv(df_merge_2170, "/home/yue1118/Bioskryb_data_analysis/NCIH2170_ssGSEA_PVT1_splice_variant_df_07262025.csv")

df_merge_2170 <- read.csv("/home/yue1118/Bioskryb_data_analysis/NCIH2170_ssGSEA_PVT1_splice_variant_df_07262025.csv")

View(df_merge_2170)

cor_vec_2170 <- apply(df_merge_2170[ , 3:52], 2, function(x) cor(x, df_merge_2170[["chr8_127794735_127890588"]], use = "complete.obs"))

cor_df_2170 <- data.frame(
  variable = names(cor_vec_2170),
  correlation = cor_vec_2170
)

View(cor_df_2170)

svglite("/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/NCIH2170_PVT1_exon1_exon2_exp_vs_MYC_targets_V1_07242025.svg", width = 8, height = 6)
plot(df_merge_2170$chr8_127794735_127890588, df_merge_2170$HALLMARK_MYC_TARGETS_V1)
abline(lm(as.vector(df_merge_2170$HALLMARK_MYC_TARGETS_V1) ~ as.vector(df_merge_2170$chr8_127794735_127890588)), col = "hotpink2", lwd = 4)
dev.off()
cor.test(df_merge_2170$HALLMARK_MYC_TARGETS_V1, df_merge_2170$chr8_127794735_127890588) #cor=0.3741366,p-value = 0.001206


svglite("/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/NCIH2170_PVT1_exon1_exon2_exp_vs_MYC_targets_V2_07242025.svg", width = 8, height = 6)
plot(df_merge_2170$chr8_127794735_127890588, df_merge_2170$HALLMARK_MYC_TARGETS_V2)
abline(lm(as.vector(df_merge_2170$HALLMARK_MYC_TARGETS_V2) ~ as.vector(df_merge_2170$chr8_127794735_127890588)), col = "hotpink2", lwd = 4)
dev.off()
cor.test(df_merge_2170$HALLMARK_MYC_TARGETS_V2, df_merge_2170$chr8_127794735_127890588) #cor=0.3792603, p-value = 0.001018

#ssGSEA for SNU16 of Harper's/the 2nd plate 
seurat_SNU16 <- readRDS('/home/yue1118/Bioskryb_data_analysis/SNU16_plate2_Bioskryb_RNA_logTPM_FGFR2_MYC_ecDNA_counts_rds_09092025.rds')
X_SNU16 <- as.matrix(seurat_SNU16@assays$RNA@data)

gsvaPar_SNU16 <- gsvaParam(X_SNU16, gene_sets_list)
es_SNU16 <- gsva(gsvaPar_SNU16, verbose=FALSE)
es_SNU16 <- as.data.frame(t(es_SNU16))
View(es_SNU16)

table(rownames(es_SNU16) == colnames(seurat_SNU16))

write.csv(es_SNU16, "/home/yue1118/Bioskryb_data_analysis/SNU16_ssGSEA_Bioskryb_2nd_harpers_plate_09092025.csv")



#ssGSEA for SNU16
seurat_SNU16 <- readRDS('/home/yue1118/Bioskryb_data_analysis/SNU16_Bioskryb_RNA_logTPM_ecDNA_counts_rds_04292025.rds')
X_SNU16 <- as.matrix(seurat_SNU16@assays$RNA@data)

gsvaPar_SNU16 <- gsvaParam(X_SNU16, gene_sets_list)
es_SNU16 <- gsva(gsvaPar_SNU16, verbose=FALSE)
es_SNU16 <- as.data.frame(t(es_SNU16))
View(es_SNU16)

rownames(es_SNU16) == colnames(seurat_SNU16)
es_SNU16$ecDNA_CN <- seurat_SNU16$chr10_copy_number
es_SNU16$PVT1_exp <- X_SNU16["PVT1",]
write.csv(es_SNU16, "/home/yue1118/Bioskryb_data_analysis/SNU16_ssGSEA_PVT1_splice_variant_df_07262025.csv")

cor_results_SNU16 <- as.data.frame(sapply(es_SNU16[, -ncol(es_SNU16)], function(x) cor(x, es_SNU16[[ncol(es_SNU16)]], use = "complete.obs")))
View(cor_results_SNU16)

es_SNU16 <- read.csv("/home/yue1118/Bioskryb_data_analysis/SNU16_ssGSEA_PVT1_splice_variant_df_07262025.csv")


cor.test(es_SNU16$HALLMARK_MYC_TARGETS_V1, es_SNU16$PVT1_exp)#cor=-0.2959858; p-value = 0.008959
cor.test(es_SNU16$HALLMARK_MYC_TARGETS_V2, es_SNU16$PVT1_exp)#cor=-0.2633604; p-value = 0.02066
cor.test(es_SNU16$HALLMARK_MYC_TARGETS_V2, es_SNU16$ecDNA_CN) #cor=0.2038662; p-value = 0.07534
cor.test(es_SNU16$HALLMARK_MYC_TARGETS_V1, es_SNU16$ecDNA_CN)#cor=0.1368838;p-value = 0.2352

svglite("/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/SNU16_PVT1_exon1_exon2_exp_vs_MYC_targets_V1_07242025.svg", width = 8, height = 6)
plot(es_SNU16$PVT1_exp, es_SNU16$HALLMARK_MYC_TARGETS_V1)
abline(lm(as.vector(es_SNU16$HALLMARK_MYC_TARGETS_V1) ~ as.vector(es_SNU16$PVT1_exp)), col = "hotpink2", lwd = 4)
dev.off()
cor.test(es_SNU16$HALLMARK_MYC_TARGETS_V1, es_SNU16$PVT1_exp)#cor=-0.2959858 ;p-value = 0.008959


svglite("/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/SNU16_PVT1_exon1_exon2_exp_vs_MYC_targets_V2_07242025.svg", width = 8, height = 6)
plot(es_SNU16$PVT1_exp, es_SNU16$HALLMARK_MYC_TARGETS_V2)
abline(lm(as.vector(es_SNU16$HALLMARK_MYC_TARGETS_V2) ~ as.vector(es_SNU16$PVT1_exp)), col = "hotpink2", lwd = 4)
dev.off()
cor.test(es_SNU16$HALLMARK_MYC_TARGETS_V2, es_SNU16$PVT1_exp)#cor=-0.2633604;p-value = 0.02066

svglite("/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/SNU16_ecDNA_CNvs_MYC_targets_V2_07242025.svg", width = 8, height = 6)
plot(es_SNU16$ecDNA_CN, es_SNU16$HALLMARK_MYC_TARGETS_V2)
abline(lm(as.vector(es_SNU16$HALLMARK_MYC_TARGETS_V2) ~ as.vector(es_SNU16$ecDNA_CN)), col = "hotpink2", lwd = 4)
dev.off()
cor.test(es_SNU16$HALLMARK_MYC_TARGETS_V2, es_SNU16$ecDNA_CN) #cor=0.2038662;p-value = 0.07534

svglite("/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/SNU16_ecDNA_CNvs_MYC_targets_V1_07242025.svg", width = 8, height = 6)
plot(es_SNU16$ecDNA_CN, es_SNU16$HALLMARK_MYC_TARGETS_V1)
abline(lm(as.vector(es_SNU16$HALLMARK_MYC_TARGETS_V1) ~ as.vector(es_SNU16$ecDNA_CN)), col = "hotpink2", lwd = 4)
dev.off()
cor.test(es_SNU16$HALLMARK_MYC_TARGETS_V1, es_SNU16$ecDNA_CN) #cor=0.1368838; p-value = 0.2352


##ssGSEA for COLO320DM
seurat_DM_HSR <- readRDS("/home/yue1118/Bioskryb_data_analysis/Bioskryb_COLO320DM_HSR_naive_nonsorted_RNA_seurat_with_DNA_counts_06042025.rds")
X_DM_HSR <- as.matrix(seurat_DM_HSR@assays$RNA@data)

gsvaPar_DM_HSR <- gsvaParam(X_DM_HSR, gene_sets_list)
es_DM_HSR <- gsva(gsvaPar_DM_HSR, verbose=FALSE)
es_DM_HSR <- as.data.frame(t(es_DM_HSR))
View(es_DM_HSR)

table(rownames(es_DM_HSR) == seurat_DM_HSR$barcodes)
es_DM_HSR$ecDNA_CN <- seurat_DM_HSR$chr8_region_CN
es_DM_HSR$Cell <- seurat_DM_HSR$Cell

es_320DM <- es_DM_HSR[grep("ColoDM", rownames(es_DM_HSR)), ]


View(es_320DM)

#####Now incorporate SJ single cell information for COLO320DM####
df_SJ_320DMhigh <- read.delim("/datastore/lbcfs/labs/brunk_lab/private/Bioskryb_COLO320DM_Naive_RNA/COLO320DM_high_RNA/STAR_alignment_for_single_cells/COLO320DM_high_PVT1_chr8_SJs.tsv")
df_SJ_320DMmed <- read.delim("/datastore/lbcfs/labs/brunk_lab/private/Bioskryb_COLO320DM_Naive_RNA/COLO320DM_med_RNA/STAR_alignment_for_single_cells/COLO320DM_med_PVT1_chr8_SJs.tsv")
df_SJ_320DMlow <- read.delim("/datastore/lbcfs/labs/brunk_lab/private/Bioskryb_COLO320DM_Naive_RNA/COLO320DM_low_RNA/STAR_alignment_for_single_cells/COLO320DM_low_PVT1_chr8_SJs.tsv")
df_SJ_320DMall <- rbind(df_SJ_320DMhigh, df_SJ_320DMmed, df_SJ_320DMlow)

df_SJ_320DMall$splice_isoform <- paste(df_SJ_320DMall$chr, df_SJ_320DMall$start, df_SJ_320DMall$end, sep = "_")
df_SJ_320DMall <- df_SJ_320DMall[, c("Sample", "numUnique", "splice_isoform")]
df_SJ_320DMwide <- pivot_wider(df_SJ_320DMall,  names_from = splice_isoform, values_from = numUnique)


df_SJ_320DMwide$Cell <- sub(".*-R-(.*)_S.*", "\\1", df_SJ_320DMwide$Sample)
df_SJ_320DMwide$Cell <- sub("_.*", "", df_SJ_320DMwide$Cell)

df_merge_320 <- merge(es_320DM, df_SJ_320DMwide, by.x="Cell", by.y="Cell", all.x=TRUE)
write.csv(df_merge_320, "/home/yue1118/Bioskryb_data_analysis/COLO320DM_ssGSEA_PVT1_splice_variant_df_07262025.csv")

df_merge_320 <- read.csv("/home/yue1118/Bioskryb_data_analysis/COLO320DM_ssGSEA_PVT1_splice_variant_df_07262025.csv")
df_merge_320_subset <- df_merge_320[df_merge_320$chr8_127794735_127890588 <30,]
plot(df_merge_320_subset$chr8_127794735_127890588, df_merge_320_subset$HALLMARK_MYC_TARGETS_V2)
cor.test(df_merge_320_subset$chr8_127794735_127890588, df_merge_320_subset$HALLMARK_MYC_TARGETS_V2)

cor_vec_320 <- apply(df_merge_320[ , 3:52], 2, function(x) cor(x, df_merge_320[["chr8_127794735_127890588"]], use = "complete.obs"))

cor_df_320 <- data.frame(
  variable = names(cor_vec_320),
  correlation = cor_vec_320
)

View(cor_df_320)


cor.test(df_merge_320$HALLMARK_MYC_TARGETS_V2, df_merge_320$chr8_127794735_127890588)

svglite("/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/COLO320DM_PVT1_exon1_exon2_exp_vs_MYC_targets_V2_07242025.svg", width = 8, height = 6)
plot(df_merge_320_subset$chr8_127794735_127890588, df_merge_320_subset$HALLMARK_MYC_TARGETS_V2)
abline(lm(as.vector(df_merge_320_subset$HALLMARK_MYC_TARGETS_V2) ~ as.vector(df_merge_320_subset$chr8_127794735_127890588)), col = "hotpink2", lwd = 4)
dev.off()
cor.test(df_merge_320_subset$HALLMARK_MYC_TARGETS_V2, df_merge_320_subset$chr8_127794735_127890588) #cor=0.286473,p-value = 0.001823







cor_results_320DM <- as.data.frame(sapply(es_320DM[, -ncol(es_320DM)], function(x) cor(x, es_320DM[[ncol(es_320DM)]], use = "complete.obs")))
View(cor_results_320DM)

cor.test(es_2170$HALLMARK_MTORC1_SIGNALING, es_2170$ecDNA_CN)#p-value = 4.831e-05, cor=0.40
cor.test(es_SNU16$HALLMARK_MTORC1_SIGNALING, es_SNU16$ecDNA_CN)#p-value = 0.002765, cor=0.3365493
cor.test(es_320DM$HALLMARK_MTORC1_SIGNALING, es_320DM$ecDNA_CN)#p-value = 0.1156, cor=0.1281745 
cor.test(es_320DM$HALLMARK_MYC_TARGETS_V2, es_320DM$ecDNA_CN)#p-value = 6.11e-05, cor=0.3192214



##ssGSEA for CHP212
seurat_CHP212 <- readRDS("/datastore/lbcfs/labs/brunk_lab/public/GSM7157156_CHP212_scRNA_seq/CHP212_scRNA_seurat_obj_06282025.rds")
X_CHP212 <- as.matrix(seurat_CHP212@assays$SCT@data)

gsvaPar_CHP212 <- gsvaParam(X_CHP212, gene_sets_list)
es_CHP212 <- gsva(gsvaPar_CHP212, verbose=FALSE)
es_CHP212 <- as.data.frame(t(es_CHP212))




plot(ssgsea_scores_SNU16['HALLMARK_MTORC1_SIGNALING',], seurat_SNU16$chr10_copy_number)
cor(ssgsea_scores_SNU16['HALLMARK_MTORC1_SIGNALING',], seurat_SNU16$chr10_copy_number)

#ssGSEA for COLO320DM
seurat_DM_HSR <- readRDS("/home/yue1118/Bioskryb_data_analysis/Bioskryb_COLO320DM_HSR_naive_nonsorted_RNA_seurat_with_DNA_counts_06042025.rds")
X_DM_HSR <- as.matrix(seurat_DM_HSR@assays$RNA@data)
ssgsea_scores_DM_HSR <- gsva(X_DM_HSR, gene_sets_list, method = "ssgsea", ssgsea.norm = TRUE)

ssgsea_scores_DM_HSR <- as.data.frame(t(ssgsea_scores_DM_HSR))

table(rownames(ssgsea_scores_DM_HSR) == seurat_DM_HSR$barcodes)

ssgsea_scores_DM_HSR$ecDNA_CN <- seurat_DM_HSR$chr8_region_CN

ssgsea_scores_320DM <- ssgsea_scores_DM_HSR[grep("ColoDM", rownames(ssgsea_scores_DM_HSR)), ]

plot(as.numeric(ssgsea_scores_320DM$HALLMARK_MTORC1_SIGNALING), as.numeric(ssgsea_scores_320DM$ecDNA_CN))
cor(as.numeric(ssgsea_scores_320DM$HALLMARK_MTORC1_SIGNALING), as.numeric(ssgsea_scores_320DM$ecDNA_CN))


#seurat_2170$E2F_score <- ssgsea_scores['HALLMARK_E2F_TARGETS',]
#seurat_2170$mTORC1_score <- ssgsea_scores['HALLMARK_MTORC1_SIGNALING',]

plot(seurat_2170$E2F_score, seurat_2170$ecDNA_cn)
cor(as.numeric(seurat_2170$E2F_score), as.numeric(seurat_2170$ecDNA_cn))
cor(seurat_2170$E2F_score, seurat_2170$ecDNA_cn, use = "complete.obs")

plot(seurat_2170$mTORC1_score, seurat_2170$ecDNA_cn)
cor(seurat_2170$mTORC1_score, seurat_2170$ecDNA_cn, use = "complete.obs")
