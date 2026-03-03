###This notebook
library(GSVA)
packageVersion("GSVA")
library(msigdbr)
library(dplyr)
library(tidyr)


msig_hallmark <- msigdbr(species = "Homo sapiens", category = "H")
msig_hallmark <- msig_hallmark[, c("gs_name", "gene_symbol")]
gene_sets_list <- split(msig_hallmark$gene_symbol, msig_hallmark$gs_name)


#ssGSEA for COLO320HSR
seurat_DM_HSR <- readRDS("/home/yue1118/Bioskryb_data_analysis/Bioskryb_COLO320DM_HSR_naive_nonsorted_RNA_seurat_with_DNA_counts_06042025.rds")
X_DM_HSR <- as.matrix(seurat_DM_HSR@assays$RNA@data)

gsvaPar_DM_HSR <- gsvaParam(X_DM_HSR, gene_sets_list)
es_DM_HSR <- gsva(gsvaPar_DM_HSR, verbose=FALSE)
es_DM_HSR <- as.data.frame(t(es_DM_HSR))
View(es_DM_HSR)

table(rownames(es_DM_HSR) == seurat_DM_HSR$barcodes)
es_DM_HSR$ecDNA_CN <- seurat_DM_HSR$chr8_region_CN
es_DM_HSR$Cell <- seurat_DM_HSR$Cell

es_320HSR <- es_DM_HSR[grep("ColoHS", rownames(es_DM_HSR)), ]

write.csv(es_320HSR, "/home/yue1118/Bioskryb_data_analysis/Bioskryb_COLO320HSR_ssGSEA_08012025.csv")

#ssGSEA for KATOIII
seurat_kato <- readRDS("/home/yue1118/Bioskryb_data_analysis/KATOIII_Bioskryb_RNA_logTPM_FGFR2_copy_number_rds_05042025.rds")
X_kato <- as.matrix(seurat_kato@assays$RNA@data)

gsvaPar_kato <- gsvaParam(X_kato, gene_sets_list)
es_kato <- gsva(gsvaPar_kato, verbose=FALSE)
es_kato <- as.data.frame(t(es_kato))
View(es_kato)

table(rownames(es_kato) == colnames(seurat_kato))
es_kato$Cell <- seurat_kato$Cell
es_kato$ecDNA_CN <- seurat_kato$chr10_copy_number

write.csv(es_kato, "/home/yue1118/Bioskryb_data_analysis/Bioskryb_KATOIII_ssGSEA_df_08012025.csv")

####ssGSEA for SKGT2 ######
seurat_SKGT <- readRDS("/home/yue1118/Bioskryb_data_analysis/SKGT2_Bioskryb_RNA_logTPM_ecDNA_counts_rds_05272025.rds")
X_SKGT <- as.matrix(seurat_SKGT@assays$RNA@data)

gsvaPar_SKGT <- gsvaParam(X_SKGT, gene_sets_list)
es_SKGT <- gsva(gsvaPar_SKGT, verbose=FALSE)
es_SKGT <- as.data.frame(t(es_SKGT))
View(es_SKGT)

table(rownames(es_SKGT) == seurat_SKGT$barcodes)
es_SKGT$chr17_CN <- seurat_SKGT$ERBB2_CN
es_SKGT$chr8_coverage <- seurat_SKGT$chr8_coverage
es_SKGT$Cell <- seurat_SKGT$Cell


write.csv(es_SKGT, "/home/yue1118/Bioskryb_data_analysis/Bioskryb_SKGT2_ssGSEA_08012025.csv")
