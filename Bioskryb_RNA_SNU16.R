.libPaths()

library(Seurat)
library(Signac)
library(dplyr)
library(tidyr)
library(svglite)

#The TPM for the three cell lines NCIH2170, SNU16, and H716 has been processed and calculated before when processing NCIH2170 data.

df1_log <- read.csv("/home/yue1118/Bioskryb_data_analysis/RNA_Salmon_log_TPM_low_cells_and_low_genes_filtered_out_03052025.csv")
View(df1_log)
rownames(df1_log) <- df1_log$X
df1_log <- df1_log[,-1]
df1_log_SNU16 <- df1_log[, grep("SNU16", colnames(df1_log))]

#Plot FGFR2 expression in high, med, low groups
df2 <- as.data.frame(t(df1_log_SNU16["FGFR2",]))
rownames_vector <- rownames(df2)
df2$phenotype <- ifelse(grepl("Low", rownames_vector, ignore.case = TRUE), "Low",
                        ifelse(grepl("High", rownames_vector, ignore.case = TRUE), "High",
                               ifelse(grepl("Med", rownames_vector, ignore.case = TRUE), "Med", NA)))
df2$phenotype <- factor(df2$phenotype, levels = c("Low", "Med", "High"))
ggplot(df2, aes(x = `phenotype`, y = `FGFR2`)) +
  geom_boxplot() +
  labs(x = "Grouping Variable", y = "log(TPM+1)", title = "Boxplot by Group") +
  theme_minimal()
ggsave('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/FGFR2_exp_SNU16_high_med_low_04292025.svg', width = 6, height = 6)

wilcox.test(df2[df2$phenotype =="High","FGFR2"], df2[df2$phenotype =="Low","FGFR2"], paired = F) #p-value = 0.0007021
wilcox.test(df2[df2$phenotype =="High","FGFR2"], df2[df2$phenotype =="Med","FGFR2"], paired = F) #p-value = 0.003733
wilcox.test(df2[df2$phenotype =="Med","FGFR2"], df2[df2$phenotype =="Low","FGFR2"], paired = F) #p-value = 0.7511



# Create a Seurat object.
seurat_SNU16 <- CreateSeuratObject(
  counts = df1_log_SNU16, 
  project = "Bioskryb_SNU16", 
  assay = "RNA"
)

str(seurat_SNU16)
seurat_SNU16$FACS_sorting_phenotype <- df2$phenotype
Idents(seurat_SNU16) <- seurat_SNU16$FACS_sorting_phenotype

table(seurat_SNU16$FACS_sorting_phenotype)
View(seurat_SNU16@meta.data)

#Add SNU16 ecDNA counts to the seurat object
#When merging RNA and DNA data, the part "-D/R-SC337-SNU16High" is the part for cell barcodes. The "SC337" part is present in both RNA and DNA data.
df_counts <- read.csv("/home/yue1118/Bioskryb_data_analysis/SNU16_ecDNA_counts_04292025.csv")
df_counts$Cell <- sapply(strsplit(df_counts$Cell, "-D-"), `[`, 2)
df_counts$Cell <- sapply(strsplit(df_counts$Cell, "-_"), `[`, 1) #Double check this one!
df_metadata <- seurat_SNU16@meta.data
df_metadata$Cell <- rownames(df_metadata)
df_metadata$Cell <- sapply(strsplit(df_metadata$Cell, "\\.R\\."), `[`, 2)
df_metadata$Cell <- gsub("\\.", "-", df_metadata$Cell)
df_metadata$barcodes <- rownames(df_metadata)
df_metadata$Cell <- sapply(strsplit(df_metadata$Cell, "-_"), `[`, 1)#Double check this one!

df_metadata_2 <- merge(df_metadata, df_counts, by.x="Cell", by.y="Cell", all.x=TRUE)
rownames(df_metadata_2) <- df_metadata_2$barcodes


df_reordered <- df_metadata_2[match(colnames(seurat_SNU16), df_metadata_2$barcodes), ]
seurat_SNU16@meta.data <- df_reordered

View(seurat_SNU16@meta.data)


# Extract the values from the Seurat object metadata
values <- seurat_SNU16@meta.data[["chr10_copy_number"]]

# Compute the quantiles for 30% and 70%
q30 <- quantile(values, 0.30, na.rm = TRUE)
q70 <- quantile(values, 0.70, na.rm = TRUE)

# Assign "high", "med", or "low" based on the quantiles
seurat_SNU16@meta.data$ecDNA_bin <- ifelse(
  values >= q70, "high",
  ifelse(values <= q30, "low", "med")
)

saveRDS(seurat_SNU16, '/home/yue1118/Bioskryb_data_analysis/SNU16_Bioskryb_RNA_logTPM_ecDNA_counts_rds_04292025.rds')


#Now to integrate MYC bin copy number to the mRNA seurat object:

#Add SNU16 ecDNA MYC species counts to the seurat object
#When merging RNA and DNA data, the part "-D/R-SC337-SNU16High" is the part for cell barcodes. The "SC337" part is present in both RNA and DNA data.
df_counts <- read.csv("/home/yue1118/Bioskryb_data_analysis/SNU16_ecDNA_counts_FGFR2_species_and_MYC_species_06182025.csv")
df_counts$Cell <- sapply(strsplit(df_counts$Cell, "-D-"), `[`, 2)
df_counts$Cell <- sapply(strsplit(df_counts$Cell, "-_"), `[`, 1) #Double check this one!
df_metadata <- seurat_SNU16@meta.data
df_metadata$Cell <- rownames(df_metadata)
df_metadata$Cell <- sapply(strsplit(df_metadata$Cell, "\\.R\\."), `[`, 2)
df_metadata$Cell <- gsub("\\.", "-", df_metadata$Cell)
df_metadata$barcodes <- rownames(df_metadata)
df_metadata$Cell <- sapply(strsplit(df_metadata$Cell, "-_"), `[`, 1)#Double check this one!

df_metadata_2 <- merge(df_metadata, df_counts, by.x="Cell", by.y="Cell", all.x=TRUE)
rownames(df_metadata_2) <- df_metadata_2$barcodes


df_reordered <- df_metadata_2[match(colnames(seurat_SNU16), df_metadata_2$barcodes), ]
seurat_SNU16@meta.data <- df_reordered

View(seurat_SNU16@meta.data)




#To plot the correlation between FGFR2 expression and ecDNA counts in SNU16
seurat_SNU16_subset <- subset(seurat_SNU16, subset=FGFR2>1) #remove the one outlier
svglite("/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/inferred_ecDNA_counts_vs_FGFR2_SNU16_expression_correlation_04292025.svg", width = 8, height = 6)
plot(as.vector(seurat_SNU16_subset$chr10_copy_number), as.vector(seurat_SNU16_subset@assays$RNA@data['FGFR2',]))
model <- lm(as.vector(seurat_SNU16_subset@assays$RNA@data['FGFR2',]) ~ as.vector(seurat_SNU16_subset$chr10_copy_number))
summary_model <- summary(model)
summary_model$r.squared #0.1759619
summary_model$coefficients #P-value=1.462016e-04
abline(lm(as.vector(seurat_SNU16_subset@assays$RNA@data['FGFR2',]) ~ as.vector(seurat_SNU16_subset$chr10_copy_number)), col = "hotpink2", lwd = 4)
dev.off()
cor_test <- cor.test(as.vector(seurat_SNU16_subset@assays$RNA@data['FGFR2',]), as.vector(seurat_SNU16_subset$chr10_copy_number))
cor_test$estimate #pearson correlation =0.4194782 
cor_test$p.value #p-value=0.0001462016


#To plot the correlation between WDR11 expression and ecDNA counts in SNU16
seurat_SNU16_subset <- subset(seurat_SNU16, subset=WDR11>1) #remove the one outlier
svglite("/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/inferred_FGFR2_ecDNA_species_counts_vs_WDR11_SNU16_expression_correlation_04292025.svg", width = 8, height = 6)
plot(as.vector(seurat_SNU16_subset$chr10_copy_number), as.vector(seurat_SNU16_subset@assays$RNA@data['WDR11',]))
model <- lm(as.vector(seurat_SNU16_subset@assays$RNA@data['WDR11',]) ~ as.vector(seurat_SNU16_subset$chr10_copy_number))
summary_model <- summary(model)
summary_model$r.squared #0.09114327
summary_model$coefficients #P-value=1.533127e-02
abline(lm(as.vector(seurat_SNU16_subset@assays$RNA@data['WDR11',]) ~ as.vector(seurat_SNU16_subset$chr10_copy_number)), col = "hotpink2", lwd = 4)
dev.off()
cor_test <- cor.test(as.vector(seurat_SNU16_subset@assays$RNA@data['WDR11',]), as.vector(seurat_SNU16_subset$chr10_copy_number))
cor_test$estimate #pearson correlation =0.3018994  
cor_test$p.value #p-value=0.01533127


#To plot the correlation between MYC expression and ecDNA counts in SNU16
seurat_SNU16_subset <- subset(seurat_SNU16, subset=MYC>1) #remove the one outlier
svglite("/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/inferred_MYC_ecDNA_species_counts_vs_MYC_SNU16_expression_correlation_06182025.svg", width = 8, height = 6)
plot(as.vector(seurat_SNU16_subset$chr8_copy_number), as.vector(seurat_SNU16_subset@assays$RNA@data['MYC',]))
model <- lm(as.vector(seurat_SNU16_subset@assays$RNA@data['MYC',]) ~ as.vector(seurat_SNU16_subset$chr8_copy_number))
summary_model <- summary(model)
summary_model$r.squared #0.001481496
summary_model$coefficients #P-value=7.482079e-01
abline(lm(as.vector(seurat_SNU16_subset@assays$RNA@data['MYC',]) ~ as.vector(seurat_SNU16_subset$chr8_copy_number)), col = "hotpink2", lwd = 4)
dev.off()
cor_test <- cor.test(as.vector(seurat_SNU16_subset@assays$RNA@data['MYC',]), as.vector(seurat_SNU16_subset$chr8_copy_number))
cor_test$estimate #pearson correlation =0.03849021 
cor_test$p.value #p-value=0.7482079


#To plot the correlation between PVT1 expression and ecDNA counts in SNU16
seurat_SNU16_subset <- subset(seurat_SNU16, subset=PVT1>0) #remove the one outlier
svglite("/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/inferred_MYC_ecDNA_species_counts_vs_PVT1_SNU16_expression_correlation_06182025.svg", width = 8, height = 6)
plot(as.vector(seurat_SNU16_subset$chr8_copy_number), as.vector(seurat_SNU16_subset@assays$RNA@data['PVT1',]))
model <- lm(as.vector(seurat_SNU16_subset@assays$RNA@data['PVT1',]) ~ as.vector(seurat_SNU16_subset$chr8_copy_number))
summary_model <- summary(model)
summary_model$r.squared #0.001481496
summary_model$coefficients #P-value=7.482079e-01
abline(lm(as.vector(seurat_SNU16_subset@assays$RNA@data['PVT1',]) ~ as.vector(seurat_SNU16_subset$chr8_copy_number)), col = "hotpink2", lwd = 4)
dev.off()
cor_test <- cor.test(as.vector(seurat_SNU16_subset@assays$RNA@data['PVT1',]), as.vector(seurat_SNU16_subset$chr8_copy_number))
cor_test$estimate #pearson correlation =0.03640234  
cor_test$p.value #p-value=0.845851





df_expr_SNU16 <- as.data.frame(seurat_SNU16@assays$RNA@data)
df_expr_SNU16 <- as.data.frame(t(df_expr_SNU16))
rownames(df_expr_SNU16) <- gsub("\\.", "-", rownames(df_expr_SNU16))
rownames(df_expr_SNU16) <- gsub("-R-", "-D-", rownames(df_expr_SNU16))
rownames(df_expr_SNU16) <- sub(".*-(SC[0-9]{1,3}-SNU16(?:High|Low|Med))-.*", "\\1", rownames(df_expr_SNU16))
View(df_expr_SNU16)
df_expr_SNU16$Cell <- rownames(df_expr_SNU16)

df_SNU16_PVT1_MYC_FGFR2_wider_subset$Cell <- sub(".*-(SC[0-9]{1,3}-SNU16(?:High|Low|Med))-.*", "\\1", df_SNU16_PVT1_MYC_FGFR2_wider_subset$Cell)
df_merge_all <- merge(df_SNU16_PVT1_MYC_FGFR2_wider_subset, df_expr_SNU16, by.x="Cell", by.y="Cell", all.x=TRUE)
View(df_merge_all)
plot(df_merge_all$normalized_PVT1, df_merge_all$PVT1.y)


seurat_SNU16 <- readRDS('/home/yue1118/Bioskryb_data_analysis/SNU16_Bioskryb_RNA_logTPM_ecDNA_counts_rds_04292025.rds')

expr_matrix_snu16 <- as.data.frame(as.matrix(GetAssayData(seurat_SNU16, slot = "data")))

# Initialize results
results_snu16 <- data.frame(
  gene = character(),
  pcc = numeric(),
  p_val = numeric(),
  stringsAsFactors = FALSE
)

# Get PVT1 expression vector
PVT1_expr_snu16 <- as.numeric(expr_matrix_snu16["PVT1", ])

FGFR2_expr_snu16 <- as.numeric(expr_matrix_snu16["FGFR2", ])

SNU16_ecDNA_count_signature <- as.numeric(seurat_SNU16$chr10_copy_number)

# Loop through each gene
for (gene in rownames(expr_matrix_snu16)) {
  gene_expr <- as.numeric(expr_matrix_snu16[gene, ])
  
  cor_result <- cor.test(gene_expr, PVT1_expr_snu16, method = "pearson")
  
  results_snu16 <- rbind(results_snu16, data.frame(
    gene = gene,
    pcc = cor_result$estimate,
    p_val = cor_result$p.value
  ))
}

# Adjust p-values for multiple testing
results_snu16$p_val_adj <- p.adjust(results_snu16$p_val, method = "BH")
results_snu16 <- results_snu16[order(results_snu16$p_val_adj), ]
View(results_snu16)

msig_hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, entrez_gene)
gene_df <- bitr(results_snu16$gene, fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)
df_merged <- merge(results_snu16, gene_df, by.x = "gene", by.y = "SYMBOL")

gene_list <- df_merged$pcc
names(gene_list) <- df_merged$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

# 3. Run GSEA
gsea_result_snu16 <- GSEA(geneList = gene_list,
                    TERM2GENE = msig_hallmark,
                    pvalueCutoff = 0.05,
                    verbose = FALSE)
dotplot(gsea_result, showCategory = 10)

View(as.data.frame(gsea_result_snu16))
gsea_result_snu16 <- as.data.frame(gsea_result_snu16)
gsea_result_snu16$Description <- factor(gsea_result_snu16$Description, levels = gsea_result_snu16$Description[order(gsea_result_snu16$NES)])


ggplot(gsea_result_snu16, aes(x = Description, y = NES)) +
  geom_segment(aes(x = Description, xend = Description, y = 0, yend = NES), color = "grey50") +
  geom_point(aes(color = p.adjust), size = 5) +
  scale_color_gradient(low = "red", high = "blue", name = "Adjusted p-value") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Lollipop Plot of Pathways",
       x = NULL, y = "Normalized Enrichment Score (NES)")
ggsave('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/lollipop_GSEA_result_for_genes_correlate_with_ecDNA_counts_in_SNU16_bioskryb_06262025.svg', width = 6, height = 6)


seurat_SNU16_subset <- subset(seurat_SNU16, chr10_copy_number < 150) #26 cells
seurat_KATOIII_subset <- subset(seurat_KATOIII, chr10_copy_number < 150) #61 cells

wilcox.test(as.numeric(seurat_SNU16_subset$chr10_copy_number), as.numeric(seurat_KATOIII_subset$chr10_copy_number)) #p-value = 0.4498

mean(seurat_SNU16_subset$chr10_copy_number)

mean(seurat_KATOIII_subset$chr10_copy_number)


expr_matrix_SNU16 <- as.data.frame(GetAssayData(seurat_SNU16_subset, slot = "data")) 
expr_matrix_KATOIII <- as.data.frame(GetAssayData(seurat_KATOIII_subset, slot = "data")) 

list_common_genes <- intersect(rownames(expr_matrix_SNU16), rownames(expr_matrix_KATOIII))
# Initialize result storage
results <- data.frame(
  gene = character(),
  p_val = numeric(),
  mean_value_diff = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each gene
for (gene in list_common_genes) {
  # Get expression values for each group
  expr1 <- expr_matrix_SNU16[gene,]
  expr2 <- expr_matrix_KATOIII[gene, ]
  expr1 <- as.numeric(expr1)
  expr2 <- as.numeric(expr2)
  
  # Do Wilcoxon test
  test_result <- wilcox.test(expr1, expr2)
  
  # Calculate log2 fold change
  mean_value_diff <- mean(expr1) - mean(expr2)
  
  # Store results
  results <- rbind(results, data.frame(
    gene = gene,
    p_val = test_result$p.value,
    mean_value_diff = mean_value_diff
  ))
}

# Adjust p-values
results$p_val_adj <- p.adjust(results$p_val, method = "bonferroni")

# Sort by adjusted p-value
results <- results[order(results$p_val_adj), ]

View(results[(results$p_val_adj < 0.05) & (results$mean_value_diff >0),])



seurat_SNU16 <- readRDS('/home/yue1118/Bioskryb_data_analysis/SNU16_Bioskryb_RNA_logTPM_ecDNA_counts_rds_04292025.rds')

expr_matrix_snu16 <- as.data.frame(t(as.matrix(GetAssayData(seurat_SNU16, slot = "data"))))
table(rownames(expr_matrix_snu16) == colnames(seurat_SNU16))
expr_matrix_snu16$FGFR2_CN <- seurat_SNU16$chr10_copy_number

write.csv(expr_matrix_snu16, "/home/yue1118/Bioskryb_data_analysis/SNU16_bioskryb_RNA_FGFR2_CN_07012025.csv")

