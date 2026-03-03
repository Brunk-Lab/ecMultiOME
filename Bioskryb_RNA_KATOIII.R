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
df1_log_KATOIII <- df1_log[, grep("KATO", colnames(df1_log))]

#Plot FGFR2 expression in high, med, low groups
df2 <- as.data.frame(t(df1_log_KATOIII["FGFR2",]))
rownames_vector <- rownames(df2)
df2$phenotype <- ifelse(grepl("Low", rownames_vector, ignore.case = TRUE), "Low",
                        ifelse(grepl("High", rownames_vector, ignore.case = TRUE), "High",
                               ifelse(grepl("Med", rownames_vector, ignore.case = TRUE), "Med", NA)))
df2$phenotype <- factor(df2$phenotype, levels = c("Low", "Med", "High"))
ggplot(df2, aes(x = `phenotype`, y = `FGFR2`)) +
  geom_boxplot() +
  labs(x = "Grouping Variable", y = "log(TPM+1)", title = "Boxplot by Group") +
  theme_minimal()
ggsave('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/', width = 6, height = 6)

wilcox.test(df2[df2$phenotype =="High","FGFR2"], df2[df2$phenotype =="Low","FGFR2"], paired = F) #2.548e-05
wilcox.test(df2[df2$phenotype =="High","FGFR2"], df2[df2$phenotype =="Med","FGFR2"], paired = F) #p-value = 0.04214
wilcox.test(df2[df2$phenotype =="Med","FGFR2"], df2[df2$phenotype =="Low","FGFR2"], paired = F) #p-value = 0.1655

View(df2)

# Create a Seurat object.
seurat_KATOIII <- CreateSeuratObject(
  counts = df1_log_KATOIII, 
  project = "Bioskryb_KATOIII", 
  assay = "RNA"
)

str(seurat_KATOIII)
seurat_KATOIII$FACS_sorting_phenotype <- df2$phenotype
Idents(seurat_KATOIII) <- seurat_KATOIII$FACS_sorting_phenotype

table(seurat_KATOIII$FACS_sorting_phenotype)
View(seurat_KATOIII@meta.data)

#Add KATOIII FGFR2 copy number to the seurat object
#When merging RNA and DNA data, the part "-D/R-SC337-KATO3High" is the part for cell barcodes. The "SC337" part is present in both RNA and DNA data.
df_counts <- read.csv("/home/yue1118/Bioskryb_data_analysis/KATOIII_FGFR2_copy_number_Bioskryb_05042025.csv")
df_counts$Cell <- sapply(strsplit(df_counts$Cell, "-D-"), `[`, 2)
df_counts$Cell <- sapply(strsplit(df_counts$Cell, "-_"), `[`, 1) #Double check this one!
df_metadata <- seurat_KATOIII@meta.data
df_metadata$Cell <- rownames(df_metadata)
df_metadata$Cell <- sapply(strsplit(df_metadata$Cell, "\\.R\\."), `[`, 2)
df_metadata$Cell <- gsub("\\.", "-", df_metadata$Cell)
df_metadata$barcodes <- rownames(df_metadata)
df_metadata$Cell <- sapply(strsplit(df_metadata$Cell, "-_"), `[`, 1)#Double check this one!

df_metadata_2 <- merge(df_metadata, df_counts, by.x="Cell", by.y="Cell", all.x=TRUE)
rownames(df_metadata_2) <- df_metadata_2$barcodes


df_reordered <- df_metadata_2[match(colnames(seurat_KATOIII), df_metadata_2$barcodes), ]
seurat_KATOIII@meta.data <- df_reordered

View(seurat_KATOIII@meta.data)


# Extract the values from the Seurat object metadata
values <- seurat_KATOIII@meta.data[["chr10_copy_number"]]

# Compute the quantiles for 30% and 70%
q30 <- quantile(values, 0.30, na.rm = TRUE)
q70 <- quantile(values, 0.70, na.rm = TRUE)

# Assign "high", "med", or "low" based on the quantiles
seurat_KATOIII@meta.data$ecDNA_bin <- ifelse(
  values >= q70, "high",
  ifelse(values <= q30, "low", "med")
)

saveRDS(seurat_KATOIII, '/home/yue1118/Bioskryb_data_analysis/KATOIII_Bioskryb_RNA_logTPM_FGFR2_copy_number_rds_05042025.rds')

seurat_KATOIII <- readRDS('/home/yue1118/Bioskryb_data_analysis/KATOIII_Bioskryb_RNA_logTPM_FGFR2_copy_number_rds_05042025.rds')

#To plot the correlation between FGFR2 expression and copy number in KATOIII
seurat_KATOIII_subset <- subset(seurat_KATOIII, subset=FGFR2>1) #remove the one outlier
svglite("/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/inferred_ecDNA_counts_vs_FGFR2_KATOIII_expression_correlation_05042025.svg", width = 8, height = 6)
plot(as.vector(seurat_KATOIII_subset$chr10_copy_number), as.vector(seurat_KATOIII_subset@assays$RNA@data['FGFR2',]))
model <- lm(as.vector(seurat_KATOIII_subset@assays$RNA@data['FGFR2',]) ~ as.vector(seurat_KATOIII_subset$chr10_copy_number))
summary_model <- summary(model)
summary_model$r.squared#0.2722505
summary_model$coefficients #p-value = 5.877901e-07
abline(lm(as.vector(seurat_KATOIII_subset@assays$RNA@data['FGFR2',]) ~ as.vector(seurat_KATOIII_subset$chr10_copy_number)), col = "hotpink2", lwd = 4)
dev.off()
cor_test <- cor.test(as.vector(seurat_KATOIII_subset@assays$RNA@data['FGFR2',]), as.vector(seurat_KATOIII_subset$chr10_copy_number))
cor_test$estimate #pearson correlation =0.5217763 
cor_test$p.value #p-value=5.877901e-07

ggplot(df1, aes(x = chr10_copy_number)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  labs(title = "KDE Plot", x = "Value", y = "Density")


df_meta <- seurat_KATOIII@meta.data
View(df_meta)

df_meta$FGFR2_expression <- as.vector(seurat_KATOIII@assays$RNA@data['FGFR2',])
df_meta$ecDNA_bin <- factor(df_meta$ecDNA_bin, levels = c("low", "med", "high"))
ggplot(df_meta, aes(x = `ecDNA_bin`, y = `FGFR2_expression`)) +
  geom_boxplot() +
  labs(x = "Grouping Variable", y = "log(TPM+1)", title = "Boxplot by Group") +
  theme_minimal()
ggsave('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/FGFR2_exp_KATOIII_ecDNA_inferred_high_med_low_05052025.svg', width = 6, height = 6)

wilcox.test(df_meta[df_meta$ecDNA_bin =="high","FGFR2_expression"], df_meta[df_meta$ecDNA_bin =="low","FGFR2_expression"], paired = F) #p-value = 0.0001294
wilcox.test(df_meta[df_meta$ecDNA_bin =="high","FGFR2_expression"], df_meta[df_meta$ecDNA_bin =="med","FGFR2_expression"], paired = F) #p-value = 7.83e-05
wilcox.test(df_meta[df_meta$ecDNA_bin =="med","FGFR2_expression"], df_meta[df_meta$ecDNA_bin =="low","FGFR2_expression"], paired = F) #p-value = 0.8573


df_meta$ecDNA_bin <- factor(df_meta$ecDNA_bin, levels = c("low", "med", "high"))
ggplot(df_meta, aes(x = `ecDNA_bin`, y = `chr10_copy_number`)) +
  geom_boxplot() +
  labs(x = "Grouping Variable", y = "count", title = "Boxplot by Group") +
  theme_minimal()
ggsave('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/KATOIII_ecDNA_counts_plate1_across_inferred_high_med_low_bins_05052025.svg', width = 6, height = 6)


wilcox.test(df_meta[df_meta$ecDNA_bin =="high","chr10_copy_number"], df_meta[df_meta$ecDNA_bin =="low","chr10_copy_number"], paired = F)  #p-value = 1.582e-14
wilcox.test(df_meta[df_meta$ecDNA_bin =="high","chr10_copy_number"], df_meta[df_meta$ecDNA_bin =="med","chr10_copy_number"], paired = F) #p-value = 3.588e-16
wilcox.test(df_meta[df_meta$ecDNA_bin =="med","chr10_copy_number"], df_meta[df_meta$ecDNA_bin =="low","chr10_copy_number"], paired = F)  #p-value = 3.588e-16



seurat_KATOIII <- readRDS('/home/yue1118/Bioskryb_data_analysis/KATOIII_Bioskryb_RNA_logTPM_FGFR2_copy_number_rds_05042025.rds')

expr_matrix_KATOIII <- as.data.frame(GetAssayData(seurat_KATOIII, slot = "data")) 


# Initialize results
results_KATOIII <- data.frame(
  gene = character(),
  pcc = numeric(),
  p_val = numeric(),
  stringsAsFactors = FALSE
)

# Get PVT1 expression vector

KATOIII_amp_signature <- as.numeric(seurat_KATOIII$chr10_copy_number)
# Loop through each gene
for (gene in rownames(expr_matrix_KATOIII)) {
  
  gene_expr <- as.numeric(expr_matrix_KATOIII[gene, ])
  
  cor_result <- cor.test(gene_expr, KATOIII_amp_signature, method = "pearson")
  
  results_KATOIII <- rbind(results_KATOIII, data.frame(
    gene = gene,
    pcc = cor_result$estimate,
    p_val = cor_result$p.value
  ))
}

# Adjust p-values for multiple testing
results_KATOIII$p_val_adj <- p.adjust(results_KATOIII$p_val, method = "BH")

msig_hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, entrez_gene)
gene_df <- bitr(results_KATOIII$gene, fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)
df_merged <- merge(results_KATOIII, gene_df, by.x = "gene", by.y = "SYMBOL")

gene_list <- df_merged$pcc
names(gene_list) <- df_merged$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

# 3. Run GSEA
gsea_result_KATOIII <- GSEA(geneList = gene_list,
                    TERM2GENE = msig_hallmark,
                    pvalueCutoff = 0.05,
                    verbose = FALSE)


View(as.data.frame(gsea_result_KATOIII))

gsea_result_KATOIII <- as.data.frame(gsea_result_KATOIII)
gsea_result_KATOIII$Description <- factor(gsea_result_KATOIII$Description, levels = gsea_result_KATOIII$Description[order(gsea_result_KATOIII$NES)])


ggplot(gsea_result_KATOIII, aes(x = Description, y = NES)) +
  geom_segment(aes(x = Description, xend = Description, y = 0, yend = NES), color = "grey50") +
  geom_point(aes(color = p.adjust), size = 5) +
  scale_color_gradient(low = "red", high = "blue", name = "Adjusted p-value") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Lollipop Plot of Pathways",
       x = NULL, y = "Normalized Enrichment Score (NES)")
ggsave('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/lollipop_GSEA_result_for_genes_correlate_with_ecDNA_counts_in_KATOIII_bioskryb_06262025.svg', width = 6, height = 6)



seurat_KATOIII <- readRDS('/home/yue1118/Bioskryb_data_analysis/KATOIII_Bioskryb_RNA_logTPM_FGFR2_copy_number_rds_05042025.rds')
expr_matrix_KATOIII <- as.data.frame(t(GetAssayData(seurat_KATOIII, slot = "data"))) 
table(rownames(expr_matrix_KATOIII) == colnames(seurat_KATOIII))

expr_matrix_KATOIII$FGFR2_CN <- seurat_KATOIII$chr10_copy_number
write.csv(expr_matrix_KATOIII, "/home/yue1118/Bioskryb_data_analysis/KATOIII_bioskryb_RNA_FGFR2_CN_07012025.csv")

