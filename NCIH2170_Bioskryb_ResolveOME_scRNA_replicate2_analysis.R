.libPaths()

library(Seurat)
library(Signac)
library(dplyr)
library(tidyr)
#to load the read count matrix from the second plate produced by Salmon
df_count <-read.csv("/datastore/lbcfs/labs/brunk_lab/unc-external-downloads/bioskryb_UNC_BrunkGroup/AAG7YVWM5-OMEV2384-Beta---Brunk-Plate-2/2024_12_10T15_31_18/bj-expression_1.8.3_241119_91152/secondary_analyses/quantification_salmon/matrix_gene_counts_salmon.tsv",
                    sep = '\t')
df_count_2 <- df_count
rownames(df_count_2) <- df_count_2$gene_id
df_count_2 <- df_count_2[,3:ncol(df_count_2)]

#To filter out samples/cells with no more than 1000 genes
n_genes_per_cell <- apply(df_count_2, 2, function(col) sum(as.logical(col)))
df_count_2 <- df_count_2[,names(n_genes_per_cell[n_genes_per_cell>1000])]

#To filter out genes that are present in less than 50% cells
n_cells_per_gene <- apply(df_count_2, 1, function(row) sum(as.logical(row)))
df_count_2 <- df_count_2[names(n_cells_per_gene[n_cells_per_gene > dim(df_count_2)[2]/2]),]


#Now we calculate TPM for the filtered gene count matrix

#Read tsv file for transcript length
df_tl <- read.csv("/home/yue1118/Bioskryb_data_analysis/biomart_transcript_length.txt", sep = ",")
View(df_tl)
# Remove version numbers (everything after the dot)
df_tl$GeneID <- df_tl$Gene.stable.ID


df_count_2$gene_id <- rownames(df_count_2)
df_count_2$GeneID <- sub("\\..*", "", df_count_2$gene_id)
df_count_2 <- df_count_2[df_count_2$GeneID %in% df_tl$GeneID,]

df_tl <- df_tl[df_tl$GeneID %in% df_count_2$GeneID,]
#Combine multiple transcript lengths for one genes by taking the mean value
df_tl <- df_tl %>%
  group_by(GeneID) %>%
  summarise(
    transcript_length= mean(Transcript.length..including.UTRs.and.CDS.)  # Sum raw counts
  ) %>%
  ungroup()
df_tl <- df_tl[match(df_count_2$GeneID, df_tl$GeneID), ]

#Now normalize for the factor gene length
gene_lengths <- as.vector(df_tl$transcript_length)
df_count_3 <- df_count_2[,1:(ncol(df_count_2)-2)]
df_gene_length_normalized <- apply(df_count_3, 2, function(x) x / gene_lengths)
View(df_gene_length_normalized)

#Now normalize for sequencing depth to calculte TPM
column_sums <- colSums(df_gene_length_normalized)
scaling_factors <- column_sums / 1e6 #per million
df_tpm <- sweep(df_gene_length_normalized, 2, scaling_factors, "/")
df_tpm <- as.data.frame(df_tpm)

colSums(df_tpm)
View(df_tpm)

#Now match the gene IDs back to their gene symbols
df_gene_symbol <- df_count[, c(1,2)]
df_tpm$gene_id <- rownames(df_tpm)
df_tpm <- merge(df_tpm, df_gene_symbol, by.x="gene_id", by.y="gene_id", all.x=TRUE)

#To aggregate TPMs for ENSG IDs with the same gene symbol
num_rows <- length(unique(df_tpm$gene_symbol))
df1 <- data.frame(matrix(ncol = 0, nrow = num_rows))
for (i in colnames(df_tpm)[2:(dim(df_tpm)[2]-1)]) {
  df_test <- df_tpm[,c("gene_symbol", i)]
  colnames(df_test) <- c("gene_symbol", "TPM")
  df_agg <- df_test %>%
    group_by(gene_symbol) %>%
    summarise(
      TPM = sum(TPM)  # Sum raw counts
    ) %>%
    ungroup()
  colnames(df_agg) <- c("gene_symbol", i)
  df1 <- cbind(df1, df_agg)
}
rownames(df1) <- df1$gene_symbol
df1 <- df1[,seq(2, ncol(df1), by = 2)]

#To save the generated TPM files
View(df1)
write.csv(df1, "/home/yue1118/Bioskryb_data_analysis/Harper_2nd_plate_RNA_Salmon_TPM_low_cells_and_low_genes_filtered_out_03252025.csv", row.names=TRUE)
df1_log <- log(df1+1)
write.csv(df1_log, "/home/yue1118/Bioskryb_data_analysis/Harper_2nd_plate_RNA_Salmon_log_TPM_low_cells_and_low_genes_filtered_out_03252025.csv", row.names = TRUE)

df1_log <- read.csv("/home/yue1118/Bioskryb_data_analysis/Harper_2nd_plate_RNA_Salmon_log_TPM_low_cells_and_low_genes_filtered_out_03252025.csv")
View(df1_log)
rownames(df1_log) <- df1_log$X
df1_log <- df1_log[,-1]
df1_log_2170 <- df1_log[, grep("2170", colnames(df1_log))]

list_ecDNA_ATAC = c('MYC', 'ERBB2', 'CDC6', 'IKZF3', 'PVT1', 'CASC11', 'WIPF2',  'LINC00824', 'PNMT', 'PCAT1', 'GRB7', 'RAPGEFL1', 'MIEN1', 'PGAP3', 'CASC21', 'CASC3', 'ZPBP2', 'TCAP',
                    'POU5F1B')

#Plot ERBB2 expression in high, med, low groups
df2 <- as.data.frame(t(df1_log_2170["ERBB2",]))
rownames_vector <- rownames(df2)
df2$phenotype <- ifelse(grepl("Low", rownames_vector, ignore.case = TRUE), "Low",
                        ifelse(grepl("High", rownames_vector, ignore.case = TRUE), "High",
                               ifelse(grepl("Med", rownames_vector, ignore.case = TRUE), "Med", NA)))
df2$phenotype <- factor(df2$phenotype, levels = c("Low", "Med", "High"))
ggplot(df2, aes(x = `phenotype`, y = `ERBB2`)) +
  geom_boxplot() +
  labs(x = "Grouping Variable", y = "log(TPM+1)", title = "Boxplot by Group") +
  theme_minimal()
ggsave('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/ERBB2_exp_logTPM_2nd_plate_NCIH2170_high_med_low_03252025.svg', width = 6, height = 6)

wilcox.test(df2[df2$phenotype =="High","ERBB2"], df2[df2$phenotype =="Low","ERBB2"], paired = F) #p-value = 0.001294
wilcox.test(df2[df2$phenotype =="High","ERBB2"], df2[df2$phenotype =="Med","ERBB2"], paired = F) #p-value = 0.5268
wilcox.test(df2[df2$phenotype =="Med","ERBB2"], df2[df2$phenotype =="Low","ERBB2"], paired = F) #p-value = 0.114

#Plot CDC6 expression in high, med, low groups
df2 <- as.data.frame(t(df1_log_2170["CDC6",]))
rownames_vector <- rownames(df2)
df2$phenotype <- ifelse(grepl("Low", rownames_vector, ignore.case = TRUE), "Low",
                        ifelse(grepl("High", rownames_vector, ignore.case = TRUE), "High",
                               ifelse(grepl("Med", rownames_vector, ignore.case = TRUE), "Med", NA)))
df2$phenotype <- factor(df2$phenotype, levels = c("Low", "Med", "High"))
ggplot(df2, aes(x = `phenotype`, y = `CDC6`)) +
  geom_boxplot() +
  labs(x = "Grouping Variable", y = "log(TPM+1)", title = "Boxplot by Group") +
  theme_minimal()
ggsave('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/CDC6_exp_logTPM_2nd_plate_NCIH2170_high_med_low_03252025.svg', width = 6, height = 6)


wilcox.test(df2[df2$phenotype =="High","CDC6"], df2[df2$phenotype =="Low","CDC6"], paired = F) #p-value = 1.838e-06
wilcox.test(df2[df2$phenotype =="High","CDC6"], df2[df2$phenotype =="Med","CDC6"], paired = F) #p-value = 0.05114
wilcox.test(df2[df2$phenotype =="Med","CDC6"], df2[df2$phenotype =="Low","CDC6"], paired = F) #p-value = 0.3293


#Plot MYC expression in high, med, low groups
df2 <- as.data.frame(t(df1_log_2170["MYC",]))
rownames_vector <- rownames(df2)
df2$phenotype <- ifelse(grepl("Low", rownames_vector, ignore.case = TRUE), "Low",
                        ifelse(grepl("High", rownames_vector, ignore.case = TRUE), "High",
                               ifelse(grepl("Med", rownames_vector, ignore.case = TRUE), "Med", NA)))
df2$phenotype <- factor(df2$phenotype, levels = c("Low", "Med", "High"))
ggplot(df2, aes(x = `phenotype`, y = `MYC`)) +
  geom_boxplot() +
  labs(x = "Grouping Variable", y = "log(TPM+1)", title = "Boxplot by Group") +
  theme_minimal()
ggsave('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/MYC_exp_NCIH2170_high_med_low_03082025.svg', width = 6, height = 6)

wilcox.test(df2[df2$phenotype =="High","MYC"], df2[df2$phenotype =="Low","MYC"], paired = F) #p-value = 0.4965
wilcox.test(df2[df2$phenotype =="High","MYC"], df2[df2$phenotype =="Med","MYC"], paired = F) #p-value = 0.7234
wilcox.test(df2[df2$phenotype =="Med","MYC"], df2[df2$phenotype =="Low","MYC"], paired = F) #p-value = 0.3655



# Create a Seurat object.
seurat_2170 <- CreateSeuratObject(
  counts = df1_log_2170, 
  project = "Bioskryb_NCIH2170", 
  assay = "RNA"
)

str(seurat_2170)
seurat_2170$FACS_sorting_phenotype <- df2$phenotype
Idents(seurat_2170) <- seurat_2170$FACS_sorting_phenotype

table(seurat_2170$FACS_sorting_phenotype)
View(seurat_2170@meta.data)


#Add NCIH2170 ecDNA counts to the seurat object
#When merging RNA and DNA data, the part "-D/R-SC337-2170High" is the part for cell barcodes. The "SC337" part is present in both RNA and DNA data.
df_counts <- read.csv("/home/yue1118/Bioskryb_data_analysis/NCIH2170_plate2_ecDNA_counts_03252025.csv")
df_counts$Cell <- sapply(strsplit(df_counts$Cell, "-D-"), `[`, 2)
df_counts$Cell <- sapply(strsplit(df_counts$Cell, "-_"), `[`, 1) #Double check this one!
df_metadata <- seurat_2170@meta.data
df_metadata$Cell <- rownames(df_metadata)
df_metadata$Cell <- sapply(strsplit(df_metadata$Cell, "\\.R\\."), `[`, 2)
df_metadata$Cell <- gsub("\\.", "-", df_metadata$Cell)
df_metadata$barcodes <- rownames(df_metadata)
df_metadata$Cell <- sapply(strsplit(df_metadata$Cell, "-_"), `[`, 1)#Double check this one!

df_metadata_2 <- merge(df_metadata, df_counts, by.x="Cell", by.y="Cell", all.x=TRUE)
rownames(df_metadata_2) <- df_metadata_2$barcodes


df_reordered <- df_metadata_2[match(colnames(seurat_2170), df_metadata_2$barcodes), ]
seurat_2170@meta.data <- df_reordered

View(seurat_2170@meta.data)


# Extract the values from the Seurat object metadata
values <- seurat_2170@meta.data[["ecDNA_cn"]]

# Compute the quantiles for 30% and 70%
q30 <- quantile(values, 0.30, na.rm = TRUE)
q70 <- quantile(values, 0.70, na.rm = TRUE)

# Assign "high", "med", or "low" based on the quantiles
seurat_2170@meta.data$ecDNA_bin <- ifelse(
  values >= q70, "high",
  ifelse(values <= q30, "low", "med")
)

saveRDS(seurat_2170, '/home/yue1118/Bioskryb_data_analysis/NCIH2170_plate2_Bioskryb_RNA_logTPM_ecDNA_counts_rds_03252025.rds')

seurat_2170 <- readRDS('/home/yue1118/Bioskryb_data_analysis/NCIH2170_plate2_Bioskryb_RNA_logTPM_ecDNA_counts_rds_03252025.rds')

#To perform DEG analysis between FACS sorting labeled groups:
Idents(seurat_2170) <- seurat_2170$FACS_sorting_phenotype
df_deg <- FindMarkers(seurat_2170, ident.1 = "High", ident.2 = "Low")
View(df_deg[df_deg$p_val_adj < 0.05,]) 

volcano_plot_function <- function(input_df) {
  # Convert p-value to -log10 p-value for the volcano plot
  input_df$neg_log10_pvalue <- -log10(input_df$p_val_adj)
  
  # Add a column to classify significance (adjust thresholds as needed)
  input_df$significant <- ifelse(input_df$avg_log2FC > 0 & input_df$p_val_adj < 0.05, "Upregulated",
                                 ifelse(input_df$avg_log2FC < 0 & input_df$p_val_adj < 0.05, "Downregulated", "Not Significant"))
  
  # Create the volcano plot
  ggplot(input_df, aes(x = avg_log2FC, y = neg_log10_pvalue)) +
    geom_point(aes(color = significant), size= 4) + # Color points by significance
    scale_color_manual(values = c("blue", "gray", "red")) + # Custom colors
    labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-log10(adj p-value)") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") + # Add dashed lines for log2FC thresholds
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") + 
    theme(
      panel.grid.major.x = element_blank(),               # Remove major grid lines on x-axis if needed
      panel.grid.major.y = element_blank()                # Remove major grid lines on y-axis if needed
    )# Horizontal line for p-value threshold
}

volcano_plot_function(df_deg)
ggsave('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/DEGs_volcano_plot_NCIH2170_plate2_FACS_sorting_high_vs_low_03252025.svg', width = 6, height = 6)

#To perform DEG analysis between inferred ecDNA count labeled groups:
Idents(seurat_2170) <- seurat_2170$ecDNA_bin
df_deg_ecDNA <- FindMarkers(seurat_2170, ident.1 = "high", ident.2 = "low")
volcano_plot_function(df_deg_ecDNA)
ggsave('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/DEGs_volcano_plot_NCIH2170_plate2_ecDNA_inferred_counts_high_vs_low_03252025.svg', width = 6, height = 6)

View(df_deg_ecDNA)


list_gene <- c()
list_pcc <- c()
list_pvalue_pcc <- c()
df_2170 <- as.data.frame(t(seurat_2170@assays$RNA@data))
df_2170$ecDNA_counts <- seurat_2170$ecDNA_cn
for (i in rownames(seurat_2170)) {
  df_2170_subset <- df_2170[df_2170[i]>1,]
  result <- try(cor_test <- cor.test(as.vector(df_2170_subset[, i]), as.vector(df_2170_subset[, "ecDNA_counts"])), silent = TRUE) 
  # Check if the command failed
  if (inherits(result, "try-error")) {
    list_pcc <- c(list_pcc, "doesnt_exist")
    list_pvalue_pcc <- c(list_pvalue_pcc, "doesnt_exist")
  } else {
    list_pcc <- c(list_pcc, cor_test$estimate)
    list_pvalue_pcc <- c(list_pvalue_pcc, cor_test$p.value)
  }
  list_gene <- c(list_gene, i)
}

df_cor_all <- data.frame(gene = list_gene, pcc = list_pcc, pvalue_pcc = list_pvalue_pcc)
df_cor_all_filtered <- df_cor_all[!df_cor_all$gene %in% c("PRSS3"), ]
df_cor_all_filtered$pvalue_pcc <- as.numeric(df_cor_all_filtered$pvalue_pcc)
#df_cor_all_filtered <- df_cor_all_filtered[df_cor_all_filtered$pvalue_pcc < 0.01,]
df_cor_all_filtered$pcc <- as.numeric(df_cor_all_filtered$pcc)
View(df_cor_all_filtered)
View(df_cor_all)

hist(df_cor_all_filtered$pcc)

#perform Bonferroni correction on the p-values
df_cor_all_filtered$p_val_adj <- df_cor_all_filtered$pvalue_pcc * dim(seurat_2170)[1]
View(df_cor_all_filtered[df_cor_all_filtered$p_val_adj < 0.05,])

write.csv(df_cor_all_filtered, '/home/yue1118/Bioskryb_data_analysis/Bioskryb_plate2_gene_expression_TPM_correlation_with_ecDNA_counts_NCIH2170_03252025.csv')


#To plot the correlation between ERBB2 expression and ecDNA counts in NCIH2170
seurat_2170_subset <- subset(seurat_2170, subset=ERBB2>1) #remove the one outlier
svglite("/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/plate2_inferred_ecDNA_counts_vs_ERBB2_expression_correlation_03252025.svg", width = 8, height = 6)
plot(as.vector(seurat_2170_subset$ecDNA_cn), as.vector(seurat_2170_subset@assays$RNA@data['ERBB2',]))
model <- lm(as.vector(seurat_2170_subset@assays$RNA@data['ERBB2',]) ~ as.vector(seurat_2170_subset$ecDNA_cn))
summary_model <- summary(model)
summary_model$r.squared #0.1848542
summary_model$coefficients #P-value=1.180463e-06
abline(lm(as.vector(seurat_2170_subset@assays$RNA@data['ERBB2',]) ~ as.vector(seurat_2170_subset$ecDNA_cn)), col = "hotpink2", lwd = 4)
dev.off()
cor_test <- cor.test(as.vector(seurat_2170_subset@assays$RNA@data['ERBB2',]), as.vector(seurat_2170_subset$ecDNA_cn)) 
cor_test$estimate #pearson correlation =0.4299467 
cor_test$p.value #1.180463e-06


#To plot the correlation between WIPF2 expression and ecDNA counts in NCIH2170
seurat_2170_subset <- subset(seurat_2170, subset=WIPF2>1) #remove the one outlier
svglite("/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/plate2_inferred_ecDNA_counts_vs_WIPF2_expression_correlation_03252025.svg", width = 8, height = 6)
plot(seurat_2170_subset$ecDNA_cn, seurat_2170_subset@assays$RNA@data['WIPF2',])
model <- lm(as.vector(seurat_2170_subset@assays$RNA@data['WIPF2',]) ~ as.vector(seurat_2170_subset$ecDNA_cn))
summary_model <- summary(model)
summary_model$r.squared #0.4733022
summary_model$coefficients #P-value=1.445104e-17
abline(lm(as.vector(seurat_2170_subset@assays$RNA@data['WIPF2',]) ~ as.vector(seurat_2170_subset$ecDNA_cn)), col = "hotpink2", lwd = 4)
dev.off()
cor_test <- cor.test(as.vector(seurat_2170_subset@assays$RNA@data['WIPF2',]), as.vector(seurat_2170_subset$ecDNA_cn))
cor_test$estimate #pearson correlation =0.6879696 
cor_test$p.value #p-value=1.445104e-17


#To plot the correlation between MIEN1 expression and ecDNA counts in NCIH2170
seurat_2170_subset <- subset(seurat_2170, subset=MIEN1>1) #remove the one outlier
svglite("/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/plate2_inferred_ecDNA_counts_vs_MIEN1_expression_correlation_03252025.svg", width = 8, height = 6)
plot(seurat_2170_subset$ecDNA_cn, seurat_2170_subset@assays$RNA@data['MIEN1',])
model <- lm(as.vector(seurat_2170_subset@assays$RNA@data['MIEN1',]) ~ as.vector(seurat_2170_subset$ecDNA_cn))
summary_model <- summary(model)
summary_model$r.squared #0.7447722
summary_model$coefficients #P-value=3.414583e-36
abline(lm(as.vector(seurat_2170_subset@assays$RNA@data['MIEN1',]) ~ as.vector(seurat_2170_subset$ecDNA_cn)), col = "hotpink2", lwd = 4)
dev.off()
cor_test <- cor.test(as.vector(seurat_2170_subset@assays$RNA@data['MIEN1',]), as.vector(seurat_2170_subset$ecDNA_cn)) 
cor_test$estimate #pearson correlation =0.8630018 
cor_test$p.value #p-value=3.414583e-36


#To plot the correlation between PVT1 expression and ecDNA counts in NCIH2170
seurat_2170_subset <- subset(seurat_2170, subset=PVT1>1) #remove the one outlier
svglite("/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/plate2_inferred_ecDNA_counts_vs_PVT1_expression_correlation_03252025.svg", width = 8, height = 6)
plot(seurat_2170_subset$ecDNA_cn, seurat_2170_subset@assays$RNA@data['PVT1',])
model <- lm(as.vector(seurat_2170_subset@assays$RNA@data['PVT1',]) ~ as.vector(seurat_2170_subset$ecDNA_cn))
summary_model <- summary(model)
summary_model$r.squared #0.4715761
summary_model$coefficients #P-value=9.14008e-18
abline(lm(as.vector(seurat_2170_subset@assays$RNA@data['PVT1',]) ~ as.vector(seurat_2170_subset$ecDNA_cn)), col = "hotpink2", lwd = 4)
dev.off()
cor_test <- cor.test(as.vector(seurat_2170_subset@assays$RNA@data['PVT1',]), as.vector(seurat_2170_subset$ecDNA_cn))
cor_test$estimate #pearson correlation =0.686714
cor_test$p.value #p-value=9.14008e-18
