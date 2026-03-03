.libPaths()

library(Seurat)
library(Signac)
library(dplyr)
library(tidyr)
#to load the read count matrix produced by Salmon
df_count <-read.csv("/home/yue1118/Bioskryb_data_analysis/matrix_gene_counts_salmon_SKGT2_COLO320.tsv", sep = '\t')
df_count_2 <- df_count
rownames(df_count_2) <- df_count_2$gene_id
df_count_2 <- df_count_2[,3:ncol(df_count_2)]

#To filter out samples/cells with no more than 1000 genes
n_genes_per_cell <- apply(df_count_2, 2, function(col) sum(as.logical(col)))
df_count_2 <- df_count_2[,names(n_genes_per_cell[n_genes_per_cell>1000])]

#To filter out genes that are present in less than 50% cells
n_cells_per_gene <- apply(df_count_2, 1, function(row) sum(as.logical(row)))
df_count_2 <- df_count_2[names(n_cells_per_gene[n_cells_per_gene > dim(df_count_2)[2]/2]),]

dim(df_count_2[,grepl("Colo", colnames(df_count_2))]) #None of the COLO320 cells passed this filtering step due to permeabilization to stain intracellular MYC protein.

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
write.csv(df1, "/home/yue1118/Bioskryb_data_analysis/SKGT2_Bioskryb_RNA_Salmon_TPM_low_cells_and_low_genes_filtered_out_05272025.csv", row.names=TRUE)
df1_log <- log(df1+1)
write.csv(df1_log, "/home/yue1118/Bioskryb_data_analysis/SKGT2_Bioskryb_RNA_Salmon_log_TPM_low_cells_and_low_genes_filtered_out_05272025.csv", row.names = TRUE)

df1_log <- read.csv("/home/yue1118/Bioskryb_data_analysis/SKGT2_Bioskryb_RNA_Salmon_log_TPM_low_cells_and_low_genes_filtered_out_05272025.csv")
View(df1_log)
rownames(df1_log) <- df1_log$X
df1_log <- df1_log[,-1]
df1_log_SKGT2 <- df1_log[, grep("SKGT", colnames(df1_log))]

list_ecDNA_ATAC = c('MYC', 'ERBB2', 'CDC6', 'IKZF3', 'PVT1', 'CASC11', 'WIPF2',  'LINC00824', 'PNMT', 'PCAT1', 'GRB7', 'RAPGEFL1', 'MIEN1', 'PGAP3', 'CASC21', 'CASC3', 'ZPBP2', 'TCAP',
                    'POU5F1B')

#Plot ERBB2 expression in high, med, low groups
df2 <- as.data.frame(t(df1_log_SKGT2["ERBB2",]))
rownames_vector <- rownames(df2)
df2$phenotype <- ifelse(grepl("\\.L\\.", rownames_vector, ignore.case = TRUE), "Low",
                        ifelse(grepl("\\.H\\.", rownames_vector, ignore.case = TRUE), "High",
                               ifelse(grepl("\\.M\\.", rownames_vector, ignore.case = TRUE), "Med", NA)))
df2$phenotype <- factor(df2$phenotype, levels = c("Low", "Med", "High"))
ggplot(df2, aes(x = `phenotype`, y = `ERBB2`)) +
  geom_boxplot() +
  labs(x = "Grouping Variable", y = "log(TPM+1)", title = "Boxplot by Group") +
  theme_minimal()
ggsave('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/ERBB2_exp_SKGT2_high_med_low_05272025.svg', width = 6, height = 6)

wilcox.test(df2[df2$phenotype =="High","ERBB2"], df2[df2$phenotype =="Low","ERBB2"], paired = F) #p-value = 0.2363
wilcox.test(df2[df2$phenotype =="High","ERBB2"], df2[df2$phenotype =="Med","ERBB2"], paired = F) #p-value = 0.2855
wilcox.test(df2[df2$phenotype =="Med","ERBB2"], df2[df2$phenotype =="Low","ERBB2"], paired = F) #p-value = 0.8379


#Plot MIEN1 expression in high, med, low groups
df2 <- as.data.frame(t(df1_log_SKGT2["MIEN1",]))
rownames_vector <- rownames(df2)
df2$phenotype <- ifelse(grepl("\\.L\\.", rownames_vector, ignore.case = TRUE), "Low",
                        ifelse(grepl("\\.H\\.", rownames_vector, ignore.case = TRUE), "High",
                               ifelse(grepl("\\.M\\.", rownames_vector, ignore.case = TRUE), "Med", NA)))
df2$phenotype <- factor(df2$phenotype, levels = c("Low", "Med", "High"))
ggplot(df2, aes(x = `phenotype`, y = `MIEN1`)) +
  geom_boxplot() +
  labs(x = "Grouping Variable", y = "log(TPM+1)", title = "Boxplot by Group") +
  theme_minimal()
ggsave('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/MIEN1_exp_SKGT2_high_med_low_05272025.svg', width = 6, height = 6)


wilcox.test(df2[df2$phenotype =="High","MIEN1"], df2[df2$phenotype =="Low","MIEN1"], paired = F) #p-value = 0.1525
wilcox.test(df2[df2$phenotype =="High","MIEN1"], df2[df2$phenotype =="Med","MIEN1"], paired = F) #p-value = 0.434
wilcox.test(df2[df2$phenotype =="Med","MIEN1"], df2[df2$phenotype =="Low","MIEN1"], paired = F) #p-value = 0.513


#Plot MIEN1 expression in high, med, low groups
df2 <- as.data.frame(t(df1_log_SKGT2["IKZF3",]))
rownames_vector <- rownames(df2)
df2$phenotype <- ifelse(grepl("\\.L\\.", rownames_vector, ignore.case = TRUE), "Low",
                        ifelse(grepl("\\.H\\.", rownames_vector, ignore.case = TRUE), "High",
                               ifelse(grepl("\\.M\\.", rownames_vector, ignore.case = TRUE), "Med", NA)))
df2$phenotype <- factor(df2$phenotype, levels = c("Low", "Med", "High"))
ggplot(df2, aes(x = `phenotype`, y = `IKZF3`)) +
  geom_boxplot() +
  labs(x = "Grouping Variable", y = "log(TPM+1)", title = "Boxplot by Group") +
  theme_minimal()
ggsave('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/IKZF3_exp_SKGT2_high_med_low_05272025.svg', width = 6, height = 6)


# Create a Seurat object.
seurat_SKGT <- CreateSeuratObject(
  counts = df1_log_SKGT2, 
  project = "Bioskryb_SKGT2", 
  assay = "RNA"
)

str(seurat_SKGT)
seurat_SKGT$FACS_sorting_phenotype <- df2$phenotype
Idents(seurat_SKGT) <- seurat_SKGT$FACS_sorting_phenotype

table(seurat_SKGT$FACS_sorting_phenotype)
View(seurat_SKGT@meta.data)


#Add SKGT2 ecDNA counts to the seurat object
#When merging RNA and DNA data, the part "-D/R-SC50-SKGT-M" is the part for cell barcodes. The "SC50" part is present in both RNA and DNA data.
df_counts <- read.csv("/home/yue1118/Bioskryb_data_analysis/SKGT2_Bioskryb_ecDNA_counts_05272025.csv")
df_counts$Cell <- sapply(strsplit(df_counts$Cell, "-D-"), `[`, 2)
df_counts$Cell <- sapply(strsplit(df_counts$Cell, "-_"), `[`, 1) #Double check this one!
df_metadata <- seurat_SKGT@meta.data
df_metadata$Cell <- rownames(df_metadata)
df_metadata$Cell <- sapply(strsplit(df_metadata$Cell, "\\.R\\."), `[`, 2)
df_metadata$Cell <- gsub("\\.", "-", df_metadata$Cell)
df_metadata$barcodes <- rownames(df_metadata)
#df_metadata$Cell <- sapply(strsplit(df_metadata$Cell, "-_"), `[`, 1)#Double check this one!

df_metadata_2 <- merge(df_metadata, df_counts, by.x="Cell", by.y="Cell", all.x=TRUE)
rownames(df_metadata_2) <- df_metadata_2$barcodes


df_reordered <- df_metadata_2[match(colnames(seurat_SKGT), df_metadata_2$barcodes), ]
seurat_SKGT@meta.data <- df_reordered

View(seurat_SKGT@meta.data)

# Extract the values from the Seurat object metadata
values <- seurat_SKGT@meta.data[["ERBB2_CN"]]

# Compute the quantiles for 30% and 70%
q30 <- quantile(values, 0.30, na.rm = TRUE)
q70 <- quantile(values, 0.70, na.rm = TRUE)

# Assign "high", "med", or "low" based on the quantiles
seurat_SKGT@meta.data$ERBB2_bin <- ifelse(
  values >= q70, "high",
  ifelse(values <= q30, "low", "med")
)

saveRDS(seurat_SKGT, '/home/yue1118/Bioskryb_data_analysis/SKGT2_Bioskryb_RNA_logTPM_ecDNA_counts_rds_05272025.rds')

seurat_SKGT <- readRDS('/home/yue1118/Bioskryb_data_analysis/SKGT2_Bioskryb_RNA_logTPM_ecDNA_counts_rds_05272025.rds')

#To plot the correlation between ERBB2 expression and ecDNA counts in SKGT2
seurat_SKGT_subset <- subset(seurat_SKGT, subset=ERBB2>1) #remove the one outlier
svglite("/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/SKGT2_inferred_ecDNA_counts_vs_ERBB2_expression_correlation_05272025.svg", width = 8, height = 6)
plot(as.vector(seurat_SKGT_subset$ERBB2_CN), as.vector(seurat_SKGT_subset@assays$RNA@data['ERBB2',]))
model <- lm(as.vector(seurat_SKGT_subset@assays$RNA@data['ERBB2',]) ~ as.vector(seurat_SKGT_subset$ERBB2_CN))
summary_model <- summary(model)
summary_model$r.squared #0.02118823
summary_model$coefficients #P-value=1.038822e-01
abline(lm(as.vector(seurat_SKGT_subset@assays$RNA@data['ERBB2',]) ~ as.vector(seurat_SKGT_subset$ERBB2_CN)), col = "hotpink2", lwd = 4)
dev.off()
cor_test <- cor.test(as.vector(seurat_SKGT_subset@assays$RNA@data['ERBB2',]), as.vector(seurat_SKGT_subset$ERBB2_CN)) 
cor_test$estimate #pearson correlation =0.1455618 
cor_test$p.value #p-value=0.1038822


#To plot the correlation between MIEN1 expression and ecDNA counts in SKGT2
seurat_SKGT_subset <- subset(seurat_SKGT, subset=MIEN1>1) #remove the one outlier
svglite("/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/SKGT2_inferred_ecDNA_counts_vs_MIEN1_expression_correlation_05272025.svg", width = 8, height = 6)
plot(as.vector(seurat_SKGT_subset$ERBB2_CN), as.vector(seurat_SKGT_subset@assays$RNA@data['MIEN1',]))
model <- lm(as.vector(seurat_SKGT_subset@assays$RNA@data['MIEN1',]) ~ as.vector(seurat_SKGT_subset$ERBB2_CN))
summary_model <- summary(model)
summary_model$r.squared #0.01205574
summary_model$coefficients #P-value=2.209871e-01
abline(lm(as.vector(seurat_SKGT_subset@assays$RNA@data['MIEN1',]) ~ as.vector(seurat_SKGT_subset$ERBB2_CN)), col = "hotpink2", lwd = 4)
dev.off()
cor_test <- cor.test(as.vector(seurat_SKGT_subset@assays$RNA@data['MIEN1',]), as.vector(seurat_SKGT_subset$ERBB2_CN)) 
cor_test$estimate #pearson correlation =0.1097986 
cor_test$p.value #p-value=0.2209871


#To plot the correlation between PVT1 expression and ecDNA counts in SKGT2
seurat_SKGT$MYC_CN <- 2.5*seurat_SKGT$chr8_coverage/seurat_SKGT$mean_coverage
seurat_SKGT_subset <- subset(seurat_SKGT, subset=PVT1>1) #remove the one outlier
svglite("/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/SKGT2_inferred_chr8_ecDNA_counts_vs_PVT1_expression_correlation_05272025.svg", width = 8, height = 6)
plot(as.vector(seurat_SKGT_subset$MYC_CN), as.vector(seurat_SKGT_subset@assays$RNA@data['PVT1',]))
model <- lm(as.vector(seurat_SKGT_subset@assays$RNA@data['PVT1',]) ~ as.vector(seurat_SKGT_subset$MYC_CN))
summary_model <- summary(model)
summary_model$r.squared #0.02087403
summary_model$coefficients #P-value=1.093854e-01
abline(lm(as.vector(seurat_SKGT_subset@assays$RNA@data['PVT1',]) ~ as.vector(seurat_SKGT_subset$MYC_CN)), col = "hotpink2", lwd = 4)
dev.off()
cor_test <- cor.test(as.vector(seurat_SKGT_subset@assays$RNA@data['PVT1',]), as.vector(seurat_SKGT_subset$MYC_CN)) 
cor_test$estimate #pearson correlation =0.1444785  
cor_test$p.value #p-value=0.1093854




seurat_SKGT <- readRDS('/home/yue1118/Bioskryb_data_analysis/SKGT2_Bioskryb_RNA_logTPM_ecDNA_counts_rds_05272025.rds')

df1 <- as.data.frame(as.matrix(seurat_SKGT@assays$RNA@data))
View(df1)
df1 <- as.data.frame(t(df1))
df1$HSR_count <- seurat_SKGT$ERBB2_CN
table(rownames(df1) == colnames(seurat_SKGT))
write.csv(df1, "/home/yue1118/Bioskryb_data_analysis/SKGT2_Bioskryb_RNA_logTPM_ecDNA_counts_dataframe_06162025.csv")






#To compute gene expression correlation with PVT1
seurat_SKGT2 <- readRDS('/home/yue1118/Bioskryb_data_analysis/SKGT2_Bioskryb_RNA_logTPM_ecDNA_counts_rds_05272025.rds')

expr_matrix_SKGT2 <- as.data.frame(GetAssayData(seurat_SKGT2, slot = "data")) 


# Initialize results
results_SKGT2 <- data.frame(
  gene = character(),
  pcc = numeric(),
  p_val = numeric(),
  stringsAsFactors = FALSE
)

# Get PVT1 expression vector
pvt1_expr_SKGT2 <- as.numeric(expr_matrix_SKGT2["PVT1", ])
SKGT2_amp_signature <- as.numeric(seurat_SKGT2$ERBB2_CN)
# Loop through each gene
for (gene in rownames(expr_matrix_SKGT2)) {
  
  gene_expr <- as.numeric(expr_matrix_SKGT2[gene, ])
  
  cor_result <- cor.test(gene_expr, SKGT2_amp_signature, method = "pearson")
  
  results_SKGT2 <- rbind(results_SKGT2, data.frame(
    gene = gene,
    pcc = cor_result$estimate,
    p_val = cor_result$p.value
  ))
}

# Adjust p-values for multiple testing
results_SKGT2$p_val_adj <- p.adjust(results_SKGT2$p_val, method = "BH")

msig_hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, entrez_gene)
gene_df <- bitr(results_SKGT2$gene, fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)
df_merged <- merge(results_SKGT2, gene_df, by.x = "gene", by.y = "SYMBOL")

gene_list <- df_merged$pcc
names(gene_list) <- df_merged$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

# 3. Run GSEA
gsea_result_SKGT2 <- GSEA(geneList = gene_list,
                    TERM2GENE = msig_hallmark,
                    pvalueCutoff = 0.05,
                    verbose = FALSE)
dotplot(gsea_result_SKGT2, showCategory = 15)

View(as.data.frame(gsea_result_SKGT2))

View(df_merged)



View(results_SKGT2[(results_SKGT2$p_val_adj < 0.05) & (results_SKGT2$pcc >0),])



# Get expression matrix (can use data slot or counts depending on normalization)
expr_matrix <- as.data.frame(GetAssayData(seurat_DM_HSR, slot = "data")) # log-normalized expression
df_metadata <- seurat_DM_HSR@meta.data
colnames(expr_matrix) <- rownames(df_metadata)

# Get cells in each group
cells_group1 <- rownames(df_metadata[df_metadata$cell_line_ecDNA_bin == "COLO320DMlow",])
cells_group2 <- rownames(df_metadata[df_metadata$cell_line_ecDNA_bin == "COLO320HSRlow",])

# Initialize result storage
results <- data.frame(
  gene = character(),
  p_val = numeric(),
  mean_value_diff = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each gene
for (gene in genes) {
  # Get expression values for each group
  expr1 <- expr_matrix[gene, cells_group1]
  expr2 <- expr_matrix[gene, cells_group2]
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
results$p_val_adj <- p.adjust(results$p_val, method = "BH")

# Sort by adjusted p-value
results <- results[order(results$p_val_adj), ]

View(results[(results$p_val_adj < 0.05) & (results$mean_value_diff >0),])





seurat_SKGT2 <- readRDS('/home/yue1118/Bioskryb_data_analysis/SKGT2_Bioskryb_RNA_logTPM_ecDNA_counts_rds_05272025.rds')

expr_matrix_SKGT2 <- as.data.frame(t(GetAssayData(seurat_SKGT2, slot = "data")))
expr_matrix_SKGT2$ERBB2_CN <- seurat_SKGT2$ERBB2_CN

write.csv(expr_matrix_SKGT2, "/home/yue1118/Bioskryb_data_analysis/SKGT2_bioskryb_RNA_gene_amp_CN_07012025.csv")


View(expr_matrix_SKGT2)

table(rownames(expr_matrix_SKGT2) == colnames(seurat_SKGT2))


seurat_SKGT2$ERBB2_CN

