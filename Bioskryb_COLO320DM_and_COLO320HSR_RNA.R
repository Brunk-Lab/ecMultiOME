.libPaths()

library(Seurat)
library(Signac)
library(dplyr)
library(tidyr)
library(ggplot2)
library(svglite)
#to load the read count matrix produced by Salmon

df_count <-read.csv("/home/yue1118/Bioskryb_data_analysis/Bioskryb_320DM_HSR_Naive_nonsorted_matrix_gene_counts_salmon_JD_with_cutadapt_05302025.tsv", sep = '\t')
df_count_2 <- df_count
rownames(df_count_2) <- df_count_2$gene_id
df_count_2 <- df_count_2[,3:ncol(df_count_2)]

#To filter out samples/cells with no more than 1000 genes
n_genes_per_cell <- apply(df_count_2, 2, function(col) sum(as.logical(col)))
df_count_2 <- df_count_2[,names(n_genes_per_cell[n_genes_per_cell>1000])]

#To filter out genes that are present in less than 50% cells
n_cells_per_gene <- apply(df_count_2, 1, function(row) sum(as.logical(row)))
df_count_2 <- df_count_2[names(n_cells_per_gene[n_cells_per_gene > dim(df_count_2)[2]/2]),]
View(df_count_2)


#Read tsv file for transcript length
df_tl <- read.csv("/home/yue1118/Bioskryb_data_analysis/biomart_transcript_length.txt", sep = ",")
View(df_tl)
# Remove version numbers (everything after the dot)
df_tl$GeneID <- df_tl$Gene.stable.ID


df_count_2$gene_id <- rownames(df_count_2)
df_count_2$GeneID <- sub("\\..*", "", df_count_2$gene_id)
df_count_2 <- df_count_2[df_count_2$GeneID %in% df_tl$GeneID,]

df_tl <- df_tl[df_tl$GeneID %in% df_count_2$GeneID,]
#Combine multiple transcript lengths for one gene by taking the mean value
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
write.csv(df1, "/home/yue1118/Bioskryb_data_analysis/Bioskryb_320DM_HSR_naive_nonsorted_RNA_Salmon_TPM_filtered_out_06042025.csv", row.names=TRUE)
df1_log <- log(df1+1)
write.csv(df1_log, "/home/yue1118/Bioskryb_data_analysis/Bioskryb_320DM_HSR_naive_nonsorted_RNA_Salmon_TPM_filtered_out_06042025.csv", row.names = TRUE)


# Create a Seurat object.
seurat_DM_HSR <- CreateSeuratObject(
  counts = df1_log, 
  project = "Naive_320DM_HSR", 
  assay = "RNA"
)

cell_names <- Cells(seurat_DM_HSR)

# Create a new metadata column based on pattern in cell names
source_label <- ifelse(grepl("DM", cell_names), "COLO320DM",
                       ifelse(grepl("HS", cell_names), "COLO320HSR", NA))

# Add to Seurat object's metadata
seurat_DM_HSR$DM_or_HSR <- source_label
seurat_DM_HSR$Cell <- colnames(seurat_DM_HSR)
seurat_DM_HSR$barcodes <- colnames(seurat_DM_HSR)
seurat_DM_HSR$Cell <- gsub("\\.", "-", seurat_DM_HSR$Cell)
seurat_DM_HSR$Cell <- sapply(strsplit(seurat_DM_HSR$Cell, "-R-"), `[`, 2)
seurat_DM_HSR$Cell <- sapply(strsplit(seurat_DM_HSR$Cell, "-_"), `[`, 1)
View(seurat_DM_HSR@meta.data)

#Add chr8 ecDNA region count to the metadata
df_counts_DM <- read.csv("/home/yue1118/Bioskryb_data_analysis/COLO320DM_Naive_nonsorted_Bioskryb_ecDNA_chr8-126598157-127691734_bin_counts_06042025.csv")
df_counts_HSR <- read.csv("/home/yue1118/Bioskryb_data_analysis/COLO320HSR_Naive_nonsorted_Bioskryb_ecDNA_chr8-126598157-127691734_bin_counts_06042025.csv")
df_DNA_counts <- rbind(df_counts_DM, df_counts_HSR)
df_DNA_counts$Cell <- sapply(strsplit(df_DNA_counts$Cell, "-D-"), `[`, 2)
df_DNA_counts$Cell <- sapply(strsplit(df_DNA_counts$Cell, "-_"), `[`, 1) #Double check this one!

df_metadata <- seurat_DM_HSR@meta.data
df_metadata$barcodes <- rownames(df_metadata)
df_metadata_2 <- merge(df_metadata, df_DNA_counts, by.x="Cell", by.y="Cell", all.x=TRUE)
View(df_metadata_2)

df_reordered <- df_metadata_2[match(colnames(seurat_DM_HSR), df_metadata_2$barcodes), ]
View(df_reordered)

table(df_reordered$barcodes==colnames(seurat_DM_HSR))
seurat_DM_HSR@meta.data <- df_reordered

#saveRDS(seurat_DM_HSR, "/home/yue1118/Bioskryb_data_analysis/Bioskryb_COLO320DM_HSR_naive_nonsorted_RNA_seurat_with_DNA_counts_06042025.rds")

seurat_DM_HSR <- readRDS("/home/yue1118/Bioskryb_data_analysis/Bioskryb_COLO320DM_HSR_naive_nonsorted_RNA_seurat_with_DNA_counts_06042025.rds")


#### To compare gene expression between HSR and DM  #####

list_ecDNA_genes <- c("MYC", "PVT1", "LINC02912", "PCAT1", "CASC19", "CASC11", "MIR1204", "POU5F1B", "CASC8","PRNCR1", "CCAT2",
                      "CDX2", "FAM84B", "RNU6-869P", "PLUT", "LINC00543", "PDX1", "ATP5EP2","CCAT1", "CASC21","TMEM75")

genes_in_data <- intersect(list_ecDNA_genes, rownames(seurat_DM_HSR))

# Step 2: Extract the expression matrix and subset to these genes
expr_mat <- as.data.frame(seurat_DM_HSR@assays$RNA@data[genes_in_data, ])

# Step 3: Calculate mean expression per cell (column-wise mean)
mean_expr_per_cell <- as.data.frame(colMeans(expr_mat))
mean_expr_per_cell$cell_name <- rownames(mean_expr_per_cell)
df_PVT1_MYC <- as.data.frame(t(expr_mat[c('PVT1','MYC'),]))
mean_expr_per_cell$PVT1 <- df_PVT1_MYC$PVT1
mean_expr_per_cell$MYC <- df_PVT1_MYC$MYC

# Step 4: Add gene expressions to metadata
seurat_DM_HSR$ecDNA_signature <- mean_expr_per_cell$`colMeans(expr_mat)`
seurat_DM_HSR$PVT1 <- mean_expr_per_cell$PVT1
seurat_DM_HSR$MYC <- mean_expr_per_cell$MYC
#seurat_DM_HSR$MYC <- seurat_DM_HSR@assays$RNA@data['MYC',]
#seurat_DM_HSR$PCAT1 <- seurat_DM_HSR@assays$RNA@data['PCAT1',]
#seurat_DM_HSR$CASC19 <- seurat_DM_HSR@assays$RNA@data['CASC19',]
#seurat_DM_HSR$CASC11 <- seurat_DM_HSR@assays$RNA@data['CASC11',]
#seurat_DM_HSR$POU5F1B <- seurat_DM_HSR@assays$RNA@data['POU5F1B',]
#seurat_DM_HSR$CDX2 <- seurat_DM_HSR@assays$RNA@data['CDX2',]

df_metadata <- seurat_DM_HSR@meta.data
View(df_metadata)


#To compare chr8_region_CN copy number between DM and HSR:
ggplot(df_metadata, aes(x = `DM_or_HSR`, y = `chr8_region_CN`)) +
  geom_boxplot() +
  labs(x = "Grouping Variable", y = "Value", title = "Boxplot by Group") +
  theme_minimal()
ggsave('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/COLO320DM_vs_HSR_chr8_ecDNA_bin_CN_06052025.svg', width = 6, height = 6)

wilcox.test(df_metadata[df_metadata$DM_or_HSR =="COLO320DM","chr8_region_CN"], df_metadata[df_metadata$DM_or_HSR =="COLO320HSR","chr8_region_CN"], paired = F) #p-value < 2.2e-16


#To compare mean ecDNA gene expression between DM and HSR:
ggplot(df_metadata, aes(x = `DM_or_HSR`, y = `ecDNA_signature`)) +
  geom_boxplot() +
  labs(x = "Grouping Variable", y = "Value", title = "Boxplot by Group") +
  theme_minimal()
ggsave('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/COLO320DM_vs_HSR_chr8_ecDNA_mean_expression_06052025.svg', width = 6, height = 6)

wilcox.test(df_metadata[df_metadata$DM_or_HSR =="COLO320DM","ecDNA_signature"], df_metadata[df_metadata$DM_or_HSR =="COLO320HSR","ecDNA_signature"], paired = F) #p-value < 2.2e-16

#To compare PVT1 gene expression between DM and HSR:
ggplot(df_metadata, aes(x = `DM_or_HSR`, y = `PVT1`)) +
  geom_boxplot() +
  labs(x = "Grouping Variable", y = "Value", title = "Boxplot by Group") +
  theme_minimal()
ggsave('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/COLO320DM_vs_HSR_chr8_PVT1_expression_06162025.svg', width = 6, height = 6)

wilcox.test(df_metadata[df_metadata$DM_or_HSR =="COLO320DM","PVT1"], df_metadata[df_metadata$DM_or_HSR =="COLO320HSR","PVT1"], paired = F) #p-value < 2.2e-16


#To plot scatterplot for chr8_region_CN and ecDNA mean expression between DM and HSR
ggplot(df_metadata, aes(x = chr8_region_CN, y = ecDNA_signature, color = DM_or_HSR)) +
  geom_point(size=3) +
  labs(x = "chr8_region_CN", y = "ecDNA mean expression", color = "Group", title = "Scatterplot") +
  theme_minimal()
ggsave('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/COLO320DM_vs_HSR_scatterplot_chr8_bin_CN_vs_ecDNA_mean_expression_06052025.svg', width = 10, height = 6)

#To plot scatterplot for chr8_region_CN and PVT1 expression between DM and HSR
ggplot(df_metadata, aes(x = chr8_region_CN, y = PVT1, color = DM_or_HSR)) +
  geom_point(size=3) +
  labs(x = "chr8_region_CN", y = "PVT1 expression", color = "Group", title = "Scatterplot") +
  theme_minimal()
ggsave('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/COLO320DM_vs_HSR_scatterplot_chr8_bin_CN_vs_PVT1_expression_06052025.svg', width = 10, height = 6)


#To compare gene expression between ecDNA and HSR on a comparable gene dosage, we notice that most of CN of HSR is below, so:
df_metadata_125 <- df_metadata[!is.na(df_metadata$chr8_region_CN) & (df_metadata$chr8_region_CN < 125),]


#We first make sure the gene dosage is comparable between HSR and DM with CN below 120.
ggplot(df_metadata_125, aes(x = `DM_or_HSR`, y = `chr8_region_CN`)) +
  geom_boxplot() +
  labs(x = "Grouping Variable", y = "Value", title = "Boxplot by Group") +
  theme_minimal()
ggsave('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/COLO320DM_vs_HSR_chr8_ecDNA_bin_CN_below_125_06052025.svg', width =6, height = 6)

wilcox.test(df_metadata_125[df_metadata_125$DM_or_HSR =="COLO320DM","chr8_region_CN"], df_metadata_125[df_metadata_125$DM_or_HSR =="COLO320HSR","chr8_region_CN"], paired = F) #p-value = 0.928


#We first compare the mean level:
ggplot(df_metadata_125, aes(x = `DM_or_HSR`, y = `ecDNA_signature`)) +
  geom_boxplot() +
  labs(x = "Grouping Variable", y = "log(TPM+1)", title = "Boxplot by Group") +
  theme_minimal()
ggsave('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/COLO320DM_vs_HSR_chr8_ecDNA_mean_expression_for_cells_with_CN_below_125_06052025.svg', width =6, height = 6)

wilcox.test(df_metadata_125[df_metadata_125$DM_or_HSR =="COLO320DM","ecDNA_signature"], df_metadata_125[df_metadata_125$DM_or_HSR =="COLO320HSR","ecDNA_signature"], paired = F) #p-value = 0.3009
t.test(df_metadata_125[df_metadata_125$DM_or_HSR =="COLO320DM","ecDNA_signature"], df_metadata_125[df_metadata_125$DM_or_HSR =="COLO320HSR","ecDNA_signature"], paired = F) #p-value = 0.3637

#We noticed that on the mean level, no difference in gene expression from a comparable gene dosage. Then we can compare individual ecDNA gene expression:

ggplot(df_metadata_125, aes(x = `DM_or_HSR`, y = `MYC`)) +
  geom_boxplot() +
  labs(x = "Grouping Variable", y = "log(TPM+1)", title = "Boxplot by Group") +
  theme_minimal()
ggsave('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/COLO320DM_vs_HSR_MYC_expression_for_cells_with_CN_below_125_06052025.svg', width =6, height = 6)

wilcox.test(df_metadata_125[df_metadata_125$DM_or_HSR =="COLO320DM","MYC"], df_metadata_125[df_metadata_125$DM_or_HSR =="COLO320HSR","MYC"], paired = F) #p-value = 0.5008


ggplot(df_metadata_125, aes(x = `DM_or_HSR`, y = `PVT1`)) +
  geom_boxplot() +
  labs(x = "Grouping Variable", y = "log(TPM+1)", title = "Boxplot by Group") +
  theme_minimal()
wilcox.test(df_metadata_125[df_metadata_125$DM_or_HSR =="COLO320DM","PVT1"], df_metadata_125[df_metadata_125$DM_or_HSR =="COLO320HSR","PVT1"], paired = F) #p-value = 3.649e-07
ggsave('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/COLO320DM_vs_HSR_PVT1_expression_for_cells_with_CN_below_125_06052025.svg', width =6, height = 6)


ggplot(df_metadata_125, aes(x = `DM_or_HSR`, y = `CASC19`)) +
  geom_boxplot() +
  labs(x = "Grouping Variable", y = "log(TPM+1)", title = "Boxplot by Group") +
  theme_minimal()
wilcox.test(df_metadata_125[df_metadata_125$DM_or_HSR =="COLO320DM","CASC19"], df_metadata_125[df_metadata_125$DM_or_HSR =="COLO320HSR","CASC19"], paired = F) #p-value = 4.67e-05
ggsave('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/COLO320DM_vs_HSR_CASC19_expression_for_cells_with_CN_below_125_06052025.svg', width =6, height = 6)


ggplot(df_metadata_125, aes(x = `DM_or_HSR`, y = `CASC11`)) +
  geom_boxplot() +
  labs(x = "Grouping Variable", y = "log(TPM+1)", title = "Boxplot by Group") +
  theme_minimal()
wilcox.test(df_metadata_125[df_metadata_125$DM_or_HSR =="COLO320DM","CASC11"], df_metadata_125[df_metadata_125$DM_or_HSR =="COLO320HSR","CASC11"], paired = F) #p-value = 0.03532
ggsave('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/COLO320DM_vs_HSR_CASC11_expression_for_cells_with_CN_below_125_06052025.svg', width =6, height = 6)

ggplot(df_metadata_125, aes(x = `DM_or_HSR`, y = `POU5F1B`)) +
  geom_boxplot() +
  labs(x = "Grouping Variable", y = "log(TPM+1)", title = "Boxplot by Group") +
  theme_minimal()
wilcox.test(df_metadata_125[df_metadata_125$DM_or_HSR =="COLO320DM","POU5F1B"], df_metadata_125[df_metadata_125$DM_or_HSR =="COLO320HSR","POU5F1B"], paired = F) #p-value = 0.0007611
ggsave('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/COLO320DM_vs_HSR_POU5F1B_expression_for_cells_with_CN_below_125_06052025.svg', width =6, height = 6)




###Now to explore the linear relationship of gene dosage and gene expression between HSR and DM
#For HSR
df_metadata_HSR <- df_metadata[df_metadata$DM_or_HSR == "COLO320HSR",]
df_metadata_HSR <- df_metadata_HSR[order(-df_metadata_HSR$ecDNA_signature), ]
svglite("/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/COLO320HSR_inferred_ecDNA_bin_counts_vs_mean_ecDNA_gene_exp_06052025.svg", width = 8, height = 6)
plot(df_metadata_HSR$chr8_region_CN, df_metadata_HSR$ecDNA_signature)
model <- lm(as.vector(df_metadata_HSR$ecDNA_signature) ~ as.vector(df_metadata_HSR$chr8_region_CN))
summary_model <- summary(model)
summary_model$r.squared #0.03257967
summary_model$coefficients #P-value=2.606288e-02
abline(lm(as.vector(df_metadata_HSR$ecDNA_signature) ~ as.vector(df_metadata_HSR$chr8_region_CN)), col = "hotpink2", lwd = 4)
dev.off()
cor_test <- cor.test(as.vector(df_metadata_HSR$ecDNA_signature), as.vector(df_metadata_HSR$chr8_region_CN)) 
cor_test$estimate #pearson correlation =0.1804984  
cor_test$p.value #p-value=0.02606288

#For HSR PVT1 gene
df_metadata_HSR <- df_metadata[df_metadata$DM_or_HSR == "COLO320HSR",]
df_metadata_HSR <- df_metadata_HSR[order(-df_metadata_HSR$PVT1), ]
svglite("/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/COLO320HSR_inferred_ecDNA_bin_counts_vs_PVT1_gene_exp_06162025.svg", width = 8, height = 6)
plot(df_metadata_HSR$chr8_region_CN, df_metadata_HSR$PVT1)
model <- lm(as.vector(df_metadata_HSR$PVT1) ~ as.vector(df_metadata_HSR$chr8_region_CN))
summary_model <- summary(model)
summary_model$r.squared #0.004450119
summary_model$coefficients #P-value=4.141758e-01
abline(lm(as.vector(df_metadata_HSR$PVT1) ~ as.vector(df_metadata_HSR$chr8_region_CN)), col = "hotpink2", lwd = 4)
dev.off()
cor_test <- cor.test(as.vector(df_metadata_HSR$PVT1), as.vector(df_metadata_HSR$chr8_region_CN)) 
cor_test$estimate #pearson correlation =0.06670921  
cor_test$p.value #p-value=0.4141758

#For MYC in HSR
df_metadata_HSR <- df_metadata[df_metadata$DM_or_HSR == "COLO320HSR",]
df_metadata_HSR <- df_metadata_HSR[order(-df_metadata_HSR$MYC), ]
svglite("/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/COLO320HSR_inferred_ecDNA_bin_counts_vs_MYC_gene_exp_06182025.svg", width = 8, height = 6)
plot(df_metadata_HSR$chr8_region_CN, df_metadata_HSR$MYC)
model <- lm(as.vector(df_metadata_HSR$MYC) ~ as.vector(df_metadata_HSR$chr8_region_CN))
summary_model <- summary(model)
summary_model$r.squared #0.03603071
summary_model$coefficients #1.916809e-02
abline(lm(as.vector(df_metadata_HSR$MYC) ~ as.vector(df_metadata_HSR$chr8_region_CN)), col = "hotpink2", lwd = 4)
dev.off()
cor_test <- cor.test(as.vector(df_metadata_HSR$MYC), as.vector(df_metadata_HSR$chr8_region_CN)) 
cor_test$estimate #pearson correlation =0.1898176 
cor_test$p.value #p-value=0.01916809


#For DM:
df_metadata_DM <- df_metadata[df_metadata$DM_or_HSR == "COLO320DM",]
df_metadata_DM <- df_metadata_DM[order(-df_metadata_DM$ecDNA_signature), ]
svglite("/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/COLO320DM_inferred_ecDNA_bin_counts_vs_mean_ecDNA_gene_exp_06052025.svg", width = 8, height = 6)
plot(df_metadata_DM$chr8_region_CN, df_metadata_DM$ecDNA_signature)
model <- lm(as.vector(df_metadata_DM$ecDNA_signature) ~ as.vector(df_metadata_DM$chr8_region_CN))
summary_model <- summary(model)
summary_model$r.squared #0.2079768
summary_model$coefficients #P-value=3.541817e-09
abline(lm(as.vector(df_metadata_DM$ecDNA_signature) ~ as.vector(df_metadata_DM$chr8_region_CN)), col = "hotpink2", lwd = 4)
dev.off()
cor_test <- cor.test(as.vector(df_metadata_DM$ecDNA_signature), as.vector(df_metadata_DM$chr8_region_CN)) 
cor_test$estimate #pearson correlation =0.4560448  
cor_test$p.value #p-value=3.541817e-09


#For PVT1 in DM
df_metadata_DM <- df_metadata[df_metadata$DM_or_HSR == "COLO320DM",]
df_metadata_DM <- df_metadata_DM[order(-df_metadata_DM$PVT1), ]
svglite("/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/COLO320DM_inferred_ecDNA_bin_counts_vs_PVT1_gene_exp_06162025.svg", width = 8, height = 6)
plot(df_metadata_DM$chr8_region_CN, df_metadata_DM$PVT1)
model <- lm(as.vector(df_metadata_DM$PVT1) ~ as.vector(df_metadata_DM$chr8_region_CN))
summary_model <- summary(model)
summary_model$r.squared #0.3368367
summary_model$coefficients #4.629581e-15
abline(lm(as.vector(df_metadata_DM$PVT1) ~ as.vector(df_metadata_DM$chr8_region_CN)), col = "hotpink2", lwd = 4)
dev.off()
cor_test <- cor.test(as.vector(df_metadata_DM$PVT1), as.vector(df_metadata_DM$chr8_region_CN)) 
cor_test$estimate #pearson correlation =0.5803763  
cor_test$p.value #p-value=4.629581e-15


#For MYC in DM
df_metadata_DM <- df_metadata[df_metadata$DM_or_HSR == "COLO320DM",]
df_metadata_DM <- df_metadata_DM[order(-df_metadata_DM$MYC), ]
svglite("/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/COLO320DM_inferred_ecDNA_bin_counts_vs_MYC_gene_exp_06182025.svg", width = 8, height = 6)
plot(df_metadata_DM$chr8_region_CN, df_metadata_DM$MYC)
model <- lm(as.vector(df_metadata_DM$MYC) ~ as.vector(df_metadata_DM$chr8_region_CN))
summary_model <- summary(model)
summary_model$r.squared #0.2031773
summary_model$coefficients #5.633692e-09
abline(lm(as.vector(df_metadata_DM$MYC) ~ as.vector(df_metadata_DM$chr8_region_CN)), col = "hotpink2", lwd = 4)
dev.off()
cor_test <- cor.test(as.vector(df_metadata_DM$MYC), as.vector(df_metadata_DM$chr8_region_CN)) 
cor_test$estimate #pearson correlation =0.4507519  
cor_test$p.value #p-value=5.633692e-09

#For MYC in DM 125
df_metadata_DM_125 <- df_metadata_125[df_metadata_125$DM_or_HSR == "COLO320DM",]
svglite("/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/COLO320DM_CN_less_than_125_inferred_ecDNA_bin_counts_vs_MYC_gene_exp_06182025.svg", width = 8, height = 6)
plot(df_metadata_DM_125$chr8_region_CN, df_metadata_DM_125$MYC)
model <- lm(as.vector(df_metadata_DM_125$MYC) ~ as.vector(df_metadata_DM_125$chr8_region_CN))
summary_model <- summary(model)
summary_model$r.squared #0.08475569
summary_model$coefficients #1.406712e-01
abline(lm(as.vector(df_metadata_DM_125$MYC) ~ as.vector(df_metadata_DM_125$chr8_region_CN)), col = "hotpink2", lwd = 4)
dev.off()
cor_test <- cor.test(as.vector(df_metadata_DM_125$MYC), as.vector(df_metadata_DM_125$chr8_region_CN)) 
cor_test$estimate #pearson correlation =0.2911283  
cor_test$p.value #p-value=0.14


#For PVT1 in DM 125
df_metadata_DM_125 <- df_metadata_125[df_metadata_125$DM_or_HSR == "COLO320DM",]
svglite("/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/COLO320DM_CN_less_than_125_inferred_ecDNA_bin_counts_vs_PVT1_gene_exp_06182025.svg", width = 8, height = 6)
plot(df_metadata_DM_125$chr8_region_CN, df_metadata_DM_125$PVT1)
model <- lm(as.vector(df_metadata_DM_125$PVT1) ~ as.vector(df_metadata_DM_125$chr8_region_CN))
summary_model <- summary(model)
summary_model$r.squared #0.1705434
summary_model$coefficients #0.0322806970
abline(lm(as.vector(df_metadata_DM_125$PVT1) ~ as.vector(df_metadata_DM_125$chr8_region_CN)), col = "hotpink2", lwd = 4)
dev.off()
cor_test <- cor.test(as.vector(df_metadata_DM_125$PVT1), as.vector(df_metadata_DM_125$chr8_region_CN)) 
cor_test$estimate #pearson correlation =0.412969   
cor_test$p.value #p-value=0.0322


#For MYC in HSR 125
df_metadata_HSR_125 <- df_metadata_125[df_metadata_125$DM_or_HSR == "COLO320HSR",]
svglite("/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/COLO320HSR_CN_less_than_125_inferred_ecDNA_bin_counts_vs_MYC_gene_exp_06182025.svg", width = 8, height = 6)
plot(df_metadata_HSR_125$chr8_region_CN, df_metadata_HSR_125$MYC)
model <- lm(as.vector(df_metadata_HSR_125$MYC) ~ as.vector(df_metadata_HSR_125$chr8_region_CN))
summary_model <- summary(model)
summary_model$r.squared #0.04485201
summary_model$coefficients #9.518467e-03
abline(lm(as.vector(df_metadata_HSR_125$MYC) ~ as.vector(df_metadata_HSR_125$chr8_region_CN)), col = "hotpink2", lwd = 4)
dev.off()
cor_test <- cor.test(as.vector(df_metadata_HSR_125$MYC), as.vector(df_metadata_HSR_125$chr8_region_CN)) 
cor_test$estimate #pearson correlation =0.2117829 
cor_test$p.value #p-value=0.01


#For PVT1 in DM 125
df_metadata_HSR_125 <- df_metadata_125[df_metadata_125$DM_or_HSR == "COLO320HSR",]
svglite("/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/COLO320HSR_CN_less_than_125_inferred_ecDNA_bin_counts_vs_PVT1_gene_exp_06182025.svg", width = 8, height = 6)
plot(df_metadata_HSR_125$chr8_region_CN, df_metadata_HSR_125$PVT1)
model <- lm(as.vector(df_metadata_HSR_125$PVT1) ~ as.vector(df_metadata_HSR_125$chr8_region_CN))
summary_model <- summary(model)
summary_model$r.squared #0.08572693
summary_model$coefficients # 0.0002905088
abline(lm(as.vector(df_metadata_HSR_125$PVT1) ~ as.vector(df_metadata_HSR_125$chr8_region_CN)), col = "hotpink2", lwd = 4)
dev.off()
cor_test <- cor.test(as.vector(df_metadata_HSR_125$PVT1), as.vector(df_metadata_HSR_125$chr8_region_CN)) 
cor_test$estimate #pearson correlation =0.2927916 
cor_test$p.value #p-value=0.0002905088


#Consistent with 10X multiome data volcano plots, the gene dosage vs expression correlation is much stronger in DM model than HSR.
#On the comparable gene dosage level, ecDNA gene expression overall on the mean level is not significantly different between ecDNA and HSR. Interestingly, some gene is the same between
#HSR and DM like MYC; other genes are higher in ecDNA, like PVT1; other genes are higher in HSR, like ...
#Overall, I think ecDNA differentiates from HSR in the way that it raise a much higher level of continuous variaion in gene dosage, that lead to higher RNA level and wider RNA distribution.




####To assign ecDNA count bins to single cells and do DEG analysis in COLO320DM cells###
df_metadata <- seurat_DM_HSR@meta.data

# Compute the quantiles for 30% and 70%
df_metadata_DM <- df_metadata[df_metadata$DM_or_HSR == "COLO320DM",]
df_metadata_HSR <- df_metadata[df_metadata$DM_or_HSR == "COLO320HSR",]

q30_dm <- quantile(df_metadata_DM$chr8_region_CN, 0.30, na.rm = TRUE)
q70_dm <- quantile(df_metadata_DM$chr8_region_CN, 0.70, na.rm = TRUE)

# Assign "high", "med", or "low" based on the quantiles
values <- seurat_DM_HSR@meta.data[["chr8_region_CN"]]
seurat_DM_HSR@meta.data$ecDNA_bin <- ifelse(
  values >= q70_dm, "high",
  ifelse(values <= q30_dm, "low", "med")
)

View(seurat_DM_HSR@meta.data)

seurat_DM_HSR@meta.data$cell_line_ecDNA_bin <- paste0(seurat_DM_HSR@meta.data$DM_or_HSR, seurat_DM_HSR@meta.data$ecDNA_bin)

Idents(seurat_DM_HSR) <- "cell_line_ecDNA_bin"

DefaultAssay(seurat_DM_HSR) <- "RNA"
df_dm_markers <- FindMarkers(seurat_DM_HSR, ident.1 = "COLO320DMhigh", ident.2 = "COLO320DMlow", min.pct = 0.25, logfc.threshold = 0.10, recorrect_umi=F)


View(df_dm_markers[df_dm_markers$p_val_adj < 0.05,])


#To create a function to make volcano plots for differentially expressed genes after treatment
volcano_plot_function <- function(input_df) {
  # Convert p-value to -log10 p-value for the volcano plot
  input_df$neg_log10_pvalue <- -log10(input_df$p_val_adj)
  
  # Add a column to classify significance (adjust thresholds as needed)
  input_df$significant <- ifelse(input_df$avg_log2FC > 0 & input_df$p_val_adj < 0.05, "Upregulated",
                                 ifelse(input_df$avg_log2FC < 0 & input_df$p_val_adj < 0.05, "Downregulated", "Not Significant"))
  
  # Create the volcano plot
  ggplot(input_df, aes(x = avg_log2FC, y = neg_log10_pvalue)) +
    geom_point(aes(color = significant), size=5) + # Color points by significance
    scale_color_manual(values = c("blue", "gray", "red")) + # Custom colors
    labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-log10(adj p-value)") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") + # Add dashed lines for log2FC thresholds
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") + 
    theme(
      panel.grid.major.x = element_blank(),               # Remove major grid lines on x-axis if needed
      panel.grid.major.y = element_blank()                # Remove major grid lines on y-axis if needed
    )# Horizontal line for p-value threshold
}

volcano_plot_function(df_dm_markers)
ggsave("/home/yue1118/ecDNA_project_data_and_figures/FACS_FISH_FULL_figures/Bioskryb_COLO320DM_naive_DEGs_between_MYC_bin_CN_high_vs_low_06062025.svg", width = 8, height = 6)



####To assign ecDNA count bins to single cells and do DEG analysis in COLO320HSR cells###
df_metadata_HSR <- df_metadata[df_metadata$DM_or_HSR == "COLO320HSR",]

q30_hsr <- quantile(df_metadata_HSR$chr8_region_CN, 0.30, na.rm = TRUE)
q70_hsr <- quantile(df_metadata_HSR$chr8_region_CN, 0.70, na.rm = TRUE)

# Assign "high", "med", or "low" based on the quantiles
values <- seurat_DM_HSR@meta.data[["chr8_region_CN"]]
seurat_DM_HSR@meta.data$ecDNA_bin <- ifelse(
  values >= q70_hsr, "high",
  ifelse(values <= q30_hsr, "low", "med")
)


seurat_DM_HSR@meta.data$cell_line_ecDNA_bin <- paste0(seurat_DM_HSR@meta.data$DM_or_HSR, seurat_DM_HSR@meta.data$ecDNA_bin)

Idents(seurat_DM_HSR) <- "cell_line_ecDNA_bin"

DefaultAssay(seurat_DM_HSR) <- "RNA"
df_hsr_markers <- FindMarkers(seurat_DM_HSR, ident.1 = "COLO320HSRhigh", ident.2 = "COLO320HSRlow", min.pct = 0.25, logfc.threshold = 0.10, recorrect_umi=F)


View(df_hsr_markers[df_hsr_markers$p_val_adj < 0.05,])



volcano_plot_function(df_hsr_markers)
ggsave("/home/yue1118/ecDNA_project_data_and_figures/FACS_FISH_FULL_figures/Bioskryb_COLO320HSR_naive_DEGs_between_MYC_bin_CN_high_vs_low_06062025.svg", width = 8, height = 6)


#To perform DEG between DM and HSR cells with ecDNA bin copy number below 125:
seurat_DM_HSR <- readRDS("/home/yue1118/Bioskryb_data_analysis/Bioskryb_COLO320DM_HSR_naive_nonsorted_RNA_seurat_with_DNA_counts_06042025.rds")

# Assign "high", "med", or "low" based on the quantiles
values <- seurat_DM_HSR@meta.data[["chr8_region_CN"]]
seurat_DM_HSR@meta.data$ecDNA_bin <- ifelse(
  values <= 125, "low", "others"
)


seurat_DM_HSR@meta.data$cell_line_ecDNA_bin <- paste0(seurat_DM_HSR@meta.data$DM_or_HSR, seurat_DM_HSR@meta.data$ecDNA_bin)

Idents(seurat_DM_HSR) <- "cell_line_ecDNA_bin"
cells_1 <- WhichCells(seurat_DM_HSR, idents = "COLO320DMlow")
cells_2 <- WhichCells(seurat_DM_HSR, idents = "COLO320HSRlow")

length(cells_1)
length(cells_2)

# Then run FindMarkers using vectors of cell names
df_marker <- FindMarkers(seurat_DM_HSR, ident.1 = cells_1, ident.2 = cells_2)
View(df_marker[df_marker$p_val_adj < 0.05,])

Idents(seurat_DM_HSR) <- "DM_or_HSR"
cells_1 <- WhichCells(seurat_DM_HSR, idents = "COLO320DM")
cells_2 <- WhichCells(seurat_DM_HSR, idents = "COLO320HSR")

df_marker_2 <- FindMarkers(seurat_DM_HSR, ident.1 = cells_1, ident.2 = cells_2)
df_marker_2$gene <- rownames(df_marker_2)
View(df_marker_2[df_marker_2$p_val_adj < 0.05,])


msig_hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, entrez_gene)
gene_df <- bitr(df_marker_2$gene, fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)
df_merged <- merge(df_marker_2, gene_df, by.x = "gene", by.y = "SYMBOL")

gene_list <- df_merged$avg_log2FC
names(gene_list) <- df_merged$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

# Run GSEA
gsea_result_dm_vs_hsr <- GSEA(geneList = gene_list,
                         TERM2GENE = msig_hallmark,
                         pvalueCutoff = 1,
                         verbose = FALSE)
View(as.data.frame(gsea_result_colo320_vs_hsr))

View(seurat_DM_HSR@meta.data)
DefaultAssay(seurat_DM_HSR) <- "RNA"
Idents(seurat_DM_HSR) <- "cell_line_ecDNA_bin"
df_dm_vs_hsr_markers <- FindMarkers(seurat_DM_HSR, ident.1 = "COLO320DMlow", ident.2 = "COLO320HSRlow", min.pct = 0.05, logfc.threshold = 0.10, recorrect_umi=F)

View(df_dm_vs_hsr_markers[df_dm_vs_hsr_markers$p_val_adj < 0.05,])



# Assume `seurat_obj` is your Seurat object
# And you've set identities (e.g., seurat_obj$group or Idents(seurat_obj))

# Specify the two groups you want to compare
group1 <- "COLO320HSRlow"
group2 <- "COLO320DMlow"

# Set identities to the grouping column if not done already
Idents(seurat_DM_HSR) <- "cell_line_ecDNA_bin"  # or your column name

# Get gene names
genes <- rownames(seurat_DM_HSR)

# Get expression matrix (can use data slot or counts depending on normalization)
expr_matrix <- as.data.frame(as.matrix(GetAssayData(seurat_DM_HSR, slot = "data"))) # log-normalized expression
df_metadata <- seurat_DM_HSR@meta.data
colnames(expr_matrix) <- df_metadata$barcodes
genes <- rownames(expr_matrix)
# Get cells in each group
cells_group1 <- df_metadata[df_metadata$DM_or_HSR == "COLO320DM",]$barcodes
cells_group2 <- df_metadata[df_metadata$DM_or_HSR == "COLO320HSR",]$barcodes

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

msig_hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, entrez_gene)
gene_df <- bitr(results$gene, fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)
df_merged <- merge(results, gene_df, by.x = "gene", by.y = "SYMBOL")

gene_list <- df_merged$mean_value_diff
names(gene_list) <- df_merged$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

# 3. Run GSEA
gsea_result_colo320_vs_hsr <- GSEA(geneList = gene_list,
                          TERM2GENE = msig_hallmark,
                          pvalueCutoff = 0.05,
                          verbose = FALSE)


View(as.data.frame(gsea_result_colo320_vs_hsr))



# To compute gene expression correlation of all the genes with PVT1 for DM cells with CN < 125 or just DM cells
seurat_DM_HSR <- readRDS("/home/yue1118/Bioskryb_data_analysis/Bioskryb_COLO320DM_HSR_naive_nonsorted_RNA_seurat_with_DNA_counts_06042025.rds")
expr_matrix <- as.data.frame(as.matrix(GetAssayData(seurat_DM_HSR, slot = "data")))
##expr_matrix_125 <- expr_matrix[, colnames(expr_matrix) %in% df_metadata_125$barcodes]
#View(expr_matrix_125)
expr_matrix_dm <- expr_matrix[,grepl("DM",colnames(expr_matrix))]
df_metadata <- seurat_DM_HSR@meta.data
View(df_metadata)


# Initialize results for PVT1 correlation
results_320_dm <- data.frame(
  gene = character(),
  pcc = numeric(),
  p_val = numeric(),
  stringsAsFactors = FALSE
)

# Get PVT1 expression vector
pvt1_expr_dm <- as.numeric(expr_matrix_dm["PVT1", ])
ecDNA_count_signature_320 <- as.numeric(df_metadata[df_metadata$DM_or_HSR == "COLO320DM",]$chr8_region_CN)

# Loop through each gene
for (gene in rownames(expr_matrix_dm)) {

  gene_expr <- as.numeric(expr_matrix_dm[gene, ])
  
  cor_result <- cor.test(gene_expr, ecDNA_count_signature_320, method = "pearson")
  
  results_320_dm <- rbind(results_320_dm, data.frame(
    gene = gene,
    pcc = cor_result$estimate,
    p_val = cor_result$p.value
  ))
}

# Adjust p-values for multiple testing
results_320_dm$p_val_adj <- p.adjust(results_320_dm$p_val, method = "BH")
results_320_dm_filtered <- results_320_dm[results_320_dm$p_val < 0.1,]

msig_hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, entrez_gene)
gene_df <- bitr(results_320_dm$gene, fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)
df_merged <- merge(results_320_dm, gene_df, by.x = "gene", by.y = "SYMBOL")

gene_list <- df_merged$pcc
names(gene_list) <- df_merged$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

# 3. Run GSEA
COLO320DM_gsea_result <- GSEA(geneList = gene_list,
                    TERM2GENE = msig_hallmark,
                    pvalueCutoff = 0.05,
                    verbose = FALSE)
dotplot(COLO32DM_gsea_result, showCategory = 15)
ggsave('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/GSEA_result_for_genes_correlate_with_PVT1_mRNA_in_COLO320DM_bioskryb_06252025.svg', width = 10, height = 10)

View(as.data.frame(COLO320DM_gsea_result))
COLO320DM_gsea_result <- as.data.frame(COLO320DM_gsea_result)
COLO320DM_gsea_result$Description <- factor(COLO320DM_gsea_result$Description, levels = COLO320DM_gsea_result$Description[order(COLO320DM_gsea_result$NES)])
COLO320DM_gsea_result <- COLO320DM_gsea_result[COLO320DM_gsea_result]

ggplot(COLO320DM_gsea_result, aes(x = Description, y = NES)) +
  geom_segment(aes(x = Description, xend = Description, y = 0, yend = NES), color = "grey50") +
  geom_point(aes(color = p.adjust), size = 5) +
  scale_color_gradient(low = "red", high = "blue", name = "Adjusted p-value") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Lollipop Plot of Pathways",
       x = NULL, y = "Normalized Enrichment Score (NES)")
ggsave('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/lollipop_GSEA_result_for_genes_correlate_with_ecDNA_counts_in_COLO320DM_bioskryb_06262025.svg', width = 6, height = 6)




# To compute gene expression correlation of all the genes with PVT1 or amp counts for HSR cells
seurat_DM_HSR <- readRDS("/home/yue1118/Bioskryb_data_analysis/Bioskryb_COLO320DM_HSR_naive_nonsorted_RNA_seurat_with_DNA_counts_06042025.rds")
expr_matrix <- as.data.frame(as.matrix(GetAssayData(seurat_DM_HSR, slot = "data")))

expr_matrix_hsr <- expr_matrix[,grepl("HS",colnames(expr_matrix))]
df_metadata <- seurat_DM_HSR@meta.data
View(df_metadata)


# Get PVT1 expression vector
pvt1_expr_hsr <- as.numeric(expr_matrix_hsr["PVT1", ])
HSR_count_signature_320 <- as.numeric(df_metadata[df_metadata$DM_or_HSR == "COLO320HSR",]$chr8_region_CN)

table(df_metadata[df_metadata$DM_or_HSR == "COLO320HSR",]$barcodes == colnames(expr_matrix_hsr))

# Initialize results for PVT1 correlation
results_hsr <- data.frame(
  gene = character(),
  pcc = numeric(),
  p_val = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each gene
for (gene in rownames(expr_matrix_hsr)) {

  gene_expr <- as.numeric(expr_matrix_hsr[gene, ])
  
  cor_result <- cor.test(gene_expr, HSR_count_signature_320, method = "pearson")
  
  results_hsr <- rbind(results_hsr, data.frame(
    gene = gene,
    pcc = cor_result$estimate,
    p_val = cor_result$p.value
  ))
}

# Adjust p-values for multiple testing
results_hsr$p_val_adj <- p.adjust(results_hsr$p_val, method = "BH")

msig_hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, entrez_gene)
gene_df <- bitr(results_hsr$gene, fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)
df_merged <- merge(results_hsr, gene_df, by.x = "gene", by.y = "SYMBOL")

gene_list <- df_merged$pcc
names(gene_list) <- df_merged$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

# 3. Run GSEA
COLO320HSR_gsea_result <- GSEA(geneList = gene_list,
                    TERM2GENE = msig_hallmark,
                    pvalueCutoff = 0.05,
                    verbose = FALSE)
dotplot(gsea_result, showCategory = 15)

View(as.data.frame(COLO320HSR_gsea_result))

COLO320HSR_gsea_result <- as.data.frame(COLO320HSR_gsea_result)
COLO320HSR_gsea_result$Description <- factor(COLO320HSR_gsea_result$Description, levels = COLO320HSR_gsea_result$Description[order(COLO320HSR_gsea_result$NES)])


ggplot(COLO320HSR_gsea_result, aes(x = Description, y = NES)) +
  geom_segment(aes(x = Description, xend = Description, y = 0, yend = NES), color = "grey50") +
  geom_point(aes(color = p.adjust), size = 5) +
  scale_color_gradient(low = "red", high = "blue", name = "Adjusted p-value") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Lollipop Plot of Pathways",
       x = NULL, y = "Normalized Enrichment Score (NES)")
ggsave('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/lollipop_GSEA_result_for_genes_correlate_with_ecDNA_counts_in_COLO320HSR_bioskryb_06262025.svg', width = 6, height = 6)


View(results_dm[(results_dm$p_val_adj < 0.05) & (results_dm$pcc >0),])
View(results_dm[(results_dm$p_val_adj < 0.05) & (results_dm$pcc <0),])

write.csv(results_dm[(results_dm$p_val_adj < 0.05) & (results_dm$pcc >0),], "/home/yue1118/Bioskryb_data_analysis/COLO320DM_bioskryb_PVT1_exp_correlation_across_genes_06252025.csv")

gene = "SERBP1"
gene_expr <- as.numeric(expr_matrix_dm_125[gene, ])
plot(gene_expr, pvt1_expr_dm)

# Initialize results for MYC mRNA correlation
results_dm <- data.frame(
  gene = character(),
  pcc = numeric(),
  p_val = numeric(),
  stringsAsFactors = FALSE
)

# Get PVT1 expression vector
myc_expr_dm <- as.numeric(expr_matrix_dm["MYC", ])

# Loop through each gene
for (gene in rownames(expr_matrix_dm)) {
  if (gene == "MYC") next  # skip self-correlation
  
  gene_expr <- as.numeric(expr_matrix_dm[gene, ])
  
  cor_result <- cor.test(gene_expr, myc_expr_dm, method = "pearson")
  
  results_dm <- rbind(results_dm, data.frame(
    gene = gene,
    pcc = cor_result$estimate,
    p_val = cor_result$p.value
  ))
}

# Adjust p-values for multiple testing
results_dm$p_val_adj <- p.adjust(results_dm$p_val, method = "bonferroni")
results_dm <- results_dm[order(results_dm$p_val_adj), ]


View(results_dm[(results_dm$p_val_adj < 0.05) & (results_dm$pcc >0),])

write.csv(results_dm[(results_dm$p_val_adj < 0.05) & (results_dm$pcc >0),], "/home/yue1118/Bioskryb_data_analysis/COLO320DM_bioskryb_MYC_exp_correlation_across_genes_06182025.csv")

#To do the same analysis with HSR
expr_matrix_hsr <- expr_matrix[,grepl("HS",colnames(expr_matrix))]

# Initialize results
results_hsr <- data.frame(
  gene = character(),
  pcc = numeric(),
  p_val = numeric(),
  stringsAsFactors = FALSE
)

# Get PVT1 expression vector
pvt1_expr_hsr <- as.numeric(expr_matrix_hsr["PVT1", ])

# Loop through each gene
for (gene in rownames(expr_matrix_hsr)) {
  if (gene == "PVT1") next  # skip self-correlation
  
  gene_expr <- as.numeric(expr_matrix_hsr[gene, ])
  
  cor_result <- cor.test(gene_expr, pvt1_expr_hsr, method = "pearson")
  
  results_hsr <- rbind(results_hsr, data.frame(
    gene = gene,
    pcc = cor_result$estimate,
    p_val = cor_result$p.value
  ))
}

# Adjust p-values for multiple testing
results_hsr$p_val_adj <- p.adjust(results_hsr$p_val, method = "BH")

# Sort by correlation strength or p-value
results_hsr <- results_hsr[order(results_hsr$p_val_adj), ]


View(results_hsr[(results_hsr$p_val_adj < 0.05) & (results_hsr$pcc >0),])











seurat_DM_HSR <- readRDS("/home/yue1118/Bioskryb_data_analysis/Bioskryb_COLO320DM_HSR_naive_nonsorted_RNA_seurat_with_DNA_counts_06042025.rds")
df1 <- as.data.frame(as.matrix(seurat_DM_HSR@assays$RNA@data))
df1 <- as.data.frame(t(df1))
df1$ecDNA_or_HSR_count <- seurat_DM_HSR$chr8_region_CN
table(rownames(df1) == colnames(seurat_DM_HSR))
View(df1)
write.csv(df1, "/home/yue1118/Bioskryb_data_analysis/COLO32DM_HSR_naive_Bioskryb_RNA_logTPM_ecDNA_counts_dataframe_06162025.csv")
