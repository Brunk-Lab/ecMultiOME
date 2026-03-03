#Bioskryb data analysis for cell lines NCIH2170
library(ggplot2)
setwd('/home/yue1118')
library(plotly)
library(pheatmap)
library(stringr)
library(dplyr)
library(tidyr)

#FGFR2 gene GC content is 0.454076
df1 <- read.table("/home/yue1118/Bioskryb_SNU16_bam_files/hg38_genome_3Mb_gc_content.txt", header = FALSE, sep = "\t")
View(df1)
colnames(df1) <- c("chr", "start", "end", "AT_content", "GC_content", "num_A", "num_C", "num_G", "num_T", "num_N", "num_oth", "seq_len")

df1$FGFR_difference <- (df1$GC_content-0.454076)^2

df1 <- df1[order(df1$FGFR_difference), ]
df1 <- head(df1, 50)

#To get rid of potential ecDNA regions:
df1 <- df1[!df1$chr %in% c("chr10", "chr8", 'chr13'),]
View(df1)

write.table(df1[, 1:3], "/home/yue1118/Bioskryb_SNU16_bam_files/FGFR2_closest_3Mb_bins_output.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


df2 <- read.table("/home/yue1118/Bioskryb_SNU16_bam_files/FGFR2_coverage_per_cell.tsv", header = T)
df3 <- read.table("/home/yue1118/Bioskryb_SNU16_bam_files/mean_coverage_per_cell.tsv", header = T)
View(df3)

df4 <- merge(df2, df3, by.x="Cell", by.y="Cell")
View(df4)
df4$FGFR2_cn <- 4*(df4$MeanCoverage.x/df4$MeanCoverage.y)
df4 <- df4[df4$FGFR2_cn < 2000,] #get rid of outlier
hist(df4$FGFR2_cn)
df4$FGFR2_cn <- as.numeric(df4$FGFR2_cn)
mean(df4$FGFR2_cn)
median(df4$FGFR2_cn)
mean(as.vector(df4$FGFR2_cn), na.rm=TRUE)
median(as.vector(df4$FGFR2_cn), na.rm=TRUE)

df_fish <- read.csv("/home/yue1118/Bioskryb_data_analysis/EcDNA_FISH_SNU16_counts.csv")
View(df_fish)
df_fish_control <- df_fish[df_fish$Desc == "CNTR_ORIG_2022",]
mean(df_fish_control$FGFR2.)
median(df_fish_control$FGFR2.)
wilcox.test(df_fish_control$FGFR2., df4$FGFR2_cn)

View(df4)

df4$phenotype <- ifelse(grepl("Low", df4$Cell, ignore.case = TRUE), "Low",
                        ifelse(grepl("High", df4$Cell, ignore.case = TRUE), "High",
                               ifelse(grepl("Med", df4$Cell, ignore.case = TRUE), "Med", NA)))
df4$phenotype <- factor(df4$phenotype, levels = c("Low", "Med", "High"))

ggplot(df4, aes(x = `phenotype`, y = `FGFR2_cn`)) +
  geom_boxplot() +
  labs(x = "Grouping Variable", y = "Value", title = "Boxplot by Group") +
  theme_minimal()

df_ccle <- read.csv("/home/yue1118/Bioskryb_SNU16_bam_files/SNU16_Copy_Number_log2_relative_to_ploidy.csv")
colnames(df_ccle) <- c("gene", "CN")
mean(as.vector(df_ccle$CN), na.rm=T) #1.03 mean value
View(df_ccle)


####Ginkgo output####
rds1 <- readRDS("/home/yue1118/ecDNA_project_data_and_figures/ginkgo_res.binsize_1000000.RDS")

#To extract the normalized coverage numbers
df_coverage <- rds1$SegNorm
rownames(df_coverage) <- paste(df_coverage[,1], df_coverage[,2], df_coverage[,3], sep = "-")
df_coverage <- df_coverage[,4:ncol(df_coverage)]
df_coverage_SNU16 <- df_coverage[, grepl("SNU16", names(df_coverage))]

View(df_coverage_SNU16)

View(df_coverage_SNU16[grepl("chr8", rownames(df_coverage_SNU16)),])


coverage_quantiles <- apply(df_coverage_SNU16, 2, quantile, probs = 1.00, na.rm = TRUE)
# Function to filter values below the 100th percentile and compute median
filtered_coverage_means <- sapply(1:ncol(df_coverage_SNU16), function(i) {
  col_values <- df_coverage_SNU16[, i]  # Extract column values
  filtered_values <- col_values[col_values < coverage_quantiles[i]]  # Keep values below 100th percentile
  mean(filtered_values, na.rm = TRUE)  #We use the mean value as coverage in "normal" chromosomal regions #I think we can use this way to identity the coverage for regions without CNV on chromosome
})                                      #in that for DepMap SNU16 CN profile (Copy number log2 (relative to ploidy +1)), it centers at 1, which means that most of genes don't have CNV events.
#This further means that the mean coverage values in the chromosomal regions are the coverage for "normal" regions without CNV events.


# Convert result to a dataframe
df1 <- data.frame(Cell = names(df_coverage_SNU16), "mean_coverage" = filtered_coverage_means)
View(df1)

#Now extract the chr10 region that is on ecDNA: chr10-121134447-122166917, which covers FGFR2
df1$chr10_coverage <- as.vector(t(df_coverage_SNU16[c("chr10-121134447-122166917"),]))
df1$chr10_copy_number <- 4*df1$chr10_coverage/df1$mean_coverage
hist(df1$chr10_copy_number)
mean(df1$chr10_copy_number)
median(df1$chr10_copy_number)

wilcox.test(df1$chr10_copy_number, df_fish_control$FGFR2.)

df1$phenotype <- ifelse(grepl("Low", df1$Cell, ignore.case = TRUE), "Low",
                        ifelse(grepl("High", df1$Cell, ignore.case = TRUE), "High",
                               ifelse(grepl("Med", df1$Cell, ignore.case = TRUE), "Med", NA)))
df1$phenotype <- factor(df1$phenotype, levels = c("Low", "Med", "High"))

ggplot(df1, aes(x = `phenotype`, y = `chr10_copy_number`)) +
  geom_boxplot() +
  labs(x = "Grouping Variable", y = "Value", title = "Boxplot by Group") +
  theme_minimal()
ggsave('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/SNU16_FGFR2_ecDNA_counts_boxplot_high_med_low_05012025.svg', width = 6, height = 6)

wilcox.test(df1[df1$phenotype =="High","chr10_copy_number"], df1[df1$phenotype =="Low","chr10_copy_number"], paired = F) #p-value = 0.0008464
wilcox.test(df1[df1$phenotype =="High","chr10_copy_number"], df1[df1$phenotype =="Med","chr10_copy_number"], paired = F) #p-value = 0.3167
wilcox.test(df1[df1$phenotype =="Med","chr10_copy_number"], df1[df1$phenotype =="Low","chr10_copy_number"], paired = F) #p-value = 0.1905

write.csv(df1, "/home/yue1118/Bioskryb_data_analysis/SNU16_ecDNA_counts_04292025.csv")

df1 <- read.csv("/home/yue1118/Bioskryb_data_analysis/SNU16_ecDNA_counts_04292025.csv")

df2 <- data.frame(Cell = names(df_coverage_SNU16), "mean_coverage" = filtered_coverage_means)
View(df2)

#Now extract the chr8 region that is on ecDNA: chr8-127412678-128444056, which covers MYC and PVT1
df2$chr8_coverage <- as.vector(t(df_coverage_SNU16[c("chr8-127412678-128444056"),]))
df2$chr8_copy_number <- 4*df2$chr8_coverage/df2$mean_coverage

df3 <- merge(df1, df2, by.x= "Cell", by.y= "Cell")
write.csv(df3, "/home/yue1118/Bioskryb_data_analysis/SNU16_ecDNA_counts_FGFR2_species_and_MYC_species_06182025.csv")

svg('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/SNU16_plate1_chr10_FGFR2_CN_vs_chr8_MYC_CN_07292025.svg',width = 6, height = 6)
plot(df3$chr10_copy_number, df3$chr8_copy_number)
cor.test(df3$chr10_copy_number, df3$chr8_copy_number) #cor=0.59, p-value=7.607e-12
dev.off()

df1_SNU16 <- read.csv("/home/yue1118/Bioskryb_data_analysis/SNU16_ecDNA_counts_04292025.csv")
df1_SNU16 <- df1_SNU16[order(-df1_SNU16$chr10_copy_number), ]
#To plot FGFR2 copy number distribution across single cells
svg('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/FGFR2_SNU16_plate1_copy_number_barplot_05052025.svg',width = 6, height = 6)
barplot(df1_SNU16$chr10_copy_number)
dev.off()

#Compare with FISH imaging ecDNA counts
df_fish <- read.csv("/home/yue1118/Bioskryb_data_analysis/EcDNA_FISH_SNU16_counts.csv")
View(df_fish)
df_fish_control <- df_fish[df_fish$Desc == "CNTR_ORIG_2022",]
df_fish_control$technique <- "FISH_imaging"

df1$technique <- "Bioskryb_scDNA"
df1 <- df1[, c("chr10_copy_number","technique")]
colnames(df1) <- c("ecDNA_cn", "technique")
df2 <- df_fish_control[, c("FGFR2.", "technique")]
colnames(df2) <- c("ecDNA_cn", "technique")
df3 <- rbind(df1, df2)

ggplot(df3, aes(x = `technique`, y = `ecDNA_cn`)) +
  geom_boxplot() +
  labs(x = "Grouping Variable", y = "Value", title = "Boxplot by Group") +
  theme_minimal()
ggsave('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/FGFR2_ecDNA_counts_boxplots_Bioskryb_vs_FISH_SNU16_05012025.svg', width = 6, height = 6)


wilcox.test(df1$ecDNA_cn, df2$ecDNA_cn)#p-value = 0.3483



#Now extract the chr8 region that is on ecDNA: chr8-127214678-128444056, which covers MYC
View(df_coverage_SNU16[grepl("chr8", rownames(df_coverage_SNU16)), ])

df1 <- data.frame(Cell = names(df_coverage_SNU16), "mean_coverage" = filtered_coverage_means)
df1$chr8_coverage <- as.vector(t(df_coverage_SNU16[c("chr8-127412678-128444056"),]))
df1$chr8_copy_number <- 4*df1$chr8_coverage/df1$mean_coverage
hist(df1$chr8_copy_number)
mean(df1$chr8_copy_number)
median(df1$chr8_copy_number)


df1$phenotype <- ifelse(grepl("Low", df1$Cell, ignore.case = TRUE), "Low",
                        ifelse(grepl("High", df1$Cell, ignore.case = TRUE), "High",
                               ifelse(grepl("Med", df1$Cell, ignore.case = TRUE), "Med", NA)))
df1$phenotype <- factor(df1$phenotype, levels = c("Low", "Med", "High"))

ggplot(df1, aes(x = `phenotype`, y = `chr8_copy_number`)) +
  geom_boxplot() +
  labs(x = "Grouping Variable", y = "Value", title = "Boxplot by Group") +
  theme_minimal()



df_SNU16_PVT1_MYC <- read.delim("/home/yue1118/Bioskryb_data_analysis/PVT1_MYC_read_counts_across_single_cells_SNU16_Bioskryb.tsv")
View(df_SNU16_PVT1_MYC)
df_SNU16_PVT1_MYC_wider <- pivot_wider(df_SNU16_PVT1_MYC,  names_from = Gene, values_from = ReadCount)
#normalized by gene length: PVT1 length: 393698; MYC length: 7517
df_SNU16_PVT1_MYC_wider$normalized_PVT1 <- df_SNU16_PVT1_MYC_wider$PVT1/393698
df_SNU16_PVT1_MYC_wider$normalized_MYC <- df_SNU16_PVT1_MYC_wider$MYC/7517
View(df_SNU16_PVT1_MYC_wider)
plot(df_SNU16_PVT1_MYC_wider$normalized_PVT1, df_SNU16_PVT1_MYC_wider$normalized_MYC)
df_SNU16_PVT1_MYC_wider_subset <- df_SNU16_PVT1_MYC_wider[df_SNU16_PVT1_MYC_wider$normalized_PVT1 <0.002,]
plot(df_SNU16_PVT1_MYC_wider_subset$normalized_PVT1, df_SNU16_PVT1_MYC_wider_subset$normalized_MYC)
cor.test(df_SNU16_PVT1_MYC_wider_subset$normalized_PVT1, df_SNU16_PVT1_MYC_wider_subset$normalized_MYC) #cor=0.3694701; p-value = 0.0002117

svg('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/PVT1_vs_MYC_DNA_coverage_correlation_Bioskryb_SNU16_10172025.svg',width = 6, height = 6)
plot(df_SNU16_PVT1_MYC_wider$normalized_PVT1, df_SNU16_PVT1_MYC_wider$normalized_MYC)
dev.off()
cor.test(df_SNU16_PVT1_MYC_wider$normalized_PVT1, df_SNU16_PVT1_MYC_wider$normalized_MYC)


df_NCIH2170_PVT1_MYC <- read.delim("/home/yue1118/Bioskryb_data_analysis/PVT1_MYC_read_counts_across_single_cells_NCIH2170_Bioskryb.tsv")
View(df_NCIH2170_PVT1_MYC)
df_NCIH2170_PVT1_MYC_wider <- pivot_wider(df_NCIH2170_PVT1_MYC,  names_from = Gene, values_from = ReadCount)
#normalized by gene length: PVT1 length: 393698; MYC length: 7517
df_NCIH2170_PVT1_MYC_wider$normalized_PVT1 <- df_NCIH2170_PVT1_MYC_wider$PVT1/393698
df_NCIH2170_PVT1_MYC_wider$normalized_MYC <- df_NCIH2170_PVT1_MYC_wider$MYC/7517

svg('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/PVT1_vs_MYC_DNA_coverage_correlation_Bioskryb_NCIH2170_08082025.svg',width = 6, height = 6)
plot(df_NCIH2170_PVT1_MYC_wider$normalized_PVT1, df_NCIH2170_PVT1_MYC_wider$normalized_MYC)
dev.off()

cor.test(df_NCIH2170_PVT1_MYC_wider$normalized_PVT1, df_NCIH2170_PVT1_MYC_wider$normalized_MYC) #cor=0.9664445; p-value < 2.2e-16



df_COLO320_PVT1_MYC <- read.delim("/home/yue1118/Bioskryb_data_analysis/PVT1_MYC_read_counts_across_single_cells_COLO320DM_Bioskryb.tsv")
df_COLO320_PVT1_MYC_wider <- pivot_wider(df_COLO320_PVT1_MYC,  names_from = Gene, values_from = ReadCount)
#normalized by gene length: PVT1 length: 393698; MYC length: 7517
df_COLO320_PVT1_MYC_wider$normalized_PVT1 <- df_COLO320_PVT1_MYC_wider$PVT1/393698
df_COLO320_PVT1_MYC_wider$normalized_MYC <- df_COLO320_PVT1_MYC_wider$MYC/7517

svg('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/PVT1_vs_MYC_DNA_coverage_correlation_Bioskryb_COLO320DM_08082025.svg',width = 6, height = 6)
plot(df_COLO320_PVT1_MYC_wider$normalized_PVT1, df_COLO320_PVT1_MYC_wider$normalized_MYC)
dev.off()

cor.test(df_COLO320_PVT1_MYC_wider$normalized_PVT1, df_COLO320_PVT1_MYC_wider$normalized_MYC) #cor=0.988123; p-value < 2.2e-16


df_SNU16_PVT1_MYC_FGFR2 <- read.delim("/home/yue1118/Bioskryb_data_analysis/PVT1_MYC_FGFR2_PDHX_read_counts_across_single_cells_SNU16_Bioskryb.tsv")

df_SNU16_PVT1_MYC_FGFR2_wider <- pivot_wider(df_SNU16_PVT1_MYC_FGFR2,  names_from = Gene, values_from = ReadCount)
#normalized by gene length: PVT1 length: 393698; MYC length: 7517; FGFR2 length: 120128
df_SNU16_PVT1_MYC_FGFR2_wider$normalized_PVT1 <- df_SNU16_PVT1_MYC_FGFR2_wider$PVT1/393698
df_SNU16_PVT1_MYC_FGFR2_wider$normalized_MYC <- df_SNU16_PVT1_MYC_FGFR2_wider$MYC/7517
df_SNU16_PVT1_MYC_FGFR2_wider$normalized_FGFR2 <- df_SNU16_PVT1_MYC_FGFR2_wider$FGFR2/120128
df_SNU16_PVT1_MYC_FGFR2_wider$normalized_PDHX <- df_SNU16_PVT1_MYC_FGFR2_wider$PDHX/104762
df_SNU16_PVT1_MYC_FGFR2_wider_subset <- df_SNU16_PVT1_MYC_FGFR2_wider[df_SNU16_PVT1_MYC_FGFR2_wider$normalized_PVT1 <0.002,]
df_SNU16_PVT1_MYC_FGFR2_wider_subset_2 <- df_SNU16_PVT1_MYC_FGFR2_wider[df_SNU16_PVT1_MYC_FGFR2_wider$normalized_PVT1 > 0.002,]
plot(df_SNU16_PVT1_MYC_FGFR2_wider$normalized_PVT1, df_SNU16_PVT1_MYC_FGFR2_wider$normalized_PDHX)
plot(df_SNU16_PVT1_MYC_FGFR2_wider$normalized_PVT1, df_SNU16_PVT1_MYC_FGFR2_wider$normalized_MYC)
plot(df_SNU16_PVT1_MYC_FGFR2_wider$normalized_PVT1, df_SNU16_PVT1_MYC_FGFR2_wider$normalized_MYC)
plot(df_SNU16_PVT1_MYC_FGFR2_wider$normalized_PVT1, df_SNU16_PVT1_MYC_FGFR2_wider$normalized_FGFR2)
plot(df_SNU16_PVT1_MYC_FGFR2_wider$normalized_MYC, df_SNU16_PVT1_MYC_FGFR2_wider$normalized_FGFR2)
plot(df_SNU16_PVT1_MYC_FGFR2_wider$normalized_MYC, df_SNU16_PVT1_MYC_FGFR2_wider$normalized_PDHX)
plot(df_SNU16_PVT1_MYC_FGFR2_wider$normalized_MYC, df_SNU16_PVT1_MYC_FGFR2_wider$normalized_FGFR2)
plot(df_SNU16_PVT1_MYC_FGFR2_wider_subset$normalized_PVT1, df_SNU16_PVT1_MYC_FGFR2_wider_subset$normalized_FGFR2)
cor.test(df_SNU16_PVT1_MYC_FGFR2_wider_subset$normalized_PVT1, df_SNU16_PVT1_MYC_FGFR2_wider_subset$normalized_FGFR2)

hist(df_SNU16_PVT1_MYC_FGFR2_wider$normalized_PDHX)

svg('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/PVT1_vs_MYC_DNA_coverage_correlation_Bioskryb_SNU16_08082025.svg',width = 6, height = 6)
plot(df_SNU16_PVT1_MYC_FGFR2_wider_subset$normalized_PVT1, df_SNU16_PVT1_MYC_FGFR2_wider_subset$normalized_MYC)
dev.off()

cor.test(df_SNU16_PVT1_MYC_FGFR2_wider_subset$normalized_PVT1, df_SNU16_PVT1_MYC_FGFR2_wider_subset$normalized_MYC) #cor=0.38; p-value=0.0002117

svg('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/PVT1_DNA_reads_histogram_Bioskryb_SNU16_08082025.svg',width = 6, height = 6)
hist(df_SNU16_PVT1_MYC_FGFR2_wider_subset$normalized_PVT1)
dev.off()

df_SNU16_PVT1_MYC_FGFR2_wider_subset <- df_SNU16_PVT1_MYC_FGFR2_wider_subset[order(-df_SNU16_PVT1_MYC_FGFR2_wider_subset$normalized_PVT1), ]
svg('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/PVT1_DNA_reads_barplot_Bioskryb_SNU16_08082025.svg',width = 6, height = 6)
barplot(df_SNU16_PVT1_MYC_FGFR2_wider_subset$normalized_PVT1)
dev.off()

df_SNU16_PVT1_MYC_FGFR2_wider_subset <- df_SNU16_PVT1_MYC_FGFR2_wider_subset[order(-df_SNU16_PVT1_MYC_FGFR2_wider_subset$normalized_MYC), ]
barplot(df_SNU16_PVT1_MYC_FGFR2_wider_subset$normalized_MYC)

df_SNU16_PVT1_MYC_FGFR2_wider_subset <- df_SNU16_PVT1_MYC_FGFR2_wider_subset[order(-df_SNU16_PVT1_MYC_FGFR2_wider_subset$normalized_PDHX), ]
barplot(df_SNU16_PVT1_MYC_FGFR2_wider_subset$normalized_PDHX)

df_SNU16_PVT1_MYC_FGFR2_wider_subset <- df_SNU16_PVT1_MYC_FGFR2_wider_subset[order(-df_SNU16_PVT1_MYC_FGFR2_wider_subset$normalized_FGFR2), ]
barplot(df_SNU16_PVT1_MYC_FGFR2_wider_subset$normalized_FGFR2)


dim(df_SNU16_PVT1_MYC_FGFR2_wider_subset)

df_cn <- df_SNU16_PVT1_MYC_FGFR2_wider_subset[, c("normalized_PVT1", "normalized_MYC", "normalized_FGFR2")]
df_scaled <- scale(df_cn)

# Compute distance matrix (Euclidean)
dist_mat <- dist(df_scaled)

# Perform hierarchical clustering
hc <- hclust(dist_mat, method = "ward.D2")
clusters <- cutree(hc, k = 3)
df_cn$cluster <- factor(clusters)

# 3. Plot dendrogram with colored labels
library(dendextend)
dend <- as.dendrogram(hc)
dend <- color_branches(dend, k = 3)


pheatmap(df_scaled,
         cluster_rows = TRUE,   # cluster genes
         cluster_cols = TRUE,   # cluster cells
         fontsize_row = 10,
         fontsize_col = 10,
         main = "ecDNA Copy Number Heatmap",
         color = colorRampPalette(c("blue", "white", "red"))(100))

# Plot dendrogram
plot(hc, main = "Hierarchical Clustering of ecDNA CN", xlab = "Cells")

pheatmap(df_cn[,1:3],
         cluster_rows = TRUE,   # cluster genes
         cluster_cols = TRUE,   # cluster cells
         scale = "none",        # <<< THIS ensures data is not scaled
         fontsize_row = 10,
         fontsize_col = 10,
         main = "Unscaled CN Heatmap",
         color = colorRampPalette(c("white", "orange", "red"))(100))



df_SNU16_vcf <- read.table("/home/yue1118/Bioskryb_data_analysis/SNU16_ecDNA_sniffle_output.vcf", comment.char = "#", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(df_SNU16_vcf)[1:8] <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
df_SNU16_vcf_bnd <- df_SNU16_vcf[grepl("SVTYPE=BND", df_SNU16_vcf$INFO), ]
df_SNU16_vcf_tra <- df_SNU16_vcf[grepl("SVTYPE=TRA", df_SNU16_vcf$INFO), ]
df_SNU16_vcf_del <- df_SNU16_vcf[grepl("SVTYPE=DEL", df_SNU16_vcf$INFO), ]
df_SNU16_vcf_inv <- df_SNU16_vcf[grepl("SVTYPE=INV", df_SNU16_vcf$INFO), ]

View(df_SNU16_vcf_bnd)
View(df_SNU16_vcf_tra)
View(df_SNU16_vcf_inv)


df_SNU16_wgs_vcf <- read.table("/home/yue1118/Bioskryb_data_analysis/SNU16_WGS_sniffle_output.vcf", comment.char = "#", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(df_SNU16_wgs_vcf)[1:8] <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
df_SNU16_vcf_wgs_bnd <- df_SNU16_wgs_vcf[grepl("SVTYPE=BND", df_SNU16_wgs_vcf$INFO), ]
df_SNU16_vcf_wgs_tra <- df_SNU16_wgs_vcf[grepl("SVTYPE=TRA", df_SNU16_wgs_vcf$INFO), ]
df_SNU16_vcf_wgs_inv <- df_SNU16_wgs_vcf[grepl("SVTYPE=INV", df_SNU16_wgs_vcf$INFO), ]
View(df_SNU16_vcf_wgs_bnd)
View(df_SNU16_vcf_wgs_inv)



df_320DM_vcf <- read.table("/home/yue1118/Bioskryb_data_analysis/COLO320DM_sniffle_output.vcf", comment.char = "#", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(df_320DM_vcf)[1:8] <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
df_320DM_vcf_bnd <- df_320DM_vcf[grepl("SVTYPE=BND", df_320DM_vcf$INFO), ]
df_320DM_vcf_tra <- df_320DM_vcf[grepl("SVTYPE=TRA", df_320DM_vcf$INFO), ]
df_320DM_vcf_inv <- df_320DM_vcf[grepl("SVTYPE=INV", df_320DM_vcf$INFO), ]
View(df_320DM_vcf_bnd)
View(df_320DM_vcf_tra)
View(df_320DM_vcf_inv)

df_2170_vcf <- read.table("/home/yue1118/Bioskryb_data_analysis/NCIH2170_sniffle_output.vcf", comment.char = "#", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(df_2170_vcf)[1:8] <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
df_2170_vcf_bnd <- df_2170_vcf[grepl("SVTYPE=BND", df_2170_vcf$INFO), ]
df_2170_vcf_tra <- df_2170_vcf[grepl("SVTYPE=TRA", df_2170_vcf$INFO), ]
df_2170_vcf_inv <- df_2170_vcf[grepl("SVTYPE=INV", df_2170_vcf$INFO), ]
View(df_2170_vcf_bnd)
View(df_2170_vcf_tra)
View(df_2170_vcf_inv)

View(df_2170_vcf)


df_SNU16_PDHX_SA <- read.delim("/home/yue1118/Bioskryb_data_analysis/SNU16_ONT_LRS_chr11_34.88M_35.34M_SA_tags.tsv", header = FALSE)
df_SNU16_PDHX_SA <- df_SNU16_PDHX_SA %>% separate(V2, into = c("chr", "pos", "strand", "CIGAR", "V6", "V7"), sep = ",")
df_SNU16_PDHX_SA_chr8 <- df_SNU16_PDHX_SA[df_SNU16_PDHX_SA$chr == "SA:Z:chr8",]
df_SNU16_PDHX_SA_chr8$pos <- as.numeric(df_SNU16_PDHX_SA_chr8$pos)
df_SNU16_PDHX_SA_chr8$matched_length <- str_extract(df_SNU16_PDHX_SA_chr8$CIGAR, "(?<=S).*?(?=M)")
df_SNU16_PDHX_SA_chr8$matched_length <- as.numeric(df_SNU16_PDHX_SA_chr8$matched_length)
df_SNU16_PDHX_SA_chr8 <- df_SNU16_PDHX_SA_chr8[(df_SNU16_PDHX_SA_chr8$pos > 127000000) & (df_SNU16_PDHX_SA_chr8$pos < 130000000),]
df_SNU16_PDHX_SA_chr8$strand <- as.character(df_SNU16_PDHX_SA_chr8$strand)
df_SNU16_PDHX_SA_chr8 <- df_SNU16_PDHX_SA_chr8 %>% mutate(new_pos = ifelse(strand == "+", pos + matched_length, pos - matched_length))
View(df_SNU16_PDHX_SA_chr8)

#To filter for the reads in exon1 and exon2
df_SNU16_PDHX_SA_chr8_exon1_2 <- df_SNU16_PDHX_SA_chr8[(df_SNU16_PDHX_SA_chr8$pos > 127794524) & (df_SNU16_PDHX_SA_chr8$pos < 127890998),]
View(df_SNU16_PDHX_SA_chr8[(df_SNU16_PDHX_SA_chr8$new_pos > 127794524) & (df_SNU16_PDHX_SA_chr8$new_pos < 127890998),])


table(df_SNU16_PDHX_SA_chr8$strand)

df_SNU16_PVT1_SA <- read.delim("/home/yue1118/Bioskryb_data_analysis/SNU16_ONT_LRS_PVT1_127.79M_128.05M.bam_SA_tags.tsv", header = FALSE)
df_SNU16_PVT1_SA <- df_SNU16_PVT1_SA %>% separate(V2, into = c("chr", "pos", "strand", "CIGAR", "V6", "V7"), sep = ",")
table(df_SNU16_PVT1_SA$chr)

df_SNU16_PVT1_SA_chr11 <- df_SNU16_PVT1_SA[df_SNU16_PVT1_SA$chr == "SA:Z:chr11",]
df_SNU16_PVT1_SA_chr11$pos <- as.numeric(df_SNU16_PVT1_SA_chr11$pos)
df_SNU16_PVT1_SA_chr11_ecDNA <- df_SNU16_PVT1_SA_chr11[(df_SNU16_PVT1_SA_chr11$pos > 34880000) & (df_SNU16_PVT1_SA_chr11$pos < 35340000) ,]
View(df_SNU16_PVT1_SA_chr11_ecDNA)
df_SNU16_PVT1_SA_chr8 <- df_SNU16_PVT1_SA[df_SNU16_PVT1_SA$chr == "SA:Z:chr8",]





######PVT1 exon1&exon2 SA tags in SNU16 ##########
df_SNU16_PVT1_exon12_SA <- read.delim("/home/yue1118/Bioskryb_data_analysis/SNU16_PVT1_exon1_exon2_SA_tags.tsv", header = FALSE)
df_SNU16_PVT1_exon12_SA <- df_SNU16_PVT1_exon12_SA %>% separate(V2, into = c("chr", "pos", "strand", "CIGAR", "V6", "V7"), sep = ",")
table(df_SNU16_PVT1_exon12_SA$chr)

df_SNU16_PVT1_exon12_SA_chr8 <- df_SNU16_PVT1_exon12_SA[df_SNU16_PVT1_exon12_SA$chr == "SA:Z:chr8",]
df_SNU16_PVT1_exon12_SA_chr8$pos <- as.numeric(df_SNU16_PVT1_exon12_SA_chr8$pos)
df_SNU16_PVT1_exon12_SA_chr8 <- df_SNU16_PVT1_exon12_SA_chr8[(df_SNU16_PVT1_exon12_SA_chr8$pos > 127000000) & (df_SNU16_PVT1_exon12_SA_chr8$pos < 130000000),]
hist(df_SNU16_PVT1_exon12_SA_chr8$pos, breaks = 300)
View(df_SNU16_PVT1_exon12_SA_chr8)
table(df_SNU16_PVT1_exon12_SA$chr)

######PVT1 exon1&exon2 SA tags in NCIH2170 ##########
df_2170_PVT1_exon12_SA <- read.delim("/home/yue1118/Bioskryb_data_analysis/NCIH2170_PVT1_exon1_exon2_SA_tags.tsv", header = FALSE)
df_2170_PVT1_exon12_SA <- df_2170_PVT1_exon12_SA %>% separate(V2, into = c("chr", "pos", "strand", "CIGAR", "V6", "V7"), sep = ",")
table(df_2170_PVT1_exon12_SA$chr)

df_2170_PVT1_exon12_SA_chr8 <- df_2170_PVT1_exon12_SA[df_2170_PVT1_exon12_SA$chr == "SA:Z:chr8",]
df_2170_PVT1_exon12_SA_chr8$pos <- as.numeric(df_2170_PVT1_exon12_SA_chr8$pos)
df_2170_PVT1_exon12_SA_chr8 <- df_2170_PVT1_exon12_SA_chr8[(df_2170_PVT1_exon12_SA_chr8$pos > 127000000) & (df_2170_PVT1_exon12_SA_chr8$pos < 130000000),]
hist(df_2170_PVT1_exon12_SA_chr8$pos, breaks = 300)
View(df_2170_PVT1_exon12_SA_chr8)

######PVT1 exon1&exon2 SA tags in COLO320DM ##########
df_320_PVT1_exon12_SA <- read.delim("/home/yue1118/Bioskryb_data_analysis/COLO320DM_PVT1_exon1_exon2_SA_tags.tsv", header = FALSE)
df_320_PVT1_exon12_SA <- df_320_PVT1_exon12_SA %>% separate(V2, into = c("chr", "pos", "strand", "CIGAR", "V6", "V7"), sep = ",")
table(df_320_PVT1_exon12_SA$chr) #only 95 reads with SA tags, indicating less rearrangement event happening within PVT1.

df_320_PVT1_exon12_SA_chr8 <- df_320_PVT1_exon12_SA[df_320_PVT1_exon12_SA$chr == "SA:Z:chr8",]
df_320_PVT1_exon12_SA_chr8$pos <- as.numeric(df_320_PVT1_exon12_SA_chr8$pos)
df_2170_PVT1_exon12_SA_chr8 <- df_320_PVT1_exon12_SA_chr8[(df_320_PVT1_exon12_SA_chr8$pos > 127000000) & (df_320_PVT1_exon12_SA_chr8$pos < 130000000),]
hist(df_320_PVT1_exon12_SA_chr8$pos, breaks = 300)
View(df_320_PVT1_exon12_SA_chr8)


df_SNU16_H3K27ac <- read.table('/home/yue1118/Bioskryb_data_analysis/GSE159972_SNU16_H3K27ac_peaks.narrowPeak', sep = '\t')
df_SNU16_H3K27ac_chr8 <- df_SNU16_H3K27ac[df_SNU16_H3K27ac$V1 == "chr8",]
View(df_SNU16_H3K27ac_chr8)


##############PDHX & APIP SA tags in SNU16#############
df_SNU16_PDHX_APIP_SA <- read.delim("/home/yue1118/Bioskryb_data_analysis/SNU16_PDHX_APIP_SA_tags.tsv", header = FALSE)
df_SNU16_PDHX_APIP_SA <- df_SNU16_PDHX_APIP_SA %>% separate(V2, into = c("chr", "pos", "strand", "CIGAR", "V6", "V7"), sep = ",")
df_SNU16_PDHX_APIP_SA$pos <- as.numeric(df_SNU16_PDHX_APIP_SA$pos)
table(df_SNU16_PDHX_APIP_SA$chr) 

df_SNU16_PDHX_APIP_SA_chr10 <- df_SNU16_PDHX_APIP_SA[df_SNU16_PDHX_APIP_SA$chr=="SA:Z:chr10",]
df_SNU16_PDHX_APIP_SA_chr1 <- df_SNU16_PDHX_APIP_SA[df_SNU16_PDHX_APIP_SA$chr=="SA:Z:chr1",]
df_SNU16_PDHX_APIP_SA_chr8 <- df_SNU16_PDHX_APIP_SA[df_SNU16_PDHX_APIP_SA$chr=="SA:Z:chr8",]

View(df_SNU16_PDHX_APIP_SA[df_SNU16_PDHX_APIP_SA$chr=="SA:Z:chr1",])

View(df_SNU16_PDHX_APIP_SA[df_SNU16_PDHX_APIP_SA$chr=="SA:Z:chr10",])

View(df_SNU16_PDHX_APIP_SA[df_SNU16_PDHX_APIP_SA$chr=="SA:Z:chr8",])


hist(df_SNU16_PDHX_APIP_SA_chr10$pos, breaks = 200)
hist(df_SNU16_PDHX_APIP_SA_chr1$pos, breaks = 200)
hist(df_SNU16_PDHX_APIP_SA_chr8$pos, breaks = 200)


df_SNU16_PDHX_APIP_SA_chr8_subset <- df_SNU16_PDHX_APIP_SA_chr8[(df_SNU16_PDHX_APIP_SA_chr8$pos >126000000) & (df_SNU16_PDHX_APIP_SA_chr8$pos <130000000),]

hist(df_SNU16_PDHX_APIP_SA_chr8_subset$pos, breaks = 200)

df_COLO320DM_H3K27ac <- read.table('/home/yue1118/Bioskryb_data_analysis/GSE159972_COLO320DM_DMSO_H3K27ac_peaks.narrowPeak', sep = '\t')
df_COLO320DM_H3K27ac_chr8 <- df_COLO320DM_H3K27ac[df_COLO320DM_H3K27ac$V1 == "chr8",]
View(df_COLO320DM_H3K27ac_chr8)


df_SNU16_BRD4 <- read.table('/home/yue1118/Bioskryb_data_analysis/GSE159972_SNU16_Brd4_peaks.narrowPeak', sep = '\t')
df_SNU16_BRD4_chr8 <- df_SNU16_BRD4[df_SNU16_BRD4$V1 == "chr8",]
View(df_SNU16_BRD4_chr8)


df_COLO320DM_BRD4 <- read.table('/home/yue1118/Bioskryb_data_analysis/GSE159972_COLO320DM_DMSO_Brd4_peaks.narrowPeak', sep = '\t')
df_COLO320DM_BRD4_chr8 <- df_COLO320DM_BRD4[df_COLO320DM_BRD4$V1 == "chr8",]
View(df_COLO320DM_BRD4_chr8)



####################################################################################################
df_SNU16_PVT1_read <- read.delim("/datastore/lbcfs/labs/brunk_lab/private/Bioskryb_SNU16_DNA/SNU16_PVT1_region_read_counts.tsv", header = TRUE)
df_SNU16_PVT1_read$sample <- sub(".*(SC[0-9]+-SNU16(?:Low|High|Med)).*", "\\1", df_SNU16_PVT1_read$sample)
df_SNU16_PVT1_read_wide <- pivot_wider(df_SNU16_PVT1_read,  names_from = region, values_from = read_count)
colnames(df_SNU16_PVT1_read_wide)[1] <- "Cell"
View(df_SNU16_PVT1_read_wide)

ecDNA_samples <- c("SC81-SNU16High", "SC55-SNU16Med", "SC38-SNU16Low", "SC53-SNU16Med", "SC10-SNU16Low", "SC44-SNU16Low", "SC43-SNU16Low", "SC61-SNU16Med", "SC17-SNU16Low", "SC45-SNU16Low", "SC36-SNU16Low", "SC3-SNU16Low", "SC48-SNU16Low", "SC30-SNU16Low", "SC24-SNU16Low")

df_SNU16_PVT1_read_wide <- df_SNU16_PVT1_read_wide %>% mutate(group = ifelse(Cell %in% ecDNA_samples, "ecDNA", "HSR"))

View(df_SNU16_PVT1_read_wide)

df_DNA_read_long <- df_SNU16_PVT1_read_wide %>%
  pivot_longer(cols = 5:21, 
               names_to = "region", 
               values_to = "count")

ggplot(df_DNA_read_long, aes(x = group, y = count, fill = group)) +
  geom_boxplot() +
  facet_wrap(~ region, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#Now normalize the read numbers by genomic region length:
# Suppose your dataframe is called df
colnames_df <- colnames(df_SNU16_PVT1_read_wide)

# Function to get length from "chr:start-end"
get_length <- function(x) {
  parts <- unlist(strsplit(x, "[:-]"))
  if (length(parts) == 3) {
    start <- as.numeric(parts[2])
    end <- as.numeric(parts[3])
    return(end - start)
  } else {
    return(NA)  # for columns like "sample"
  }
}

# Compute lengths for each column
interval_lengths <- sapply(colnames_df, get_length)

# Now divide each column (except "sample") by its length
df_norm <- df_SNU16_PVT1_read_wide
for (i in 2:ncol(df_SNU16_PVT1_read_wide)) {  # start at 2 to skip "sample"
  df_norm[[i]] <- df_SNU16_PVT1_read_wide[[i]] / interval_lengths[i]
}
rownames(df_norm) <- df_norm$sample
df_norm[,2:20] <- df_norm[,2:20]*10000

View(df_norm)


####################################################################################################
df_SNU16_PVT1_RNA_read <- read.delim("/datastore/lbcfs/labs/brunk_lab/private/Bioskryb_SNU16_RNA/STAR_alignment_for_single_cells/SNU16_PVT1_region_RNA_read_counts.tsv", header = TRUE)
df_SNU16_PVT1_RNA_read$sample <- sub(".*(SC[0-9]+-SNU16(?:Low|High|Med)).*", "\\1", df_SNU16_PVT1_RNA_read$sample)
df_SNU16_PVT1_RNA_read_wide <- pivot_wider(df_SNU16_PVT1_RNA_read,  names_from = region, values_from = read_count)

View(df_SNU16_PVT1_RNA_read_wide)

ecDNA_samples <- c("SC81-SNU16High", "SC55-SNU16Med", "SC38-SNU16Low", "SC53-SNU16Med", "SC10-SNU16Low", "SC44-SNU16Low", "SC43-SNU16Low", "SC61-SNU16Med", "SC17-SNU16Low", "SC45-SNU16Low", "SC36-SNU16Low", "SC3-SNU16Low", "SC48-SNU16Low", "SC30-SNU16Low", "SC24-SNU16Low")

df_SNU16_PVT1_RNA_read_wide <- df_SNU16_PVT1_RNA_read_wide %>% mutate(group = ifelse(sample %in% ecDNA_samples, "ecDNA", "HSR"))

View(df_SNU16_PVT1_RNA_read_wide)

df_RNA_read_long <- df_SNU16_PVT1_RNA_read_wide %>%
  pivot_longer(cols = 5:21, 
               names_to = "region", 
               values_to = "count")

ggplot(df_RNA_read_long, aes(x = group, y = count, fill = group)) +
  geom_boxplot() +
  facet_wrap(~ region, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Here to import normalized PVT1 RNA reads for all PVT1 regions for replicate1
df_SNU16_PVT1_RNA_read_norm <- read.csv("/home/yue1118/Bioskryb_data_analysis/Bioskryb_SNU16_exp_TPM_ecDNA_FGFR2_and_MYC_CN_PVT1_splice_junction_reads_ssGSEA_score_PVT1_all_regions_RNA_reads_09062025.csv")
View(df_SNU16_PVT1_RNA_read_norm)


df_SNU16_PVT1_RNA_read_norm_2 <- df_SNU16_PVT1_RNA_read_norm[, c("Cell","exon1_RNA_normalized_read_count", "intron_between_exon1_exon2_RNA_normalized_read_count","exon2_RNA_normalized_read_count",
                                                               "intron_between_exon2_exon3_RNA_normalized_read_count", "exon3_RNA_normalized_read_count", "intron_between_exon3_exon4_RNA_normalized_read_count",
                                                               "exon4_RNA_normalized_read_count","intron_between_exon4_exon5_RNA_normalized_read_count", "exon5_RNA_normalized_read_count",
                                                               "intron_between_exon5_exon6_RNA_normalized_read_count", "exon6_RNA_normalized_read_count", "intron_between_exon6_exon7_RNA_normalized_read_count",
                                                               "exon7_RNA_normalized_read_count", "intron_between_exon7_exon8_RNA_normalized_read_count", "exon8_RNA_normalized_read_count", 
                                                               "intron_between_exon8_exon9_RNA_normalized_read_count", "exon9_RNA_normalized_read_count")]
View(df_SNU16_PVT1_RNA_read_norm_2)

df_SNU16_PVT1_RNA_read_norm_2$group <-ifelse(df_SNU16_PVT1_RNA_read_norm_2$Cell %in% ecDNA_samples, "ecDNA", "HSR")

df_long <- df_SNU16_PVT1_RNA_read_norm_2 %>% pivot_longer(cols = c("exon1_RNA_normalized_read_count", "intron_between_exon1_exon2_RNA_normalized_read_count","exon2_RNA_normalized_read_count",
                                        "intron_between_exon2_exon3_RNA_normalized_read_count", "exon3_RNA_normalized_read_count", "intron_between_exon3_exon4_RNA_normalized_read_count",
                                        "exon4_RNA_normalized_read_count","intron_between_exon4_exon5_RNA_normalized_read_count", "exon5_RNA_normalized_read_count",
                                        "intron_between_exon5_exon6_RNA_normalized_read_count", "exon6_RNA_normalized_read_count", "intron_between_exon6_exon7_RNA_normalized_read_count",
                                        "exon7_RNA_normalized_read_count", "intron_between_exon7_exon8_RNA_normalized_read_count", "exon8_RNA_normalized_read_count", 
                                        "intron_between_exon8_exon9_RNA_normalized_read_count", "exon9_RNA_normalized_read_count"), names_to = "variable", values_to = "value")
ggplot(df_long, aes(x = variable, y = value, fill = group)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/Bioskryb_SNU16_PVT1_all_region_exp_between_PVT1_ecDNA_and_PVT1_HSR_groups_09132025.svg",width = 8, height = 6)


#Here to import normalized PVT1 RNA reads for all PVT1 regions for replicate2
df_SNU16_PVT1_RNA_read_norm_replicate2 <- read.csv("/home/yue1118/Bioskryb_data_analysis/Bioskryb_SNU16_plate2_exp_TPM_ecDNA_FGFR2_and_MYC_CN_PVT1_splice_junction_reads_ssGSEA_score_PVT1_all_regions_RNA_reads_09092025.csv")
View(df_SNU16_PVT1_RNA_read_norm_replicate2)


df_SNU16_PVT1_RNA_read_norm_replicate2_2 <- df_SNU16_PVT1_RNA_read_norm_replicate2[, c("Cell","exon1_RNA_normalized_read_count", "intron_between_exon1_exon2_RNA_normalized_read_count","exon2_RNA_normalized_read_count",
                                                                 "intron_between_exon2_exon3_RNA_normalized_read_count", "exon3_RNA_normalized_read_count", "intron_between_exon3_exon4_RNA_normalized_read_count",
                                                                 "exon4_RNA_normalized_read_count","intron_between_exon4_exon5_RNA_normalized_read_count", "exon5_RNA_normalized_read_count",
                                                                 "intron_between_exon5_exon6_RNA_normalized_read_count", "exon6_RNA_normalized_read_count", "intron_between_exon6_exon7_RNA_normalized_read_count",
                                                                 "exon7_RNA_normalized_read_count", "intron_between_exon7_exon8_RNA_normalized_read_count", "exon8_RNA_normalized_read_count", 
                                                                 "intron_between_exon8_exon9_RNA_normalized_read_count", "exon9_RNA_normalized_read_count", "PVT1_amp_format")]

View(df_SNU16_PVT1_RNA_read_norm_replicate2_2)

df_long_replicate2 <- df_SNU16_PVT1_RNA_read_norm_replicate2_2 %>% pivot_longer(cols = c("exon1_RNA_normalized_read_count", "intron_between_exon1_exon2_RNA_normalized_read_count","exon2_RNA_normalized_read_count",
                                                                   "intron_between_exon2_exon3_RNA_normalized_read_count", "exon3_RNA_normalized_read_count", "intron_between_exon3_exon4_RNA_normalized_read_count",
                                                                   "exon4_RNA_normalized_read_count","intron_between_exon4_exon5_RNA_normalized_read_count", "exon5_RNA_normalized_read_count",
                                                                   "intron_between_exon5_exon6_RNA_normalized_read_count", "exon6_RNA_normalized_read_count", "intron_between_exon6_exon7_RNA_normalized_read_count",
                                                                   "exon7_RNA_normalized_read_count", "intron_between_exon7_exon8_RNA_normalized_read_count", "exon8_RNA_normalized_read_count", 
                                                                   "intron_between_exon8_exon9_RNA_normalized_read_count", "exon9_RNA_normalized_read_count"), names_to = "variable", values_to = "value")
View(df_long_replicate2)
#To get rid of outliers:
df_long_replicate2_subset <- df_long_replicate2[df_long_replicate2$value <0.4,]
ggplot(df_long_replicate2_subset, aes(x = variable, y = value, fill = PVT1_amp_format)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/Bioskryb_SNU16_PVT1_all_region_exp_between_PVT1_ecDNA_and_PVT1_HSR_groups_replicate2_09132025.svg",width = 8, height = 6)


df_SNU16_PVT1_MYC_FGFR2 <- read.delim("/home/yue1118/Bioskryb_data_analysis/PVT1_MYC_FGFR2_PDHX_read_counts_across_single_cells_SNU16_Bioskryb.tsv")
df_SNU16_PVT1_MYC_FGFR2$Cell <- sub(".*(SC[0-9]+-SNU16(?:Low|High|Med)).*", "\\1", df_SNU16_PVT1_MYC_FGFR2$Cell)
df_SNU16_PVT1_MYC_FGFR2_wide <- pivot_wider(df_SNU16_PVT1_MYC_FGFR2,  names_from = Gene, values_from = ReadCount)
View(df_SNU16_PVT1_MYC_FGFR2_wide)

df_SNU16_PVT1_merged <- merge(df_SNU16_PVT1_read_wide, df_SNU16_PVT1_MYC_FGFR2_wide, by.x="Cell", by.y="Cell", all.x=TRUE)
View(df_SNU16_PVT1_merged)
df_SNU16_PVT1_merged_subset1 <- df_SNU16_PVT1_merged[df_SNU16_PVT1_merged$`chr8:128027822-128049689` > 100,]

plot(df_SNU16_PVT1_merged$MYC, df_SNU16_PVT1_merged$`chr8:128027822-128049689`)


plot(df_SNU16_PVT1_merged_subset1$MYC, df_SNU16_PVT1_merged_subset1$`chr8:128027822-128049689`)
cor.test(df_SNU16_PVT1_merged_subset1$`chr8:128027822-128049689`, df_SNU16_PVT1_merged_subset1$MYC)

plot(df_SNU16_PVT1_merged$PVT1, df_SNU16_PVT1_merged$MYC)
cor.test(df_SNU16_PVT1_merged$PVT1, df_SNU16_PVT1_merged$MYC)

cor.test(df_SNU16_PVT1_pieces_MYC$`chr8:128027822-128049689`, df_SNU16_PVT1_pieces_MYC$ReadCount)

plot(df_SNU16_PVT1_pieces_MYC$`chr8:128003630-128026000`, df_SNU16_PVT1_pieces_MYC$ReadCount)
plot(df_SNU16_PVT1_pieces_MYC$`chr8:127986321-128002822`, df_SNU16_PVT1_pieces_MYC$ReadCount)

df_SNU16_ssGSEA <- read.csv("/home/yue1118/Bioskryb_data_analysis/SNU16_ssGSEA_PVT1_splice_variant_df_07262025.csv")
df_SNU16_ssGSEA <- df_SNU16_ssGSEA[, 1:51]
df_SNU16_ssGSEA$X <- gsub("\\.", "-", df_SNU16_ssGSEA$X)
df_SNU16_ssGSEA$Cell <- sub(".*-R-(.*)_S.*", "\\1", df_SNU16_ssGSEA$X)
df_SNU16_ssGSEA$Cell <- sub("-$", "", df_SNU16_ssGSEA$Cell)
View(df_SNU16_ssGSEA)

df_SNU16_PVT1_merged <- merge(df_SNU16_PVT1_merged, df_SNU16_ssGSEA,by.x="Cell", by.y="Cell", all.x=TRUE)
df_SNU16_PVT1_merged_subset1 <- df_SNU16_PVT1_merged[df_SNU16_PVT1_merged$`chr8:128003630-128026000` > 130,]
View(df_SNU16_PVT1_merged_subset1)
plot(df_SNU16_PVT1_merged_subset1$PVT1, df_SNU16_PVT1_merged_subset1$HALLMARK_MYC_TARGETS_V1)
plot(df_SNU16_PVT1_merged_subset1$PVT1, df_SNU16_PVT1_merged_subset1$HALLMARK_MYC_TARGETS_V2)
plot(df_SNU16_PVT1_merged_subset1$PVT1, df_SNU16_PVT1_merged_subset1$HALLMARK_E2F_TARGETS)
plot(df_SNU16_PVT1_merged_subset1$PVT1, df_SNU16_PVT1_merged_subset1$HALLMARK_G2M_CHECKPOINT)
plot(df_SNU16_PVT1_merged_subset1$PVT1, df_SNU16_PVT1_merged_subset1$HALLMARK_MTORC1_SIGNALING)
cor.test(df_SNU16_PVT1_merged_subset1$PVT1, df_SNU16_PVT1_merged_subset1$HALLMARK_E2F_TARGETS)


df_SNU16_DNA_plate2 <- read.delim("/home/yue1118/Bioskryb_data_analysis/PVT1_MYC_FGFR2_read_counts_across_single_cells_SNU16_Bioskryb_plate2.tsv")
df_SNU16_DNA_plate2$Cell <- sub(".*(SC[0-9]+-SNU16(?:Low|High|Med)).*", "\\1", df_SNU16_DNA_plate2$Cell)
df_SNU16_DNA_plate2_wide <- pivot_wider(df_SNU16_DNA_plate2,  names_from = Gene, values_from = ReadCount)
View(df_SNU16_DNA_plate2_wide)

df_SNU16_DNA_plate2_wide_subset <- df_SNU16_DNA_plate2_wide[df_SNU16_DNA_plate2_wide$PVT1< 30000,]
plot(df_SNU16_DNA_plate2_wide_subset$PVT1, df_SNU16_DNA_plate2_wide_subset$MYC)
