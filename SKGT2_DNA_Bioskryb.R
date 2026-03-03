
library(ggplot2)
library(svglite)


rds1 <- readRDS("/home/yue1118/Bioskryb_data_analysis/SKGT2_COLO320_ginkgo_res.binsize_1000000.RDS")

#To extract the normalized coverage numbers
df_coverage <- rds1$SegNorm
rownames(df_coverage) <- paste(df_coverage[,1], df_coverage[,2], df_coverage[,3], sep = "-")
df_coverage <- df_coverage[,4:ncol(df_coverage)]
df_coverage_SKGT2 <- df_coverage[, grepl("SKGT", names(df_coverage))]

View(df_coverage_SKGT2[grepl("chr17",rownames(df_coverage_SKGT2)),])
View(df_coverage_SKGT2[grepl("chr8",rownames(df_coverage_SKGT2)),])


coverage_quantiles <- apply(df_coverage_SKGT2, 2, quantile, probs = 1.0, na.rm = TRUE)
# Function to filter values below the 90th percentile and compute median
filtered_coverage_means <- sapply(1:ncol(df_coverage_SKGT2), function(i) {
  col_values <- df_coverage_SKGT2[, i]  # Extract column values
  filtered_values <- col_values[col_values < coverage_quantiles[i]]  # Keep values below 90th percentile
  mean(filtered_values, na.rm = TRUE)  #We use the mean value as coverage in "normal" chromosomal regions #I think we can use this way to identity the coverage for regions without CNV on chromosome
})                                      #in that for DepMap NCIH2170 CN profile (Copy number log2 (relative to ploidy +1)), it centers at 1, which means that most of genes don't have CNV events.
#This further means that the mean coverage values in the chromosomal regions are the coverage for "normal" regions without CNV events.

# Convert result to a dataframe
df1 <- data.frame(Cell = names(df_coverage_SKGT2), "mean_coverage" = filtered_coverage_means)
View(df1)


df1$chr17_coverage <- as.vector(t(df_coverage_SKGT2[c("chr17-39686887-40830641"),]))
df1$chr8_coverage <- as.vector(t(df_coverage_SKGT2[c("chr8-126598157-127691734"),])) #PCAT1, PRNCR1, CASC8, CCAT2, POU5F1B, CASC11, MYC are amplified in SKTG2, therefore this bin is selected.
df1$ERBB2_CN <- 2.5*df1$chr17_coverage/df1$mean_coverage

write.csv(df1, "/home/yue1118/Bioskryb_data_analysis/SKGT2_Bioskryb_ecDNA_counts_05272025.csv")

View(df1)
rownames_vector <- df1$Cell
df1$phenotype <- ifelse(grepl("-L-", rownames_vector, ignore.case = TRUE), "Low",
                        ifelse(grepl("-H-", rownames_vector, ignore.case = TRUE), "High",
                               ifelse(grepl("-M-", rownames_vector, ignore.case = TRUE), "Med", NA)))
df1$phenotype <- factor(df1$phenotype, levels = c("Low", "Med", "High"))
df1 <- df1[order(-df1$ERBB2_CN), ]


#To plot ERBB2 copy number distribution across single cells
svg('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/SKGT2_ERBB2_copy_number_by_scDNA_seq_05272025.svg',width = 6, height = 6)
barplot(df1$ERBB2_CN)
dev.off()

# Create a boxplot
ggplot(df1, aes(x = `phenotype`, y = `ERBB2_CN`)) +
  geom_boxplot() +
  labs(x = "Grouping Variable", y = "Value", title = "Boxplot by Group") +
  theme_minimal()
ggsave('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/SKGT2_ERBB2_ecDNA_counts_boxplot_high_med_low_05272025.svg', width = 6, height = 6)



wilcox.test(df1[df1$phenotype =="High","chr17_coverage"], df1[df1$phenotype =="Low","chr17_coverage"], paired = F) #p-value = 0.1364 
wilcox.test(df1[df1$phenotype =="High","chr17_coverage"], df1[df1$phenotype =="Med","chr17_coverage"], paired = F) #p-value = 0.7289
wilcox.test(df1[df1$phenotype =="Med","chr17_coverage"], df1[df1$phenotype =="Low","chr17_coverage"], paired = F) #p-value = 0.1717




df_cn_skgt2 <- read.csv("/home/yue1118/Bioskryb_data_analysis/SKGT2_Copy_Number_log2_relative_to_ploidy_1.csv")
rownames(df_cn_skgt2) <- df_cn_skgt2$gene
colnames(df_cn_skgt2) <- c("gene", "CN")
df_cn_skgt2$CN <- as.numeric(df_cn_skgt2$CN)
View(df_cn_skgt2)
mean(df_cn_skgt2$CN, na.rm = T) #Mean copy number is 1.024; karyotype is 2.52 (SKGT2 modal number is 58) according to https://www.sigmaaldrich.com/US/en/product/sigma/11012008?srsltid=AfmBOooeUaIF_9bKeznsUxDvgyAhNG7ldiXcX6ZJyZZtf64zzptMfdBB




#To compute the correlation between bin chr8-127412678-128444056 (where MYC is located) and bin chr17-39284904-40331026 (where ERBB2 is located)

rds1 <- readRDS("/home/yue1118/Bioskryb_data_analysis/SKGT2_COLO320_ginkgo_res.binsize_1000000.RDS")
df_coverage <- rds1$SegNorm
rownames(df_coverage) <- paste(df_coverage[,1], df_coverage[,2], df_coverage[,3], sep = "-")
df_coverage <- df_coverage[,4:ncol(df_coverage)]
df_coverage_SKGT <- df_coverage[, grepl("SKGT", names(df_coverage))]
df_coverage_SKGT_2 <- t(df_coverage_SKGT)
df_coverage_SKGT_2 <- as.data.frame(df_coverage_SKGT_2)
df_t <- df_coverage_SKGT_2[,c("chr17-39686887-40830641","chr8-126598157-127691734")]
colnames(df_t) <- c("ERBB2", "MYC")
svglite("/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/MYC_vs_ERBB2_copy_number_correlation_SKGT2_05272025.svg", width = 8, height = 6)
plot(df_t[,"ERBB2"], df_t[,"MYC"])
model <- lm(as.vector(df_t[,"MYC"]) ~ as.vector(df_t[,"ERBB2"]))
summary_model <- summary(model)
summary_model$r.squared #0.6856078
summary_model$coefficients #P-value=1.870845e-33
abline(lm(as.vector(df_t[,"MYC"]) ~ as.vector(df_t[,"ERBB2"])), col = "hotpink2", lwd = 4)
cor(df_t[,"ERBB2"], df_t[,"MYC"])#0.8280143
dev.off()



#To plot kde plots for ERBB2 gene copy number/ecDNA counts across single cells in NCIH2170 and SKGT2
df1_2170 <- read.csv("/home/yue1118/Bioskryb_data_analysis/NCIH2170_ecDNA_counts_03082025.csv")
View(df1_2170)
df1_2170$cell_line <- "NCIH2170"
df1_2170 <- df1_2170[,-1]
df1_2170 <- df1_2170[, c("ecDNA_cn","cell_line")]
colnames(df1_2170) <- c("ERBB2_CN","cell_line")
df1$cell_line <- "SKGT2"
df1_SKGT2 <- df1[, c("ERBB2_CN", "cell_line")]
df4 <- rbind(df1_SKGT2, df1_2170)
ggplot(df4, aes(x = `ERBB2_CN`, color = `cell_line`, fill = `cell_line`)) +
  geom_density(alpha = 0.3) +
  labs(title = "KDE of Two Cell Lines", x = "Value", y = "Density") +
  theme_minimal()
ggsave('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/SKGT2_NCIH2170_ERBB2_CN_KDE_plot_05272025.svg', width = 6, height = 6)
