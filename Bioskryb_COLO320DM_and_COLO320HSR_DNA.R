library(ggplot2)
library(dplyr)
rds1 <- readRDS("/home/yue1118/Bioskryb_data_analysis/ginkgo_res.binsize_1000000_COLO320DM_HSR_naive_JD_05302025.RDS")


##To process the COLO320DM data:
#To extract the normalized coverage numbers
df_coverage <- rds1$SegNorm
rownames(df_coverage) <- paste(df_coverage[,1], df_coverage[,2], df_coverage[,3], sep = "-")
df_coverage <- df_coverage[,4:ncol(df_coverage)]
df_coverage_COLO320DM <- df_coverage[, grepl("ColoDM", names(df_coverage))]

View(df_coverage_COLO320DM[grepl("chr8",rownames(df_coverage_COLO320DM)),])

coverage_quantiles <- apply(df_coverage_COLO320DM, 2, quantile, probs = 1.0, na.rm = TRUE)
# Function to filter values below the 100th percentile and compute median
filtered_coverage_means <- sapply(1:ncol(df_coverage_COLO320DM), function(i) {
  col_values <- df_coverage_COLO320DM[, i]  # Extract column values
  filtered_values <- col_values[col_values < coverage_quantiles[i]]  # Keep values below 90th percentile
  mean(filtered_values, na.rm = TRUE)  #We use the mean value as coverage in "normal" chromosomal regions #I think we can use this way to identity the coverage for regions without CNV on chromosome
})                                      #in that for DepMap COLO320DM CN profile (Copy number log2 (relative to ploidy +1)), it centers at 1, which means that most of genes don't have CNV events.
#This further means that the mean coverage values in the chromosomal regions are the coverage for "normal" regions without CNV events.

df1 <- data.frame(Cell = names(df_coverage_COLO320DM), "mean_coverage" = filtered_coverage_means)
df1$chr8_coverage <- as.vector(t(df_coverage_COLO320DM[c("chr8-126598157-127691734"),]))
df1$chr8_region_CN <- 2.4*df1$chr8_coverage/df1$mean_coverage
df1 <- df1[order(-df1$chr8_region_CN), ]

write.csv(df1, "/home/yue1118/Bioskryb_data_analysis/COLO320DM_Naive_nonsorted_Bioskryb_ecDNA_chr8-126598157-127691734_bin_counts_06042025.csv")

df1 <- read.csv("/home/yue1118/Bioskryb_data_analysis/COLO320DM_Naive_nonsorted_Bioskryb_ecDNA_chr8-126598157-127691734_bin_counts_06042025.csv")

View(df1)

df1$chr8_region_real_CN <- df1$chr8_region_CN/4 #from the ecDNA hub nature paper by Kung, we know there are 4 copies of "CCAT1+PVT1" on ecDNA

#To import FISH ground-truth imaging data
df_fish <- read.csv("/home/yue1118/Bioskryb_data_analysis/EcDNA_FISH_ COLO320DM_qPCR_for_python_analysis.csv")
View(df_fish)


wilcox.test(df1$chr8_region_real_CN, df_fish$Total_ecDNA)

svg("/home/yue1118/Bioskryb_data_analysis/COLO320DM_ecDNA_inferred_copy_number_vs_FISH_imaging.svg", width = 6, height = 4)
boxplot(df1$chr8_region_real_CN, df_fish$Total_ecDNA,
        names = c("scDNA", "FISH_imaging"),
        ylab = "Value")
dev.off()

svg("/home/yue1118/Bioskryb_data_analysis/COLO320DM_ecDNA_actual_inferred_copy_number_11302026.svg", width = 6, height = 6)
barplot(df1$chr8_region_real_CN)
dev.off()

median(df1$chr8_region_real_CN) #43.8
median(df_fish$Total_ecDNA) #41

##To process the COLO320HSR data:
df_coverage_COLO320HSR <- df_coverage[, grepl("ColoHS", names(df_coverage))]

View(df_coverage_COLO320HSR[grepl("chr8",rownames(df_coverage_COLO320HSR)),])

coverage_quantiles <- apply(df_coverage_COLO320HSR, 2, quantile, probs = 1.0, na.rm = TRUE)
# Function to filter values below the 100th percentile and compute median
filtered_coverage_means <- sapply(1:ncol(df_coverage_COLO320HSR), function(i) {
  col_values <- df_coverage_COLO320HSR[, i]  # Extract column values
  filtered_values <- col_values[col_values < coverage_quantiles[i]]  # Keep values below 90th percentile
  mean(filtered_values, na.rm = TRUE)  #We use the mean value as coverage in "normal" chromosomal regions #I think we can use this way to identity the coverage for regions without CNV on chromosome
})                                      #in that for DepMap COLO320DM CN profile (Copy number log2 (relative to ploidy +1)), it centers at 1, which means that most of genes don't have CNV events.
#This further means that the mean coverage values in the chromosomal regions are the coverage for "normal" regions without CNV events.

df1 <- data.frame(Cell = names(df_coverage_COLO320HSR), "mean_coverage" = filtered_coverage_means)
df1$chr8_coverage <- as.vector(t(df_coverage_COLO320HSR[c("chr8-126598157-127691734"),]))
df1$chr8_region_CN <- 2.4*df1$chr8_coverage/df1$mean_coverage
df1 <- df1[order(-df1$chr8_region_CN), ]
View(df1)
write.csv(df1, "/home/yue1118/Bioskryb_data_analysis/COLO320HSR_Naive_nonsorted_Bioskryb_ecDNA_chr8-126598157-127691734_bin_counts_06042025.csv")


df_counts_DM <- read.csv("/home/yue1118/Bioskryb_data_analysis/COLO320DM_Naive_nonsorted_Bioskryb_ecDNA_chr8-126598157-127691734_bin_counts_06042025.csv")
q30 <- quantile(df_counts_DM$chr8_region_CN, 0.3, na.rm = TRUE)
q70 <- quantile(df_counts_DM$chr8_region_CN, 0.7, na.rm = TRUE)

# Assign categories
df_counts_DM <- df_counts_DM %>%
  mutate(CN_category = case_when(
    chr8_region_CN < q30 ~ "low",
    chr8_region_CN > q70 ~ "high",
    TRUE ~ "med"
  ))

write.csv(df_counts_DM,"/home/yue1118/Bioskryb_data_analysis/COLO320DM_Naive_nonsorted_Bioskryb_ecDNA_chr8-126598157-127691734_bin_counts_07192025.csv")
