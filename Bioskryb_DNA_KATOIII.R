#Bioskryb data analysis for cell lines NCIH2170
library(ggplot2)
setwd('/home/yue1118')

df_ccle <- read.csv("/home/yue1118/Bioskryb_data_analysis/KATOIII_DepMap_CN_log2_relative_to_ploidy.csv")
colnames(df_ccle) <- c("gene", "CN")
mean(as.vector(df_ccle$CN), na.rm=T) #1.036262 mean value
View(df_ccle)


####Ginkgo output####
rds1 <- readRDS("/home/yue1118/ecDNA_project_data_and_figures/ginkgo_res.binsize_1000000.RDS")

#To extract the normalized coverage numbers
df_coverage <- rds1$SegNorm
rownames(df_coverage) <- paste(df_coverage[,1], df_coverage[,2], df_coverage[,3], sep = "-")
df_coverage <- df_coverage[,4:ncol(df_coverage)]
df_coverage_KATOIII <- df_coverage[, grepl("KATO", names(df_coverage))]

View()


coverage_quantiles <- apply(df_coverage_KATOIII, 2, quantile, probs = 1.00, na.rm = TRUE)
# Function to filter values below the 100th percentile and compute median
filtered_coverage_means <- sapply(1:ncol(df_coverage_KATOIII), function(i) {
  col_values <- df_coverage_KATOIII[, i]  # Extract column values
  filtered_values <- col_values[col_values < coverage_quantiles[i]]  # Keep values below 100th percentile
  mean(filtered_values, na.rm = TRUE)  #We use the mean value as coverage in "normal" chromosomal regions #I think we can use this way to identity the coverage for regions without CNV on chromosome
})                                      #in that for DepMap KATOIII CN profile (Copy number log2 (relative to ploidy +1)), it centers at 1, which means that most of genes don't have CNV events.
#This further means that the mean coverage values in the chromosomal regions are the coverage for "normal" regions without CNV events.


# Convert result to a dataframe
df1 <- data.frame(Cell = names(df_coverage_KATOIII), "mean_coverage" = filtered_coverage_means)
View(df1)


#Now extract the chr10 region that is on ecDNA: chr10-121134447-122166917, which covers FGFR2
df1$chr10_coverage <- as.vector(t(df_coverage_KATOIII[c("chr10-121134447-122166917"),]))
df1$chr10_copy_number <- 4*df1$chr10_coverage/df1$mean_coverage #KATOIII is hypotetraploid according to ATCC (https://www.atcc.org/products/htb-103)
hist(df1$chr10_copy_number)
mean(df1$chr10_copy_number)
median(df1$chr10_copy_number)


df1$phenotype <- ifelse(grepl("Low", df1$Cell, ignore.case = TRUE), "Low",
                        ifelse(grepl("High", df1$Cell, ignore.case = TRUE), "High",
                               ifelse(grepl("Med", df1$Cell, ignore.case = TRUE), "Med", NA)))
df1$phenotype <- factor(df1$phenotype, levels = c("Low", "Med", "High"))

ggplot(df1, aes(x = `phenotype`, y = `chr10_copy_number`)) +
  geom_boxplot() +
  labs(x = "Grouping Variable", y = "Value", title = "Boxplot by Group") +
  theme_minimal()
ggsave('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/inferred_ecDNA_counts_plate1_boxplot_FACS_FGFR2_high_med_low_KATOIII_05052025.svg', width = 6, height = 6)

write.csv(df1, "/home/yue1118/Bioskryb_data_analysis/KATOIII_FGFR2_copy_number_Bioskryb_05042025.csv")
View(df1)

wilcox.test(df1[df1$phenotype =="High","chr10_copy_number"], df1[df1$phenotype =="Low","chr10_copy_number"], paired = F) #p-value = 2.851e-09
wilcox.test(df1[df1$phenotype =="High","chr10_copy_number"], df1[df1$phenotype =="Med","chr10_copy_number"], paired = F) #p-value = 0.0005598
wilcox.test(df1[df1$phenotype =="Med","chr10_copy_number"], df1[df1$phenotype =="Low","chr10_copy_number"], paired = F) #p-value = 0.5038


df1 <- df1[order(-df1$chr10_copy_number), ]

#To plot FGFR2 copy number distribution across single cells
svg('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/FGFR2_KATOIII_plate1_copy_number_barplot_by_scDNA_seq_05052025.svg',width = 6, height = 6)
barplot(df1$chr10_copy_number)
dev.off()

df1_SNU16 <- read.csv("/home/yue1118/Bioskryb_data_analysis/SNU16_ecDNA_counts_04292025.csv")
View(df1_SNU16)
df1_SNU16$cell_line <- "SNU16"
df1_SNU16 <- df1_SNU16[,-1]
df1$cell_line <- "KATOIII"
df4 <- rbind(df1, df1_SNU16)
ggplot(df4, aes(x = `chr10_copy_number`, color = `cell_line`, fill = `cell_line`)) +
  geom_density(alpha = 0.3) +
  labs(title = "KDE of Three Groups", x = "Value", y = "Density") +
  theme_minimal()
ggsave('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/SNU16_KATOIII_FGFR2_CN_KDE_plot_05052025.svg', width = 6, height = 6)
