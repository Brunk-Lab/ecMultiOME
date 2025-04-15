#Bioskryb data analysis for cell lines NCIH2170
library(ggplot2)
setwd('/home/yue1118')


#After investigation, we find in the "ginkgo_res.binsize_1000000.RDS" file:
#1. SegNorm is the normalized coverage (read counts) for each bin
#2. SegFixed is the raw copy number profile (RCNP), with mean value of 1 in each cell according to Ginkgo source publication: "Interactive analysis and assessment of single-cell copy-number variations"
#I think the numbers in SegFixed can be considered as copy number of genes per chromosome. Therefore it is the CN before scaling and need proper scaling by multiply the number of chromosomes (ploidy).
#And I think mean=1 means that for "normal" genes without amplification or depletion on the HSR level, its CN on a chromosome should be 1.
#3. SegCopy has the calculated final copy number call for each bin. (SegCopy=the rounded up value of SegFixed*X)
#4. According to the numbers in SegFixed, it seems that after segmentation, chr8 regions next to but outside of ecDNA have the same number of chr8 regions on ecDNA. This combined, longer chr8 region
#is regarded with the same copy number, which is not accurate, and not proper to be used for ecDNA gene copy number calculation
#5. It seems that Ginkgo only have the ploidy number range of 1.5-6, which is not accurate in the scenario of ecDNA copy number estimation either.

#Therefore, considering these factors above, we can't use the numbers in CN final call in SegCopy to estimate copy number of regions on ecDNA.

#Alternatively, we use another way to estimate ecDNA copy number
#We extract the estimated copy number (CN1) and normalized read count (Count1) values for "normal" (without amplification or depletion) chromosomal regions.
#We then extract the normalized read count values (Count2) for ecDNA regions. The copy number of chromosomal region amplified on ecDNA would be: CN1*Count2/Count1.
#This is supported by the sentence "the same number of reads per bin should separate every sequential CN state, e.g., ~50 reads for CN 1, ~100 reads for CN 2, ~150 reads for CN 3, etc" 
#from the Ginkgo source publication
#For CN1, in the source paper, it says CN1=RCNP*X

#In the source paper, it also says "The most direct approach for determining the CN state of each cell is 
#available for users that have a priori knowledge of the ploidy of each sample."
#I think this is the optimal way of determining X, even better than the numerical optimization (SoS) method described in the paper.
#In prior, we know that the ploidy of NCIH2170 is near triploid to near hexaploid , so here we take X=4.5 and choose not to use the values in SegCopy (using the numerical optimization method, 
#X for many cells are estimated to be 2 or less than 2, which is not accurate)
#In this way, the CN of a "normal" region would be 1*4.5=4.5


rds1 <- readRDS("/home/yue1118/ecDNA_project_data_and_figures/ginkgo_res.binsize_1000000.RDS")

#To extract the normalized coverage numbers
df_coverage <- rds1$SegNorm
rownames(df_coverage) <- paste(df_coverage[,1], df_coverage[,2], df_coverage[,3], sep = "-")
df_coverage <- df_coverage[,4:ncol(df_coverage)]
df_coverage_2170 <- df_coverage[, grepl("2170", names(df_coverage))]


#We use the chr8 region to estimate copy number in the case of 1000kb bin size, since the chr17 region on ecDNA is more separated into many bins. 
#In the case of NCIH2170, for the region of chr8-127412678-128444056, 90% of it is on ecDNA, so it is a good one to estimate CN of ecDNA.
#We extract the coverage value for "normal" (without amplification or deletion) chromosome regions (by using quantile < 0.9) for each cell
coverage_quantiles <- apply(df_coverage_2170, 2, quantile, probs = 0.90, na.rm = TRUE)
# Function to filter values below the 90th percentile and compute median
filtered_coverage_means <- sapply(1:ncol(df_coverage_2170), function(i) {
  col_values <- df_coverage_2170[, i]  # Extract column values
  filtered_values <- col_values[col_values < coverage_quantiles[i]]  # Keep values below 90th percentile
  mean(filtered_values, na.rm = TRUE)  #We use the mean value as coverage in "normal" chromosomal regions #I think we can use this way to identity the coverage for regions without CNV on chromosome
})                                      #in that for DepMap NCIH2170 CN profile (Copy number log2 (relative to ploidy +1)), it centers at 1, which means that most of genes don't have CNV events.
                                        #This further means that the mean coverage values in the chromosomal regions are the coverage for "normal" regions without CNV events.

# Convert result to a dataframe
df1 <- data.frame(Cell = names(df_coverage_2170), "mean_coverage" = filtered_coverage_means)
View(df1)

#Now extract the chr8 region that is on ecDNA: chr8-127412678-128444056
df1$chr8_coverage <- as.vector(t(df_coverage_2170[c("chr8-127412678-128444056"),]))
df1$chr8_copy_number <- 4.5*df1$chr8_coverage/df1$mean_coverage
df1$ecDNA_cn <- df1$chr8_copy_number/2
write.csv(df1, "/home/yue1118/Bioskryb_data_analysis/NCIH2170_ecDNA_counts_03082025.csv")

hist(df1$chr8_copy_number)
boxplot(df1$chr8_copy_number)
boxplot(df1$ecDNA_cn)

rownames_vector <- df1$Cell
# Assign "phenotype" based on row names
df1$phenotype <- ifelse(grepl("Low", rownames_vector, ignore.case = TRUE), "Low",
                       ifelse(grepl("High", rownames_vector, ignore.case = TRUE), "High",
                              ifelse(grepl("Med", rownames_vector, ignore.case = TRUE), "Med", NA)))
df1$phenotype <- factor(df1$phenotype, levels = c("Low", "Med", "High"))
df1 <- df1[order(-df1$ecDNA_cn), ]
View(df1)

#To plot MYC copy number distribution across single cells
svg('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/MYC_copy_number_by_scDNA_seq_0307025.svg',width = 6, height = 6)
barplot(df1$chr8_copy_number)
dev.off()

#To plot ecDNA counts distribution across single cells
svg('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/ecDNA_counts_distribution_by_scDNA_seq_0307025.svg',width = 6, height = 6)
barplot(df1$ecDNA_cn)
dev.off()

# Create a boxplot
ggplot(df1, aes(x = `phenotype`, y = `ecDNA_cn`)) +
  geom_boxplot() +
  labs(x = "Grouping Variable", y = "Value", title = "Boxplot by Group") +
  theme_minimal()
ggsave('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/ecDNA_counts_boxplot_high_med_low_NCIH2170_03072025.svg', width = 6, height = 6)

wilcox.test(df1[df1$phenotype =="High","ecDNA_cn"], df1[df1$phenotype =="Low","ecDNA_cn"], paired = F) #p-value = 1.62e-12
wilcox.test(df1[df1$phenotype =="High","ecDNA_cn"], df1[df1$phenotype =="Med","ecDNA_cn"], paired = F) #p-value = 0.04983
wilcox.test(df1[df1$phenotype =="Med","ecDNA_cn"], df1[df1$phenotype =="Low","ecDNA_cn"], paired = F) #p-value = 2.142e-06




#To compare inferred ecDNA counts with FISH imaging ecDNA counts
#We include only the "control 2170" cells to compare with Bioskryb data
df_inferred_ecDNA <- read.csv("/home/yue1118/Bioskryb_data_analysis/NCIH2170_ecDNA_counts_03082025.csv")
df_FISH_ecDNA <- read.csv("/home/yue1118/ecDNA_project_data_and_figures/2025_NCIH2170_Control_FISH_counts_03142025.csv")
df_FISH_ecDNA <- df_FISH_ecDNA[df_FISH_ecDNA$Desc == "2170_control",]
df_FISH_ecDNA$MYC_and_ERBB2 <- df_FISH_ecDNA$Total_ecDNA - df_FISH_ecDNA$ERBB2. - df_FISH_ecDNA$Unlabeled - df_FISH_ecDNA$MYC. #To only compare the ecDNAs with both MYC and ERBB2

df_inferred_ecDNA$technique <- "Bioskryb_single_cell_sequencing"
df_FISH_ecDNA$technique <- "FISH_imaging"
View(df_inferred_ecDNA)
View(df_FISH_ecDNA)

df1 <- df_inferred_ecDNA[, c("ecDNA_cn","technique")]
df2 <- df_FISH_ecDNA[, c("MYC_and_ERBB2", "technique")]
colnames(df2) <- c("ecDNA_cn", "technique")
df3 <- rbind(df1, df2)

ggplot(df3, aes(x = `technique`, y = `ecDNA_cn`)) +
  geom_boxplot() +
  labs(x = "Grouping Variable", y = "Value", title = "Boxplot by Group") +
  theme_minimal()
ggsave('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/ecDNA_counts_boxplots_Bioskryb_vs_all_control_MYC_and_ERBB2_colocalized_ecDNA_NCIH2170_FISH_03142025.svg', width = 6, height = 6)

wilcox.test(df1$ecDNA_cn, df2$ecDNA_cn) #p-value = 0.1857




#To compute the correlation between bin chr8-127412678-128444056 (where MYC is located) and bin chr17-39284904-40331026 (where ERBB2 is located)

rds1 <- readRDS("/home/yue1118/ecDNA_project_data_and_figures/ginkgo_res.binsize_1000000.RDS")
df_coverage <- rds1$SegNorm
rownames(df_coverage) <- paste(df_coverage[,1], df_coverage[,2], df_coverage[,3], sep = "-")
df_coverage <- df_coverage[,4:ncol(df_coverage)]
df_coverage_2170 <- df_coverage[, grepl("2170", names(df_coverage))]
df_coverage_2170_2 <- t(df_coverage_2170)
df_coverage_2170_2 <- as.data.frame(df_coverage_2170_2)
df_t <- df_coverage_2170_2[,c('chr17-39284904-40331026','chr8-127412678-128444056')]
plot(df_t[,1], df_t[,2])
cor(df_t[,1], df_t[,2])
