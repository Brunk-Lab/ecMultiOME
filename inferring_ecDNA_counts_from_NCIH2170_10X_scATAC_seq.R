.libPaths()
library(epiAneufinder)
library(parallel)
library(doParallel)
library(data.table)
library(SummarizedExperiment)
library(ggplot2)

detectCores()
cl <- makeCluster(60)
registerDoParallel(cl)

epiAneufinder(input="/datastore/lbcfs/labs/brunk_lab/private/NCIH2170_multiome/cell_ranger_arc_results/NCIH2170_control/outs/atac_fragments.tsv.gz", #Enter path to your fragments.tsv file or the folder containing bam files
              outdir="epiAneufinder_results", #Path to the directory where results should be written 
              blacklist="/home/yue1118/epiAneufinder/hg38-blacklist.v2.bed", #Path to bed file that contains the blacklisted regions of your genome
              windowSize=100000, 
              genome="BSgenome.Hsapiens.UCSC.hg38", #Substitute with relevant BSgenome
              exclude=c('chrX','chrY','chrM'), 
              reuse.existing=TRUE,
              title_karyo="Karyogram of sample data", 
              ncores=60,
              minFrags=20000,
              minsizeCNV=0,
              k=4,
              plotKaryo=TRUE)


#100kb bin size
threshold_blacklist_bins=0.85 #(change this accordingly)
counts <- readRDS('/home/yue1118/epiAneufinder/epiAneufinder_results/epiAneufinder_results/count_summary.rds')
peaks <- as.data.table(assays(counts)$counts)
colnames(peaks) <- paste0('cell-', colnames(peaks))
rowinfo <- as.data.table(rowRanges(counts))
peaks <- cbind(rowinfo, peaks)
corrected_counts <- readRDS('/home/yue1118/epiAneufinder/epiAneufinder_results/epiAneufinder_results/counts_gc_corrected.rds')
corrected_counts <- cbind(rowinfo, corrected_counts)
#zeroes_per_bin <- peaks[, rowSums(.SD==0), .SDcols = patterns("cell-")]
#ncells <- length(grep("cell-", colnames(peaks)))
#peaks <- peaks[zeroes_per_bin<(threshold_blacklist_bins*ncells)]
corrected_counts$combined_position <- paste(corrected_counts$seqnames, corrected_counts$start, corrected_counts$end, sep = '_')
corrected_counts <- as.data.frame(corrected_counts)
rownames(corrected_counts) <- corrected_counts$combined_position
corrected_counts <- corrected_counts[, 13:(ncol(corrected_counts)-1)]
View(corrected_counts[grepl("chr8", rownames(corrected_counts)), ])  



#MYC(chr8) region vs ERBB2(chr17) region (100kb bin size)
df_chr8_region <- corrected_counts[rownames(corrected_counts) %in% c('chr8_127300001_127400000','chr8_127400001_127500000','chr8_127500001_127600000','chr8_127600001_127700000','chr8_127700001_127800000', 'chr8_127800001_127900000', 'chr8_127900001_128000000', 'chr8_128000001_128100000', 'chr8_128100001_128200000', 'chr8_128200001_128300000', 'chr8_128300001_128400000', 'chr8_128400001_128500000', 'chr8_128500001_128600000', 'chr8_128600001_128700000','chr8_128700001_128800000','chr8_128800001_128900000'),]
df_chr8_region_mean <- as.data.frame(apply(df_chr8_region, 2, mean))
colnames(df_chr8_region_mean) <- c('chr8_region_mean')
View(df_chr8_region_mean)

df_1 <- as.data.frame(colMeans(corrected_counts))

df_chr8_region_mean$MYC_ratio <- df_chr8_region_mean$chr8_region_mean/df_1$`colMeans(corrected_counts)`
df_chr8_region_mean$MYC_copy_number <- df_chr8_region_mean$MYC_ratio*4.5
df_chr8_region_mean$ecDNA_cn <- df_chr8_region_mean$MYC_ratio*4.5/2
df_chr8_region_mean <- df_chr8_region_mean[df_chr8_region_mean$ecDNA_cn >0,]
hist(df_chr8_region_mean$ecDNA_cn)
median(df_chr8_region_mean$ecDNA_cn) #244.2693
mean(df_chr8_region_mean$ecDNA_cn) #259.9585
wilcox.test(df_FISH_ecDNA$MYC_and_ERBB2, df_chr8_region_mean$ecDNA_cn)#p-value = 0.3821


df_FISH_ecDNA <- read.csv("/home/yue1118/ecDNA_project_data_and_figures/2025_NCIH2170_Control_FISH_counts_03142025.csv")
df_FISH_ecDNA <- df_FISH_ecDNA[df_FISH_ecDNA$Desc == "2170_control",]
df_FISH_ecDNA$MYC_and_ERBB2 <- df_FISH_ecDNA$Total_ecDNA - df_FISH_ecDNA$ERBB2. - df_FISH_ecDNA$Unlabeled - df_FISH_ecDNA$MYC. #To only compare the ecDNAs with both MYC and ERBB2
View(df_FISH_ecDNA)

df_inferred_ecDNA_plate1 <- read.csv("/home/yue1118/Bioskryb_data_analysis/NCIH2170_ecDNA_counts_03082025.csv")
df_inferred_ecDNA_plate1$technique <- "Bioskryb_scDNAseq"
df3 <- df_inferred_ecDNA_plate1[, c("ecDNA_cn","technique")]

df_chr8_region_mean$technique <- "10X_scATAC"
df_FISH_ecDNA$technique <- "FISH_imaging"
df1 <- df_FISH_ecDNA[,c("technique", "MYC_and_ERBB2")]
colnames(df1) <- c("technique", "ecDNA_cn")
df2 <- df_chr8_region_mean[,c("technique", "ecDNA_cn")]
df4 <- rbind(df1, df2, df3)


ggplot(df4, aes(x = `technique`, y = `ecDNA_cn`)) +
  geom_boxplot() +
  labs(x = "Grouping Variable", y = "Value", title = "Boxplot by Group") +
  theme_minimal()
ggsave('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/ecDNA_counts_inferred_by_10X_multiome_scATAC_data_04232025.svg', width = 6, height = 6)

ggplot(df4, aes(x = `ecDNA_cn`, color = `technique`, fill = `technique`)) +
  geom_density(alpha = 0.3) +
  labs(title = "KDE of Three Groups", x = "Value", y = "Density") +
  theme_minimal()
ggsave('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/inferred_ecDNA_counts_kde_10X_scATAC_Bioskryb_scDNA_FISH_imaging_04252025.svg', width = 6, height = 6)


df2 <- df2[order(-df2$ecDNA_cn), ]
#To plot ecDNA counts distribution across single cells
svg('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/ecDNA_counts_distribution_by_10X_scATAC_04232025.svg',width = 6, height = 6)
barplot(df2$ecDNA_cn)
dev.off()

df_chr8_region_mean <- df_chr8_region_mean[order(-df_chr8_region_mean$MYC_copy_number), ]
svg('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/MYC_copy_number_counts_distribution_by_10X_scATAC_04232025.svg',width = 6, height = 6)
barplot(df_chr8_region_mean$MYC_copy_number)
dev.off()
