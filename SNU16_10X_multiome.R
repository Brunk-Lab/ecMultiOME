#Set working directory
setwd("/home/yue1118/TNBC_scRNA_analysis_workflow/")
.libPaths("/home/yue1118/R/x86_64-pc-linux-gnu-library/4.2")
.libPaths("/datastore/clusterapps/rocky9/R-4.2.1/lib64/R/library")

library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(EnsDb.Hsapiens.v86)
library(MACSr)
library(dplyr)
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg38)
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(EnhancedVolcano)
library(circlize)
library(reshape)
library(reshape2)
require(plyr)
library(ComplexHeatmap)
library(gridBase)
library(svglite)
library(sets)
library(clusterProfiler)
library(stringr)
library(pheatmap)
library(viridis)
library(DoubletFinder)
library(GSVA)
library(msigdbr)
library(tidyr)

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

# load the RNA and ATAC data
control.multiome <- Read10X_h5("/datastore/lbcfs/labs/brunk_lab/public/SNU16_10X_multiome/cell_ranger_arc_results/SNU16_multiome/outs/filtered_feature_bc_matrix.h5")
fragpath <- "/datastore/lbcfs/labs/brunk_lab/public/SNU16_10X_multiome/cell_ranger_arc_results/SNU16_multiome/outs/atac_fragments.tsv.gz"

# create a Seurat object containing the RNA data
SNU16_multiome <- CreateSeuratObject(counts = control.multiome$`Gene Expression`,assay = "RNA")

# create ATAC assay and add it to the object
SNU16_multiome[["ATAC"]] <- CreateChromatinAssay(counts = control.multiome$Peaks, sep = c(":", "-"), fragments = fragpath, annotation = annotation)

DefaultAssay(SNU16_multiome) <- 'RNA'
SNU16_multiome[["percent.mt"]] <- PercentageFeatureSet(SNU16_multiome, pattern = "^MT-")

DefaultAssay(SNU16_multiome) <- "ATAC"

SNU16_multiome <- NucleosomeSignal(SNU16_multiome)
SNU16_multiome <- TSSEnrichment(SNU16_multiome)

p <- VlnPlot(SNU16_multiome, features = c("nCount_ATAC", "nCount_RNA", "nFeature_RNA", "TSS.enrichment", "nucleosome_signal", "percent.mt"), ncol = 6,
             log = TRUE, pt.size = 0) + NoLegend()

p

SNU16_multiome_2 <- subset(
  x = SNU16_multiome,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 1000 &
    nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    nFeature_RNA > 500  &
    nucleosome_signal < 2 &
    TSS.enrichment > 1 &
    percent.mt < 20
)

dim(SNU16_multiome_2)

# ATAC analysis
DefaultAssay(SNU16_multiome_2) <- "ATAC"
SNU16_multiome_peaks <- CallPeaks(SNU16_multiome_2, macs2.path = "~/miniconda3/envs/macs3_env/bin/macs3")
SNU16_multiome_peaks <- keepStandardChromosomes(SNU16_multiome_peaks, pruning.mode = "coarse")
SNU16_multiome_peaks <- subsetByOverlaps(x = SNU16_multiome_peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
SNU16_macs2_counts <- FeatureMatrix(
  fragments = Fragments(SNU16_multiome_2),
  features = SNU16_multiome_peaks,
  cells = colnames(SNU16_multiome_2)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
SNU16_multiome_2[["peaks"]] <- CreateChromatinAssay(
  counts = SNU16_macs2_counts,
  fragments = fragpath,
  annotation = annotation
)


DefaultAssay(SNU16_multiome_2) <- "peaks"
SNU16_multiome_2 <- FindTopFeatures(SNU16_multiome_2, min.cutoff = 5)
SNU16_multiome_2 <- RunTFIDF(SNU16_multiome_2)
SNU16_multiome_2 <- RunSVD(SNU16_multiome_2)
SNU16_multiome_2 <- RunUMAP(SNU16_multiome_2, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

SNU16.gene.activities <- GeneActivity(SNU16_multiome_2, biotypes = NULL)# to include not only protein coding genes but also long non coding genes
SNU16_multiome_2[['gene_activity']] <- CreateAssayObject(counts = SNU16.gene.activities)
SNU16_multiome_2 <- NormalizeData(
  object = SNU16_multiome_2,
  assay = 'gene_activity',
  normalization.method = 'LogNormalize',
  scale.factor = median(SNU16_multiome_2$nCount_RNA)
)

#saveRDS(SNU16_multiome_2, "/home/yue1118/ecDNA_project_data_and_figures/SNU16_multiome_intermediate_object_08022025.rds")
SNU16_multiome_2 <- readRDS("/home/yue1118/ecDNA_project_data_and_figures/SNU16_multiome_intermediate_object_08022025.rds")

# RNA analysis
DefaultAssay(SNU16_multiome_2) <- "RNA"
SNU16_multiome_2 <- SCTransform(SNU16_multiome_2, vars.to.regress = "percent.mt", verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_', seed.use = 7)



#Now run DoubletFinder to remove potential doublets from the dataset
suppressMessages(require(DoubletFinder))
#To look for the best pK value
sweep.res.list_control <- paramSweep(SNU16_multiome_2, PCs = 1:20, sct=TRUE)
sweep.stats_control <- summarizeSweep(sweep.res.list_control, GT = FALSE)
bcmvn_control <- find.pK(sweep.stats_control)
pK <- bcmvn_control %>% filter(BCmetric == max(BCmetric)) %>% select(pK) #pK=0.04
#Homotypic Doublet Proportion Estimate
cluster_annotation <- SNU16_multiome_2@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(cluster_annotation)
#To estimate the number of doublets
nExp_poi <- round(0.012*nrow(SNU16_multiome_2@meta.data)) #assuming there are 1.2% doublets
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
SNU16_multiome_2 <- DoubletFinder::doubletFinder(SNU16_multiome_2, pN = 0.25, pK = 0.04, nExp = nExp_poi.adj, PCs = 1:10, sct = TRUE)

table(SNU16_multiome_2$DF.classifications_0.25_0.04_1)

#saveRDS(SNU16_multiome_2, "/home/yue1118/ecDNA_project_data_and_figures/SNU16_multiome_intermediate_object_08022025.rds")
SNU16_multiome_2 <- readRDS("/home/yue1118/ecDNA_project_data_and_figures/SNU16_multiome_intermediate_object_08022025.rds")

##ChromVAR assay
DefaultAssay(SNU16_multiome_2) <- "ATAC"
main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
keep.peaks <- which(as.character(seqnames(granges(SNU16_multiome_2))) %in% main.chroms)
SNU16_multiome_2[["ATAC"]] <- subset(SNU16_multiome_2[["ATAC"]], features = rownames(SNU16_multiome_2[["ATAC"]])[keep.peaks])

DefaultAssay(SNU16_multiome_2) <- "ATAC"
pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE))
motif.matrix <- CreateMotifMatrix(features = granges(SNU16_multiome_2), pwm = pwm_set, genome = 'hg38', use.counts = FALSE)
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)
SNU16_multiome_2 <- SetAssayData(SNU16_multiome_2, assay = 'ATAC', slot ='motifs', new.data = motif.object)


# Note that this step can take 30-60 minutes 
SNU16_multiome_2 <- RunChromVAR(
  object = SNU16_multiome_2,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

saveRDS(SNU16_multiome_2, "/home/yue1118/ecDNA_project_data_and_figures/SNU16_multiome_doublet_correctly_removed_with_ChromVAR_08022025.rds")



# Volcano plots for DEG and differential ATAC peaks (gene activity) between ERBB2 mRNA high and ERBB2 mRNA low cells in NCIH2170
group.func <- function(x){
  if (x %in% rownames(df3)) {
    "high"
  } else if (x %in% rownames(df4)) {
    "low"
  } else {
    "medium"
  }
}

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

i = "FGFR2"
df_sct <- as.data.frame(SNU16_multiome_2@assays$gene_activity@data[i,])
colnames(df_sct) <- 'gene'
df_sct$cells <- rownames(df_sct)
df3 <- df_sct[df_sct$gene >= quantile(as.numeric(SNU16_multiome_2@assays$gene_activity@data[i,]), .90, na.rm=TRUE), ]
df3$group <- 'higher'
df4 <- df_sct[df_sct$gene <= quantile(as.numeric(SNU16_multiome_2@assays$gene_activity@data[i,]), .10, na.rm=TRUE), ]
df4$group <- 'lower'
SNU16_multiome_2@meta.data$group <- lapply(colnames(SNU16_multiome_2), group.func)
Idents(SNU16_multiome_2) <- SNU16_multiome_2$group
DefaultAssay(SNU16_multiome_2) <- "SCT"
df_mRNA_markers<- FindMarkers(SNU16_multiome_2, ident.1 = "high", ident.2 = "low", min.pct = 0.25, logfc.threshold = 0.10, recorrect_umi=F)
DefaultAssay(SNU16_multiome_2) <- "gene_activity"
df_atac_markers <- FindMarkers(SNU16_multiome_2, ident.1 = "high", ident.2 = "low", min.pct = 0.25, logfc.threshold = 0.10, recorrect_umi=F)


i = "MYC"
df_sct <- as.data.frame(SNU16_multiome_2@assays$gene_activity@data[i,])
colnames(df_sct) <- 'gene'
df_sct$cells <- rownames(df_sct)
df3 <- df_sct[df_sct$gene >= quantile(as.numeric(SNU16_multiome_2@assays$gene_activity@data[i,]), .90, na.rm=TRUE), ]
df3$group <- 'higher'
df4 <- df_sct[df_sct$gene <= quantile(as.numeric(SNU16_multiome_2@assays$gene_activity@data[i,]), .10, na.rm=TRUE), ]
df4$group <- 'lower'
SNU16_multiome_2@meta.data$group <- lapply(colnames(SNU16_multiome_2), group.func)
Idents(SNU16_multiome_2) <- SNU16_multiome_2$group
DefaultAssay(SNU16_multiome_2) <- "SCT"
df_mRNA_markers_MYC<- FindMarkers(SNU16_multiome_2, ident.1 = "high", ident.2 = "low", min.pct = 0.25, logfc.threshold = 0.10, recorrect_umi=F)
DefaultAssay(SNU16_multiome_2) <- "gene_activity"
df_atac_markers_MYC <- FindMarkers(SNU16_multiome_2, ident.1 = "high", ident.2 = "low", min.pct = 0.25, logfc.threshold = 0.10, recorrect_umi=F)


View(df_mRNA_markers_MYC[df_mRNA_markers_MYC$p_val_adj < 0.05,])
View(df_atac_markers_MYC[df_atac_markers_MYC$p_val_adj < 0.05,])


#SNU16 pseudobulk

Idents(SNU16_multiome_2) <- 'WT'
df_sc_aggre <- AggregateExpression(SNU16_multiome_2, assays = 'SCT', slot = 'data')
df_sc_aggre <- data.frame(df_sc_aggre)
colnames(df_sc_aggre) <- c('aggregated_value')
df_sc_aggre$log2_value <- log2(df_sc_aggre$aggregated_value)
df_sc_aggre <- subset(df_sc_aggre, log2_value >0)
df_sc_aggre <- df_sc_aggre[order(df_sc_aggre$log2_value,decreasing=TRUE),]
df_sc_aggre$gene <- rownames(df_sc_aggre)
df_sc_aggre <- df_sc_aggre[1:100,]
View(df_sc_aggre)

##Now run ssGSEA for SNU16 single cells
msig_hallmark <- msigdbr(species = "Homo sapiens", category = "H")
msig_hallmark <- msig_hallmark[, c("gs_name", "gene_symbol")]
gene_sets_list <- split(msig_hallmark$gene_symbol, msig_hallmark$gs_name)

#ssGSEA for SNU16
SNU16_multiome_2 <- readRDS('/home/yue1118/ecDNA_project_data_and_figures/SNU16_multiome_doublet_correctly_removed_with_ChromVAR_08022025.rds')
X_SNU16_multiome <- as.matrix(SNU16_multiome_2@assays$SCT@data)

gsvaPar_SNU16_multiome <- gsvaParam(X_SNU16_multiome, gene_sets_list)
es_SNU16_multiome <- gsva(gsvaPar_SNU16_multiome, verbose=FALSE)
es_SNU16_multiome <- as.data.frame(t(es_SNU16_multiome))
View(es_SNU16_multiome)

#write.csv(es_SNU16_multiome, '/home/yue1118/ecDNA_project_data_and_figures/SNU16_10X_multiome_ssGSEA_on_MsigDB_hallmark_pathways_08022025.csv')
es_SNU16_multiome <- read.csv("/home/yue1118/ecDNA_project_data_and_figures/SNU16_10X_multiome_ssGSEA_on_MsigDB_hallmark_pathways_08022025.csv")
table(es_SNU16_multiome$X == colnames(SNU16_multiome_2))
es_SNU16_multiome$PVT1_exp <- SNU16_multiome_2@assays$SCT@data['PVT1',]
es_SNU16_multiome$PVT1_ATAC <- SNU16_multiome_2@assays$gene_activity@data['PVT1',]
es_SNU16_multiome$FGFR2_ATAC <- SNU16_multiome_2@assays$gene_activity@data['FGFR2',]
es_SNU16_multiome$FGFR2_exp <- SNU16_multiome_2@assays$SCT@data['FGFR2',]
es_SNU16_multiome$MYC_ATAC <- SNU16_multiome_2@assays$gene_activity@data['MYC',]


plot(es_SNU16_multiome$HALLMARK_MYC_TARGETS_V1, es_SNU16_multiome$PVT1_exp)
plot(es_SNU16_multiome$HALLMARK_MYC_TARGETS_V1, es_SNU16_multiome$PVT1_ATAC)
plot(es_SNU16_multiome$HALLMARK_MYC_TARGETS_V1, es_SNU16_multiome$FGFR2_ATAC)
plot(es_SNU16_multiome$HALLMARK_MYC_TARGETS_V1, es_SNU16_multiome$FGFR2_exp)
plot(es_SNU16_multiome$FGFR2_ATAC, es_SNU16_multiome$FGFR2_exp)
cor.test(es_SNU16_multiome$HALLMARK_MYC_TARGETS_V1, es_SNU16_multiome$PVT1_exp)
cor.test(es_SNU16_multiome$HALLMARK_MYC_TARGETS_V2, es_SNU16_multiome$PVT1_exp)
cor.test(es_SNU16_multiome$HALLMARK_MYC_TARGETS_V1, es_SNU16_multiome$FGFR2_ATAC)
cor.test(es_SNU16_multiome$HALLMARK_MYC_TARGETS_V1, es_SNU16_multiome$FGFR2_exp)

SNU16_multiome_2 <- readRDS("/home/yue1118/ecDNA_project_data_and_figures/SNU16_multiome_doublet_correctly_removed_with_ChromVAR_08022025.rds")

svg('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/PVT1_vs_MYC_ATAC_correlation_10X_multiome_SNU16_08082025.svg',width = 6, height = 6)
plot(SNU16_multiome_2@assays$gene_activity@data['PVT1',], SNU16_multiome_2@assays$gene_activity@data['MYC',])
dev.off()

cor.test(SNU16_multiome_2@assays$gene_activity@data['PVT1',], SNU16_multiome_2@assays$gene_activity@data['MYC',]) #cor=0.197; p-value = 5.791e-05

dim(SNU16_multiome_2@meta.data) #413 cells


######################To compute correlation of each peak accessibility with PVT1 mRNA across single cells######################
SNU16_multiome_2 <- readRDS("/home/yue1118/ecDNA_project_data_and_figures/SNU16_multiome_doublet_correctly_removed_with_ChromVAR_08022025.rds")

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene) 
library(org.Hs.eg.db)
library(Matrix)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

#To calculate PVT1 mRNA-peak correlation
SNU16_pvt1_expr <- GetAssayData(SNU16_multiome_2, assay = "SCT", slot = "data")["PVT1", ]
SNU16_peak_mat <- GetAssayData(SNU16_multiome_2, assay = "peaks", slot = "data")


# Compute correlations
SNU16_cor_vals <- apply(SNU16_peak_mat, 1, function(x) cor(x, SNU16_pvt1_expr, use = "complete.obs", method = "pearson"))

# Optionally compute p-values
SNU16_pvals <- apply(SNU16_peak_mat, 1, function(x) cor.test(x, SNU16_pvt1_expr)$p.value)

# Put in a dataframe
SNU16_cor_df <- data.frame(
  peak = rownames(SNU16_peak_mat),
  correlation = SNU16_cor_vals,
  p_value = SNU16_pvals
)

SNU16_cor_df$padj <- p.adjust(SNU16_cor_df$p_value, method = "BH")
SNU16_cor_df <- SNU16_cor_df[SNU16_cor_df$padj < 0.05, ]
SNU16_cor_df <- na.omit(SNU16_cor_df)
SNU16_cor_df <- SNU16_cor_df[order(-SNU16_cor_df$correlation),]


peaks_split_SNU16 <- do.call(rbind, strsplit(SNU16_cor_df$peak, "[-]"))
colnames(peaks_split_SNU16) <- c("chr", "start", "end")
peaks_df_SNU16 <- data.frame( chr = peaks_split_SNU16[, 1],start = as.numeric(peaks_split_SNU16[, 2]),end = as.numeric(peaks_split_SNU16[, 3]))
granges_obj_SNU16 <- GRanges(seqnames = peaks_df_SNU16$chr,ranges = IRanges(start = peaks_df_SNU16$start, end = peaks_df_SNU16$end))


peakAnno_SNU16_PVT1_isoform <- annotatePeak(granges_obj_SNU16, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db") 
df_peakAnno_SNU16_PVT1_isoform <- as.data.frame(peakAnno_SNU16_PVT1_isoform)


df_peakAnno_SNU16_PVT1_isoform$peak <- paste(df_peakAnno_SNU16_PVT1_isoform$seqnames, df_peakAnno_SNU16_PVT1_isoform$start, df_peakAnno_SNU16_PVT1_isoform$end, sep = "-")
df_peakAnno_SNU16_PVT1_isoform <- merge(df_peakAnno_SNU16_PVT1_isoform, SNU16_cor_df, by.x="peak", by.y= "peak", all.x=TRUE)


df_peakAnno_SNU16_PVT1_isoform_promoter <- df_peakAnno_SNU16_PVT1_isoform[grepl("Promoter",df_peakAnno_SNU16_PVT1_isoform$annotation) ,]

df_peakAnno_SNU16_PVT1_isoform_promoter_plot <- df_peakAnno_SNU16_PVT1_isoform_promoter[(df_peakAnno_SNU16_PVT1_isoform_promoter$correlation < -0.24) | (df_peakAnno_SNU16_PVT1_isoform_promoter$correlation > 0.3),]
df_peakAnno_SNU16_PVT1_isoform_promoter_plot <- df_peakAnno_SNU16_PVT1_isoform_promoter_plot[order(-df_peakAnno_SNU16_PVT1_isoform_promoter_plot$correlation),]
View(df_peakAnno_SNU16_PVT1_isoform_promoter_plot)

mat <- matrix(df_peakAnno_SNU16_PVT1_isoform_promoter_plot$correlation, ncol = 1)
rownames(mat) <- df_peakAnno_SNU16_PVT1_isoform_promoter_plot$peak
colnames(mat) <- "Correlation"

svg('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/PVT1_exp_peak_correlation_10X_multiome_SNU16_08082025.svg',width = 10, height = 10)
pheatmap(mat, cluster_rows = FALSE, cluster_cols = FALSE, color = colorRampPalette(c("blue", "white", "red"))(50))
dev.off()

#Since PVT1, PDHX and APIP are both amplified on both ecDNA and HSR, we can see there are two populations of cells
svg('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/PVT1_ATAC_vs_EXP_correlation_10X_multiome_SNU16_08172025.svg',width = 6, height = 6)
plot(SNU16_multiome_2@assays$gene_activity@data['PVT1',], SNU16_multiome_2@assays$SCT@data['PVT1',])
dev.off()
svg('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/PDHX_ATAC_vs_EXP_correlation_10X_multiome_SNU16_08172025.svg',width = 6, height = 6)
plot(SNU16_multiome_2@assays$gene_activity@data['PDHX',], SNU16_multiome_2@assays$SCT@data['PDHX',])
dev.off()
svg('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/APIP_ATAC_vs_EXP_correlation_10X_multiome_SNU16_08172025.svg',width = 6, height = 6)
plot(SNU16_multiome_2@assays$gene_activity@data['APIP',], SNU16_multiome_2@assays$SCT@data['APIP',])
dev.off()
svg('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/MYC_ATAC_vs_EXP_correlation_10X_multiome_SNU16_08172025.svg',width = 6, height = 6)
plot(SNU16_multiome_2@assays$gene_activity@data['MYC',], SNU16_multiome_2@assays$SCT@data['MYC',])
dev.off()
svg('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/FGFR2_ATAC_vs_EXP_correlation_10X_multiome_SNU16_08172025.svg',width = 6, height = 6)
plot(SNU16_multiome_2@assays$gene_activity@data['FGFR2',], SNU16_multiome_2@assays$SCT@data['FGFR2',])
dev.off()

#In this subpopulation maybe PVT1 is on ecDNA (subset cells based on PVT1 ATAC)
DefaultAssay(SNU16_multiome_2) <- "gene_activity"
SNU16_multiome_2_subset_1 <- subset(x = SNU16_multiome_2, cells = colnames(SNU16_multiome_2)[FetchData(SNU16_multiome_2, vars = "PVT1")[,1] > 3.2])
cor.test(SNU16_multiome_2_subset_1@assays$gene_activity@data['MYC',], SNU16_multiome_2_subset_1@assays$gene_activity@data['PVT1',])#cor=0.4950367, p-value=2.2e-16
cor.test(SNU16_multiome_2_subset_1@assays$gene_activity@data['CASC11',], SNU16_multiome_2_subset_1@assays$gene_activity@data['PVT1',])#cor=0.4711789, p-value=2.2e-16
#APIP and PDHX forms fusion with PVT1 potentially on both ecDNA and HSR so we don't use these two genes to test the correlation.
svg('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/MYC_vs_PVT1_ATAC_correlation_subpopulation1_10X_multiome_SNU16_08172025.svg',width = 6, height = 6)
plot(SNU16_multiome_2_subset_1@assays$gene_activity@data['MYC',], SNU16_multiome_2_subset_1@assays$gene_activity@data['PVT1',])
dev.off()

svg('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/CASC11_vs_PVT1_ATAC_correlation_subpopulation1_10X_multiome_SNU16_08172025.svg',width = 6, height = 6)
plot(SNU16_multiome_2_subset_1@assays$gene_activity@data['CASC11',], SNU16_multiome_2_subset_1@assays$gene_activity@data['PVT1',])
dev.off()


#In this subpopulation maybe PVT1 is on HSR (subset cells based on PVT1 ATAC)
SNU16_multiome_2_subset_2 <- subset(x = SNU16_multiome_2, cells = colnames(SNU16_multiome_2)[FetchData(SNU16_multiome_2, vars = "PVT1")[,1] < 3.2])
cor.test(SNU16_multiome_2_subset_2@assays$gene_activity@data['MYC',], SNU16_multiome_2_subset_2@assays$gene_activity@data['PVT1',]) #cor=-0.07835157, p-value= 0.4655
cor.test(SNU16_multiome_2_subset_2@assays$gene_activity@data['CASC11',], SNU16_multiome_2_subset_2@assays$gene_activity@data['PVT1',]) #cor=-0.01239567, p-value = 0.9082
#APIP and PDHX forms fusion with PVT1 potentially on both ecDNA and HSR so we don't use these two genes to test the correlation.

svg('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/MYC_vs_PVT1_ATAC_correlation_subpopulation2_10X_multiome_SNU16_08172025.svg',width = 6, height = 6)
plot(SNU16_multiome_2_subset_2@assays$gene_activity@data['MYC',], SNU16_multiome_2_subset_2@assays$gene_activity@data['PVT1',])
dev.off()

svg('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/CASC11_vs_PVT1_ATAC_correlation_subpopulation2_10X_multiome_SNU16_08172025.svg',width = 6, height = 6)
plot(SNU16_multiome_2_subset_2@assays$gene_activity@data['CASC11',], SNU16_multiome_2_subset_2@assays$gene_activity@data['PVT1',])
dev.off()


cor.test(SNU16_multiome_2@assays$gene_activity@data['FGFR2',], SNU16_multiome_2@assays$gene_activity@data['PVT1',])
plot(SNU16_multiome_2@assays$gene_activity@data['FGFR2',], SNU16_multiome_2@assays$gene_activity@data['PVT1',])

#To compare FGFR2 gene expression between PVT1 HSR group and PVT1 ecDNA group
DefaultAssay(SNU16_multiome_2_subset_1) <- "SCT"
DefaultAssay(SNU16_multiome_2_subset_2) <- "SCT"
expr1 <- FetchData(SNU16_multiome_2_subset_1, vars = "FGFR2")[,1]
expr2 <- FetchData(SNU16_multiome_2_subset_2, vars = "FGFR2")[,1]
df <- data.frame(expression = c(expr1, expr2), group = c(rep("ecDNA", length(expr1)), rep("HSR", length(expr2))))
ggplot(df, aes(x = group, y = expression, fill = group)) +
  geom_boxplot() +
  theme_classic()
ggsave('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/FGFR2_exp_subpop1_subpop2_10X_multiome_SNU16_08172025.svg',width = 6, height = 6)
wilcox.test(expr1, expr2)


#To compare FGFR2 ATAC between PVT1 HSR group and PVT1 ecDNA group
DefaultAssay(SNU16_multiome_2_subset_1) <- "gene_activity"
DefaultAssay(SNU16_multiome_2_subset_2) <- "gene_activity"
expr1 <- FetchData(SNU16_multiome_2_subset_1, vars = "FGFR2")[,1]
expr2 <- FetchData(SNU16_multiome_2_subset_2, vars = "FGFR2")[,1]
df <- data.frame(ATAC = c(expr1, expr2), group = c(rep("ecDNA", length(expr1)), rep("HSR", length(expr2))))
ggplot(df, aes(x = group, y = ATAC, fill = group)) +
  geom_boxplot() +
  theme_classic()
ggsave('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/FGFR2_ATAC_subpop1_subpop2_10X_multiome_SNU16_08172025.svg',width = 6, height = 6)
wilcox.test(expr1, expr2)

#To compare MYC gene expression between PVT1 HSR group and PVT1 ecDNA group
DefaultAssay(SNU16_multiome_2_subset_1) <- "SCT"
DefaultAssay(SNU16_multiome_2_subset_2) <- "SCT"
expr1 <- FetchData(SNU16_multiome_2_subset_1, vars = "MYC")[,1]
expr2 <- FetchData(SNU16_multiome_2_subset_2, vars = "MYC")[,1]
df <- data.frame(expression = c(expr1, expr2), group = c(rep("ecDNA", length(expr1)), rep("HSR", length(expr2))))
ggplot(df, aes(x = group, y = expression, fill = group)) +
  geom_boxplot() +
  theme_classic()

wilcox.test(expr1, expr2)


#To perform DEG between HSR and ecDNA subpopulations in SNU16
DefaultAssay(SNU16_multiome_2) <- "gene_activity"

# pull per-cell PVT1 values (log/normalized "data" slot is typical)
pvt1 <- FetchData(SNU16_multiome_2, vars = "PVT1")[,1]   # or: GetAssayData(obj, slot="data")["PVT1", ]

# make labels
state <- ifelse(pvt1 > 3.2, "ecDNA", "HSR")
state[is.na(pvt1)] <- NA  # keep NAs if any

# add to metadata (as an ordered factor)
SNU16_multiome_2$PVT1_amp_state <-state
table(SNU16_multiome_2$PVT1_amp_state)

Idents(SNU16_multiome_2) <- "PVT1_amp_state"
DefaultAssay(SNU16_multiome_2) <- "SCT"
df_mRNA_markers <- FindMarkers(SNU16_multiome_2, ident.1 = "HSR", ident.2 = "ecDNA", min.pct = 0.25, logfc.threshold = 0.10, recorrect_umi=F)
View(df_mRNA_markers[df_mRNA_markers$p_val_adj < 0.05,])


#To see which program/phenotypes are associated with PVT1-ecDNA species
es_SNU16_multiome <- read.csv("/home/yue1118/ecDNA_project_data_and_figures/SNU16_10X_multiome_ssGSEA_on_MsigDB_hallmark_pathways_08022025.csv")
table(es_SNU16_multiome$X == colnames(SNU16_multiome_2))
es_SNU16_multiome$PVT1_exp <- SNU16_multiome_2@assays$SCT@data['PVT1',]
es_SNU16_multiome$PVT1_ATAC <- SNU16_multiome_2@assays$gene_activity@data['PVT1',]
es_SNU16_multiome$MYC_ATAC <- SNU16_multiome_2@assays$gene_activity@data['MYC',]
es_SNU16_multiome$FGFR2_ATAC <- SNU16_multiome_2@assays$gene_activity@data['FGFR2',]
es_SNU16_multiome$FGFR2_exp <- SNU16_multiome_2@assays$SCT@data['FGFR2',]

es_SNU16_multiome_subset_1 <- es_SNU16_multiome[es_SNU16_multiome$X %in% colnames(SNU16_multiome_2_subset_1),]


svg('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/MYC_target_V1_vs_PVT1_exp_correlation_subpopulation1_10X_multiome_SNU16_08172025.svg',width = 6, height = 6)
plot(es_SNU16_multiome_subset_1$HALLMARK_MYC_TARGETS_V1, es_SNU16_multiome_subset_1$PVT1_exp)
dev.off()
cor.test(es_SNU16_multiome_subset_1$HALLMARK_MYC_TARGETS_V1, es_SNU16_multiome_subset_1$PVT1_exp) #cor=-0.62,p-value < 2.2e-16

svg('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/mTORC1_signaling_vs_PVT1_exp_correlation_subpopulation1_10X_multiome_SNU16_08172025.svg',width = 6, height = 6)
plot(es_SNU16_multiome_subset_1$HALLMARK_MTORC1_SIGNALING, es_SNU16_multiome_subset_1$PVT1_exp)
dev.off()
cor.test(es_SNU16_multiome_subset_1$HALLMARK_MTORC1_SIGNALING, es_SNU16_multiome_subset_1$PVT1_exp) #cor=-0.4822505, p-value < 2.2e-16


svg('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/E2F_targets_vs_PVT1_exp_correlation_subpopulation1_10X_multiome_SNU16_08172025.svg',width = 6, height = 6)
plot(es_SNU16_multiome_subset_1$HALLMARK_E2F_TARGETS, es_SNU16_multiome_subset_1$PVT1_exp)
dev.off()
cor.test(es_SNU16_multiome_subset_1$HALLMARK_E2F_TARGETS, es_SNU16_multiome_subset_1$PVT1_exp)#cor=-0.4043984, p-value = 3.533e-14

#To test the same correlations in the whole population:
cor.test(es_SNU16_multiome$HALLMARK_MYC_TARGETS_V1, es_SNU16_multiome$PVT1_exp)#cor=-0.4630459,p-value < 2.2e-16
cor.test(es_SNU16_multiome$HALLMARK_MTORC1_SIGNALING, es_SNU16_multiome$PVT1_exp)#cor=-0.402333,p-value < 2.2e-16
cor.test(es_SNU16_multiome$HALLMARK_E2F_TARGETS, es_SNU16_multiome$PVT1_exp)#cor=-0.3910192,p-value < 2.2e-16

plot(es_SNU16_multiome_subset_1$HALLMARK_MYC_TARGETS_V1, es_SNU16_multiome_subset_1$FGFR2_exp)
plot(es_SNU16_multiome_subset_1$HALLMARK_MTORC1_SIGNALING, es_SNU16_multiome_subset_1$FGFR2_exp)
plot(es_SNU16_multiome_subset_1$HALLMARK_BILE_ACID_METABOLISM, es_SNU16_multiome_subset_1$FGFR2_exp)

