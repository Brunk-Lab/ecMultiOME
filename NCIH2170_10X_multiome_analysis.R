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


# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

# load the RNA and ATAC data
control.multiome <- Read10X_h5("/datastore/lbcfs/labs/brunk_lab/private/NCIH2170_multiome/cell_ranger_arc_results/NCIH2170_control/outs/filtered_feature_bc_matrix.h5")
fragpath <- "/datastore/lbcfs/labs/brunk_lab/private/NCIH2170_multiome/cell_ranger_arc_results/NCIH2170_control/outs/atac_fragments.tsv.gz"

# create a Seurat object containing the RNA data
NCIH2170_multiome <- CreateSeuratObject(counts = control.multiome$`Gene Expression`,assay = "RNA")

# create ATAC assay and add it to the object
NCIH2170_multiome[["ATAC"]] <- CreateChromatinAssay(counts = control.multiome$Peaks, sep = c(":", "-"), fragments = fragpath, annotation = annotation)

DefaultAssay(NCIH2170_multiome) <- 'RNA'
NCIH2170_multiome[["percent.mt"]] <- PercentageFeatureSet(NCIH2170_multiome, pattern = "^MT-")

DefaultAssay(NCIH2170_multiome) <- "ATAC"

NCIH2170_multiome <- NucleosomeSignal(NCIH2170_multiome)
NCIH2170_multiome <- TSSEnrichment(NCIH2170_multiome)

p <- VlnPlot(NCIH2170_multiome, features = c("nCount_ATAC", "nCount_RNA", "nFeature_RNA", "TSS.enrichment", "nucleosome_signal", "percent.mt"), ncol = 6,
             log = TRUE, pt.size = 0) + NoLegend()

p

NCIH2170_multiome_2 <- subset(
  x = NCIH2170_multiome,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 1000 &
    nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    nFeature_RNA > 500  &
    nucleosome_signal < 2 &
    TSS.enrichment > 1 &
    percent.mt < 20
)


# ATAC analysis
DefaultAssay(NCIH2170_multiome_2) <- "ATAC"
NCIH2170_multiome_peaks <- CallPeaks(NCIH2170_multiome_2)
NCIH2170_multiome_peaks <- keepStandardChromosomes(NCIH2170_multiome_peaks, pruning.mode = "coarse")
NCIH2170_multiome_peaks <- subsetByOverlaps(x = NCIH2170_multiome_peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
NCIH2170_macs2_counts <- FeatureMatrix(
  fragments = Fragments(NCIH2170_multiome_2),
  features = NCIH2170_multiome_peaks,
  cells = colnames(NCIH2170_multiome_2)
)
# create a new assay using the MACS2 peak set and add it to the Seurat object
NCIH2170_multiome_2[["peaks"]] <- CreateChromatinAssay(
  counts = NCIH2170_macs2_counts,
  fragments = fragpath,
  annotation = annotation
)

DefaultAssay(NCIH2170_multiome_2) <- "peaks"
NCIH2170_multiome_2 <- FindTopFeatures(NCIH2170_multiome_2, min.cutoff = 5)
NCIH2170_multiome_2 <- RunTFIDF(NCIH2170_multiome_2)
NCIH2170_multiome_2 <- RunSVD(NCIH2170_multiome_2)
NCIH2170_multiome_2 <- RunUMAP(NCIH2170_multiome_2, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

NCIH2170.gene.activities <- GeneActivity(NCIH2170_multiome_2, biotypes = NULL)# to include not only protein coding genes but also long non coding genes
NCIH2170_multiome_2[['gene_activity']] <- CreateAssayObject(counts = NCIH2170.gene.activities)
NCIH2170_multiome_2 <- NormalizeData(
  object = NCIH2170_multiome_2,
  assay = 'gene_activity',
  normalization.method = 'LogNormalize',
  scale.factor = median(NCIH2170_multiome_2$nCount_RNA)
)

#saveRDS(NCIH2170_multiome_2, "/home/yue1118/ecDNA_project_data_and_figures/NCIH2170_multiome_intermediate_object_02192025.rds")
#NCIH2170_multiome_2 <- readRDS("/home/yue1118/ecDNA_project_data_and_figures/NCIH2170_multiome_intermediate_object_02192025.rds")

# RNA analysis
DefaultAssay(NCIH2170_multiome_2) <- "RNA"
NCIH2170_multiome_2 <- SCTransform(NCIH2170_multiome_2, vars.to.regress = "percent.mt", verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_', seed.use = 7)
#saveRDS(NCIH2170_multiome_2, "/home/yue1118/ecDNA_project_data_and_figures/NCIH2170_multiome_intermediate_object_02192025.rds")

#Now run DoubletFinder to remove potential doublets from the dataset
suppressMessages(require(DoubletFinder))
#To look for the best pK value
sweep.res.list_control <- paramSweep(NCIH2170_multiome_2, PCs = 1:20, sct=TRUE)
sweep.stats_control <- summarizeSweep(sweep.res.list_control, GT = FALSE)
bcmvn_control <- find.pK(sweep.stats_control)
pK <- bcmvn_control %>% filter(BCmetric == max(BCmetric)) %>% select(pK) #pK=0.14
#Homotypic Doublet Proportion Estimate
cluster_annotation <- NCIH2170_multiome_2@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(cluster_annotation)
#To estimate the number of doublets
nExp_poi <- round(0.06*nrow(NCIH2170_multiome_2@meta.data)) #assuming there are 6% doublets
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
NCIH2170_multiome_2 <- DoubletFinder::doubletFinder(NCIH2170_multiome_2, pN = 0.25, pK = 0.14, nExp = nExp_poi.adj, PCs = 1:10, sct = TRUE)


##ChromVAR assay
DefaultAssay(NCIH2170_multiome_2) <- "ATAC"
main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
keep.peaks <- which(as.character(seqnames(granges(NCIH2170_multiome_2))) %in% main.chroms)
NCIH2170_multiome_2[["ATAC"]] <- subset(NCIH2170_multiome_2[["ATAC"]], features = rownames(NCIH2170_multiome_2[["ATAC"]])[keep.peaks])

DefaultAssay(NCIH2170_multiome_2) <- "ATAC"
pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE))
motif.matrix <- CreateMotifMatrix(features = granges(NCIH2170_multiome_2), pwm = pwm_set, genome = 'hg38', use.counts = FALSE)
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)
NCIH2170_multiome_2 <- SetAssayData(NCIH2170_multiome_2, assay = 'ATAC', slot ='motifs', new.data = motif.object)


# Note that this step can take 30-60 minutes 
NCIH2170_multiome_2 <- RunChromVAR(
  object = NCIH2170_multiome_2,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

saveRDS(NCIH2170_multiome_2, "/home/yue1118/ecDNA_project_data_and_figures/NCIH2170_multiome_doublet_correctly_removed_with_ChromVAR_02192025.rds")

NCIH2170_multiome_2 <- readRDS("/home/yue1118/ecDNA_project_data_and_figures/NCIH2170_multiome_doublet_correctly_removed_with_ChromVAR_02192025.rds")

#To only select the predicted singlets for downstream analysis
Idents(NCIH2170_multiome_2) <- "DF.classifications_0.25_0.14_477"
NCIH2170_multiome_2 <- subset(NCIH2170_multiome_2, idents = "Singlet")

NCIH2170_multiome_2$condition <- "WT"

#Codes for Extended Data Fig.5b
list_ecDNA_ATAC = c('MYC', 'ERBB2', 'CDC6', 'IKZF3', 'PVT1', 'CASC11', 'WIPF2',  'LINC00824', 'PNMT', 'PCAT1', 'GRB7', 'RAPGEFL1', 'MIEN1', 'PGAP3', 'CASC21', 'CASC3', 'ZPBP2', 'TCAP',
                    'POU5F1B')

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

i = "ERBB2"
df_sct <- as.data.frame(NCIH2170_multiome_2@assays$SCT@data[i,])
colnames(df_sct) <- 'gene'
df_sct$cells <- rownames(df_sct)
df3 <- df_sct[df_sct$gene >= quantile(as.numeric(NCIH2170_multiome_2@assays$SCT@data[i,]), .90, na.rm=TRUE), ]
df3$group <- 'higher'
df4 <- df_sct[df_sct$gene <= quantile(as.numeric(NCIH2170_multiome_2@assays$SCT@data[i,]), .10, na.rm=TRUE), ]
df4$group <- 'lower'
NCIH2170_multiome_2@meta.data$group <- lapply(colnames(NCIH2170_multiome_2), group.func)
Idents(NCIH2170_multiome_2) <- NCIH2170_multiome_2$group
DefaultAssay(NCIH2170_multiome_2) <- "SCT"
df_mRNA_markers<- FindMarkers(NCIH2170_multiome_2, ident.1 = "high", ident.2 = "low", min.pct = 0.25, logfc.threshold = 0.10, recorrect_umi=F)
DefaultAssay(NCIH2170_multiome_2) <- "gene_activity"
df_atac_markers <- FindMarkers(NCIH2170_multiome_2, ident.1 = "high", ident.2 = "low", min.pct = 0.25, logfc.threshold = 0.10, recorrect_umi=F)

View(df_mRNA_markers[df_mRNA_markers$p_val_adj < 0.05,])
volcano_plot_function(df_mRNA_markers)
ggsave('/home/yue1118/ecDNA_project_data_and_figures/FACS_FISH_FULL_figures/volcanoplot_DEG_ERBB2_mRNA_high_vs_low_cells_NCIH2170_doublet_correctly_removed_02202025.svg', width = 8, height = 6)

#Codes for Extended Data Fig.5c
View(df_atac_markers[(df_atac_markers$p_val_adj < 0.05),])
volcano_plot_function(df_atac_markers)
ggsave('/home/yue1118/ecDNA_project_data_and_figures/FACS_FISH_FULL_figures/volcanoplot_Differential_ATAC_gene_activities_ERBB2_mRNA_high_vs_low_cells_NCIH2170_doublet_correctly_removed_02202025.svg', width = 8, height = 6)


#To perform genome-wide searching to look for genes whose expression or chromatin accessibility is associated with ecDNA genes.
list_ecDNA_ATAC = c('MYC', 'ERBB2', 'CDC6', 'IKZF3', 'PVT1', 'CASC11', 'WIPF2',  'LINC00824', 'PNMT', 'PCAT1', 'GRB7', 'RAPGEFL1', 'MIEN1', 'PGAP3', 'CASC21', 'CASC3', 'ZPBP2', 'TCAP',
                    'POU5F1B')

group.func <- function(x){
  if (x %in% rownames(df3)) {
    "high"
  } else if (x %in% rownames(df4)) {
    "low"
  } else {
    "medium"
  }
}
columns= c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")
df_mRNA_markers = data.frame(matrix(nrow = 0, ncol = length(columns)))
df_atac_markers = data.frame(matrix(nrow = 0, ncol = length(columns)))
for (i in list_ecDNA_ATAC){
  df_atac <- as.data.frame(NCIH2170_multiome_2@assays$gene_activity@data[i,])
  colnames(df_atac) <- 'gene'
  df_atac$cells <- rownames(df_atac)
  df3 <- df_atac[df_atac$gene >= quantile(as.numeric(NCIH2170_multiome_2@assays$gene_activity@data[i,]), .90, na.rm=TRUE), ]
  df3$group <- 'higher'
  df4 <- df_atac[df_atac$gene <= quantile(as.numeric(NCIH2170_multiome_2@assays$gene_activity@data[i,]), .10, na.rm=TRUE), ]
  df4$group <- 'lower'
  NCIH2170_multiome_2@meta.data$group <- lapply(colnames(NCIH2170_multiome_2), group.func)
  Idents(NCIH2170_multiome_2) <- NCIH2170_multiome_2$group
  DefaultAssay(NCIH2170_multiome_2) <- "SCT"
  df_mRNA_markers_tmp <- FindMarkers(NCIH2170_multiome_2, ident.1 = "high", ident.2 = "low", min.pct = 0.25, logfc.threshold = 0.10, recorrect_umi=F)
  df_mRNA_markers_tmp$gene <- rownames(df_mRNA_markers_tmp)
  rownames(df_mRNA_markers_tmp)<-NULL
  df_mRNA_markers_tmp$ecDNA_gene <- i
  df_mRNA_markers <- rbind(df_mRNA_markers, df_mRNA_markers_tmp)
  DefaultAssay(NCIH2170_multiome_2) <- "gene_activity"
  df_atac_markers_tmp <- FindMarkers(NCIH2170_multiome_2, ident.1 = "high", ident.2 = "low", min.pct = 0.25, logfc.threshold = 0.10, recorrect_umi=F)
  df_atac_markers_tmp$gene <- rownames(df_atac_markers_tmp)
  rownames(df_atac_markers_tmp)<-NULL
  df_atac_markers_tmp$ecDNA_gene <- i
  df_atac_markers <- rbind(df_atac_markers, df_atac_markers_tmp)
}

write.csv(df_atac_markers, '/home/yue1118/ecDNA_project_data_and_figures/FACS_FISH_FULL_figures/find_atac_markers_between_atac_high_vs_low_ecDNA_genes_NCIH2170_multiome_WT_wilcoxon_rank_sum_02212025.csv')
write.csv(df_mRNA_markers, '/home/yue1118/ecDNA_project_data_and_figures/FACS_FISH_FULL_figures/find_mRNA_markers_between_atac_high_vs_low_ecDNA_genes_NCIH2170_multiome_WT_wilcoxon_rank_sum_02212025.csv')


#Codes for Extended Data Fig.5d
# Plot circular heatmap with links between ATAC and gene expressions

#From the "df_atac_markers" (ATAC-ATAC association) dataframe obtained above, we choose the most frequent genes showing significant correlation with ecDNA genes and below are all of them:
list_gene_to_plot <- c("MYC", "ERBB2", "CDC6", "IKZF3", "PVT1", "CASC11", "WIPF2", "LINC00824", "PNMT", "PCAT1", "GRB7", "RAPGEFL1", "MIEN1", "PGAP3", "CASC21",  
                       "CASC3", "ZPBP2", "TCAP", "CCAT2", "POU5F1B","LINC00976", "STARD3", "RARA", "DENND2C", "CSDE1", "TRIM33", "ZBTB7C", "MYH14", "MARK4", "RAB27B", "ATP7B", "NTSR1", "BRMS1L", "SH3RF3", "CDH26")

Idents(NCIH2170_multiome_2) <- NCIH2170_multiome_2$condition
df_atac_aggre <- AggregateExpression(NCIH2170_multiome_2, assays = 'gene_activity', slot = 'counts')
df_atac_aggre <- data.frame(df_atac_aggre)
colnames(df_atac_aggre) <- c('aggregated_value')
write.csv(df_atac_aggre, '/home/yue1118/ecDNA_project_data_and_figures/FACS_FISH_FULL_figures/NCIH2170_multiome_WT_ATAC_gene_activity_aggregated_02202025.csv')


df_scRNA_aggre <- AggregateExpression(NCIH2170_multiome_2, assays = 'SCT', slot = 'counts')
df_scRNA_aggre <- data.frame(df_scRNA_aggre)
colnames(df_scRNA_aggre) <- c('aggregated_mRNA_value')
write.csv(df_scRNA_aggre, '/home/yue1118/ecDNA_project_data_and_figures/FACS_FISH_FULL_figures/NCIH2170_multiome_WT_gene_expression_aggregated_02202025.csv')


group.func <- function(x){
  if (x %in% rownames(df3)) {
    "high"
  } else if (x %in% rownames(df4)) {
    "low"
  } else {
    "medium"
  }
}
columns= c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")
df_mRNA_markers_circular_heatmap = data.frame(matrix(nrow = 0, ncol = length(columns)))
for (i in list_gene_to_plot){
  DefaultAssay(NCIH2170_multiome_2) <- "gene_activity"
  if (i %in% rownames(NCIH2170_multiome_2)) {
    df_atac <- as.data.frame(NCIH2170_multiome_2@assays$gene_activity@data[i,])
    colnames(df_atac) <- 'gene'
    df_atac$cells <- rownames(df_atac)
    df3 <- df_atac[df_atac$gene >= quantile(as.numeric(NCIH2170_multiome_2@assays$gene_activity@data[i,]), .90, na.rm=TRUE), ]
    df3$group <- 'higher'
    df4 <- df_atac[df_atac$gene <= quantile(as.numeric(NCIH2170_multiome_2@assays$gene_activity@data[i,]), .10, na.rm=TRUE), ]
    df4$group <- 'lower'
    NCIH2170_multiome_2@meta.data$group <- lapply(colnames(NCIH2170_multiome_2), group.func)
    Idents(NCIH2170_multiome_2) <- NCIH2170_multiome_2$group
    DefaultAssay(NCIH2170_multiome_2) <- "SCT"
    df_mRNA_markers_tmp <- FindMarkers(NCIH2170_multiome_2, ident.1 = "high", ident.2 = "low", min.pct = 0.25, logfc.threshold = 0.10, recorrect_umi=F)
    df_mRNA_markers_tmp$gene <- rownames(df_mRNA_markers_tmp)
    rownames(df_mRNA_markers_tmp)<-NULL
    df_mRNA_markers_tmp$potential_ecDNA_gene <- i
    df_mRNA_markers_circular_heatmap <- rbind(df_mRNA_markers_circular_heatmap, df_mRNA_markers_tmp)
  }
}
df_mRNA_markers_circular_heatmap <- subset(df_mRNA_markers_circular_heatmap, p_val_adj < 0.05 & gene %in% list_gene_to_plot)


df_scge_bind <- read.csv('/home/yue1118/ecDNA_project_data_and_figures/FACS_FISH_FULL_figures/NCIH2170_multiome_WT_aggregated_single_cell_expression_bind_02202025.csv')
df_scge_bind$aggregated_expression <- as.numeric(df_scge_bind$aggregated_expression)
df_scge_bind$aggregated_expression <- log2(df_scge_bind$aggregated_expression+1)
rownames(df_scge_bind) <- df_scge_bind$X
df_scge_bind <- df_scge_bind[,-1]
split_scge = factor(as.character(df_scge_bind[,2]), levels = c('chr1','chr2', 'chr8', 'chr13', 'chr14', 'chr17','chr18','chr19','chr20'))
col_scge = colorRamp2(c(0, 0, 21), c("white", "white", "orange2"))

df_ge_bind <- read.csv('/home/yue1118/ecDNA_project_data_and_figures/FACS_FISH_FULL_figures/NCIH2170_multiome_WT_GE_bind_02202025.csv')
df_ge_bind$gene_expression <- as.numeric(df_ge_bind$gene_expression)
rownames(df_ge_bind) <- df_ge_bind$X
df_ge_bind <- df_ge_bind[,-1]
split_ge = factor(as.character(df_ge_bind[,2]), levels = c('chr1','chr2', 'chr8', 'chr13', 'chr14', 'chr17','chr18','chr19','chr20'))
col_ge = colorRamp2(c(0, 0, 13), c("white", "white", "brown"))


df_cn_bind <- read.csv('/home/yue1118/ecDNA_project_data_and_figures/FACS_FISH_FULL_figures/NCIH2170_multiome_WT_CNV_bind_02202025.csv')
df_cn_bind$copy_number <- as.numeric(df_cn_bind$copy_number)
rownames(df_cn_bind) <- df_cn_bind$X
df_cn_bind <- df_cn_bind[,-1]
split_cn = factor(as.character(df_cn_bind[,2]), levels = c('chr1','chr2', 'chr8', 'chr13', 'chr14', 'chr17','chr18','chr19','chr20'))
col_cn = colorRamp2(c(0, 0, 7), c("white", "white", "darkgreen"))
#set_track_gap(gap = 0.02)


df_atac_bind <- read.csv('/home/yue1118/ecDNA_project_data_and_figures/FACS_FISH_FULL_figures/NCIH2170_multiome_WT_ATAC_bind_02202025.csv')
df_atac_bind$atac_value <- as.numeric(df_atac_bind$atac_value)
df_atac_bind$atac_value <- log2(df_atac_bind$atac_value)
rownames(df_atac_bind) <- df_atac_bind$X
df_atac_bind <- df_atac_bind[,-1]
split_atac = factor(as.character(df_atac_bind[,2]), levels = c('chr1','chr2', 'chr8', 'chr13', 'chr14', 'chr17','chr18','chr19','chr20'))
col_atac = colorRamp2(c(0, 0, 21), c("white", "white", "blue4"))

df_link <- data.frame(matrix(ncol = 0, nrow = nrow(df_mRNA_markers_circular_heatmap)))
df_link$from_index <- mapvalues(df_mRNA_markers_circular_heatmap$potential_ecDNA_gene, 
                                from= as.vector(rownames(df_ge_bind)), 
                                to=as.vector(seq_len(nrow(df_ge_bind))))
df_link$to_index <- mapvalues(df_mRNA_markers_circular_heatmap$gene, 
                              from= as.vector(rownames(df_ge_bind)), 
                              to=as.vector(seq_len(nrow(df_ge_bind))))
df_link$from_index <- as.numeric(df_link$from_index)
df_link$to_index <- as.numeric(df_link$to_index)


colors30 = c("#466791","#60bf37","#953ada","#4fbe6c","#ce49d3","#a7b43d","#5a51dc","#d49f36","#552095","#507f2d","#db37aa","#84b67c","#a06fda","#df462a","#5b83db","#c76c2d","#4f49a3","#82702d","#dd6bbb","#334c22","#d83979","#55baad","#dc4555","#62aad3","#8c3025","#417d61","#862977","#bba672","#403367","#da8a6d")#,"#a79cd4","#71482c","#c689d0","#6b2940","#d593a7","#895c8b","#bd5975")

df_link$color <- mapvalues(df_link$from_index, 
                           from= unique(df_link$from_index), 
                           to=colors30)

View(df_link)
View(df_mRNA_markers_circular_heatmap)
# To plot the heatmap and save it as a svg file

svglite("/home/yue1118/ecDNA_project_data_and_figures/FACS_FISH_FULL_figures/circle_heatmap_NCIH2170_multiome_WT_02202025.svg", width = 10, height = 10)
circos.par(gap.degree = 6)  # Increase gap between sectors if needed
circos.heatmap(as.matrix(df_ge_bind[1]), split = split_ge, col = col_ge, na.col = "grey", cluster = FALSE, track.height = 0.01, rownames.side = "outside")
circos.heatmap(as.matrix(df_scge_bind[1]), split = split_scge, col = col_scge, na.col = "grey", cluster = FALSE, track.height = 0.02)
circos.heatmap(as.matrix(df_cn_bind[1]), split = split_cn, col = col_cn, na.col = "grey", cluster = FALSE, track.height = 0.02)
circos.heatmap(as.matrix(df_atac_bind[1]), split = split_atac, col = col_atac, na.col = "grey", cluster = FALSE, track.height = 0.02)
for(i in seq_len(nrow(df_link))) {
  circos.heatmap.link(df_link$from_index[i],
                      df_link$to_index[i],
                      col = df_link$color[i])
}
dev.off()
circos.clear()

#To plot the heatmap legends
lgd_expr = Legend(title = "Expression", col_fun = col_ge)
lgd_cn = Legend(title = "Copy Number", col_fun = col_cn)
lgd_atac = Legend(title = "ATAC signal", col_fun = col_atac)
lgd_scge = Legend(title = 'aggregated scRNA signal', col_fun = col_scge)
svglite("/home/yue1118/ecDNA_project_data_and_figures/FACS_FISH_FULL_figures/circle_heatmap_legends_NCIH2170_ATAC_02202025.svg", width = 10, height = 10)
draw(lgd_atac, just = "center")
dev.off()

svglite("/home/yue1118/ecDNA_project_data_and_figures/FACS_FISH_FULL_figures/circle_heatmap_legends_NCIH2170_scRNA_02202025.svg", width = 10, height = 10)
draw(lgd_scge, just = "center")
dev.off()

svglite("/home/yue1118/ecDNA_project_data_and_figures/FACS_FISH_FULL_figures/circle_heatmap_legends_NCIH2170_bulk_exp_02202025.svg", width = 10, height = 10)
draw(lgd_expr, just = "center")
dev.off()

svglite("/home/yue1118/ecDNA_project_data_and_figures/FACS_FISH_FULL_figures/circle_heatmap_legends_NCIH2170_bulk_CN_02202025.svg", width = 10, height = 10)
draw(lgd_cn, just = "center")
dev.off()
