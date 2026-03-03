.libPaths("/home/yue1118/R/x86_64-pc-linux-gnu-library/4.2")


library(Seurat)
library(Signac)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(MACSr)
library(dplyr)
library(ggplot2)
library(EnsDb.Hsapiens.v75)
library(Rsamtools)

#To compress the fragment.tsv.gz file first
# Step 1: bgzip
#bgzip("/datastore/lbcfs/labs/brunk_lab/public/COLO320_multiome_hg19/DM_fragments/DM_fragments_chr.sorted.tsv", overwrite = TRUE)

# Step 2: take the output file and do tabix index
#indexTabix("/datastore/lbcfs/labs/brunk_lab/public/COLO320_multiome_hg19/DM_fragments/DM_fragments_chr.sorted.tsv.bgz", format = "bed")
#Step 3: rename the files to "tsv.gz" and "tsv.gz.tbi" instead of "bgz"

# get gene annotations for hg19
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
genome(annotation) <- "hg19"

#Now create the ChromatinAssayObject
genome <- BSgenome.Hsapiens.UCSC.hg19
seqinfo <- seqinfo(genome)
tiles <- tileGenome(seqlengths = seqinfo, tilewidth = 5000, cut.last.tile.in.chrom = TRUE)

frags <- CreateFragmentObject(path = "/datastore/lbcfs/labs/brunk_lab/public/COLO320_multiome_hg19/DM_fragments/DM_fragments_chr.sorted.tsv.gz")
counts <- FeatureMatrix(fragments = frags,features = tiles,cells = NULL)
dim(counts)
DM_chrom_assay <-  CreateChromatinAssay(counts = counts, fragments = frags, annotation = annotation, genome = "hg19")

#saveRDS(DM_chrom_assay, file = "/datastore/lbcfs/labs/brunk_lab/public/COLO320_multiome_hg19/DM_fragments/DM_chrom_assay_temporary.rds")
DM_chrom_assay <- readRDS("/datastore/lbcfs/labs/brunk_lab/public/COLO320_multiome_hg19/DM_fragments/DM_chrom_assay_temporary.rds")

DM_chrom_obj <- CreateSeuratObject(counts = DM_chrom_assay, assay = "ATAC")
DefaultAssay(DM_chrom_obj) <- "ATAC"
DM_chrom_obj <- NucleosomeSignal(DM_chrom_obj)
DM_chrom_obj <- TSSEnrichment(DM_chrom_obj)

p <- VlnPlot(DM_chrom_obj, features = c("nCount_ATAC", "TSS.enrichment", "nucleosome_signal"), ncol = 3, log = TRUE, pt.size = 0) + NoLegend()

p

DM_chrom_obj <- subset(x = DM_chrom_obj, subset = nCount_ATAC < 7e4 & nCount_ATAC > 1000 & nucleosome_signal < 2 & TSS.enrichment > 1)

dim(DM_chrom_obj)


# calling peaks from ATAC
DefaultAssay(DM_chrom_obj) <- "ATAC"
DM_chrom_peaks <- CallPeaks(DM_chrom_obj, macs2.path = "~/miniconda3/envs/macs3_env/bin/macs3")
DM_chrom_peaks <- keepStandardChromosomes(DM_chrom_peaks, pruning.mode = "coarse")
DM_chrom_peaks <- subsetByOverlaps(x = DM_chrom_peaks, ranges = blacklist_hg19, invert = TRUE)

# quantify counts in each peak
DM_chrom_peak_counts <- FeatureMatrix(
  fragments = Fragments(DM_chrom_obj),
  features = DM_chrom_peaks,
  cells = colnames(DM_chrom_obj)
)
# create a new assay using the MACS2 peak set and add it to the Seurat object
DM_chrom_obj[["peaks"]] <- CreateChromatinAssay(
  counts = DM_chrom_peak_counts,
  fragments = "/datastore/lbcfs/labs/brunk_lab/public/COLO320_multiome_hg19/DM_fragments/DM_fragments_chr.sorted.tsv.gz",
  annotation = annotation
)

DefaultAssay(DM_chrom_obj) <- "peaks"
DM_chrom_obj <- FindTopFeatures(DM_chrom_obj, min.cutoff = 5)
DM_chrom_obj <- RunTFIDF(DM_chrom_obj)
DM_chrom_obj <- RunSVD(DM_chrom_obj)
DM_chrom_obj <- RunUMAP(DM_chrom_obj, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")


DM_chrom_obj.gene.activities <- GeneActivity(DM_chrom_obj, biotypes = NULL)# to include not only protein coding genes but also long non coding genes
DM_chrom_obj[['gene_activity']] <- CreateAssayObject(counts = DM_chrom_obj.gene.activities)
DM_chrom_obj <- NormalizeData(
  object = DM_chrom_obj,
  assay = 'gene_activity',
  normalization.method = 'LogNormalize',
  scale.factor = 10000
)

saveRDS(DM_chrom_obj, "/datastore/lbcfs/labs/brunk_lab/public/COLO320_multiome_hg19/DM_fragments/COLO320DM_scATAC_object_05162025.rds")



#Now to examine the pseudobulk gene_activity of ecDNA genes:
df_gene_activity <- as.data.frame(AggregateExpression(DM_chrom_obj, assays = "gene_activity"))
View(df_gene_activity)
#We can see that ecDNA genes have the most chromatin accessibility.


#Now integrate scRNA data into the DM_chrom_obj seurat object:
COLO320_data <- Read10X(data.dir = "/datastore/lbcfs/labs/brunk_lab/public/COLO320_multiome_hg19/scRNA/")
COLO320_RNA_seurat_obj <- CreateSeuratObject(counts = COLO320_data, project = "COLO320")
COLO320_RNA_seurat_obj$cell_barcode <- colnames(COLO320_RNA_seurat_obj)
COLO320_RNA_seurat_obj$cleaned_barcodes <- sub(".*_rep[1-7]_", "", COLO320_RNA_seurat_obj$cell_barcode)
COLO320_RNA_seurat_obj_subset <- subset(COLO320_RNA_seurat_obj, subset = orig.ident == "COLO320DM")

cell_names <- COLO320_RNA_seurat_obj_subset$cleaned_barcodes
cells_to_keep <- cell_names[!duplicated(cell_names)]
COLO320_RNA_seurat_obj_subset <- subset(COLO320_RNA_seurat_obj, subset = cleaned_barcodes %in% cells_to_keep)

new_barcodes <- as.character(COLO320_RNA_seurat_obj_subset$cleaned_barcodes)
length(new_barcodes) == length(Cells(COLO320_RNA_seurat_obj_subset))
new_barcodes_unique <- make.unique(new_barcodes)
COLO320_RNA_seurat_obj_subset <- RenameCells(COLO320_RNA_seurat_obj_subset, new.names = new_barcodes_unique)

cells_to_keep <- as.character(intersect(COLO320_RNA_seurat_obj_subset$cleaned_barcodes, colnames(DM_chrom_obj)))
DM_chrom_obj_subset <- subset(DM_chrom_obj, cells = cells_to_keep)
COLO320_RNA_seurat_obj_subset <- subset(COLO320_RNA_seurat_obj_subset, cells = cells_to_keep)

table(colnames(COLO320_RNA_seurat_obj_subset) == colnames(DM_chrom_obj_subset)) #TRUE

DM_chrom_obj_subset[['RNA']] <- COLO320_RNA_seurat_obj_subset[['RNA']]

DefaultAssay(DM_chrom_obj_subset) <- 'RNA'
DM_chrom_obj_subset[["percent.mt"]] <- PercentageFeatureSet(DM_chrom_obj_subset, pattern = "^MT-")


p <- VlnPlot(DM_chrom_obj_subset, features = c("nCount_ATAC", "nCount_RNA", "nFeature_RNA", "TSS.enrichment", "nucleosome_signal", "percent.mt"), ncol = 6,
             log = TRUE, pt.size = 0) + NoLegend()

p

DM_chrom_obj_subset <- subset(x = DM_chrom_obj_subset, subset = nCount_RNA < 25000 & nCount_RNA > 1000 & nFeature_RNA > 500  & percent.mt < 10)

#saveRDS(DM_chrom_obj_subset, "/datastore/lbcfs/labs/brunk_lab/public/COLO320_multiome_hg19/DM_fragments/COLO320DM_scATAC_scRNA_multiome_reintegrated_object_05172025.rds")
#DM_chrom_obj <- readRDS("/datastore/lbcfs/labs/brunk_lab/public/COLO320_multiome_hg19/DM_fragments/COLO320DM_scATAC_scRNA_multiome_reintegrated_object_05172025.rds")

#DefaultAssay(DM_chrom_obj) <- "RNA"
#DM_chrom_obj <- SCTransform(DM_chrom_obj, vars.to.regress = "percent.mt", verbose = FALSE)
#saveRDS(DM_chrom_obj, "/datastore/lbcfs/labs/brunk_lab/public/COLO320_multiome_hg19/DM_fragments/COLO320DM_scATAC_scRNA_multiome_reintegrated_object_05202025.rds")


DM_chrom_obj <- readRDS("/datastore/lbcfs/labs/brunk_lab/public/COLO320_multiome_hg19/DM_fragments/COLO320DM_scATAC_scRNA_multiome_reintegrated_object_05202025.rds")

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

group.func <- function(x){
  if (x %in% rownames(df3)) {
    "high"
  } else if (x %in% rownames(df4)) {
    "low"
  } else {
    "medium"
  }
}


i = "MYC"
df_sct <- as.data.frame(DM_chrom_obj@assays$gene_activity@data[i,])
colnames(df_sct) <- 'gene'
df_sct$cells <- rownames(df_sct)
df3 <- df_sct[df_sct$gene >= quantile(as.numeric(DM_chrom_obj@assays$gene_activity@data[i,]), .90, na.rm=TRUE), ]
df3$group <- 'higher'
df4 <- df_sct[df_sct$gene <= quantile(as.numeric(DM_chrom_obj@assays$gene_activity@data[i,]), .10, na.rm=TRUE), ]
df4$group <- 'lower'
DM_chrom_obj@meta.data$group <- lapply(colnames(DM_chrom_obj), group.func)
Idents(DM_chrom_obj) <- DM_chrom_obj$group

DefaultAssay(DM_chrom_obj) <- "gene_activity"
df_atac_markers <- FindMarkers(DM_chrom_obj, ident.1 = "high", ident.2 = "low", min.pct = 0.25, logfc.threshold = 0.10, recorrect_umi=F)

DefaultAssay(DM_chrom_obj) <- "SCT"
df_sct_markers <- FindMarkers(DM_chrom_obj, ident.1 = "high", ident.2 = "low", min.pct = 0.25, logfc.threshold = 0.10, recorrect_umi=F)


list_ecDNA_genes <- c("MYC", "PVT1", "LINC02912", "PCAT1", "CASC19", "CASC11", "MIR1204", "POU5F1B", "CASC8","PRNCR1", "CCAT2",
                      "CDX2", "FAM84B", "RNU6-869P", "PLUT", "LINC00543", "PDX1", "ATP5EP2","CCAT1", "CASC21","TMEM75")

View(df_atac_markers[(df_atac_markers$p_val_adj < 0.05) & (rownames(df_atac_markers) %in% list_ecDNA_genes),])

View(df_atac_markers[(df_atac_markers$p_val_adj < 0.05),])

View(df_sct_markers[(df_sct_markers$p_val_adj < 0.05) & (rownames(df_sct_markers) %in% list_ecDNA_genes),])
View(df_sct_markers[df_sct_markers$p_val_adj < 0.05,])
View(df_sct_markers)

volcano_plot_function(df_atac_markers)
ggsave("/home/yue1118/ecDNA_project_data_and_figures/FACS_FISH_FULL_figures/COLO320DM_differential_peaks_MYC_gene_activity_high_vs_low_05292025.svg", width = 8, height = 6)


volcano_plot_function(df_sct_markers)
ggsave("/home/yue1118/ecDNA_project_data_and_figures/FACS_FISH_FULL_figures/COLO320DM_differential_expressed_gene_between_MYC_gene_activity_high_vs_low_05292025.svg", width = 8, height = 6)





###To plot pseudobulk gene expression of all the genes in the genome and highlight PVT1
DM_chrom_obj <- readRDS("/datastore/lbcfs/labs/brunk_lab/public/COLO320_multiome_hg19/DM_fragments/COLO320DM_scATAC_scRNA_multiome_reintegrated_object_05202025.rds")
Idents(DM_chrom_obj) <- 'WT'
df_sc_aggre <- AggregateExpression(DM_chrom_obj, assays = 'SCT', slot = 'data')
df_sc_aggre <- data.frame(df_sc_aggre)
colnames(df_sc_aggre) <- c('aggregated_value')
df_sc_aggre$log2_value <- log2(df_sc_aggre$aggregated_value)
df_sc_aggre <- subset(df_sc_aggre, log2_value >0)
df_sc_aggre <- df_sc_aggre[order(df_sc_aggre$log2_value,decreasing=TRUE),]
df_sc_aggre$gene <- rownames(df_sc_aggre)
df_sc_aggre <- df_sc_aggre[1:100,]

df_sc_aggre$highlight <- NA
for (i in (1: nrow(df_sc_aggre))) {
  if (df_sc_aggre[i, 'gene'] %in% c("PVT1")){
    df_sc_aggre[i, 'highlight'] <- "yes"
  } else {df_sc_aggre[i, 'highlight'] <- "no"}
}

df_sc_aggre$gene <- factor(df_sc_aggre$gene, levels = rownames(df_sc_aggre))

ggplot(df_sc_aggre, aes( x = gene, y = log2_value, fill = highlight ) ) +
  geom_bar( stat = "identity",width = 1, position = "stack") +
  scale_fill_manual( values = c( "yes"="tomato", "no"="beige" )) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave('/home/yue1118/Bioskryb_data_analysis/10x_multiome_gene_expression_across_the_genome_COLO320HSR_muiltiome_WT_06182025.svg', height = 5, width = 5)


#To explore different ATAC peaks associated with PVT1 mRNA

group.func <- function(x){
  if (x %in% rownames(df3)) {
    "high"
  } else if (x %in% rownames(df4)) {
    "low"
  } else {
    "medium"
  }
}

i = "PVT1"
df_sct <- as.data.frame(DM_chrom_obj@assays$SCT@data[i,])
colnames(df_sct) <- 'gene'
df_sct$cells <- rownames(df_sct)
df3 <- df_sct[df_sct$gene >= quantile(as.numeric(DM_chrom_obj@assays$SCT@data[i,]), .90, na.rm=TRUE), ]
df3$group <- 'higher'
df4 <- df_sct[df_sct$gene <= quantile(as.numeric(DM_chrom_obj@assays$SCT@data[i,]), .10, na.rm=TRUE), ]
df4$group <- 'lower'
DM_chrom_obj@meta.data$group <- lapply(colnames(DM_chrom_obj), group.func)
Idents(DM_chrom_obj) <- DM_chrom_obj$group
DefaultAssay(DM_chrom_obj) <- "SCT"
df_mRNA_markers<- FindMarkers(DM_chrom_obj, ident.1 = "high", ident.2 = "low", min.pct = 0.25, logfc.threshold = 0.10, recorrect_umi=F)
DefaultAssay(DM_chrom_obj) <- "gene_activity"
df_atac_markers <- FindMarkers(DM_chrom_obj, ident.1 = "high", ident.2 = "low", min.pct = 0.25, logfc.threshold = 0.10, recorrect_umi=F)

View(df_atac_markers[(df_atac_markers$p_val_adj < 0.05) & (df_atac_markers$avg_log2FC >0.2),])


df_atac_marker_test <- df_atac_markers[(df_atac_markers$p_val_adj < 0.05) & (df_atac_markers$avg_log2FC >0.2),]


write.csv(df_atac_marker_test, "/home/yue1118/TNBC_scRNA_analysis_workflow/COLO320DM_ATAC_associated_with_PVT1_mRNA_07292025.csv")

######################To compute correlation of each peak accessibility with PVT1 mRNA across single cells with subsampling######################

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)  #COLO320DM sample is aligned using hg19
library(org.Hs.eg.db)
library(Matrix)

txdb_hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene

set.seed(7)
COLO320DM_selected_cells <- sample(colnames(DM_chrom_obj), 400)
COLO320DM_seurat_downsampled <- subset(DM_chrom_obj, cells = COLO320DM_selected_cells)

#To calculate PVT1 mRNA-peak correlation
COLO320DM_pvt1_expr <- GetAssayData(COLO320DM_seurat_downsampled, assay = "SCT", slot = "data")["PVT1", ]
COLO320DM_peak_mat <- GetAssayData(COLO320DM_seurat_downsampled, assay = "peaks", slot = "data")


# Ensure the cell order matches between assays
COLO320DM_peak_mat <- COLO320DM_peak_mat[, colnames(COLO320DM_seurat_downsampled)]
COLO320DM_pvt1_expr <- COLO320DM_pvt1_expr[colnames(COLO320DM_seurat_downsampled)]

# Compute correlations
COLO320DM_cor_vals <- apply(COLO320DM_peak_mat, 1, function(x) cor(x, COLO320DM_pvt1_expr, use = "complete.obs", method = "pearson"))

# Optionally compute p-values
COLO320DM_pvals <- apply(COLO320DM_peak_mat, 1, function(x) cor.test(x, COLO320DM_pvt1_expr)$p.value)

# Put in a dataframe
COLO320DM_cor_df <- data.frame(
  peak = rownames(COLO320DM_peak_mat),
  correlation = COLO320DM_cor_vals,
  p_value = COLO320DM_pvals
)

COLO320DM_cor_df$padj <- p.adjust(COLO320DM_cor_df$p_value, method = "BH")
COLO320DM_cor_df <- COLO320DM_cor_df[COLO320DM_cor_df$padj < 0.05, ]
COLO320DM_cor_df <- na.omit(COLO320DM_cor_df)
COLO320DM_cor_df <- COLO320DM_cor_df[order(-COLO320DM_cor_df$correlation),]

View(COLO320DM_cor_df)


peaks_split_COLO320DM <- do.call(rbind, strsplit(COLO320DM_cor_df$peak, "[-]"))
colnames(peaks_split_COLO320DM) <- c("chr", "start", "end")
peaks_df_COLO320DM <- data.frame( chr = peaks_split_COLO320DM[, 1],start = as.numeric(peaks_split_COLO320DM[, 2]),end = as.numeric(peaks_split_COLO320DM[, 3]))
granges_obj_COLO320DM <- GRanges(seqnames = peaks_df_COLO320DM$chr,ranges = IRanges(start = peaks_df_COLO320DM$start, end = peaks_df_COLO320DM$end))


peakAnno_COLO320DM_PVT1_isoform <- annotatePeak(granges_obj_COLO320DM, tssRegion = c(-3000, 3000), TxDb = txdb_hg19, annoDb = "org.Hs.eg.db") 
df_peakAnno_COLO320DM_PVT1_isoform <- as.data.frame(peakAnno_COLO320DM_PVT1_isoform)

df_peakAnno_COLO320DM_PVT1_isoform$peak <- paste(df_peakAnno_COLO320DM_PVT1_isoform$seqnames, df_peakAnno_COLO320DM_PVT1_isoform$start, df_peakAnno_COLO320DM_PVT1_isoform$end, sep = "-")
df_peakAnno_COLO320DM_PVT1_isoform <- merge(df_peakAnno_COLO320DM_PVT1_isoform, COLO320DM_cor_df, by.x="peak", by.y= "peak", all.x=TRUE)

write.csv(df_peakAnno_COLO320DM_PVT1_isoform, "/home/yue1118/Bioskryb_data_analysis/peaks_associated_with_PVT1_mRNA_in_COLO320DM_08152025.csv")

View(df_peakAnno_COLO320DM_PVT1_isoform)
View(df_peakAnno_COLO320DM_PVT1_isoform[grepl("Promoter",df_peakAnno_COLO320DM_PVT1_isoform$annotation),])







#Now we perform differential expressed genes/peaks using inferred ecDNA counts:

#Now import inferred ecDNA counts for rep7:
df_ecDNA <- read.csv("/home/yue1118/epiAneufinder/COLO320DM_rep7_inferred_ecDNA_counts_from_scATAC_05142025.csv")
df_ecDNA$cell_barcode <- sub("cell-", "", df_ecDNA$X)
View(df_ecDNA)


DM_chrom_obj_subset <- subset(DM_chrom_obj_subset, cells = df_ecDNA$cell_barcode)
df_ecDNA <- df_ecDNA[df_ecDNA$cell_barcode %in% colnames(DM_chrom_obj_subset),]

table(df_ecDNA$cell_barcode == colnames(DM_chrom_obj_subset))

DM_chrom_obj_subset$ecDNA_count <- df_ecDNA$MYC_true_CN


# Extract the values from the Seurat object metadata
values <- DM_chrom_obj_subset@meta.data[["ecDNA_count"]]

# Compute the quantiles for 30% and 70%
q30 <- quantile(values, 0.30, na.rm = TRUE)
q70 <- quantile(values, 0.70, na.rm = TRUE)

# Assign "high", "med", or "low" based on the quantiles
DM_chrom_obj_subset@meta.data$ecDNA_bin <- ifelse(
  values >= q70, "high",
  ifelse(values <= q30, "low", "med")
)

Idents(DM_chrom_obj_subset) <- "ecDNA_bin"
DefaultAssay(DM_chrom_obj_subset) <- "gene_activity"
df_atac_markers_ecDNA <- FindMarkers(DM_chrom_obj_subset, ident.1 = "high", ident.2 = "low", min.pct = 0.25, logfc.threshold = 0.10, recorrect_umi=F)

View(df_atac_markers_ecDNA[df_atac_markers_ecDNA$p_val_adj < 0.05,])


Idents(DM_chrom_obj_subset) <- "ecDNA_bin"
DefaultAssay(DM_chrom_obj_subset) <- "RNA"
df_RNA_markers_ecDNA <- FindMarkers(DM_chrom_obj_subset, ident.1 = "high", ident.2 = "low", min.pct = 0.25, logfc.threshold = 0.10, recorrect_umi=F)
View(df_RNA_markers_ecDNA[df_RNA_markers_ecDNA$p_val_adj < 0.05,])





















View(DM_chrom_obj_subset@meta.data)

COLO320_RNA_seurat_obj <- RenameCells(COLO320_RNA_seurat_obj, new.names = new_barcodes)
COLO320_RNA_seurat_obj_subset <- subset(COLO320_RNA_seurat_obj, cells = cells_to_keep)
dim(DM_chrom_obj_subset)

View(COLO320_RNA_seurat_obj@meta.data)



#To test if colnames (cell names) in RNA object are the same as teh colnames (cell names) in ATAC object
length(intersect(colnames(COLO320_ATAC_seurat_obj), colnames(COLO320_RNA_seurat_obj)))
setequal(colnames(COLO320_ATAC_seurat_obj), colnames(COLO320_RNA_seurat_obj)) 


df_ecDNA <- read.csv("/home/yue1118/epiAneufinder/COLO320DM_rep7_inferred_ecDNA_counts_from_scATAC_05142025.csv")
df_ecDNA$cell_barcode <- sub(".*cell-", "", df_ecDNA$X)
df_ecDNA$cell_barcode <- paste0("COLO320DM_5K_rep7_", df_ecDNA$cell_barcode)
View(df_ecDNA)

length(intersect(df_ecDNA$cell_barcode, COLO320_RNA_seurat_obj@meta.data$cell_barcode))

length(intersect(df_ecDNA$cell_barcode, COLO320_ATAC_seurat_obj@meta.data$cell_barcode))

#########COLO320 RNA#########
COLO320_RNA_seurat_obj_rep7 <- subset(COLO320_RNA_seurat_obj, cells = df_ecDNA$cell_barcode)
df_ecDNA <- df_ecDNA[df_ecDNA$cell_barcode %in% colnames(COLO320_RNA_seurat_obj_rep7),]

table(colnames(COLO320_RNA_seurat_obj_rep7) == df_ecDNA$cell_barcode)

COLO320_RNA_seurat_obj_rep7$ecDNA_count <- df_ecDNA$MYC_true_CN
DefaultAssay(COLO320_RNA_seurat_obj_rep7) <- "RNA"

# Extract the values from the Seurat object metadata
values <- as.vector(COLO320_RNA_seurat_obj_rep7$ecDNA_count)

# Compute the quantiles for 30% and 70%
q30 <- quantile(values, 0.10, na.rm = TRUE)
q70 <- quantile(values, 0.90, na.rm = TRUE)

# Assign "high", "med", or "low" based on the quantiles
COLO320_RNA_seurat_obj_rep7@meta.data$ecDNA_bin <- ifelse(
  values >= q70, "high",
  ifelse(values <= q30, "low", "med")
)

Idents(COLO320_RNA_seurat_obj_rep7) <- "ecDNA_bin"
df_ecDNA_markers<- FindMarkers(COLO320_RNA_seurat_obj_rep7, ident.1 = "high", ident.2 = "low", min.pct = 0.25, logfc.threshold = 0.10, recorrect_umi=F)
View(df_ecDNA_markers[df_ecDNA_markers$p_val_adj < 0.05,])


#########COLO320 ATAC#########
COLO320_ATAC_seurat_obj_rep7 <- subset(COLO320_ATAC_seurat_obj, cells = df_ecDNA$cell_barcode)
df_ecDNA <- df_ecDNA[df_ecDNA$cell_barcode %in% colnames(COLO320_ATAC_seurat_obj_rep7),]

table(colnames(COLO320_ATAC_seurat_obj_rep7) == df_ecDNA$cell_barcode)

COLO320_ATAC_seurat_obj_rep7$ecDNA_count <- df_ecDNA$MYC_true_CN
DefaultAssay(COLO320_ATAC_seurat_obj_rep7) <- "RNA"

# Extract the values from the Seurat object metadata
values <- as.vector(COLO320_ATAC_seurat_obj_rep7$ecDNA_count)

# Compute the quantiles for 30% and 70%
q30 <- quantile(values, 0.10, na.rm = TRUE)
q70 <- quantile(values, 0.90, na.rm = TRUE)

# Assign "high", "med", or "low" based on the quantiles
COLO320_ATAC_seurat_obj_rep7@meta.data$ecDNA_bin <- ifelse(
  values >= q70, "high",
  ifelse(values <= q30, "low", "med")
)

Idents(COLO320_ATAC_seurat_obj_rep7) <- "ecDNA_bin"
df_ecDNA_ATAC_markers<- FindMarkers(COLO320_ATAC_seurat_obj_rep7, ident.1 = "high", ident.2 = "low", test.use = "fisher")
View(df_ecDNA_ATAC_markers[df_ecDNA_ATAC_markers$p_val_adj < 0.05,])


svg('/home/yue1118/Bioskryb_data_analysis/Bioskryb_plots_for_FACS_FISH_paper/PVT1_vs_MYC_ATAC_correlation_10X_multiome_COLO320DM_08082025.svg',width = 6, height = 6)
plot(DM_chrom_obj@assays$gene_activity@data['PVT1',], DM_chrom_obj@assays$gene_activity@data['MYC',])
dev.off()

cor.test(DM_chrom_obj@assays$gene_activity@data['PVT1',], DM_chrom_obj@assays$gene_activity@data['MYC',]) #0.4607193; p-value < 2.2e-16

dim(DM_chrom_obj@meta.data) #12263 cells
