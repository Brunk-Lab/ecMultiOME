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
bgzip("/datastore/lbcfs/labs/brunk_lab/public/COLO320_multiome_hg19/HSR_fragments/HSR_fragments_chr.sorted.tsv", overwrite = TRUE)

# Step 2: take the output file and do tabix index
indexTabix("/datastore/lbcfs/labs/brunk_lab/public/COLO320_multiome_hg19/HSR_fragments/HSR_fragments_chr.sorted.tsv.bgz", format = "bed")
#Step 3: rename the files to "tsv.gz" and "tsv.gz.tbi" instead of "bgz"


# get gene annotations for hg19
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
genome(annotation) <- "hg19"

#Now create the ChromatinAssayObject
genome <- BSgenome.Hsapiens.UCSC.hg19
seqinfo <- seqinfo(genome)
tiles <- tileGenome(seqlengths = seqinfo, tilewidth = 5000, cut.last.tile.in.chrom = TRUE)

frags <- CreateFragmentObject(path = "/datastore/lbcfs/labs/brunk_lab/public/COLO320_multiome_hg19/HSR_fragments/HSR_fragments_chr.sorted.tsv.gz")
counts <- FeatureMatrix(fragments = frags,features = tiles,cells = NULL)
dim(counts)
HSR_chrom_assay <-  CreateChromatinAssay(counts = counts, fragments = frags, annotation = annotation, genome = "hg19")

saveRDS(HSR_chrom_assay, file = "/datastore/lbcfs/labs/brunk_lab/public/COLO320_multiome_hg19/HSR_fragments/HSR_chrom_assay_temporary.rds")


HSR_chrom_obj <- CreateSeuratObject(counts = HSR_chrom_assay, assay = "ATAC")
DefaultAssay(HSR_chrom_obj) <- "ATAC"
HSR_chrom_obj <- NucleosomeSignal(HSR_chrom_obj)
HSR_chrom_obj <- TSSEnrichment(HSR_chrom_obj)

p <- VlnPlot(HSR_chrom_obj, features = c("nCount_ATAC", "TSS.enrichment", "nucleosome_signal"), ncol = 3, log = TRUE, pt.size = 0) + NoLegend()

p

HSR_chrom_obj <- subset(x = HSR_chrom_obj, subset = nCount_ATAC < 7e4 & nCount_ATAC > 1000 & nucleosome_signal < 2 & TSS.enrichment > 1)

# calling peaks from ATAC
DefaultAssay(HSR_chrom_obj) <- "ATAC"
HSR_chrom_peaks <- CallPeaks(HSR_chrom_obj, macs2.path = "~/miniconda3/envs/macs3_env/bin/macs3")
HSR_chrom_peaks <- keepStandardChromosomes(HSR_chrom_peaks, pruning.mode = "coarse")
HSR_chrom_peaks <- subsetByOverlaps(x = HSR_chrom_peaks, ranges = blacklist_hg19, invert = TRUE)

# quantify counts in each peak
HSR_chrom_peak_counts <- FeatureMatrix(
  fragments = Fragments(HSR_chrom_obj),
  features = HSR_chrom_peaks,
  cells = colnames(HSR_chrom_obj)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
HSR_chrom_obj[["peaks"]] <- CreateChromatinAssay(
  counts = HSR_chrom_peak_counts,
  fragments = "/datastore/lbcfs/labs/brunk_lab/public/COLO320_multiome_hg19/HSR_fragments/HSR_fragments_chr.sorted.tsv.gz",
  annotation = annotation
)

DefaultAssay(HSR_chrom_obj) <- "peaks"
HSR_chrom_obj <- FindTopFeatures(HSR_chrom_obj, min.cutoff = 5)
HSR_chrom_obj <- RunTFIDF(HSR_chrom_obj)
HSR_chrom_obj <- RunSVD(HSR_chrom_obj)
HSR_chrom_obj <- RunUMAP(HSR_chrom_obj, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

HSR_chrom_obj.gene.activities <- GeneActivity(HSR_chrom_obj, biotypes = NULL)# to include not only protein coding genes but also long non coding genes
HSR_chrom_obj[['gene_activity']] <- CreateAssayObject(counts = HSR_chrom_obj.gene.activities)
HSR_chrom_obj <- NormalizeData(
  object = HSR_chrom_obj,
  assay = 'gene_activity',
  normalization.method = 'LogNormalize',
  scale.factor = 10000
)

#saveRDS(HSR_chrom_obj, "/datastore/lbcfs/labs/brunk_lab/public/COLO320_multiome_hg19/HSR_fragments/COLO320HSR_scATAC_object_05212025.rds")

#Now to examine the pseudobulk gene_activity of ecDNA genes:
df_gene_activity <- as.data.frame(AggregateExpression(HSR_chrom_obj, assays = "gene_activity"))
View(df_gene_activity)


#Now integrate scRNA data into the DM_chrom_obj seurat object:
COLO320_data <- Read10X(data.dir = "/datastore/lbcfs/labs/brunk_lab/public/COLO320_multiome_hg19/scRNA/")
COLO320_RNA_seurat_obj <- CreateSeuratObject(counts = COLO320_data, project = "COLO320")
COLO320_RNA_seurat_obj$cell_barcode <- colnames(COLO320_RNA_seurat_obj)
COLO320_RNA_seurat_obj$cleaned_barcodes <- sub(".*_rep[1-7]_", "", COLO320_RNA_seurat_obj$cell_barcode)
COLO320_RNA_seurat_obj_subset <- subset(COLO320_RNA_seurat_obj, subset = orig.ident == "COLO320HSR")

cell_names <- COLO320_RNA_seurat_obj_subset$cleaned_barcodes
cells_to_keep <- cell_names[!duplicated(cell_names)]
COLO320_RNA_seurat_obj_subset <- subset(COLO320_RNA_seurat_obj_subset, subset = cleaned_barcodes %in% cells_to_keep)

new_barcodes <- as.character(COLO320_RNA_seurat_obj_subset$cleaned_barcodes)
length(new_barcodes) == length(Cells(COLO320_RNA_seurat_obj_subset))
new_barcodes_unique <- make.unique(new_barcodes)
COLO320_RNA_seurat_obj_subset <- RenameCells(COLO320_RNA_seurat_obj_subset, new.names = new_barcodes_unique)


cells_to_keep <- as.character(intersect(COLO320_RNA_seurat_obj_subset$cleaned_barcodes, colnames(HSR_chrom_obj)))
HSR_chrom_obj_subset <- subset(HSR_chrom_obj, cells = cells_to_keep)
COLO320_RNA_seurat_obj_subset <- subset(COLO320_RNA_seurat_obj_subset, cells = cells_to_keep)

table(colnames(COLO320_RNA_seurat_obj_subset) == colnames(HSR_chrom_obj_subset)) #TRUE

HSR_chrom_obj_subset[['RNA']] <- COLO320_RNA_seurat_obj_subset[['RNA']]

DefaultAssay(HSR_chrom_obj_subset) <- 'RNA'
HSR_chrom_obj_subset[["percent.mt"]] <- PercentageFeatureSet(HSR_chrom_obj_subset, pattern = "^MT-")

p <- VlnPlot(HSR_chrom_obj_subset, features = c("nCount_ATAC", "nCount_RNA", "nFeature_RNA", "TSS.enrichment", "nucleosome_signal", "percent.mt"), ncol = 6,
             log = TRUE, pt.size = 0) + NoLegend()

p

HSR_chrom_obj_subset <- subset(x = HSR_chrom_obj_subset, subset = nCount_RNA < 25000 & nCount_RNA > 1000 & nFeature_RNA > 500  & percent.mt < 10)

saveRDS(HSR_chrom_obj_subset, "/datastore/lbcfs/labs/brunk_lab/public/COLO320_multiome_hg19/HSR_fragments/COLO320HSR_scATAC_scRNA_multiome_reintegrated_object_05212025.rds")

HSR_chrom_obj <- readRDS("/datastore/lbcfs/labs/brunk_lab/public/COLO320_multiome_hg19/HSR_fragments/COLO320HSR_scATAC_scRNA_multiome_reintegrated_object_05212025.rds")

DefaultAssay(HSR_chrom_obj) <- "RNA"
HSR_chrom_obj <- SCTransform(HSR_chrom_obj, vars.to.regress = "percent.mt", verbose = FALSE)
saveRDS(HSR_chrom_obj, "/datastore/lbcfs/labs/brunk_lab/public/COLO320_multiome_hg19/HSR_fragments/COLO320HSR_scATAC_scRNA_multiome_reintegrated_object_05212025.rds")





HSR_chrom_obj <- readRDS("/datastore/lbcfs/labs/brunk_lab/public/COLO320_multiome_hg19/HSR_fragments/COLO320HSR_scATAC_scRNA_multiome_reintegrated_object_05212025.rds")


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
df_sct <- as.data.frame(HSR_chrom_obj@assays$gene_activity@data[i,])
colnames(df_sct) <- 'gene'
df_sct$cells <- rownames(df_sct)
df3 <- df_sct[df_sct$gene >= quantile(as.numeric(HSR_chrom_obj@assays$gene_activity@data[i,]), .90, na.rm=TRUE), ]
df3$group <- 'higher'
df4 <- df_sct[df_sct$gene <= quantile(as.numeric(HSR_chrom_obj@assays$gene_activity@data[i,]), .10, na.rm=TRUE), ]
df4$group <- 'lower'
HSR_chrom_obj@meta.data$group <- lapply(colnames(HSR_chrom_obj), group.func)
Idents(HSR_chrom_obj) <- HSR_chrom_obj$group

DefaultAssay(HSR_chrom_obj) <- "gene_activity"
df_hsr_atac_markers <- FindMarkers(HSR_chrom_obj, ident.1 = "high", ident.2 = "low", min.pct = 0.25, logfc.threshold = 0.10, recorrect_umi=F)

DefaultAssay(HSR_chrom_obj) <- "SCT"
df_hsr_sct_markers <- FindMarkers(HSR_chrom_obj, ident.1 = "high", ident.2 = "low", min.pct = 0.25, logfc.threshold = 0.10, recorrect_umi=F)


list_ecDNA_genes <- c("MYC", "PVT1", "LINC02912", "PCAT1", "CASC19", "CASC11", "MIR1204", "POU5F1B", "CASC8","PRNCR1", "CCAT2",
                      "CDX2", "FAM84B", "RNU6-869P", "PLUT", "LINC00543", "PDX1", "ATP5EP2","CCAT1", "CASC21","TMEM75")

View(df_hsr_atac_markers[(df_hsr_atac_markers$p_val_adj < 0.05) & (rownames(df_hsr_atac_markers) %in% list_ecDNA_genes),])
View(df_hsr_sct_markers[(df_hsr_sct_markers$p_val_adj < 0.05) & (rownames(df_hsr_sct_markers) %in% list_ecDNA_genes),])

volcano_plot_function(df_hsr_atac_markers)
ggsave("/home/yue1118/ecDNA_project_data_and_figures/FACS_FISH_FULL_figures/COLO320HSR_differential_peaks_MYC_gene_activity_high_vs_low_05292025.svg", width = 8, height = 6)


volcano_plot_function(df_hsr_sct_markers)
ggsave("/home/yue1118/ecDNA_project_data_and_figures/FACS_FISH_FULL_figures/COLO320HSR_differential_expressed_genes_between_MYC_gene_activity_high_vs_low_05292025.svg", width = 8, height = 6)







###To plot pseudobulk gene expression of all the genes in the genome and highlight PVT1
HSR_chrom_obj <- readRDS("/datastore/lbcfs/labs/brunk_lab/public/COLO320_multiome_hg19/HSR_fragments/COLO320HSR_scATAC_scRNA_multiome_reintegrated_object_05212025.rds")
Idents(HSR_chrom_obj) <- 'WT'
df_sc_aggre <- AggregateExpression(HSR_chrom_obj, assays = 'SCT', slot = 'data')
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
ggsave('/home/yue1118/Bioskryb_data_analysis/10x_multiome_gene_expression_across_the_genome_COLO320DM_muiltiome_WT_06182025.svg', height = 5, width = 5)







