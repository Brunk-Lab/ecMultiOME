#### Link genes to peaks using NCIH2170 10X Multiome data (scRNA-seq + scATAC-seq)

# Load Seurat object

NCIH2170 <- readRDS("path/to/NCIH2170.RDS")

######################################  LINK PEAKS TO GENES FUNCTION

library(biomaRt)
library(GenomicRanges)
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)

DefaultAssay(NCIH2170) <- "RNA"


Reductions(NCIH2170)
Assays(NCIH2170)



# parameters
seurat_object <- NCIH2170  # your Seurat object
peak_assay <- "peaks"   # your peaks assay, for better results use MACS2 peak calling
expression_assay <- "SCT"  ## Normalized RNA
reduction_method <- "harmony"  # or "umap", lsi", "pca"

generate_metacells <- function(seurat_object, peak_assay, expression_assay, reduction = "pca", dims_to_use = 1:30, k = 50, knn_iterations = 500, overlap_cutoff = 0.8, verbose = TRUE) {
  # Ensure the specified reduction is computed
  if (!(reduction %in% names(seurat_object@reductions))) {
    stop(paste("Reduction method", reduction, "not found in the Seurat object."))
  }
  
  # Create a kNN graph for the cells using the specified reduction
  seurat_object <- FindNeighbors(seurat_object, reduction = reduction, dims = dims_to_use, verbose = verbose)
  
  # Identify non-overlapping kNN groups
  kNN_groups <- list()
  all_cells <- Cells(seurat_object)
  for (i in 1:knn_iterations) {
    group <- sample(all_cells, size = k)
    if (length(intersect(unlist(kNN_groups), group)) / length(group) < overlap_cutoff) {
      kNN_groups[[length(kNN_groups) + 1]] <- group
    }
  }
  
  # Debug: Print the number of kNN groups
  if (verbose) {
    cat("Number of kNN groups identified:", length(kNN_groups), "\n")
  }
  
  # Aggregate peak counts and gene expression counts within each group
  peak_data <- GetAssayData(seurat_object, assay = peak_assay, slot = "counts")
  expression_data <- GetAssayData(seurat_object, assay = expression_assay, slot = "data")
  
  aggregated_peaks <- do.call(cbind, lapply(kNN_groups, function(group) {
    rowSums(peak_data[, group, drop = FALSE])
  }))
  
  aggregated_expression <- do.call(cbind, lapply(kNN_groups, function(group) {
    rowMeans(expression_data[, group, drop = FALSE])
  }))
  
  # Debug: Print the dimensions of the aggregated matrices
  if (verbose) {
    cat("Dimensions of aggregated peaks matrix:", dim(aggregated_peaks), "\n")
    cat("Dimensions of aggregated expression matrix:", dim(aggregated_expression), "\n")
  }
  
  return(list(metacells_peaks = aggregated_peaks, metacells_expression = aggregated_expression))
}



metacells_data <- generate_metacells(seurat_object, peak_assay, expression_assay, reduction = reduction_method)

# Inspect the results
print(dim(metacells_data$metacells_peaks))
print(dim(metacells_data$metacells_expression))

metacells_data
str(metacells_data)


# Function to extract genomic ranges from a Seurat object
extract_genomic_ranges <- function(seurat_object, assay) {
  peaks_gr <- rownames(seurat_object[[assay]])
  
  # Print the first few entries to understand the structure
  print(peaks_gr[1:10])
  
  peaks_df <- do.call(rbind, strsplit(peaks_gr, "-"))
  
  # Ensure we have three columns after splitting the peak identifiers
  if (ncol(peaks_df) != 3) {
    stop("Expected 3 columns after splitting the peak identifiers")
  }
  
  gr <- GRanges(seqnames = peaks_df[, 1], ranges = IRanges(start = as.numeric(peaks_df[, 2]), end = as.numeric(peaks_df[, 3])))
  names(gr) <- peaks_gr
  return(gr)
}

# Extract genomic ranges for peaks
peaks_gr <- extract_genomic_ranges(seurat_object, peak_assay)

# Inspect the GRanges object
print(peaks_gr)



#################



# Function to fetch genomic coordinates for genes
fetch_gene_coordinates <- function(genes) {
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  gene_coords <- getBM(
    attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
    filters = "hgnc_symbol",
    values = genes,
    mart = ensembl
  )
  gene_coords <- gene_coords[gene_coords$chromosome_name %in% c(1:22, "X", "Y"), ]
  colnames(gene_coords) <- c("gene", "chromosome", "start", "end")
  gene_coords$chromosome <- paste0("chr", gene_coords$chromosome)
  return(gene_coords)
}

# Function to create a GRanges object for genes from gene_coords data frame
create_gene_granges <- function(gene_coords) {
  gr <- GRanges(
    seqnames = gene_coords$chromosome,
    ranges = IRanges(start = gene_coords$start, end = gene_coords$end),
    names = gene_coords$gene
  )
  return(gr)
}

# Function to extract genomic ranges from a Seurat object
extract_genomic_ranges <- function(seurat_object, assay) {
  peaks_gr <- rownames(seurat_object[[assay]])
  peaks_df <- do.call(rbind, strsplit(peaks_gr, "-"))
  if (ncol(peaks_df) != 3) {
    stop("Expected 3 columns after splitting the peak identifiers")
  }
  gr <- GRanges(seqnames = peaks_df[, 1], ranges = IRanges(start = as.numeric(peaks_df[, 2]), end = as.numeric(peaks_df[, 3])))
  names(gr) <- peaks_gr
  return(gr)
}




# Extract gene names from the SCT assay
genes <- rownames(seurat_object[[expression_assay]])
# Fetch gene coordinates
gene_coords <- fetch_gene_coordinates(genes)
# Create GRanges object for genes
genes_gr <- create_gene_granges(gene_coords)

# Extract genomic ranges for peaks
peaks_gr <- extract_genomic_ranges(seurat_object, peak_assay)

# Convert seqnames to characters to avoid factor level issues
seqlevelsStyle(genes_gr) <- "UCSC"
seqlevelsStyle(peaks_gr) <- "UCSC"



# Function to find nearby peaks for each gene
find_nearby_peaks <- function(gene, peaks_gr, max_distance = 50000) {
  # Find peaks that are on the same chromosome and within the max_distance
  nearby_peaks <- peaks_gr[seqnames(peaks_gr) == as.character(seqnames(gene)) & 
                             (start(peaks_gr) >= start(gene) - max_distance) &
                             (start(peaks_gr) <= end(gene) + max_distance)]
  return(nearby_peaks)
}

# Initialize a list to store the results
gene_peak_pairs <- list()

# Loop through each gene and find nearby peaks
for (i in seq_along(genes_gr)) {
  gene <- genes_gr[i]
  nearby_peaks <- find_nearby_peaks(gene, peaks_gr)
  if (length(nearby_peaks) > 0) {
    gene_peak_pairs[[as.character(genes_gr$names[i])]] <- nearby_peaks
  }
}

# Function to calculate correlations between gene expression and nearby peak accessibility
calculate_gene_peak_correlations <- function(gene_peak_pairs, metacells_data) {
  gene_peak_correlations <- list()
  
  for (gene in names(gene_peak_pairs)) {
    peaks <- gene_peak_pairs[[gene]]
    
    # Check if the gene is in the metacells expression data
    if (gene %in% rownames(metacells_data$metacells_expression)) {
      gene_expression <- metacells_data$metacells_expression[gene, ]
      
      peak_correlations <- list()
      
      for (i in seq_along(peaks)) {
        peak <- peaks[i]
        peak_id <- paste0(seqnames(peak), "-", start(peak), "-", end(peak))
        
        # Check if the peak is in the metacells peaks data
        if (peak_id %in% rownames(metacells_data$metacells_peaks)) {
          peak_accessibility <- metacells_data$metacells_peaks[peak_id, ]
          
          # Calculate the correlation
          correlation <- cor(gene_expression, peak_accessibility, method = "pearson")
          peak_correlations[[peak_id]] <- correlation
        }
      }
      
      if (length(peak_correlations) > 0) {
        gene_peak_correlations[[gene]] <- peak_correlations
      }
    }
  }
  
  return(gene_peak_correlations)
}

# Calculate correlations for all gene-peak pairs
gene_peak_correlations <- calculate_gene_peak_correlations(gene_peak_pairs, metacells_data)



# Inspect the results
head(gene_peak_correlations)



# Check if metacells_expression has proper column names
if (is.null(colnames(metacells_data$metacells_expression))) {
  cat("metacells_expression has no column names. Setting dummy column names.\n")
  colnames(metacells_data$metacells_expression) <- paste0("metacell", 1:ncol(metacells_data$metacells_expression))
}

# Verify the structure again
cat("Column names of metacells_expression after setting dummy names:\n")
print(colnames(metacells_data$metacells_expression))

##########################################


# Function to calculate null correlations for gene-peak pairs with correlation > 0.2 , set correlation value threshold desired  

calculate_null_correlations_filtered <- function(
    gene_peak_pairs,
    metacells_data,
    gene_peak_correlations,
    threshold = 0.2,
    num_permutations = 100
) {
  null_correlations <- list()
  
  for (gene in names(gene_peak_pairs)) {
    cat("Processing gene:", gene, "\n")
    
    # Skip if gene not in expression
    if (!gene %in% rownames(metacells_data$metacells_expression)) {
      cat("Gene not found in metacells_expression. Skipping.\n")
      next
    }
    
    # Extract gene expression vector
    gene_expression <- as.numeric(metacells_data$metacells_expression[gene, ])
    peaks <- gene_peak_pairs[[gene]]
    
    for (i in seq_along(peaks)) {
      peak <- peaks[i]
      peak_id <- paste0(seqnames(peak), "-", start(peak), "-", end(peak))
      
      # Skip if peak not in peaks
      if (!peak_id %in% rownames(metacells_data$metacells_peaks)) {
        cat("Peak ID not found in metacells_peaks:", peak_id, "Skipping.\n")
        next
      }
      
      peak_accessibility <- as.numeric(metacells_data$metacells_peaks[peak_id, ])
      gene_peak_key <- paste(gene, peak_id, sep = "_")
      
      # Get the observed correlation from your precomputed list
      observed_correlation <- gene_peak_correlations[[gene]][[peak_id]]
      
      # Skip dimension mismatch
      if (length(gene_expression) != length(peak_accessibility)) {
        cat("Dimension mismatch for gene:", gene, "and peak:", peak_id, "Skipping.\n")
        next
      }
      
      # Skip if correlation is NA
      if (is.na(observed_correlation)) {
        cat("Correlation is NA for gene:", gene, "peak:", peak_id, "Skipping.\n")
        next
      }
      
      # Only permute if correlation exceeds threshold
      if (abs(observed_correlation) > threshold) {
        cat("Processing gene-peak pair:", gene, "-", peak_id,
            "with observed correlation:", observed_correlation, "\n")
        null_cor_values <- numeric(num_permutations)
        
        for (j in seq_len(num_permutations)) {
          permuted_gene_expression <- sample(gene_expression)
          correlation_val <- cor(permuted_gene_expression, peak_accessibility, method = "pearson")
          null_cor_values[j] <- correlation_val
        }
        
        null_correlations[[gene_peak_key]] <- null_cor_values
        
      } else {
        cat("Skipping gene-peak pair:", gene, "-", peak_id,
            "with low observed correlation:", observed_correlation, "\n")
      }
    }
  }
  
  return(null_correlations)
}


# Example usage:
null_correlations_filtered <- calculate_null_correlations_filtered(gene_peak_pairs, metacells_data, gene_peak_correlations, threshold = 0.2)

# Print summary of filtered null correlations
summary(null_correlations_filtered)






# Function to calculate null correlations for random gene-peak pairs
calculate_random_null_correlations <- function(gene_peak_pairs, metacells_data, num_random_pairs = 10000, num_permutations = 100) {
  all_genes <- rownames(metacells_data$metacells_expression)
  all_peaks <- rownames(metacells_data$metacells_peaks)
  
  null_correlations <- list()
  
  set.seed(123) # For reproducibility
  for (i in 1:num_random_pairs) {
    gene <- sample(all_genes, 1)
    peak_id <- sample(all_peaks, 1)
    
    gene_expression <- as.numeric(metacells_data$metacells_expression[gene, ])
    peak_accessibility <- as.numeric(metacells_data$metacells_peaks[peak_id, ])
    
    gene_peak_key <- paste(gene, peak_id, sep = "_")
    
    null_correlations[[gene_peak_key]] <- numeric(num_permutations)
    
    for (j in 1:num_permutations) {
      permuted_gene_expression <- sample(gene_expression)
      correlation <- cor(permuted_gene_expression, peak_accessibility, method = "pearson")
      null_correlations[[gene_peak_key]][j] <- correlation
    }
    
    if (i %% 100 == 0) {
      cat("Processed", i, "random gene-peak pairs\n")
    }
  }
  
  return(null_correlations)
}

# Execute the function with your data
null_correlations_random <- calculate_random_null_correlations(gene_peak_pairs, metacells_data, num_random_pairs = 10000, num_permutations = 100)


combined_null_correlations <- c(null_correlations_filtered, null_correlations_random)



############## 


# Function to calculate p-value for a single gene-peak pair using Gaussian distribution
calculate_p_value <- function(null_correlations, actual_correlation) {
  # Fit a Gaussian distribution to the null correlations
  mean_null <- mean(null_correlations)
  sd_null <- sd(null_correlations)
  
  # Calculate the z-score for the actual correlation
  z_score <- (actual_correlation - mean_null) / sd_null
  
  # Calculate the p-value
  p_value <- pnorm(-abs(z_score))
  return(p_value)
}

# Extract all gene-peak pairs from null_correlations_filtered
gene_peak_pairs <- names(null_correlations_filtered)

# Initialize a list to store results
results <- list()

# Iterate through each gene-peak pair
for (gene_peak in gene_peak_pairs) {
  # Extract the gene and peak information
  gene_peak_split <- strsplit(gene_peak, "_")[[1]]
  gene <- gene_peak_split[1]
  peak <- paste(gene_peak_split[-1], collapse = "_")
  
  # Extract the null correlations for the current gene-peak pair
  null_correlations <- null_correlations_filtered[[gene_peak]]
  
  # Extract the actual correlation from gene_peak_correlations
  if (!is.null(gene_peak_correlations[[gene]])) {
    actual_correlation <- gene_peak_correlations[[gene]][[peak]]
    
    if (!is.null(actual_correlation)) {
      # Calculate the p-value for the current gene-peak pair
      p_value <- calculate_p_value(null_correlations, actual_correlation)
      
      # Store the results
      results[[gene_peak]] <- list(
        gene_peak_pair = gene_peak,
        actual_correlation = actual_correlation,
        p_value = p_value
      )
    }
  }
}

# Convert the list of results to a data frame
results_df <- do.call(rbind, lapply(results, as.data.frame))

# Adjust p-values for multiple testing using the Benjamini-Hochberg method
results_df$fdr <- p.adjust(results_df$p_value, method = "BH")

# Print the head of the results data frame
print(head(results_df))

# View the complete results data frame
results_df


write.csv(results_df, "gene_peak_correlations.csv")







# Function to process the data and prepare for making the density plot:
process_peak_data <- function(data) {
  # Filter for adjusted p-values (FDR) less than 0.1
  filtered_data <- dplyr::filter(data, fdr < 0.1)
  
  # Extract the gene_peak column
  peaks <- dplyr::select(filtered_data, gene_peak_pair)
  
  # Split the gene_peak into gene, chromosome, start, and end
  peaks <- peaks %>%
    separate(gene_peak_pair, into = c("gene", "peak"), sep = "_") %>%
    separate(peak, into = c("chr", "start", "end"), sep = "-") %>%
    mutate(start = as.numeric(start), end = as.numeric(end))
  
  # Calculate peak density per 50 kb window
  bin_size <- 50000
  peaks <- peaks %>%
    mutate(bin_start = floor(start / bin_size) * bin_size,
           bin_end = bin_start + bin_size)
  
  peak_density <- peaks %>%
    group_by(chr, bin_start, bin_end) %>%
    summarize(peak_count = n(), .groups = 'drop')
  
  return(peak_density)
}

# Function to plot peak density heatmap and save to a file
plot_peak_density_heatmap <- function(peak_density, output_file) {
  # Debug: Print the structure of peak_density
  print(head(peak_density))
  
  # Ensure chr is a factor for proper facet ordering
  peak_density$chr <- factor(peak_density$chr, levels = unique(peak_density$chr))
  
  # Create the heatmap plot
  p <- ggplot(peak_density, aes(x = bin_start, y = chr, fill = peak_count)) +
    geom_tile() +
    scale_fill_gradient2(low = "white", high = "red", mid = "yellow", midpoint = median(peak_density$peak_count), limits = c(0, max(peak_density$peak_count)), space = "Lab", name="Peak Count") +
    labs(title = "Peak Density per 50 kb Window",
         x = "Genomic Position (bp)",
         y = "Chromosome") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  # Save the plot to a file
  ggsave(output_file, plot = p, width = 12, height = 8)
}

# Function to annotate bins with the highest peak density
annotate_high_density_bins <- function(peak_density) {
  # Find bins with the highest peak density for each chromosome
  high_density_bins <- peak_density %>%
    group_by(chr) %>%
    dplyr::filter(peak_count == max(peak_count)) %>%
    ungroup()
  
  return(high_density_bins)
}

# Function to get genes in bins with peak density > 3
get_genes_in_high_density_bins <- function(data, peak_density) {
  high_density_bins <- peak_density %>%
    dplyr::filter(peak_count > 3)
  
  genes_in_bins <- data %>%
    separate(gene_peak_pair, into = c("gene", "peak"), sep = "_") %>%
    separate(peak, into = c("chr", "start", "end"), sep = "-") %>%
    mutate(start = as.numeric(start), end = as.numeric(end)) %>%
    mutate(bin_start = floor(start / 50000) * 50000,
           bin_end = bin_start + 50000) %>%
    inner_join(high_density_bins, by = c("chr", "bin_start", "bin_end")) %>%
    dplyr::select(gene) %>%
    distinct()
  
  return(genes_in_bins)
}

#Execute with the result file
peak_density <- process_peak_data(results_df)
high_density_bins <- annotate_high_density_bins(peak_density)

# Plot the peak density heatmap and save it to a PDF file
output_file <- "peak_density_heatmap.pdf"
plot_peak_density_heatmap(peak_density, output_file)

# Get the genes in bins with peak density > 3
genes_in_high_density_bins <- get_genes_in_high_density_bins(results_df, peak_density)

# Print the genes in high-density bins
print(genes_in_high_density_bins)




#######


#To prepare for generating a density heatmap.

# Load your metacells data
metacells_peaks <- metacells_data$metacells_peaks
metacells_expression <- metacells_data$metacells_expression

# Function to process the data
process_peak_data <- function(data) {
  filtered_data <- dplyr::filter(data, fdr < 0.1)
  peaks <- dplyr::select(filtered_data, gene_peak_pair)
  peaks <- peaks %>%
    separate(gene_peak_pair, into = c("gene", "peak"), sep = "_") %>%
    separate(peak, into = c("chr", "start", "end"), sep = "-") %>%
    mutate(start = as.numeric(start), end = as.numeric(end))
  bin_size <- 50000
  peaks <- peaks %>%
    mutate(bin_start = floor(start / bin_size) * bin_size,
           bin_end = bin_start + bin_size)
  peak_density <- peaks %>%
    group_by(chr, bin_start, bin_end) %>%
    summarize(peak_count = n(), .groups = 'drop')
  return(peak_density)
}

# Function to get genes in bins with peak density > 3
get_genes_in_high_density_bins <- function(data, peak_density) {
  high_density_bins <- peak_density %>%
    dplyr::filter(peak_count > 3)
  genes_in_bins <- data %>%
    separate(gene_peak_pair, into = c("gene", "peak"), sep = "_") %>%
    separate(peak, into = c("chr", "start", "end"), sep = "-") %>%
    mutate(start = as.numeric(start), end = as.numeric(end)) %>%
    mutate(bin_start = floor(start / 50000) * 50000,
           bin_end = bin_start + 50000) %>%
    inner_join(high_density_bins, by = c("chr", "bin_start", "bin_end")) %>%
    dplyr::select(gene) %>%
    distinct()
  return(genes_in_bins)
}

# Function to calculate gene-peak correlations
calculate_gene_peak_correlations <- function(gene, peaks_data, expression_data) {
  gene_expression <- expression_data[gene, ]
  correlations <- numeric(nrow(peaks_data))
  for (i in 1:nrow(peaks_data)) {
    peak_counts <- peaks_data[i, ]
    correlations[i] <- cor(peak_counts, gene_expression, method = "pearson")
  }
  filtered_correlations <- correlations[correlations > 0.2]
  filtered_peaks <- rownames(peaks_data)[correlations > 0.2]
  result <- data.frame(
    peak = filtered_peaks,
    correlation = filtered_correlations
  )
  return(result)
}

# Main script
peak_density <- process_peak_data(results_df)
genes_in_high_density_bins <- get_genes_in_high_density_bins(results_df, peak_density)

# Initialize a list to store results for each gene
all_results <- list()

# Iterate over each gene and perform the correlation analysis
for (gene in genes_in_high_density_bins$gene) {
  result <- calculate_gene_peak_correlations(gene, metacells_peaks, metacells_expression)
  all_results[[gene]] <- result
  print(paste("Results for gene:", gene))
  print(result)
}

# Combine all results into a single data frame
final_results <- bind_rows(all_results, .id = "gene")

# Save the final results to a CSV file
write.csv(final_results, "gene_peak_correlations.csv", row.names = FALSE)

# Function to calculate peak density per 50 kb window for each gene
calculate_peak_density_per_gene <- function(correlated_peaks) {
  peaks <- correlated_peaks %>%
    separate(peak, into = c("chr", "start", "end"), sep = "-") %>%
    mutate(start = as.numeric(start), end = as.numeric(end))
  bin_size <- 200000
  peaks <- peaks %>%
    mutate(bin_start = floor(start / bin_size) * bin_size,
           bin_end = bin_start + bin_size)
  peak_density <- peaks %>%
    group_by(chr, bin_start, bin_end) %>%
    summarize(peak_count = n(), .groups = 'drop')
  return(peak_density)
}

# Function to find high density areas
find_high_density_areas <- function(peak_density) {
  high_density_areas <- peak_density %>%
    dplyr::filter(peak_count > 3)
  return(high_density_areas)
}

# Initialize a dataframe to store results for each gene
all_high_density_areas <- data.frame()



# Iterate over each gene and perform the correlation analysis and density calculation
for (gene in genes_in_high_density_bins$gene) {
  correlated_peaks <- calculate_gene_peak_correlations(gene, metacells_peaks_NCIH2170DM, metacells_expression_NCIH2170DM)
  peak_density <- calculate_peak_density_per_gene(correlated_peaks)
  high_density_areas <- find_high_density_areas(peak_density)
  
  high_density_areas$gene <- gene
  all_high_density_areas <- rbind(all_high_density_areas, high_density_areas)
}

# Print the dataframe with all high density areas per gene
print(all_high_density_areas)

write.csv(all_high_density_areas, "high_density_areas_per_gene.csv", row.names = FALSE)



