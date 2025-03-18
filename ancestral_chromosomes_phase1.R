#!/usr/bin/env Rscript

# Ancestral Chromosome Reconstruction Phase 1: Data Integration Framework
# Date: 2025-03-18
# Author: ShrivatsaMehan
# Description: Integrates phylogenetic tree, chromosome counts, and bidirectional mapping data,
#              processing ONLY species that are common to all datasets.

# Required packages
suppressPackageStartupMessages({
  library(ape)         # For phylogenetic tree handling
  library(data.table)  # For efficient data manipulation
})

# ==============================
# Data Loading Functions
# ==============================

load_phylogenetic_tree <- function(tree_file) {
  # Load phylogenetic tree from Newick file
  message("Loading phylogenetic tree...")
  tree <- try(read.tree(tree_file), silent = TRUE)
  
  if(inherits(tree, "try-error")) {
    stop("Error reading phylogenetic tree file: ", tree_file)
  }
  
  # Basic validation
  if(is.null(tree) || length(tree$tip.label) < 3) {
    stop("Invalid tree file or insufficient taxa (< 3)")
  }
  
  message(paste("Loaded tree with", length(tree$tip.label), "species and", tree$Nnode, "internal nodes"))
  return(tree)
}

load_chromosome_counts <- function(counts_file) {
  # Load chromosome counts from CSV file
  message("Loading chromosome counts...")
  
  if(!file.exists(counts_file)) {
    stop("Chromosome counts file not found: ", counts_file)
  }
  
  counts_data <- try(read.csv(counts_file, stringsAsFactors = FALSE), silent = TRUE)
  
  if(inherits(counts_data, "try-error")) {
    stop("Error reading chromosome counts file")
  }
  
  # Check required columns
  if(!all(c("species", "count") %in% colnames(counts_data))) {
    stop("Chromosome counts file must contain 'species' and 'count' columns")
  }
  
  # Convert to named vector for easier access
  chr_counts <- counts_data$count
  names(chr_counts) <- counts_data$species
  
  message(paste("Loaded chromosome counts for", length(chr_counts), "species"))
  return(list(
    data_frame = counts_data,
    counts = chr_counts
  ))
}

load_bidirectional_maps <- function(maps_dir) {
  # Load bidirectional mapping data from TSV files
  message("Loading bidirectional mapping data...")
  
  if(!dir.exists(maps_dir)) {
    stop("Bidirectional maps directory not found: ", maps_dir)
  }
  
  # Find all mapping files
  map_files <- list.files(maps_dir, pattern = "*_bidirectional\\.tsv$", full.names = TRUE)
  
  if(length(map_files) == 0) {
    # Try any TSV files if specific pattern not found
    map_files <- list.files(maps_dir, pattern = "\\.tsv$", full.names = TRUE)
    
    if(length(map_files) == 0) {
      stop("No mapping files found in directory: ", maps_dir)
    }
  }
  
  message(paste("Found", length(map_files), "bidirectional mapping files"))
  
  # Read all files
  all_maps <- lapply(map_files, function(file) {
    tryCatch({
      map_data <- fread(file)
      
      # Check required columns
      req_cols <- c("species_A", "chromosome_A", "species_B", "chromosome_B")
      if(!all(req_cols %in% names(map_data))) {
        warning("File missing required columns, skipping: ", file)
        return(NULL)
      }
      return(map_data)
    }, error = function(e) {
      warning("Error reading file: ", file, " - ", e$message)
      return(NULL)
    })
  })
  
  # Remove NULL entries (failed reads)
  all_maps <- all_maps[!sapply(all_maps, is.null)]
  
  if(length(all_maps) == 0) {
    stop("No valid mapping data could be loaded")
  }
  
  # Combine all mapping data
  combined_maps <- rbindlist(all_maps, fill = TRUE)
  message(paste("Combined mapping data has", nrow(combined_maps), "rows"))
  
  # Extract all species from the mapping data
  mapping_species <- unique(c(combined_maps$species_A, combined_maps$species_B))
  message(paste("Mapping data contains", length(mapping_species), "species"))
  
  return(list(
    combined = combined_maps,
    species = mapping_species
  ))
}

# ==============================
# Data Integration Functions
# ==============================

find_common_species <- function(tree_species, chr_count_species, mapping_species) {
  # Find species that are common to all three datasets
  message("Finding species common to all datasets...")
  
  # Find species present in all datasets
  common_species <- Reduce(intersect, list(tree_species, chr_count_species, mapping_species))
  
  if(length(common_species) < 3) {
    stop("Less than 3 species common to all datasets. Cannot proceed with analysis.")
  }
  
  # Calculate overlap statistics
  tree_coverage <- length(common_species) / length(tree_species) * 100
  counts_coverage <- length(common_species) / length(chr_count_species) * 100
  mapping_coverage <- length(common_species) / length(mapping_species) * 100
  
  message(paste("Found", length(common_species), "species common to all datasets"))
  message(paste("Common species represent:",
                round(tree_coverage, 1), "% of tree species,",
                round(counts_coverage, 1), "% of chromosome count species,",
                round(mapping_coverage, 1), "% of mapping species"))
  
  return(common_species)
}

validate_species_consistency <- function(tree_species, chr_count_species, mapping_species, common_species) {
  # Report on species missing from different datasets
  message("Generating species consistency report...")
  
  # Find species missing in one or more datasets
  missing_in_tree <- setdiff(chr_count_species, tree_species)
  missing_in_counts <- setdiff(tree_species, chr_count_species)
  missing_in_mapping <- intersect(tree_species, chr_count_species)
  missing_in_mapping <- setdiff(missing_in_mapping, mapping_species)
  
  # Report missing species
  if(length(missing_in_tree) > 0) {
    warning(paste("Species in chromosome counts but missing from tree:", 
                  paste(missing_in_tree[1:min(5, length(missing_in_tree))], collapse=", "),
                  if(length(missing_in_tree) > 5) "... and others"))
  }
  
  if(length(missing_in_counts) > 0) {
    warning(paste("Species in tree but missing chromosome counts:", 
                  paste(missing_in_counts[1:min(5, length(missing_in_counts))], collapse=", "),
                  if(length(missing_in_counts) > 5) "... and others"))
  }
  
  if(length(missing_in_mapping) > 0) {
    warning(paste("Species in tree and counts but missing from mapping data:", 
                  paste(missing_in_mapping[1:min(5, length(missing_in_mapping))], collapse=", "),
                  if(length(missing_in_mapping) > 5) "... and others"))
  }
  
  # Verify that analysis will proceed with common species only
  message(paste("Analysis will proceed with ONLY the", length(common_species), 
                "species present in ALL datasets"))
  
  return(list(
    common_species = common_species,
    missing_in_tree = missing_in_tree,
    missing_in_counts = missing_in_counts,
    missing_in_mapping = missing_in_mapping
  ))
}

prune_phylogenetic_tree <- function(tree, common_species) {
  # Prune phylogenetic tree to keep only species present in all datasets
  message("Pruning phylogenetic tree to common species only...")
  
  pruned_tree <- try(keep.tip(tree, common_species), silent = TRUE)
  
  if(inherits(pruned_tree, "try-error")) {
    stop("Error pruning phylogenetic tree")
  }
  
  # Verify that pruning worked correctly
  if(length(pruned_tree$tip.label) != length(common_species)) {
    stop("Error: Pruned tree does not contain exactly the common species")
  }
  
  # Check if all tips in pruned tree are in common_species
  if(!all(pruned_tree$tip.label %in% common_species)) {
    stop("Error: Pruned tree contains species not in common species list")
  }
  
  message(paste("Pruned tree has", length(pruned_tree$tip.label), "species and", 
                pruned_tree$Nnode, "internal nodes"))
  
  return(pruned_tree)
}

filter_chromosome_counts <- function(chr_counts, common_species) {
  # Filter chromosome counts to keep only species present in all datasets
  message("Filtering chromosome counts to common species only...")
  
  filtered_counts <- chr_counts[common_species]
  
  # Verify filtering worked correctly
  if(length(filtered_counts) != length(common_species)) {
    stop("Error: Filtered counts do not contain exactly the common species")
  }
  
  # Check if all names in filtered_counts are in common_species
  if(!all(names(filtered_counts) %in% common_species)) {
    stop("Error: Filtered counts contain species not in common species list")
  }
  
  message(paste("Filtered to", length(filtered_counts), "species with chromosome counts"))
  
  return(filtered_counts)
}

filter_bidirectional_maps <- function(maps_data, common_species) {
  # Filter mapping data to keep only species present in all datasets
  message("Filtering bidirectional mapping data to common species only...")
  
  combined_maps <- maps_data$combined
  filtered_maps <- combined_maps[species_A %in% common_species & species_B %in% common_species]
  
  # Verify that filtering worked correctly
  if(!all(filtered_maps$species_A %in% common_species) || 
     !all(filtered_maps$species_B %in% common_species)) {
    stop("Error: Filtered maps contain species not in common species list")
  }
  
  message(paste("Filtered mapping data has", nrow(filtered_maps), "rows"))
  
  # Calculate species coverage in mapping data
  species_in_mapping <- unique(c(filtered_maps$species_A, filtered_maps$species_B))
  
  # Check if all common species appear in mapping data
  if(!all(common_species %in% species_in_mapping)) {
    missing_in_mapping <- setdiff(common_species, species_in_mapping)
    warning(paste(length(missing_in_mapping), "common species do not appear in any mapping relationship"))
  }
  
  # Calculate species coverage statistics
  species_coverage <- table(c(filtered_maps$species_A, filtered_maps$species_B))
  
  min_coverage <- min(species_coverage)
  max_coverage <- max(species_coverage)
  mean_coverage <- mean(species_coverage)
  
  message(paste("Species mapping coverage - min:", min_coverage, 
                "max:", max_coverage, 
                "mean:", round(mean_coverage, 1)))
  
  return(filtered_maps)
}

integrate_data <- function(tree_file, chr_counts_file, maps_dir) {
  # Main function to integrate all data sources
  message("Starting data integration process...")
  message("NOTE: Analysis will only include species present in ALL three datasets")
  
  # Load data
  tree <- load_phylogenetic_tree(tree_file)
  chr_counts_data <- load_chromosome_counts(chr_counts_file)
  maps_data <- load_bidirectional_maps(maps_dir)
  
  # Extract species lists
  tree_species <- tree$tip.label
  chr_count_species <- names(chr_counts_data$counts)
  mapping_species <- maps_data$species
  
  # Find common species across all datasets
  common_species <- find_common_species(tree_species, chr_count_species, mapping_species)
  
  # Validate species consistency and generate report
  consistency <- validate_species_consistency(tree_species, chr_count_species, 
                                             mapping_species, common_species)
  
  # Filter data to common species only
  pruned_tree <- prune_phylogenetic_tree(tree, common_species)
  filtered_counts <- filter_chromosome_counts(chr_counts_data$counts, common_species)
  filtered_maps <- filter_bidirectional_maps(maps_data, common_species)
  
  # Create integrated data structure
  integrated_data <- list(
    tree = pruned_tree,
    chr_counts = filtered_counts,
    bidirectional_maps = filtered_maps,
    common_species = common_species,
    consistency_report = consistency,
    original = list(
      tree_species_count = length(tree_species),
      chr_count_species_count = length(chr_count_species),
      mapping_species_count = length(mapping_species)
    )
  )
  
  message("Data integration complete - all datasets filtered to common species only")
  return(integrated_data)
}

# ==============================
# Output Functions
# ==============================

save_integrated_data <- function(integrated_data, output_dir) {
  # Save integrated data to intermediate files
  message("Saving integrated data...")
  
  # Create output directory if it doesn't exist
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message(paste("Created output directory:", output_dir))
  }
  
  # Save pruned tree
  pruned_tree_file <- file.path(output_dir, "common_species_tree.nwk")
  write.tree(integrated_data$tree, file = pruned_tree_file)
  
  # Save filtered chromosome counts
  filtered_counts_file <- file.path(output_dir, "common_species_chromosome_counts.csv")
  counts_df <- data.frame(
    species = names(integrated_data$chr_counts),
    count = integrated_data$chr_counts,
    stringsAsFactors = FALSE
  )
  write.csv(counts_df, file = filtered_counts_file, row.names = FALSE)
  
  # Save filtered mapping data
  filtered_maps_file <- file.path(output_dir, "common_species_bidirectional_maps.tsv")
  fwrite(integrated_data$bidirectional_maps, file = filtered_maps_file, sep = "\t")
  
  # Save intermediate R object for next phases
  integrated_rdata_file <- file.path(output_dir, "integrated_data.RData")
  saveRDS(integrated_data, file = integrated_rdata_file)
  
  # Save species consistency report
  consistency_file <- file.path(output_dir, "species_consistency_report.txt")
  
  sink(consistency_file)
  cat("Species Consistency Report\n")
  cat("=========================\n\n")
  
  cat("IMPORTANT: Analysis includes ONLY species present in ALL three datasets\n\n")
  
  cat(paste("Total species in original tree:", integrated_data$original$tree_species_count, "\n"))
  cat(paste("Total species in original chromosome counts:", integrated_data$original$chr_count_species_count, "\n"))
  cat(paste("Total species in original mapping data:", integrated_data$original$mapping_species_count, "\n"))
  cat(paste("Common species across all datasets:", length(integrated_data$common_species), "\n\n"))
  
  cat("Common species included in analysis:\n")
  cat(paste(" -", integrated_data$common_species, collapse = "\n"))
  cat("\n\n")
  
  if(length(integrated_data$consistency_report$missing_in_tree) > 0) {
    cat("Species in chromosome counts but missing from tree:\n")
    cat(paste(" -", integrated_data$consistency_report$missing_in_tree, collapse = "\n"))
    cat("\n\n")
  }
  
  if(length(integrated_data$consistency_report$missing_in_counts) > 0) {
    cat("Species in tree but missing chromosome counts:\n")
    cat(paste(" -", integrated_data$consistency_report$missing_in_counts, collapse = "\n"))
    cat("\n\n")
  }
  
  if(length(integrated_data$consistency_report$missing_in_mapping) > 0) {
    cat("Species in tree and counts but missing from mapping data:\n")
    cat(paste(" -", integrated_data$consistency_report$missing_in_mapping, collapse = "\n"))
    cat("\n")
  }
  
  sink()
  
  # Save list of common species for reference
  common_species_file <- file.path(output_dir, "common_species_list.txt")
  writeLines(integrated_data$common_species, common_species_file)
  
  message(paste("Saved integrated data to:", output_dir))
  
  # Return file paths
  return(list(
    common_species_tree = pruned_tree_file,
    common_species_counts = filtered_counts_file,
    common_species_maps = filtered_maps_file,
    integrated_data = integrated_rdata_file,
    consistency_report = consistency_file,
    common_species_list = common_species_file
  ))
}

# ==============================
# Main Function
# ==============================

main <- function(args) {
  # Check arguments
  if(length(args) < 4) {
    cat("Usage: Rscript ancestral_chromosomes_phase1.R tree_file chr_counts_file maps_dir output_dir\n")
    cat("Example: Rscript ancestral_chromosomes_phase1.R phylogeny.nwk chromosome_counts.csv maps_dir/ results/\n")
    return(1)
  }
  
  tree_file <- args[1]
  chr_counts_file <- args[2]
  maps_dir <- args[3]
  output_dir <- args[4]
  
  # Run data integration process
  integrated_data <- integrate_data(tree_file, chr_counts_file, maps_dir)
  
  # Save results
  output_files <- save_integrated_data(integrated_data, output_dir)
  
  # Print basic summary
  cat("\n===== Data Integration Summary =====\n")
  cat("IMPORTANT: Analysis includes ONLY species present in ALL three datasets\n")
  cat(paste("Common species count:", length(integrated_data$common_species), "\n"))
  cat(paste("Chromosome count range:", min(integrated_data$chr_counts), "to", max(integrated_data$chr_counts), "\n"))
  cat(paste("Total mapping relationships:", nrow(integrated_data$bidirectional_maps), "\n"))
  cat(paste("Results saved to:", output_dir, "\n"))
  cat("===================================\n")
  
  # Return success
  return(0)
}

# Run the script if executed directly
if(!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  exit_code <- main(args)
  quit(status = exit_code)
}
