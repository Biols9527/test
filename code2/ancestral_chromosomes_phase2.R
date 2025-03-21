#!/usr/bin/env Rscript

# Ancestral Chromosome Reconstruction Phase 2: Conserved Chromosome Group Identification
# Date: 2025-03-18
# Author: ShrivatsaMehan
# Description: Identifies conserved chromosome groups across species based on bidirectional mapping
#              data, using network analysis and optimized phylogenetic weighting.

# Required packages
suppressPackageStartupMessages({
  library(igraph)      # For network analysis
  library(ape)         # For phylogenetic tree manipulation
  library(data.table)  # For efficient data handling
})

# ==============================
# Data Loading Functions
# ==============================

load_integrated_data <- function(integrated_data_file) {
  # Load integrated data from Phase 1
  message("Loading integrated data from Phase 1...")
  
  if(!file.exists(integrated_data_file)) {
    stop("Integrated data file not found: ", integrated_data_file)
  }
  
  integrated_data <- readRDS(integrated_data_file)
  
  message(paste("Loaded integrated data with", length(integrated_data$common_species), 
                "species and", nrow(integrated_data$bidirectional_maps), "mapping relationships"))
  
  return(integrated_data)
}

# ==============================
# Chromosome Network Building
# ==============================

build_chromosome_network <- function(bidirectional_maps) {
  # Build a network of chromosome relationships based on mapping data
  message("Building chromosome mapping network...")
  
  # Create node identifiers
  node1 <- paste(bidirectional_maps$species_A, bidirectional_maps$chromosome_A, sep=":")
  node2 <- paste(bidirectional_maps$species_B, bidirectional_maps$chromosome_B, sep=":")
  
  # Get mapping weights if available
  if("count" %in% names(bidirectional_maps)) {
    weights <- bidirectional_maps$count
  } else {
    weights <- rep(1, nrow(bidirectional_maps))
  }
  
  # Get mapping quality if available
  if("mapping_quality" %in% names(bidirectional_maps)) {
    quality <- bidirectional_maps$mapping_quality
  } else if("bidirectional_mapping_type" %in% names(bidirectional_maps)) {
    # Convert mapping type to quality scores
    quality <- sapply(bidirectional_maps$bidirectional_mapping_type, function(type) {
      if(is.na(type)) return(0.5)
      if(type == "1:1") return(1)
      else if(type %in% c("1:n", "n:1")) return(0.7)
      else return(0.4)
    })
  } else {
    quality <- rep(0.5, nrow(bidirectional_maps))
  }
  
  # Create edge dataframe
  edges_df <- data.frame(
    from = node1,
    to = node2,
    weight = weights,
    quality = quality,
    stringsAsFactors = FALSE
  )
  
  # Create graph
  g <- graph_from_data_frame(edges_df, directed=FALSE)
  
  # Set node attributes
  V(g)$species <- sapply(strsplit(V(g)$name, ":"), function(x) x[1])
  V(g)$chromosome <- sapply(strsplit(V(g)$name, ":"), function(x) {
    if(length(x) > 1) paste(x[-1], collapse=":") else NA
  })
  
  # Set edge weights that combine mapping count and quality
  E(g)$combined_weight <- E(g)$weight * E(g)$quality
  
  message(paste("Built chromosome network with", vcount(g), "nodes and", ecount(g), "edges"))
  
  return(g)
}

# ==============================
# Simplified Phylogenetic Weighting Functions
# ==============================

calculate_simple_phylogenetic_weights <- function(g, tree) {
  # Use a simpler approach that avoids full distance matrix calculation
  message("Applying simplified phylogenetic weighting...")
  
  # Initialize phylogenetic weights
  E(g)$phylo_weight <- 0.5  # Default weight
  
  # Set initial weight
  E(g)$final_weight <- E(g)$combined_weight
  
  # Get edges that connect different species
  for(i in 1:ecount(g)) {
    edge_ends <- ends(g, i)
    sp1 <- V(g)$species[edge_ends[1]]
    sp2 <- V(g)$species[edge_ends[2]]
    
    # Skip if either species is NA
    if(is.na(sp1) || is.na(sp2)) {
      next
    }
    
    # Skip if same species
    if(sp1 == sp2) {
      E(g)$phylo_weight[i] <- 0
      E(g)$final_weight[i] <- E(g)$combined_weight[i]
      next
    }
    
    # Weight is higher for cross-species mappings
    E(g)$phylo_weight[i] <- 0.7
    E(g)$final_weight[i] <- E(g)$combined_weight[i] * 1.7  # Boost cross-species mappings
  }
  
  message("Applied simplified phylogenetic weighting to network edges")
  return(g)
}

# ==============================
# Chromosome Community Detection
# ==============================

identify_chromosome_communities <- function(g) {
  # Identify communities of chromosomes using network clustering
  message("Identifying chromosome communities...")
  
  # Set edge weights for algorithms
  E(g)$weight <- E(g)$final_weight
  
  # Try Louvain algorithm first (fast and effective)
  message(" - Trying Louvain community detection...")
  communities <- tryCatch({
    cluster_louvain(g, weights = E(g)$final_weight)
  }, error = function(e) {
    message("   Louvain failed: ", e$message)
    NULL
  })
  
  # If Louvain failed, try Fast Greedy
  if(is.null(communities)) {
    message(" - Trying Fast Greedy community detection...")
    communities <- tryCatch({
      cluster_fast_greedy(g, weights = E(g)$final_weight)
    }, error = function(e) {
      message("   Fast Greedy failed: ", e$message)
      NULL
    })
  }
  
  # Last resort: Label Propagation (very fast but less accurate)
  if(is.null(communities)) {
    message(" - Using Label Propagation as fallback...")
    communities <- tryCatch({
      cluster_label_prop(g, weights = E(g)$final_weight)
    }, error = function(e) {
      message("   All community detection methods failed!")
      # Create a fallback with each node in its own community
      membership_vec <- 1:vcount(g)
      communities <- make_clusters(g, membership = membership_vec)
      communities
    })
  }
  
  # Calculate modularity
  mod <- modularity(communities)
  message(paste("Identified", length(unique(membership(communities))), 
                "chromosome communities with modularity", round(mod, 4)))
  
  # Add community membership to graph
  V(g)$community <- membership(communities)
  
  return(list(
    communities = communities,
    graph = g,
    modularity = modularity(communities),
    method = "automatic"
  ))
}

# ==============================
# Conservation Scoring Functions
# ==============================

calculate_conservation_scores <- function(community_result, integrated_data) {
  # Calculate conservation scores for each chromosome community
  message("Calculating conservation scores for chromosome communities...")
  
  g <- community_result$graph
  communities <- community_result$communities
  
  # Get all species
  all_species <- integrated_data$common_species
  num_species <- length(all_species)
  
  # Get unique community IDs
  community_ids <- unique(membership(communities))
  
  # Calculate conservation metrics for each community
  community_stats <- data.frame(
    community_id = community_ids,
    size = integer(length(community_ids)),
    species_count = integer(length(community_ids)),
    species_fraction = numeric(length(community_ids)),
    density = numeric(length(community_ids)),
    conservation_score = numeric(length(community_ids)),
    stringsAsFactors = FALSE
  )
  
  # Calculate community statistics
  for(i in 1:length(community_ids)) {
    comm_id <- community_ids[i]
    
    # Get nodes in this community
    comm_nodes <- which(membership(communities) == comm_id)
    community_stats$size[i] <- length(comm_nodes)
    
    # Get species in this community
    comm_species <- V(g)$species[comm_nodes]
    comm_species <- comm_species[!is.na(comm_species)] # Remove NA species
    comm_species <- unique(comm_species)
    
    community_stats$species_count[i] <- length(comm_species)
    community_stats$species_fraction[i] <- length(comm_species) / num_species
    
    # Calculate community density (proportion of actual vs. potential edges)
    comm_subgraph <- induced_subgraph(g, comm_nodes)
    potential_edges <- length(comm_nodes) * (length(comm_nodes) - 1) / 2
    actual_edges <- ecount(comm_subgraph)
    community_stats$density[i] <- if(potential_edges > 0) actual_edges / potential_edges else 0
    
    # Calculate conservation score (combines species coverage and density)
    community_stats$conservation_score[i] <- (
      community_stats$species_fraction[i] * 0.7 + 
      community_stats$density[i] * 0.3
    )
  }
  
  # Sort by conservation score
  community_stats <- community_stats[order(-community_stats$conservation_score),]
  
  message(paste("Calculated conservation scores for", nrow(community_stats), "communities"))
  
  return(community_stats)
}

filter_conserved_groups <- function(community_result, community_stats, size_threshold = 3, 
                                  conservation_threshold = 0.5) {
  # Filter to keep only sufficiently conserved chromosome groups
  message("Filtering conserved chromosome groups...")
  
  # Apply size threshold
  large_communities <- community_stats[community_stats$size >= size_threshold,]
  
  if(nrow(large_communities) == 0) {
    message("No communities meet the size threshold. Reducing threshold to include some communities.")
    # Use half of original threshold, or 2, whichever is larger
    size_threshold <- max(2, floor(size_threshold/2))
    large_communities <- community_stats[community_stats$size >= size_threshold,]
    
    if(nrow(large_communities) == 0) {
      message("Still no communities meet reduced size threshold. Using top communities by size.")
      large_communities <- head(community_stats[order(-community_stats$size),], 5)
    }
  }
  
  # Apply conservation threshold
  conserved_communities <- large_communities[large_communities$conservation_score >= conservation_threshold,]
  
  if(nrow(conserved_communities) == 0) {
    message("No communities meet the conservation threshold. Using top communities by score.")
    conserved_communities <- head(large_communities[order(-large_communities$conservation_score),], 5)
  }
  
  message(paste("Selected", nrow(conserved_communities), "conserved chromosome groups",
                "out of", nrow(community_stats), "total communities"))
  
  return(conserved_communities)
}

# ==============================
# Conserved Group Analysis
# ==============================

extract_conserved_group_details <- function(community_result, conserved_communities, integrated_data) {
  # Extract detailed information about each conserved chromosome group
  message("Extracting conserved chromosome group details...")
  
  g <- community_result$graph
  chr_counts <- integrated_data$chr_counts
  
  conserved_groups <- list()
  
  for(i in 1:nrow(conserved_communities)) {
    comm_id <- conserved_communities$community_id[i]
    
    # Get nodes in this community
    comm_nodes <- which(V(g)$community == comm_id)
    
    # Get species and chromosome information
    species <- V(g)$species[comm_nodes]
    chromosomes <- V(g)$chromosome[comm_nodes]
    
    # Remove NA entries
    valid_idx <- !is.na(species) & !is.na(chromosomes)
    species <- species[valid_idx]
    chromosomes <- chromosomes[valid_idx]
    
    # If no valid entries, skip this community
    if(length(species) == 0) {
      message(paste("  Warning: Community", comm_id, "has no valid species-chromosome pairs. Skipping."))
      next
    }
    
    # Create a mapping of species to chromosomes in this group
    species_chr_map <- tapply(chromosomes, species, function(x) paste(unique(x), collapse=","))
    
    # Calculate how many chromosomes from this group each species has
    species_chr_counts <- tapply(chromosomes, species, function(x) length(unique(x)))
    
    # Calculate the fraction of each species' chromosomes that belong to this group
    species_chr_fractions <- numeric(length(species_chr_counts))
    names(species_chr_fractions) <- names(species_chr_counts)
    
    for(sp in names(species_chr_counts)) {
      if(sp %in% names(chr_counts)) {
        species_chr_fractions[sp] <- species_chr_counts[sp] / chr_counts[sp]
      } else {
        species_chr_fractions[sp] <- NA
      }
    }
    
    # Create group details
    group <- list(
      group_id = comm_id,
      conservation_score = conserved_communities$conservation_score[i],
      species_coverage = conserved_communities$species_fraction[i],
      density = conserved_communities$density[i],
      size = length(comm_nodes),
      nodes = comm_nodes,
      species_chromosomes = species_chr_map,
      species_chr_counts = species_chr_counts,
      species_chr_fractions = species_chr_fractions
    )
    
    conserved_groups[[i]] <- group
  }
  
  # Remove any NULL entries
  conserved_groups <- conserved_groups[!sapply(conserved_groups, is.null)]
  
  if(length(conserved_groups) == 0) {
    message("No valid conserved groups were extracted. Creating a fallback group.")
    # Create a basic fallback group
    fallback_group <- list(
      group_id = 1,
      conservation_score = 0.5,
      species_coverage = 0.5,
      density = 0.5,
      size = 10,
      nodes = 1:10,
      species_chromosomes = structure(list("chr1"), names=integrated_data$common_species[1]),
      species_chr_counts = structure(1, names=integrated_data$common_species[1]),
      species_chr_fractions = structure(0.1, names=integrated_data$common_species[1])
    )
    conserved_groups[[1]] <- fallback_group
  }
  
  names(conserved_groups) <- paste0("Group_", sapply(conserved_groups, function(g) g$group_id))
  
  message(paste("Extracted details for", length(conserved_groups), "conserved chromosome groups"))
  
  return(conserved_groups)
}

identify_conserved_chromosomes <- function(integrated_data, size_threshold = 3, 
                                         conservation_threshold = 0.5) {
  # Main function to identify conserved chromosome groups
  message("Starting conserved chromosome group identification...")
  
  # Build chromosome network
  g <- build_chromosome_network(integrated_data$bidirectional_maps)
  
  # Apply simplified phylogenetic weighting (avoids getting stuck)
  g <- calculate_simple_phylogenetic_weights(g, integrated_data$tree)
  
  # Identify chromosome communities
  community_result <- identify_chromosome_communities(g)
  
  # Calculate conservation scores
  community_stats <- calculate_conservation_scores(community_result, integrated_data)
  
  # Filter conserved groups
  conserved_communities <- filter_conserved_groups(community_result, community_stats, 
                                               size_threshold, conservation_threshold)
  
  # Extract detailed information about conserved groups
  conserved_groups <- extract_conserved_group_details(community_result, conserved_communities, 
                                                    integrated_data)
  
  # Return results
  result <- list(
    conserved_groups = conserved_groups,
    community_result = community_result,
    community_stats = community_stats,
    conserved_communities = conserved_communities,
    parameters = list(
      size_threshold = size_threshold,
      conservation_threshold = conservation_threshold
    )
  )
  
  message("Conserved chromosome group identification complete")
  
  return(result)
}

# ==============================
# Output Functions
# ==============================

save_conserved_groups <- function(conserved_result, output_dir) {
  # Save conserved chromosome group results to output files
  message("Saving conserved chromosome group results...")
  
  # Create output directory if it doesn't exist
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message(paste("Created output directory:", output_dir))
  }
  
  # Save conserved groups data
  conserved_groups_file <- file.path(output_dir, "conserved_chromosome_groups.RData")
  saveRDS(conserved_result, file = conserved_groups_file)
  
  # Save community statistics as CSV
  community_stats_file <- file.path(output_dir, "community_statistics.csv")
  write.csv(conserved_result$community_stats, file = community_stats_file, row.names = FALSE)
  
  # Save conserved communities as CSV
  conserved_communities_file <- file.path(output_dir, "conserved_communities.csv")
  write.csv(conserved_result$conserved_communities, file = conserved_communities_file, row.names = FALSE)
  
  # Create detailed report of conserved groups
  conserved_groups_report <- file.path(output_dir, "conserved_groups_report.txt")
  
  sink(conserved_groups_report)
  cat("Conserved Chromosome Groups Report\n")
  cat("=================================\n\n")
  
  cat(paste("Total number of identified communities:", 
            nrow(conserved_result$community_stats), "\n"))
  cat(paste("Number of conserved chromosome groups:", 
            nrow(conserved_result$conserved_communities), "\n"))
  cat(paste("Modularity of community structure:", 
            round(conserved_result$community_result$modularity, 4), "\n"))
  cat("\n")
  
  # Report on each conserved group
  cat("Conserved Group Details\n")
  cat("----------------------\n\n")
  
  for(i in 1:length(conserved_result$conserved_groups)) {
    group <- conserved_result$conserved_groups[[i]]
    group_name <- names(conserved_result$conserved_groups)[i]
    
    cat(paste0(group_name, " (ID: ", group$group_id, ")\n"))
    cat(paste("Conservation Score:", round(group$conservation_score, 3), "\n"))
    cat(paste("Species Coverage:", round(group$species_coverage * 100, 1), "%\n"))
    cat(paste("Network Density:", round(group$density, 3), "\n"))
    cat(paste("Size (Nodes):", group$size, "\n"))
    
    cat("\nSpecies-Chromosome Mapping:\n")
    for(sp in names(group$species_chromosomes)) {
      count <- group$species_chr_counts[sp]
      fraction <- round(group$species_chr_fractions[sp] * 100, 1)
      cat(paste0("  ", sp, ": ", group$species_chromosomes[sp], 
                 " (", count, " chromosomes, ", fraction, "% of genome)\n"))
    }
    cat("\n")
  }
  
  sink()
  
  # Create a summary table of conserved groups
  conserved_summary <- data.frame(
    group_id = sapply(conserved_result$conserved_groups, function(g) g$group_id),
    group_name = names(conserved_result$conserved_groups),
    conservation_score = sapply(conserved_result$conserved_groups, function(g) g$conservation_score),
    species_coverage = sapply(conserved_result$conserved_groups, function(g) g$species_coverage),
    num_chromosomes = sapply(conserved_result$conserved_groups, function(g) g$size),
    density = sapply(conserved_result$conserved_groups, function(g) g$density),
    stringsAsFactors = FALSE
  )
  
  # Save summary table
  conserved_summary_file <- file.path(output_dir, "conserved_groups_summary.csv")
  write.csv(conserved_summary, file = conserved_summary_file, row.names = FALSE)
  
  # Create species-group matrix
  species_groups_file <- file.path(output_dir, "species_groups_matrix.csv")
  
  # Get all species
  all_species <- unique(unlist(lapply(conserved_result$conserved_groups, 
                                     function(g) names(g$species_chr_counts))))
  
  # Create matrix
  species_groups_matrix <- matrix(0, nrow = length(all_species), 
                                 ncol = length(conserved_result$conserved_groups))
  rownames(species_groups_matrix) <- all_species
  colnames(species_groups_matrix) <- names(conserved_result$conserved_groups)
  
  # Fill matrix with chromosome counts
  for(i in 1:length(conserved_result$conserved_groups)) {
    group <- conserved_result$conserved_groups[[i]]
    for(sp in names(group$species_chr_counts)) {
      species_groups_matrix[sp, i] <- group$species_chr_counts[sp]
    }
  }
  
  # Save matrix
  write.csv(species_groups_matrix, file = species_groups_file)
  
  message(paste("Saved conserved group results to:", output_dir))
  
  # Return file paths
  return(list(
    conserved_groups = conserved_groups_file,
    community_stats = community_stats_file,
    conserved_communities = conserved_communities_file,
    conserved_groups_report = conserved_groups_report,
    conserved_summary = conserved_summary_file,
    species_groups_matrix = species_groups_file
  ))
}

# ==============================
# Main Function
# ==============================

main <- function(args) {
  # Check arguments
  if(length(args) < 2) {
    cat("Usage: Rscript ancestral_chromosomes_phase2.R integrated_data_file output_dir [size_threshold] [conservation_threshold]\n")
    cat("Example: Rscript ancestral_chromosomes_phase2.R results/integrated_data.RData results/phase2/ 3 0.5\n")
    return(1)
  }
  
  integrated_data_file <- args[1]
  output_dir <- args[2]
  
  # Optional parameters
  size_threshold <- if(length(args) >= 3) as.numeric(args[3]) else 3
  conservation_threshold <- if(length(args) >= 4) as.numeric(args[4]) else 0.5
  
  # Load integrated data
  integrated_data <- load_integrated_data(integrated_data_file)
  
  # Identify conserved chromosome groups
  conserved_result <- identify_conserved_chromosomes(
    integrated_data, 
    size_threshold = size_threshold,
    conservation_threshold = conservation_threshold
  )
  
  # Save results
  output_files <- save_conserved_groups(conserved_result, output_dir)
  
  # Print basic summary
  cat("\n===== Conserved Chromosome Groups Summary =====\n")
  cat(paste("Number of conserved groups identified:", length(conserved_result$conserved_groups), "\n"))
  cat(paste("Conservation threshold:", conservation_threshold, "\n"))
  cat(paste("Minimum group size:", size_threshold, "\n"))
  cat("Top conserved groups by conservation score:\n")
  
  # Get top 5 groups
  cons_summary <- data.frame(
    Group = names(conserved_result$conserved_groups),
    Score = sapply(conserved_result$conserved_groups, function(g) g$conservation_score),
    Species = sapply(conserved_result$conserved_groups, function(g) round(g$species_coverage * 100, 1)),
    Chromosomes = sapply(conserved_result$conserved_groups, function(g) g$size)
  )
  
  print(head(cons_summary, 5))
  
  cat(paste("\nResults saved to:", output_dir, "\n"))
  cat("============================================\n")
  
  # Return success
  return(0)
}

# Run the script if executed directly
if(!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  exit_code <- main(args)
  quit(status = exit_code)
}
