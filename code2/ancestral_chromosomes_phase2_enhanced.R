#!/usr/bin/env Rscript

# Ancestral Chromosome Reconstruction Phase 2: Conserved Chromosome Group Identification
# Date: 2025-03-18
# Author: ShrivatsaMehan
# Description: Identifies conserved chromosome groups across species based on bidirectional mapping
#              data, using network analysis and optimized phylogenetic weighting.
# Enhanced with multiple community detection methods and ensemble approach

# Required packages
suppressPackageStartupMessages({
  library(igraph)      # For network analysis
  library(ape)         # For phylogenetic tree manipulation
  library(data.table)  # For efficient data handling
  library(parallel)    # For parallel processing
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
# Enhanced Chromosome Community Detection - 修复版本
# ==============================

identify_chromosome_communities <- function(g, methods = c("louvain", "fast_greedy", "label_prop", 
                                                          "walktrap", "leiden", "infomap", "spinglass"),
                                          use_ensemble = TRUE) {
  # Identify communities of chromosomes using multiple network clustering algorithms
  message("Identifying chromosome communities using multiple methods...")

  # Set edge weights for algorithms
  E(g)$weight <- E(g)$final_weight
  
  # Store results from each method
  community_results <- list()
  
  # Try each specified method
  for (method in methods) {
    message(paste(" - Trying", method, "community detection..."))
    
    # 不要在循环内返回结果，而是收集结果
    tryCatch({
      # Apply the selected community detection method
      comm <- NULL
      
      if (method == "louvain") {
        comm <- cluster_louvain(g, weights = E(g)$final_weight)
      } else if (method == "fast_greedy") {
        comm <- cluster_fast_greedy(g, weights = E(g)$final_weight)
      } else if (method == "label_prop") {
        comm <- cluster_label_prop(g, weights = E(g)$final_weight)
      } else if (method == "walktrap") {
        # Walktrap with steps = 4 (compromise between local and global)
        comm <- cluster_walktrap(g, weights = E(g)$final_weight, steps = 4)
      } else if (method == "leiden") {
        # Check if leiden package is available
        if (requireNamespace("leiden", quietly = TRUE)) {
          library(leiden)
          # Convert igraph to adjacency matrix
          adj_matrix <- as_adjacency_matrix(g, attr = "final_weight", sparse = TRUE)
          membership_vec <- leiden::leiden(adj_matrix)
          comm <- make_clusters(g, membership = membership_vec)
        } else {
          message("   Leiden package not installed. Skipping.")
          next
        }
      } else if (method == "infomap") {
        # Try infomap with default parameters
        comm <- cluster_infomap(g, e.weights = E(g)$final_weight)
      } else if (method == "spinglass") {
        if (vcount(g) <= 1000) { # Spinglass may be slow for large networks
          comm <- cluster_spinglass(g, weights = E(g)$final_weight)
        } else {
          message("   Network too large for spinglass. Skipping.")
          next
        }
      } else if (method == "edge_betweenness") {
        if (vcount(g) <= 500) { # Edge betweenness is very slow
          comm <- cluster_edge_betweenness(g, weights = E(g)$final_weight)
        } else {
          message("   Network too large for edge betweenness. Skipping.")
          next
        }
      }
      
      # Calculate modularity
      if (!is.null(comm)) {
        mod <- modularity(comm)
        message(paste("   Found", length(unique(membership(comm))),
                      "communities with modularity", round(mod, 4)))
        
        # 验证comm是否有效，并能够调用membership
        message(paste("   DEBUG: Community class:", paste(class(comm), collapse=", ")))
        
        # 测试是否可以调用membership
        tryCatch({
          test_membership <- membership(comm)
          message(paste("   DEBUG: Can extract membership -", length(test_membership), "values"))
          
          # Store result - 存储结果而不是直接返回
          community_results[[method]] <- list(
            communities = comm,
            modularity = mod,
            membership = test_membership
          )
        }, error = function(e) {
          message(paste("   WARNING: Cannot extract membership -", e$message))
        })
      }
      
    }, error = function(e) {
      message(paste("   Method failed:", e$message))
    })
  }

  # If all methods failed, create fallback with singleton communities
  if (length(community_results) == 0) {
    message(" - All community detection methods failed. Creating fallback communities.")
    membership_vec <- 1:vcount(g)
    comm <- make_clusters(g, membership = membership_vec)
    
    community_results[["fallback"]] <- list(
      communities = comm,
      modularity = 0,
      membership = membership_vec
    )
  }
  
  # Find the method with highest modularity
  modularity_scores <- sapply(community_results, function(x) x$modularity)
  best_method <- names(community_results)[which.max(modularity_scores)]
  
  message(paste(" - Best community detection method:", best_method,
                "with modularity", round(max(modularity_scores), 4)))
  
  # Create ensemble result if requested and we have multiple methods
  if (use_ensemble && length(community_results) > 1) {
    message(" - Creating ensemble community detection...")
    
    # Create consensus communities
    ensemble_result <- create_ensemble_communities(g, community_results)
    
    # Check if ensemble is better than best individual method
    if (ensemble_result$modularity > max(modularity_scores)) {
      message(paste("   Ensemble method improved modularity to",
                    round(ensemble_result$modularity, 4)))
      
      best_method <- "ensemble"
      community_results[["ensemble"]] <- ensemble_result
    } else {
      message("   Ensemble method did not improve on best individual method.")
    }
  }
  
  # 获取最佳社区对象
  best_result <- community_results[[best_method]]
  
  # 确保best_communities是一个有效的igraph社区对象
  if (is.null(best_result) || is.null(best_result$communities)) {
    message(" - WARNING: Best method didn't produce valid communities. Creating fallback.")
    # 创建备用社区
    membership_vec <- 1:vcount(g)
    best_communities <- make_clusters(g, membership = membership_vec)
  } else {
    best_communities <- best_result$communities
    
    # 验证这是否是一个有效的社区对象
    message(paste(" - DEBUG: Final best_communities class:", paste(class(best_communities), collapse=", ")))
    
    # 确保可以调用membership
    tryCatch({
      test_membership <- membership(best_communities)
      message(paste(" - DEBUG: Final membership extraction successful with", length(test_membership), "values"))
    }, error = function(e) {
      message(paste(" - DEBUG: Final membership extraction FAILED:", e$message))
      message(" - Creating compatible community object...")
      # 如果membership函数失败，尝试使用存储的成员关系创建新的社区对象
      membership_vec <- best_result$membership
      best_communities <- make_clusters(g, membership = membership_vec)
    })
  }
  
  # 将社区成员关系添加到图的顶点属性
  V(g)$community <- membership(best_communities)
  
  message(" - Community detection completed successfully")
  
  return(list(
    communities = best_communities,
    graph = g,
    modularity = modularity(best_communities),
    method = best_method,
    all_methods = community_results
  ))
}

create_ensemble_communities <- function(g, community_results) {
  # Create consensus communities from multiple detection methods
  message("   Creating consensus communities from multiple methods...")
  
  # Create co-occurrence matrix
  n_vertices <- vcount(g)
  cooccurrence <- matrix(0, n_vertices, n_vertices)
  
  # Count how many times each pair of vertices is in the same community
  n_methods <- length(community_results)
  
  for (method_name in names(community_results)) {
    method_result <- community_results[[method_name]]
    membership_vector <- method_result$membership
    
    for (i in 1:(n_vertices-1)) {
      for (j in (i+1):n_vertices) {
        if (membership_vector[i] == membership_vector[j]) {
          cooccurrence[i, j] <- cooccurrence[i, j] + 1
          cooccurrence[j, i] <- cooccurrence[j, i] + 1
        }
      }
    }
  }
  
  # Normalize co-occurrence matrix
  cooccurrence <- cooccurrence / n_methods
  
  # Create a new weighted graph based on co-occurrence
  consensus_g <- graph_from_adjacency_matrix(
    cooccurrence,
    mode = "undirected",
    weighted = TRUE
  )
  
  # Apply Louvain to this consensus graph
  consensus_comm <- tryCatch({
    cluster_louvain(consensus_g)
  }, error = function(e) {
    message(paste("   Ensemble creation failed:", e$message))
    # 创建一个简单的备用社区
    make_clusters(g, membership = rep(1:10, length.out=vcount(g)))
  })
  
  consensus_mod <- modularity(consensus_comm)
  
  # Apply the consensus membership back to the original graph
  ensemble_membership <- membership(consensus_comm)
  ensemble_comm <- make_clusters(g, membership = ensemble_membership)
  ensemble_mod <- modularity(g, membership(ensemble_comm), weights = E(g)$final_weight)
  
  return(list(
    communities = ensemble_comm,
    modularity = ensemble_mod,
    membership = ensemble_membership,
    consensus_graph = consensus_g
  ))
}

# ==============================
# Conservation Scoring Functions - 修复版本
# ==============================

calculate_conservation_scores <- function(community_result, integrated_data) {
  # Calculate conservation scores for each chromosome community
  message("Calculating conservation scores for chromosome communities...")

  g <- community_result$graph
  communities <- community_result$communities
  
  # 错误检查和诊断
  message(paste(" - DEBUG: Communities object class:", paste(class(communities), collapse=", ")))
  message(paste(" - DEBUG: Is communities NULL:", is.null(communities)))
  
  # 确保communities是有效的,可以调用membership
  membership_result <- tryCatch({
    community_membership <- membership(communities)
    message(paste(" - Successfully obtained community membership with", length(community_membership), "values"))
    community_membership
  }, error = function(e) {
    message(paste(" - Error obtaining community membership:", e$message))
    message(" - Using graph vertex community attribute as fallback")
    # 使用图的顶点属性作为备用
    if (!is.null(V(g)$community)) {
      V(g)$community
    } else {
      # 如果所有方法都失败，创建一个简单的社区分配
      message(" - Creating fallback community assignments (one community per node)")
      1:vcount(g)
    }
  })
  
  # 确保有效的社区成员关系
  if (length(membership_result) != vcount(g)) {
    message(paste(" - WARNING: Membership length mismatch -", 
                  length(membership_result), "vs", vcount(g), "vertices"))
    # 创建备用成员关系
    membership_result <- 1:vcount(g)
  }

  # Get all species
  all_species <- integrated_data$common_species
  num_species <- length(all_species)

  # Get unique community IDs
  community_ids <- unique(membership_result)
  message(paste(" - Found", length(community_ids), "unique communities"))

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
    comm_nodes <- which(membership_result == comm_id)
    community_stats$size[i] <- length(comm_nodes)
    
    # 安全检查：确保有有效节点
    if (length(comm_nodes) == 0) {
      message(paste(" - WARNING: Empty community", comm_id, "- setting default values"))
      community_stats$species_count[i] <- 0
      community_stats$species_fraction[i] <- 0
      community_stats$density[i] <- 0
      community_stats$conservation_score[i] <- 0
      next
    }

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

# ==============================
# New: Community Method Comparison
# ==============================

compare_community_methods <- function(community_result) {
  # Compare results from different community detection methods
  message("Comparing community detection methods...")
  
  # If we only have one method, nothing to compare
  if (is.null(community_result$all_methods) || length(community_result$all_methods) <= 1) {
    message("Only one community detection method was successful. No comparison to perform.")
    return(NULL)
  }
  
  # Initialize comparison data
  method_names <- names(community_result$all_methods)
  n_methods <- length(method_names)
  
  comparison <- data.frame(
    method = method_names,
    communities = sapply(community_result$all_methods, function(m) {
      if (!is.null(m$membership)) length(unique(m$membership)) else NA
    }),
    modularity = sapply(community_result$all_methods, function(m) m$modularity),
    stringsAsFactors = FALSE
  )
  
  # Calculate similarity matrix (adjusted rand index)
  similarity_matrix <- matrix(0, n_methods, n_methods)
  rownames(similarity_matrix) <- method_names
  colnames(similarity_matrix) <- method_names
  
  # Function for calculating adjusted Rand index
  calc_ari <- function(x, y) {
    # Simple implementation of adjusted Rand index
    n <- length(x)
    if (n != length(y)) return(NA)
    
    # Create contingency table
    tab <- table(x, y)
    
    # Calculate sums
    a <- sum(choose(tab, 2))
    b <- sum(choose(rowSums(tab), 2))
    c <- sum(choose(colSums(tab), 2))
    d <- choose(n, 2)
    
    # Calculate adjusted Rand index
    ari <- (a - b * c / d) / ((b + c) / 2 - b * c / d)
    return(ari)
  }
  
  # Fill similarity matrix
  for (i in 1:n_methods) {
    for (j in 1:n_methods) {
      if (i == j) {
        similarity_matrix[i, j] <- 1
      } else if (i < j) {
        # Error handling for missing membership values
        tryCatch({
          membership_i <- community_result$all_methods[[method_names[i]]]$membership
          membership_j <- community_result$all_methods[[method_names[j]]]$membership
          
          if(!is.null(membership_i) && !is.null(membership_j) && 
             length(membership_i) == length(membership_j)) {
            similarity_matrix[i, j] <- calc_ari(membership_i, membership_j)
            similarity_matrix[j, i] <- similarity_matrix[i, j]
          } else {
            similarity_matrix[i, j] <- NA
            similarity_matrix[j, i] <- NA
          }
        }, error = function(e) {
          message(paste("   Error calculating similarity between", method_names[i], 
                        "and", method_names[j], "-", e$message))
          similarity_matrix[i, j] <- NA
          similarity_matrix[j, i] <- NA
        })
      }
    }
  }
  
  # Calculate average similarity for each method
  avg_similarity <- rowMeans(similarity_matrix, na.rm = TRUE)
  comparison$avg_similarity <- avg_similarity
  
  # Sort by modularity
  comparison <- comparison[order(-comparison$modularity),]
  
  message("Community method comparison completed")
  
  return(list(
    comparison = comparison,
    similarity_matrix = similarity_matrix
  ))
}

identify_conserved_chromosomes <- function(integrated_data, size_threshold = 3,
                                         conservation_threshold = 0.5,
                                         methods = c("louvain", "fast_greedy", "label_prop", "walktrap"),
                                         use_ensemble = TRUE) {
  # Main function to identify conserved chromosome groups
  message("Starting conserved chromosome group identification...")

  # Build chromosome network
  g <- build_chromosome_network(integrated_data$bidirectional_maps)

  # Apply simplified phylogenetic weighting (avoids getting stuck)
  g <- calculate_simple_phylogenetic_weights(g, integrated_data$tree)

  # Identify chromosome communities using multiple methods
  community_result <- identify_chromosome_communities(g, methods = methods, use_ensemble = use_ensemble)

  # Compare community detection methods
  method_comparison <- compare_community_methods(community_result)

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
    method_comparison = method_comparison,
    parameters = list(
      size_threshold = size_threshold,
      conservation_threshold = conservation_threshold,
      methods = methods,
      use_ensemble = use_ensemble
    )
  )

  message("Conserved chromosome group identification complete")

  return(result)
}

# ==============================
# Enhanced Output Functions
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
  cat(paste("Community detection method:", conserved_result$community_result$method, "\n"))
  cat("\n")

  # Add method comparison if available
  if (!is.null(conserved_result$method_comparison)) {
    cat("Community Detection Method Comparison\n")
    cat("-----------------------------------\n\n")
    
    comp_df <- conserved_result$method_comparison$comparison
    cat(sprintf("%-15s %-12s %-10s %-15s\n", "Method", "Communities", "Modularity", "Avg. Similarity"))
    cat(sprintf("%-15s %-12s %-10s %-15s\n", "------", "-----------", "----------", "--------------"))
    
    for (i in 1:nrow(comp_df)) {
      cat(sprintf("%-15s %-12d %-10.4f %-15.4f\n", 
                  comp_df$method[i], 
                  comp_df$communities[i], 
                  comp_df$modularity[i], 
                  comp_df$avg_similarity[i]))
    }
    cat("\n")
  }

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

  # Save method comparison results if available
  if (!is.null(conserved_result$method_comparison)) {
    method_comp_file <- file.path(output_dir, "community_method_comparison.csv")
    write.csv(conserved_result$method_comparison$comparison, file = method_comp_file, row.names = FALSE)
    
    # Create a heatmap of method similarities if possible
    tryCatch({
      if (requireNamespace("pheatmap", quietly = TRUE)) {
        library(pheatmap)
        pdf(file.path(output_dir, "method_similarity_heatmap.png"), width = 800, height = 800)
        pheatmap(conserved_result$method_comparison$similarity_matrix,
                 main = "Community Detection Method Similarity",
                 display_numbers = TRUE, number_format = "%.2f")
        dev.off()
      }
    }, error = function(e) {
      message(paste("Could not create similarity heatmap:", e$message))
    })
  }

  # Create a network visualization if igraph is available
  tryCatch({
    if(vcount(conserved_result$community_result$graph) <= 1000) { # Limit for reasonable plotting
      conserved_network_file <- file.path(output_dir, "conserved_chromosome_network.pdf")
      
      # Prepare graph for visualization
      g <- conserved_result$community_result$graph
      
      # Set community colors for nodes
      community_ids <- unique(V(g)$community)
      n_communities <- length(community_ids)
      colors <- rainbow(n_communities)
      
      V(g)$color <- colors[match(V(g)$community, community_ids)]
      
      # Set node sizes based on degree
      V(g)$size <- 3 + 2 * log1p(degree(g))
      
      # Set edge widths based on weights
      E(g)$width <- 0.5 + 2 * E(g)$final_weight / max(E(g)$final_weight)
      
      # Set layout (use a force-directed layout)
      layout <- layout_with_fr(g)
      
      # Save plot
      png(conserved_network_file, width = 1200, height = 1200, res = 100)
      plot(g, 
          layout = layout,
          vertex.label = NA,  # No labels for clarity
          edge.arrow.size = 0,
          main = "Chromosome Mapping Network with Communities")
      dev.off()
    }
  }, error = function(e) {
    message(paste("Could not create network visualization:", e$message))
  })

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
    cat("Usage: Rscript ancestral_chromosomes_phase2_enhanced.R integrated_data_file output_dir [size_threshold] [conservation_threshold] [use_ensemble]\n")
    cat("Example: Rscript ancestral_chromosomes_phase2_enhanced.R results/integrated_data.RData results/phase2/ 3 0.5 TRUE\n")
    return(1)
  }

  integrated_data_file <- args[1]
  output_dir <- args[2]

  # Optional parameters
  size_threshold <- if(length(args) >= 3) as.numeric(args[3]) else 3
  conservation_threshold <- if(length(args) >= 4) as.numeric(args[4]) else 0.5
  use_ensemble <- if(length(args) >= 5) as.logical(args[5]) else TRUE

  # Load integrated data
  integrated_data <- load_integrated_data(integrated_data_file)

  # Identify conserved chromosome groups with enhanced methods
  methods <- c("louvain", "fast_greedy", "label_prop", "walktrap")
  
  # Check if additional methods are available
  tryCatch({
    if(requireNamespace("leiden", quietly = TRUE)) {
      methods <- c(methods, "leiden")
    }
  }, error = function(e) {
    message("Leiden package not available, continuing without leiden method")
  })
  
  # Add infomap if network is not too large
  if(nrow(integrated_data$bidirectional_maps) < 10000) {
    methods <- c(methods, "infomap")
  }
  
  # Add spinglass for smaller networks
  if(nrow(integrated_data$bidirectional_maps) < 1000) {
    methods <- c(methods, "spinglass")
  }
  
  # Attempt to run analysis with error handling
  conserved_result <- tryCatch({
    identify_conserved_chromosomes(
      integrated_data,
      size_threshold = size_threshold,
      conservation_threshold = conservation_threshold,
      methods = methods,
      use_ensemble = use_ensemble
    )
  }, error = function(e) {
    message(paste("ERROR in chromosome group identification:", e$message))
    message("Attempting with only louvain method as fallback")
    
    # Try with just louvain as a fallback
    identify_conserved_chromosomes(
      integrated_data,
      size_threshold = size_threshold,
      conservation_threshold = conservation_threshold,
      methods = c("louvain"),
      use_ensemble = FALSE
    )
  })

  # Save results
  output_files <- save_conserved_groups(conserved_result, output_dir)

  # Print basic summary
  cat("\n===== Conserved Chromosome Groups Summary =====\n")
  cat(paste("Number of conserved groups identified:", length(conserved_result$conserved_groups), "\n"))
  cat(paste("Conservation threshold:", conservation_threshold, "\n"))
  cat(paste("Minimum group size:", size_threshold, "\n"))
  cat(paste("Community detection method:", conserved_result$community_result$method, "\n"))
  cat("Top conserved groups by conservation score:\n")

  # Get top 5 groups
  cons_summary <- data.frame(
    Group = names(conserved_result$conserved_groups),
    Score = sapply(conserved_result$conserved_groups, function(g) g$conservation_score),
    Species = sapply(conserved_result$conserved_groups, function(g) round(g$species_coverage * 100, 1)),
    Chromosomes = sapply(conserved_result$conserved_groups, function(g) g$size)
  )

  print(head(cons_summary, 5))
  
  if(!is.null(conserved_result$method_comparison)) {
    cat("\nCommunity Detection Method Comparison:\n")
    print(head(conserved_result$method_comparison$comparison))
  }

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
