#!/usr/bin/env Rscript

# Ancestral Chromosome Reconstruction Phase 3: Shared Chromosome Event Detection
# Date: 2025-03-18
# Author: ShrivatsaMehan
# Description: Identifies shared chromosome evolutionary events (fusions, fissions, etc.)
#              across the phylogenetic tree based on conserved groups and mapping data.

# Required packages
suppressPackageStartupMessages({
  library(ape)         # For phylogenetic tree manipulation
  library(data.table)  # For efficient data handling
  library(phangorn)    # For phylogenetic analysis
  library(igraph)      # For graph analysis
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

load_conserved_groups <- function(conserved_groups_file) {
  # Load conserved chromosome group data from Phase 2
  message("Loading conserved chromosome groups from Phase 2...")
  
  if(!file.exists(conserved_groups_file)) {
    stop("Conserved groups file not found: ", conserved_groups_file)
  }
  
  conserved_result <- readRDS(conserved_groups_file)
  
  message(paste("Loaded", length(conserved_result$conserved_groups), "conserved chromosome groups"))
  
  return(conserved_result)
}

# ==============================
# Tree Analysis Functions
# ==============================

extract_tree_edges <- function(tree) {
  # Extract edge information from the phylogenetic tree
  message("Extracting tree edges and nodes...")
  
  # Get number of species (tips)
  n_tips <- length(tree$tip.label)
  
  # Initialize edge data
  edge_data <- data.frame(
    edge_id = 1:nrow(tree$edge),
    parent_node_id = tree$edge[, 1],
    child_node_id = tree$edge[, 2],
    parent_is_tip = tree$edge[, 1] <= n_tips,
    child_is_tip = tree$edge[, 2] <= n_tips,
    edge_length = tree$edge.length,
    stringsAsFactors = FALSE
  )
  
  # Add node labels
  edge_data$parent_label <- sapply(edge_data$parent_node_id, function(node_id) {
    if(node_id <= n_tips) {
      tree$tip.label[node_id]
    } else {
      if(!is.null(tree$node.label)) {
        tree$node.label[node_id - n_tips]
      } else {
        paste0("Node", node_id)
      }
    }
  })
  
  edge_data$child_label <- sapply(edge_data$child_node_id, function(node_id) {
    if(node_id <= n_tips) {
      tree$tip.label[node_id]
    } else {
      if(!is.null(tree$node.label)) {
        tree$node.label[node_id - n_tips]
      } else {
        paste0("Node", node_id)
      }
    }
  })
  
  message(paste("Extracted", nrow(edge_data), "edges from phylogenetic tree"))
  
  return(edge_data)
}

identify_clades <- function(tree) {
  # Identify all clades in the phylogenetic tree
  message("Identifying clades in phylogenetic tree...")
  
  n_tips <- length(tree$tip.label)
  n_nodes <- tree$Nnode
  
  # Initialize clade data
  clade_data <- list()
  
  # For each internal node
  for(node_id in (n_tips + 1):(n_tips + n_nodes)) {
    # Get all tips descended from this node
    tips <- tryCatch({
      geiger::tips(tree, node_id)
    }, error = function(e) {
      # Alternative method if geiger is not available
      descendants <- getDescendants(tree, node_id)
      descendants <- descendants[descendants <= n_tips]
      tree$tip.label[descendants]
    })
    
    # Get node label
    if(!is.null(tree$node.label)) {
      node_label <- tree$node.label[node_id - n_tips]
    } else {
      node_label <- paste0("Node", node_id)
    }
    
    # Store clade information
    clade_data[[as.character(node_id)]] <- list(
      node_id = node_id,
      node_label = node_label,
      tips = tips,
      tip_count = length(tips)
    )
  }
  
  message(paste("Identified", length(clade_data), "clades in phylogenetic tree"))
  
  return(clade_data)
}

# Helper function to get descendants of a node
getDescendants <- function(tree, node) {
  # Recursive function to get all descendants of a node
  descendants <- c()
  edges <- which(tree$edge[, 1] == node)
  
  for (edge in edges) {
    child <- tree$edge[edge, 2]
    descendants <- c(descendants, child)
    
    if (child > length(tree$tip.label)) {
      descendants <- c(descendants, getDescendants(tree, child))
    }
  }
  
  return(unique(descendants))
}

# ==============================
# Chromosome Event Detection Functions
# ==============================

detect_count_changes <- function(edge_data, chr_counts) {
  # Detect chromosome count changes across tree edges
  message("Detecting chromosome count changes across tree edges...")
  
  # Add chromosome count data for tips
  tip_species <- unique(edge_data$child_label[edge_data$child_is_tip])
  species_with_counts <- intersect(tip_species, names(chr_counts))
  
  # Check coverage
  if(length(species_with_counts) < length(tip_species)) {
    missing_species <- setdiff(tip_species, species_with_counts)
    warning(paste("Missing chromosome counts for", length(missing_species), "species"))
  }
  
  # Initialize count change data
  count_changes <- data.frame(
    edge_id = edge_data$edge_id,
    parent_node_id = edge_data$parent_node_id,
    child_node_id = edge_data$child_node_id,
    parent_label = edge_data$parent_label,
    child_label = edge_data$child_label,
    parent_count = NA_integer_,
    child_count = NA_integer_,
    count_diff = NA_integer_,
    count_ratio = NA_real_,
    count_change_type = NA_character_,
    stringsAsFactors = FALSE
  )
  
  # Fill in known tip counts
  for(i in 1:nrow(count_changes)) {
    # If child is tip, get its chromosome count
    if(edge_data$child_is_tip[i] && edge_data$child_label[i] %in% names(chr_counts)) {
      count_changes$child_count[i] <- chr_counts[edge_data$child_label[i]]
    }
  }
  
  message("Recorded chromosome counts for tips")
  
  return(count_changes)
}

detect_mapping_patterns <- function(edge_data, integrated_data, conserved_groups) {
  # Detect mapping patterns that suggest chromosome events
  message("Analyzing chromosome mapping patterns...")
  
  # Extract bidirectional maps
  maps <- integrated_data$bidirectional_maps
  species_pairs <- data.frame(
    species_A = maps$species_A,
    species_B = maps$species_B,
    stringsAsFactors = FALSE
  )
  
  # Find unique species pairs
  unique_pairs <- unique(species_pairs)
  
  # Initialize mapping pattern data
  pattern_data <- list()
  
  # Analyze each conserved group
  for(group_name in names(conserved_groups$conserved_groups)) {
    group <- conserved_groups$conserved_groups[[group_name]]
    message(paste(" - Analyzing mapping patterns for", group_name, "..."))
    
    # Get species with chromosomes in this group
    group_species <- names(group$species_chr_counts)
    
    # Look for patterns in each species pair
    for(i in 1:nrow(unique_pairs)) {
      sp_A <- unique_pairs$species_A[i]
      sp_B <- unique_pairs$species_B[i]
      
      # Skip if either species not in group
      if(!(sp_A %in% group_species) || !(sp_B %in% group_species)) {
        next
      }
      
      # Get chromosomes for these species in this group
      chr_A <- strsplit(group$species_chromosomes[[sp_A]], ",")[[1]]
      chr_B <- strsplit(group$species_chromosomes[[sp_B]], ",")[[1]]
      
      # Look for fusion pattern (multiple chromosomes in A map to one in B)
      if(length(chr_A) > length(chr_B) && length(chr_B) == 1) {
        pattern_data[[length(pattern_data) + 1]] <- list(
          group = group_name,
          species_A = sp_A,
          species_B = sp_B,
          chromosomes_A = paste(chr_A, collapse=","),
          chromosomes_B = chr_B,
          pattern_type = "fusion_candidate",
          chromosome_ratio = length(chr_A) / length(chr_B)
        )
      }
      
      # Look for fission pattern (one chromosome in A maps to multiple in B)
      if(length(chr_A) == 1 && length(chr_B) > 1) {
        pattern_data[[length(pattern_data) + 1]] <- list(
          group = group_name,
          species_A = sp_A,
          species_B = sp_B,
          chromosomes_A = chr_A,
          chromosomes_B = paste(chr_B, collapse=","),
          pattern_type = "fission_candidate",
          chromosome_ratio = length(chr_B) / length(chr_A)
        )
      }
    }
  }
  
  # Convert to data frame
  if(length(pattern_data) > 0) {
    pattern_df <- do.call(rbind, lapply(pattern_data, function(x) {
      data.frame(
        group = x$group,
        species_A = x$species_A,
        species_B = x$species_B,
        chromosomes_A = x$chromosomes_A,
        chromosomes_B = x$chromosomes_B,
        pattern_type = x$pattern_type,
        chromosome_ratio = x$chromosome_ratio,
        stringsAsFactors = FALSE
      )
    }))
  } else {
    # Create empty dataframe if no patterns found
    pattern_df <- data.frame(
      group = character(),
      species_A = character(),
      species_B = character(),
      chromosomes_A = character(),
      chromosomes_B = character(),
      pattern_type = character(),
      chromosome_ratio = numeric(),
      stringsAsFactors = FALSE
    )
  }
  
  message(paste("Detected", nrow(pattern_df), "potential chromosome event patterns"))
  
  return(pattern_df)
}

# ==============================
# Event Inference Functions
# ==============================

infer_chromosome_events <- function(tree, edge_data, count_changes, pattern_data, chr_counts) {
  # Infer chromosome evolutionary events
  message("Inferring chromosome evolutionary events...")
  
  # Get unique pattern types
  if(nrow(pattern_data) > 0) {
    fusion_patterns <- pattern_data[pattern_data$pattern_type == "fusion_candidate",]
    fission_patterns <- pattern_data[pattern_data$pattern_type == "fission_candidate",]
  } else {
    fusion_patterns <- fission_patterns <- data.frame()
  }
  
  # Initialize events list
  events <- list()
  
  # Infer events based on chromosome count changes where child is a tip
  for(i in 1:nrow(count_changes)) {
    # Skip if no count data for child
    if(is.na(count_changes$child_count[i])) {
      next
    }
    
    # Skip if not a tip
    if(!edge_data$child_is_tip[i]) {
      next
    }
    
    child_species <- count_changes$child_label[i]
    child_count <- count_changes$child_count[i]
    
    # Look for related species (sister taxa) to compare
    sister_species <- find_related_species(tree, child_species)
    
    # Skip if no related species found
    if(length(sister_species) == 0) {
      next
    }
    
    # Get counts for related species
    sister_counts <- chr_counts[intersect(sister_species, names(chr_counts))]
    
    # Skip if no count data for related species
    if(length(sister_counts) == 0) {
      next
    }
    
    # Calculate median count of related species
    median_count <- median(sister_counts, na.rm=TRUE)
    
    # Detect potential events
    if(child_count > median_count * 1.5) {
      # Potential chromosome fission event(s)
      events[[length(events) + 1]] <- list(
        edge_id = count_changes$edge_id[i],
        parent_node = count_changes$parent_node_id[i],
        child_node = count_changes$child_node_id[i],
        parent_label = count_changes$parent_label[i],
        child_label = child_species,
        event_type = "fission",
        child_count = child_count,
        reference_count = median_count,
        count_diff = child_count - median_count,
        count_ratio = child_count / median_count,
        evidence_type = "count_comparison",
        confidence = "medium"
      )
    } else if(child_count < median_count * 0.67) {
      # Potential chromosome fusion event(s)
      events[[length(events) + 1]] <- list(
        edge_id = count_changes$edge_id[i],
        parent_node = count_changes$parent_node_id[i],
        child_node = count_changes$child_node_id[i],
        parent_label = count_changes$parent_label[i],
        child_label = child_species,
        event_type = "fusion",
        child_count = child_count,
        reference_count = median_count,
        count_diff = median_count - child_count,
        count_ratio = median_count / child_count,
        evidence_type = "count_comparison",
        confidence = "medium"
      )
    }
    
    # Check for supporting mapping patterns
    if(nrow(fusion_patterns) > 0) {
      # Look for fusion patterns where this species has fewer chromosomes
      fusion_evidence <- fusion_patterns[
        (fusion_patterns$species_B == child_species) & 
          (fusion_patterns$species_A %in% sister_species),
      ]
      
      if(nrow(fusion_evidence) > 0) {
        for(j in 1:nrow(fusion_evidence)) {
          events[[length(events) + 1]] <- list(
            edge_id = count_changes$edge_id[i],
            parent_node = count_changes$parent_node_id[i],
            child_node = count_changes$child_node_id[i],
            parent_label = count_changes$parent_label[i],
            child_label = child_species,
            event_type = "fusion",
            child_count = child_count,
            reference_count = chr_counts[fusion_evidence$species_A[j]],
            count_diff = NA,
            count_ratio = fusion_evidence$chromosome_ratio[j],
            evidence_type = "mapping_pattern",
            confidence = "high",
            group = fusion_evidence$group[j],
            chromosomes_A = fusion_evidence$chromosomes_A[j],
            chromosomes_B = fusion_evidence$chromosomes_B[j],
            reference_species = fusion_evidence$species_A[j]
          )
        }
      }
    }
    
    if(nrow(fission_patterns) > 0) {
      # Look for fission patterns where this species has more chromosomes
      fission_evidence <- fission_patterns[
        (fission_patterns$species_B == child_species) & 
          (fission_patterns$species_A %in% sister_species),
      ]
      
      if(nrow(fission_evidence) > 0) {
        for(j in 1:nrow(fission_evidence)) {
          events[[length(events) + 1]] <- list(
            edge_id = count_changes$edge_id[i],
            parent_node = count_changes$parent_node_id[i],
            child_node = count_changes$child_node_id[i],
            parent_label = count_changes$parent_label[i],
            child_label = child_species,
            event_type = "fission",
            child_count = child_count,
            reference_count = chr_counts[fission_evidence$species_A[j]],
            count_diff = NA,
            count_ratio = fission_evidence$chromosome_ratio[j],
            evidence_type = "mapping_pattern",
            confidence = "high",
            group = fission_evidence$group[j],
            chromosomes_A = fission_evidence$chromosomes_A[j],
            chromosomes_B = fission_evidence$chromosomes_B[j],
            reference_species = fission_evidence$species_A[j]
          )
        }
      }
    }
  }
  
  # Convert to data frame if events were found
  if(length(events) > 0) {
    # Handle different event structures by filling in missing columns
    events_list <- lapply(events, function(event) {
      # Create basic event data
      basic_event <- data.frame(
        edge_id = event$edge_id,
        parent_node = event$parent_node,
        child_node = event$child_node,
        parent_label = event$parent_label,
        child_label = event$child_label,
        event_type = event$event_type,
        child_count = event$child_count,
        reference_count = event$reference_count,
        count_diff = event$count_diff,
        count_ratio = event$count_ratio,
        evidence_type = event$evidence_type,
        confidence = event$confidence,
        stringsAsFactors = FALSE
      )
      
      # Add mapping pattern specific fields if present
      if(!is.null(event$group)) {
        basic_event$group <- event$group
        basic_event$chromosomes_A <- event$chromosomes_A
        basic_event$chromosomes_B <- event$chromosomes_B
        basic_event$reference_species <- event$reference_species
      } else {
        basic_event$group <- NA
        basic_event$chromosomes_A <- NA
        basic_event$chromosomes_B <- NA
        basic_event$reference_species <- NA
      }
      
      return(basic_event)
    })
    
    events_df <- do.call(rbind, events_list)
  } else {
    # Create empty dataframe if no events found
    events_df <- data.frame(
      edge_id = integer(),
      parent_node = integer(),
      child_node = integer(),
      parent_label = character(),
      child_label = character(),
      event_type = character(),
      child_count = integer(),
      reference_count = integer(),
      count_diff = integer(),
      count_ratio = numeric(),
      evidence_type = character(),
      confidence = character(),
      group = character(),
      chromosomes_A = character(),
      chromosomes_B = character(),
      reference_species = character(),
      stringsAsFactors = FALSE
    )
  }
  
  message(paste("Inferred", nrow(events_df), "chromosome evolutionary events"))
  
  return(events_df)
}

identify_shared_events <- function(tree, events_df, clades) {
  # Identify shared chromosome events across branches
  message("Identifying shared chromosome events...")
  
  # Initialize results
  shared_events <- list()
  
  # If no events were found, return empty result
  if(nrow(events_df) == 0) {
    message("No events to analyze for sharing")
    return(data.frame(
      event_type = character(),
      clade_node_id = integer(),
      clade_label = character(),
      tip_count = integer(),
      event_count = integer(),
      shared_status = character(),
      confidence = character(),
      stringsAsFactors = FALSE
    ))
  }
  
  # Get unique event types
  event_types <- unique(events_df$event_type)
  
  # For each event type and clade, check if the event is shared
  for(event_type in event_types) {
    # Filter events of this type
    type_events <- events_df[events_df$event_type == event_type,]
    
    # Check each clade
    for(clade_id in names(clades)) {
      clade <- clades[[clade_id]]
      clade_tips <- clade$tips
      
      # Count events in this clade
      clade_events <- type_events[type_events$child_label %in% clade_tips,]
      event_count <- nrow(clade_events)
      
      # Skip clades with no events
      if(event_count == 0) {
        next
      }
      
      # Calculate event prevalence in clade
      prevalence <- event_count / clade$tip_count
      
      # Determine if the event is shared across the clade
      if(prevalence >= 0.8) {
        shared_status <- "strongly_shared"
        confidence <- "high"
      } else if(prevalence >= 0.5) {
        shared_status <- "likely_shared"
        confidence <- "medium"
      } else if(prevalence >= 0.3) {
        shared_status <- "partially_shared"
        confidence <- "low"
      } else {
        # Skip clades with low prevalence
        next
      }
      
      # Record shared event
      shared_events[[length(shared_events) + 1]] <- list(
        event_type = event_type,
        clade_node_id = as.integer(clade_id),
        clade_label = clade$node_label,
        tip_count = clade$tip_count,
        event_count = event_count,
        prevalence = prevalence,
        shared_status = shared_status,
        confidence = confidence,
        tip_events = list(clade_events)
      )
    }
  }
  
  # Convert to data frame if shared events were found
  if(length(shared_events) > 0) {
    shared_df <- do.call(rbind, lapply(shared_events, function(e) {
      data.frame(
        event_type = e$event_type,
        clade_node_id = e$clade_node_id,
        clade_label = e$clade_label,
        tip_count = e$tip_count,
        event_count = e$event_count,
        prevalence = e$prevalence,
        shared_status = e$shared_status,
        confidence = e$confidence,
        stringsAsFactors = FALSE
      )
    }))
    
    # Store tip events separately because they don't fit in a data frame
    attr(shared_df, "tip_events") <- lapply(shared_events, function(e) e$tip_events[[1]])
  } else {
    # Create empty dataframe if no shared events found
    shared_df <- data.frame(
      event_type = character(),
      clade_node_id = integer(),
      clade_label = character(),
      tip_count = integer(),
      event_count = integer(),
      prevalence = numeric(),
      shared_status = character(),
      confidence = character(),
      stringsAsFactors = FALSE
    )
    attr(shared_df, "tip_events") <- list()
  }
  
  message(paste("Identified", nrow(shared_df), "shared chromosome events"))
  
  return(shared_df)
}

find_related_species <- function(tree, species) {
  # Find related species (sister taxa) for a given species
  # First get the node that this species belongs to
  tip_idx <- which(tree$tip.label == species)
  
  # If species not in tree, return empty vector
  if(length(tip_idx) == 0) {
    return(character(0))
  }
  
  # Get parent node of this species
  edges <- tree$edge
  parent_edge <- which(edges[, 2] == tip_idx)
  
  if(length(parent_edge) == 0) {
    return(character(0))
  }
  
  parent_node <- edges[parent_edge, 1]
  
  # Get all children of this parent (sister taxa)
  sister_tips <- tree$tip.label[tree$edge[tree$edge[, 1] == parent_node, 2]]
  sister_tips <- sister_tips[sister_tips != species]
  
  # If no sisters, try going up one more level
  if(length(sister_tips) == 0) {
    grandparent_edge <- which(edges[, 2] == parent_node)
    
    if(length(grandparent_edge) == 0) {
      return(character(0))
    }
    
    grandparent_node <- edges[grandparent_edge, 1]
    
    # Get all descendants of the grandparent
    descendant_nodes <- tree$edge[tree$edge[, 1] == grandparent_node, 2]
    
    # For each node, get its tip descendants
    for(node in descendant_nodes) {
      if(node <= length(tree$tip.label)) {
        # This is a tip
        sister_tips <- c(sister_tips, tree$tip.label[node])
      } else {
        # This is an internal node, get all its descendants
        descendants <- tryCatch({
          geiger::tips(tree, node)
        }, error = function(e) {
          # Alternative method if geiger is not available
          desc_nodes <- getDescendants(tree, node)
          desc_nodes <- desc_nodes[desc_nodes <= length(tree$tip.label)]
          tree$tip.label[desc_nodes]
        })
        sister_tips <- c(sister_tips, descendants)
      }
    }
    
    # Remove the original species
    sister_tips <- sister_tips[sister_tips != species]
  }
  
  return(sister_tips)
}

detect_shared_chromosome_events <- function(integrated_data, conserved_groups) {
  # Main function to detect shared chromosome events
  message("Starting shared chromosome event detection...")
  
  # Get tree and chromosome counts
  tree <- integrated_data$tree
  chr_counts <- integrated_data$chr_counts
  
  # Extract tree structure information
  edge_data <- extract_tree_edges(tree)
  clades <- identify_clades(tree)
  
  # Detect chromosome count changes
  count_changes <- detect_count_changes(edge_data, chr_counts)
  
  # Detect mapping patterns
  pattern_data <- detect_mapping_patterns(edge_data, integrated_data, conserved_groups)
  
  # Infer chromosome events
  events_df <- infer_chromosome_events(tree, edge_data, count_changes, pattern_data, chr_counts)
  
  # Identify shared events
  shared_events_df <- identify_shared_events(tree, events_df, clades)
  
  # Return results
  result <- list(
    tree = tree,
    edge_data = edge_data,
    clades = clades,
    count_changes = count_changes,
    pattern_data = pattern_data,
    events = events_df,
    shared_events = shared_events_df
  )
  
  message("Shared chromosome event detection complete")
  
  return(result)
}

# ==============================
# Output Functions
# ==============================

save_event_detection_results <- function(event_result, output_dir) {
  # Save event detection results to output files
  message("Saving event detection results...")
  
  # Create output directory if it doesn't exist
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message(paste("Created output directory:", output_dir))
  }
  
  # Save full results data
  events_rdata_file <- file.path(output_dir, "chromosome_events.RData")
  saveRDS(event_result, file = events_rdata_file)
  
  # Save individual data files
  
  # Edge data
  edge_data_file <- file.path(output_dir, "tree_edges.csv")
  write.csv(event_result$edge_data, file = edge_data_file, row.names = FALSE)
  
  # Count changes
  count_changes_file <- file.path(output_dir, "chromosome_count_changes.csv")
  write.csv(event_result$count_changes, file = count_changes_file, row.names = FALSE)
  
  # Mapping patterns
  pattern_data_file <- file.path(output_dir, "mapping_patterns.csv")
  write.csv(event_result$pattern_data, file = pattern_data_file, row.names = FALSE)
  
  # Events
  events_file <- file.path(output_dir, "inferred_events.csv")
  write.csv(event_result$events, file = events_file, row.names = FALSE)
  
  # Shared events
  shared_events_file <- file.path(output_dir, "shared_events.csv")
  write.csv(event_result$shared_events, file = shared_events_file, row.names = FALSE)
  
  # Create detailed events report
  events_report <- file.path(output_dir, "chromosome_events_report.txt")
  
  sink(events_report)
  cat("Chromosome Evolutionary Events Report\n")
  cat("====================================\n\n")
  
  cat("1. Summary of Detected Events\n")
  cat("---------------------------\n")
  if(nrow(event_result$events) > 0) {
    event_summary <- table(event_result$events$event_type)
    for(evt_type in names(event_summary)) {
      cat(paste0("  - ", evt_type, ": ", event_summary[evt_type], " events\n"))
    }
    
    evidence_summary <- table(event_result$events$evidence_type)
    cat("\nEvidence types:\n")
    for(evd_type in names(evidence_summary)) {
      cat(paste0("  - ", evd_type, ": ", evidence_summary[evd_type], " instances\n"))
    }
    
    confidence_summary <- table(event_result$events$confidence)
    cat("\nConfidence levels:\n")
    for(conf in names(confidence_summary)) {
      cat(paste0("  - ", conf, ": ", confidence_summary[conf], " events\n"))
    }
  } else {
    cat("  No chromosome events detected\n")
  }
  cat("\n")
  
  cat("2. Shared Chromosome Events\n")
  cat("-------------------------\n")
  if(nrow(event_result$shared_events) > 0) {
    for(i in 1:nrow(event_result$shared_events)) {
      se <- event_result$shared_events[i,]
      cat(paste0("  [", i, "] ", se$event_type, " in clade ", se$clade_label, "\n"))
      cat(paste0("      - Prevalence: ", round(se$prevalence * 100), "% (", 
                 se$event_count, " of ", se$tip_count, " species)\n"))
      cat(paste0("      - Status: ", se$shared_status, "\n"))
      cat(paste0("      - Confidence: ", se$confidence, "\n"))
      cat("\n")
    }
  } else {
    cat("  No shared chromosome events detected\n")
  }
  cat("\n")
  
  cat("3. Detailed Event Listing\n")
  cat("----------------------\n")
  if(nrow(event_result$events) > 0) {
    high_conf_events <- event_result$events[event_result$events$confidence == "high",]
    
    if(nrow(high_conf_events) > 0) {
      cat("High confidence events:\n")
      for(i in 1:min(10, nrow(high_conf_events))) {
        evt <- high_conf_events[i,]
        cat(paste0("  [", i, "] ", evt$event_type, " in ", evt$child_label, 
                   " (", evt$evidence_type, ")\n"))
        if(evt$evidence_type == "mapping_pattern") {
          cat(paste0("      - Group: ", evt$group, "\n"))
          cat(paste0("      - Reference species: ", evt$reference_species, "\n"))
          cat(paste0("      - Chromosomes: ", 
                     if(evt$event_type == "fusion") {
                       paste0(evt$chromosomes_A, " → ", evt$chromosomes_B)
                     } else {
                       paste0(evt$chromosomes_A, " → ", evt$chromosomes_B)
                     }, "\n"))
        }
      }
      if(nrow(high_conf_events) > 10) {
        cat(paste0("  ...and ", nrow(high_conf_events) - 10, " more high confidence events\n"))
      }
      cat("\n")
    }
  } else {
    cat("  No events detected\n")
  }
  
  sink()
  
  message(paste("Saved event detection results to:", output_dir))
  
  # Return file paths
  return(list(
    events_rdata = events_rdata_file,
    edge_data = edge_data_file,
    count_changes = count_changes_file,
    pattern_data = pattern_data_file,
    events = events_file,
    shared_events = shared_events_file,
    events_report = events_report
  ))
}

# ==============================
# Main Function
# ==============================

main <- function(args) {
  # Check arguments
  if(length(args) < 3) {
    cat("Usage: Rscript ancestral_chromosomes_phase3.R integrated_data_file conserved_groups_file output_dir\n")
    cat("Example: Rscript ancestral_chromosomes_phase3.R results/integrated_data.RData results/phase2/conserved_chromosome_groups.RData results/phase3/\n")
    return(1)
  }
  
  integrated_data_file <- args[1]
  conserved_groups_file <- args[2]
  output_dir <- args[3]
  
  # Load data
  integrated_data <- load_integrated_data(integrated_data_file)
  conserved_groups <- load_conserved_groups(conserved_groups_file)
  
  # Detect shared chromosome events
  event_result <- detect_shared_chromosome_events(integrated_data, conserved_groups)
  
  # Save results
  output_files <- save_event_detection_results(event_result, output_dir)
  
  # Print basic summary
  cat("\n===== Chromosome Evolutionary Events Summary =====\n")
  
  if(nrow(event_result$events) > 0) {
    # Summarize events by type
    event_types <- table(event_result$events$event_type)
    cat("Detected events by type:\n")
    for(evt_type in names(event_types)) {
      cat(paste0("  - ", evt_type, ": ", event_types[evt_type], "\n"))
    }
    
    # Summarize by evidence type
    evidence_types <- table(event_result$events$evidence_type)
    cat("\nEvidence types:\n")
    for(evd_type in names(evidence_types)) {
      cat(paste0("  - ", evd_type, ": ", evidence_types[evd_type], "\n"))
    }
  } else {
    cat("No chromosome events detected\n")
  }
  
  # Summarize shared events
  if(nrow(event_result$shared_events) > 0) {
    cat("\nShared events:\n")
    shared_by_type <- table(event_result$shared_events$event_type)
    for(evt_type in names(shared_by_type)) {
      cat(paste0("  - ", evt_type, ": ", shared_by_type[evt_type], "\n"))
    }
    
    # List top shared events
    cat("\nTop shared events by prevalence:\n")
    top_shared <- event_result$shared_events[order(-event_result$shared_events$prevalence),]
    for(i in 1:min(3, nrow(top_shared))) {
      se <- top_shared[i,]
      cat(paste0("  [", i, "] ", se$event_type, " in clade ", se$clade_label,
                 " (", round(se$prevalence * 100), "% prevalence)\n"))
    }
  } else {
    cat("\nNo shared chromosome events detected\n")
  }
  
  cat(paste("\nResults saved to:", output_dir, "\n"))
  cat("==============================================\n")
  
  # Return success
  return(0)
}

# ==============================
# Helper Functions
# ==============================

# If geiger package is not available, provide a fallback implementation for tips function
if(!requireNamespace("geiger", quietly = TRUE)) {
  message("Note: geiger package not found, using built-in methods for clade analysis")
  
  # Get all tip descendants of a node
  getTips <- function(tree, node) {
    if(node <= length(tree$tip.label)) {
      # This is already a tip
      return(tree$tip.label[node])
    }
    
    # Get all descendants of this node
    desc <- getDescendants(tree, node)
    
    # Filter to just tips
    tip_desc <- desc[desc <= length(tree$tip.label)]
    
    # Return tip labels
    return(tree$tip.label[tip_desc])
  }
  
  # Assign our function to the geiger::tips namespace
  geiger <- list(tips = getTips)
}

# Run the script if executed directly
if(!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  exit_code <- main(args)
  quit(status = exit_code)
}
