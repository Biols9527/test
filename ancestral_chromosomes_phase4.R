#!/usr/bin/env Rscript

# Ancestral Chromosome Reconstruction Phase 4: Ancestral Linkage Group Reconstruction
# Date: 2025-03-18
# Author: ShrivatsaMehan
# Description: Reconstructs ancestral linkage groups at each internal node of the phylogenetic tree
#              by integrating previous results and using both bottom-up and top-down approaches.

# Required packages
suppressPackageStartupMessages({
  library(ape)         # For phylogenetic tree manipulation
  library(data.table)  # For efficient data handling
  library(igraph)      # For network analysis
  library(phangorn)    # For ancestral state reconstruction
  library(ggplot2)     # For visualization
  library(ggtree)      # For phylogenetic tree visualization
  library(reshape2)    # For data reshaping
})

# ==============================
# Default Configuration Parameters
# ==============================

# Create default configuration
default_config <- list(
  # 简约法参数
  parsimony_coverage_threshold = 0.5,  # 简约法中判定祖先存在的覆盖率阈值
  
  # 网络分析参数
  community_detection_methods = c("louvain", "fast_greedy", "label_prop", "walktrap"),
  small_community_threshold = 3,  # 小社区过滤阈值
  
  # 置信度计算参数
  conservation_weight = 0.6,       # 保守性在置信度中的权重
  species_coverage_weight = 0.4,   # 物种覆盖在置信度中的权重
  
  # 方法整合参数
  common_group_boost = 1.2,        # 共同群组的置信度提升系数
  unique_group_penalty = 0.9,      # 单一方法群组的置信度惩罚系数
  
  # 事件处理参数
  fusion_confidence_penalty = 0.9, # 融合事件应用后的置信度衰减
  fission_confidence_penalty = 0.85, # 裂变事件应用后的置信度衰减
  
  # 进化模型参数
  allow_wgd_events = TRUE,         # 是否考虑全基因组复制事件
  max_chromosome_ratio = 2.5,      # 判定为WGD的染色体比例阈值
  
  # 分类群特异参数
  taxon_specific_rates = list(     # 各分类群相对变异率
    "Mammalia" = 0.8,
    "Aves" = 0.5,
    "Teleostei" = 1.2,
    "default" = 1.0
  )
)

# Initialize with default config
config <- default_config

# ==============================
# Data Loading Functions
# ==============================

load_previous_phases <- function(phase1_file, phase2_file, phase3_file) {
  # Load data from previous phases
  message("Loading data from previous phases...")
  
  if(!all(file.exists(phase1_file, phase2_file, phase3_file))) {
    missing_files <- c(phase1_file, phase2_file, phase3_file)[!file.exists(c(phase1_file, phase2_file, phase3_file))]
    stop("The following required files are missing: ", paste(missing_files, collapse=", "))
  }
  
  integrated_data <- readRDS(phase1_file)
  message("  - Loaded integrated data with ", length(integrated_data$common_species), " species")
  
  conserved_groups <- readRDS(phase2_file)
  message("  - Loaded ", length(conserved_groups$conserved_groups), " conserved chromosome groups")
  
  event_results <- readRDS(phase3_file)
  message("  - Loaded chromosome evolutionary events")
  
  return(list(
    integrated_data = integrated_data,
    conserved_groups = conserved_groups,
    event_results = event_results
  ))
}

# ==============================
# Tree Analysis Functions
# ==============================

prepare_node_reconstruction <- function(tree, event_results) {
  # Prepare for node-by-node reconstruction
  message("Preparing node-specific reconstruction context...")
  
  # Get number of species (tips)
  n_tips <- length(tree$tip.label)
  
  # Get all internal nodes
  internal_nodes <- (n_tips + 1):(n_tips + tree$Nnode)
  
  # Create context for each node
  node_contexts <- lapply(internal_nodes, function(node_id) {
    # Get node label
    if(!is.null(tree$node.label) && length(tree$node.label) >= (node_id - n_tips)) {
      node_label <- tree$node.label[node_id - n_tips]
      if(is.na(node_label) || node_label == "") {
        node_label <- paste0("Node", node_id)
      }
    } else {
      node_label <- paste0("Node", node_id)
    }
    
    # Get immediate descendants
    descendants <- get_descendants(tree, node_id)
    
    # Get tip descendants (extant species)
    tip_descendants <- get_tip_descendants(tree, node_id)
    
    # Get node's parent
    parent_id <- get_parent_node(tree, node_id)
    
    # Get node's estimated chromosome count
    node_count <- NA
    if(!is.null(event_results$ancestral_counts) && !is.null(event_results$ancestral_counts$node_counts)) {
      node_count <- event_results$ancestral_counts$node_counts[as.character(node_id)]
    }
    
    # Get events where this node is involved
    node_events <- event_results$events[event_results$events$parent_node == node_id | 
                                     event_results$events$child_node == node_id, ]
    
    # Get shared events associated with this node
    shared_events <- event_results$shared_events[event_results$shared_events$clade_node_id == node_id, ]
    
    # Try to identify taxon for applying taxon-specific rates
    # This is a simplified approach - in a real analysis, would need more sophisticated taxon assignment
    taxon <- "default"
    if(length(tip_descendants) > 0) {
      # Use the first species name to guess taxon (simplistic)
      species_name <- tree$tip.label[tip_descendants[1]]
      if(grepl("Mammalia", species_name, ignore.case = TRUE)) {
        taxon <- "Mammalia"
      } else if(grepl("Aves", species_name, ignore.case = TRUE)) {
        taxon <- "Aves"
      } else if(grepl("fish|teleost", species_name, ignore.case = TRUE)) {
        taxon <- "Teleostei"
      }
    }
    
    list(
      node_id = node_id,
      node_label = node_label,
      parent_id = parent_id,
      descendants = descendants,
      tip_descendants = tip_descendants,
      tip_species = tree$tip.label[tip_descendants],
      chr_count = node_count,
      node_events = if(is.data.frame(node_events)) node_events else data.frame(),
      shared_events = if(is.data.frame(shared_events)) shared_events else data.frame(),
      taxon = taxon
    )
  })
  
  message(paste("  - Created reconstruction context for", length(node_contexts), "internal nodes"))
  
  # Sort nodes by depth (for bottom-up and top-down approaches)
  node_depths <- get_node_depths(tree)
  
  # Add depth information to contexts
  for(i in seq_along(node_contexts)) {
    node_contexts[[i]]$depth <- node_depths[node_contexts[[i]]$node_id]
  }
  
  # Return sorted contexts for bottom-up and top-down approaches
  bottom_up_order <- order(sapply(node_contexts, function(x) x$depth), decreasing = TRUE)
  top_down_order <- order(sapply(node_contexts, function(x) x$depth))
  
  return(list(
    contexts = node_contexts,
    bottom_up_order = bottom_up_order,
    top_down_order = top_down_order
  ))
}

# Helper function to get node depths
get_node_depths <- function(tree) {
  # Calculate depth of each node from the root
  n_tips <- length(tree$tip.label)
  n_nodes <- tree$Nnode
  
  # Initialize depths (root is depth 0)
  depths <- numeric(n_tips + n_nodes)
  root_node <- n_tips + 1
  depths[root_node] <- 0
  
  # Process nodes in topological order
  for(i in 1:nrow(tree$edge)) {
    parent <- tree$edge[i, 1]
    child <- tree$edge[i, 2]
    depths[child] <- depths[parent] + 1
  }
  
  names(depths) <- 1:(n_tips + n_nodes)
  return(depths)
}

# Helper function to get descendants of a node
get_descendants <- function(tree, node) {
  # Get immediate child nodes
  children <- tree$edge[tree$edge[, 1] == node, 2]
  return(children)
}

# Helper function to get all tip descendants
get_tip_descendants <- function(tree, node) {
  # Get all tip descendants of a node
  n_tips <- length(tree$tip.label)
  
  # Handle special case where node is already a tip
  if(node <= n_tips) {
    return(node)
  }
  
  # Get all descendants
  all_desc <- c()
  to_process <- node
  
  while(length(to_process) > 0) {
    current <- to_process[1]
    to_process <- to_process[-1]
    
    children <- tree$edge[tree$edge[, 1] == current, 2]
    all_desc <- c(all_desc, children)
    
    # Add non-tip children to processing queue
    to_process <- c(to_process, children[children > n_tips])
  }
  
  # Return only tips
  return(all_desc[all_desc <= n_tips])
}

# Helper function to get parent node
get_parent_node <- function(tree, node) {
  # Find the parent of a given node
  # Root node has no parent
  if(node == length(tree$tip.label) + 1) {
    return(NA)
  }
  
  parent_edge <- which(tree$edge[, 2] == node)
  
  if(length(parent_edge) == 0) {
    return(NA)
  }
  
  return(tree$edge[parent_edge, 1])
}

# ==============================
# Linkage Marker Functions
# ==============================

create_linkage_markers <- function(conserved_groups, integrated_data) {
  # Create linkage markers from conserved groups
  message("Creating linkage markers from conserved groups...")
  
  # Extract chromosome mappings from conserved groups
  linkage_markers <- list()
  
  for(i in seq_along(conserved_groups$conserved_groups)) {
    group <- conserved_groups$conserved_groups[[i]]
    group_name <- names(conserved_groups$conserved_groups)[i]
    
    # Get all species-chromosome associations for this group
    species_chroms <- lapply(names(group$species_chromosomes), function(sp) {
      chrs <- tryCatch({
        strsplit(group$species_chromosomes[[sp]], ",")[[1]]
      }, error = function(e) {
        character(0)
      })
      
      if(length(chrs) == 0) {
        return(NULL)
      }
      
      data.frame(
        species = sp,
        chromosome = chrs,
        group = group_name,
        conservation_score = group$conservation_score,
        stringsAsFactors = FALSE
      )
    })
    
    # Remove NULL entries
    species_chroms <- species_chroms[!sapply(species_chroms, is.null)]
    
    if(length(species_chroms) > 0) {
      linkage_markers[[group_name]] <- do.call(rbind, species_chroms)
    }
  }
  
  # Combine all markers
  if(length(linkage_markers) > 0) {
    all_markers <- do.call(rbind, linkage_markers)
  } else {
    all_markers <- data.frame(
      species = character(),
      chromosome = character(),
      group = character(),
      conservation_score = numeric(),
      stringsAsFactors = FALSE
    )
  }
  
  # Enhance with bidirectional mapping data
  enhanced_markers <- enhance_with_mapping_data(all_markers, integrated_data$bidirectional_maps)
  
  message(paste("  - Created", nrow(enhanced_markers), "linkage markers across", 
                length(unique(enhanced_markers$species)), "species"))
  
  return(enhanced_markers)
}

enhance_with_mapping_data <- function(markers, bidirectional_maps) {
  # Enhance linkage markers with additional mapping information
  message("  - Enhancing markers with bidirectional mapping data...")
  
  # If no markers or no mapping data, return as is
  if(nrow(markers) == 0 || nrow(bidirectional_maps) == 0) {
    return(markers)
  }
  
  # Add columns for edge weights based on mapping quality
  markers$mapping_weight <- 1
  
  # Create a mapping key for quick lookup
  mapping_keys_A <- paste(bidirectional_maps$species_A, bidirectional_maps$chromosome_A, sep=":")
  mapping_keys_B <- paste(bidirectional_maps$species_B, bidirectional_maps$chromosome_B, sep=":")
  
  # Create keys for markers
  marker_keys <- paste(markers$species, markers$chromosome, sep=":")
  
  # Find mappings that match markers
  for(i in 1:nrow(markers)) {
    key <- marker_keys[i]
    
    # Check if this marker is in mapping data
    match_indices_A <- which(mapping_keys_A == key)
    match_indices_B <- which(mapping_keys_B == key)
    
    if(length(match_indices_A) > 0 || length(match_indices_B) > 0) {
      # Calculate average mapping weight
      weights <- c()
      
      if(length(match_indices_A) > 0) {
        if("mapping_quality" %in% names(bidirectional_maps)) {
          weights <- c(weights, bidirectional_maps$mapping_quality[match_indices_A])
        } else if("count" %in% names(bidirectional_maps)) {
          weights <- c(weights, bidirectional_maps$count[match_indices_A])
        }
      }
      
      if(length(match_indices_B) > 0) {
        if("mapping_quality" %in% names(bidirectional_maps)) {
          weights <- c(weights, bidirectional_maps$mapping_quality[match_indices_B])
        } else if("count" %in% names(bidirectional_maps)) {
          weights <- c(weights, bidirectional_maps$count[match_indices_B])
        }
      }
      
      if(length(weights) > 0) {
        markers$mapping_weight[i] <- mean(weights, na.rm = TRUE)
      }
    }
  }
  
  # Add homology links based on mapping data
  markers$homology_links <- sapply(marker_keys, function(key) {
    # Find mappings for this marker
    homologs <- c()
    
    # Check as species_A
    match_indices_A <- which(mapping_keys_A == key)
    if(length(match_indices_A) > 0) {
      homologs <- c(homologs, mapping_keys_B[match_indices_A])
    }
    
    # Check as species_B
    match_indices_B <- which(mapping_keys_B == key)
    if(length(match_indices_B) > 0) {
      homologs <- c(homologs, mapping_keys_A[match_indices_B])
    }
    
    # Return as comma-separated string
    if(length(homologs) > 0) {
      return(paste(unique(homologs), collapse=","))
    } else {
      return(NA)
    }
  })
  
  # Add gene content information if available
  if("gene_content" %in% names(bidirectional_maps)) {
    markers$gene_content <- sapply(marker_keys, function(key) {
      content <- c()
      
      # Check as species_A
      match_indices_A <- which(mapping_keys_A == key)
      if(length(match_indices_A) > 0) {
        content <- c(content, bidirectional_maps$gene_content[match_indices_A])
      }
      
      # Check as species_B
      match_indices_B <- which(mapping_keys_B == key)
      if(length(match_indices_B) > 0) {
        content <- c(content, bidirectional_maps$gene_content[match_indices_B])
      }
      
      if(length(content) > 0) {
        return(mean(content, na.rm = TRUE))
      } else {
        return(NA)
      }
    })
  } else {
    markers$gene_content <- NA
  }
  
  return(markers)
}

# ==============================
# Reconstruction Networks
# ==============================

create_marker_network <- function(markers) {
  # Create a network of markers based on homology and group relationships
  message("Creating marker relationship network...")
  
  # If no markers, return empty graph
  if(nrow(markers) == 0) {
    return(graph.empty())
  }
  
  # Create node IDs from species and chromosome
  node_ids <- paste(markers$species, markers$chromosome, sep=":")
  
  # Initialize empty edge list
  edges <- data.frame(from = character(), to = character(), weight = numeric(), stringsAsFactors = FALSE)
  
  # Add edges based on conserved groups
  for(group in unique(markers$group)) {
    group_markers <- markers[markers$group == group, ]
    group_ids <- paste(group_markers$species, group_markers$chromosome, sep=":")
    
    if(length(group_ids) > 1) {
      # Create edges between all markers in the same group
      for(i in 1:(length(group_ids)-1)) {
        for(j in (i+1):length(group_ids)) {
          # Weight is based on conservation score
          weight <- group_markers$conservation_score[1]
          
          edges <- rbind(edges, data.frame(
            from = group_ids[i],
            to = group_ids[j],
            weight = weight,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  
  # Add edges based on homology links
  for(i in 1:nrow(markers)) {
    if(!is.na(markers$homology_links[i])) {
      # Get homology links
      homologs <- strsplit(markers$homology_links[i], ",")[[1]]
      
      # Add edges to homologs
      for(homolog in homologs) {
        # Weight is based on mapping weight
        weight <- markers$mapping_weight[i]
        
        edges <- rbind(edges, data.frame(
          from = node_ids[i],
          to = homolog,
          weight = weight,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  # If no edges, return empty graph
  if(nrow(edges) == 0) {
    message("  - No edges found in marker network")
    return(graph.empty())
  }
  
  # Aggregate duplicate edges by summing weights
  edges_agg <- aggregate(weight ~ from + to, data = edges, sum)
  
  # Create graph
  g <- graph_from_data_frame(edges_agg, directed = FALSE)
  
  # 安全地为图的顶点设置属性
  # 首先获取图中所有顶点名称
  vertex_names <- V(g)$name
  
  # 准备属性数据框
  vertex_attrs <- data.frame(
    name = vertex_names,
    species = NA_character_,
    chromosome = NA_character_,
    group = NA_character_,
    stringsAsFactors = FALSE
  )
  
  # 解析顶点名称以提取物种和染色体
  split_names <- strsplit(vertex_names, ":")
  vertex_attrs$species <- sapply(split_names, function(x) x[1])
  vertex_attrs$chromosome <- sapply(split_names, function(x) {
    if(length(x) > 1) paste(x[-1], collapse=":") else NA_character_
  })
  
  # 从markers数据补充group信息
  for(i in 1:nrow(markers)) {
    marker_id <- node_ids[i]
    if(marker_id %in% vertex_names) {
      idx <- which(vertex_attrs$name == marker_id)
      if(length(idx) > 0) {
        vertex_attrs$group[idx] <- as.character(markers$group[i])
      }
    }
  }
  
  # 使用传递整个数据框的方式安全设置顶点属性
  g <- set_vertex_attr(g, "species", value = vertex_attrs$species)
  g <- set_vertex_attr(g, "chromosome", value = vertex_attrs$chromosome)
  g <- set_vertex_attr(g, "group", value = vertex_attrs$group)
  
  message(paste("  - Created network with", vcount(g), "nodes and", ecount(g), "edges"))
  
  return(g)
}

detect_communities <- function(g) {
  # Detect communities in the marker network using multiple algorithms
  message("  - Detecting communities in marker network...")
  
  # If empty graph, return empty list
  if(vcount(g) == 0) {
    return(list())
  }
  
  # List of algorithms to try (in order of preference)
  methods <- config$community_detection_methods
  
  # Try each method and calculate modularity
  results <- list()
  
  for(method in methods) {
    tryCatch({
      comm <- NULL
      
      if(method == "louvain") {
        comm <- cluster_louvain(g, weights = E(g)$weight)
      } else if(method == "fast_greedy") {
        comm <- cluster_fast_greedy(g, weights = E(g)$weight)
      } else if(method == "label_prop") {
        comm <- cluster_label_prop(g, weights = E(g)$weight)
      } else if(method == "walktrap") {
        comm <- cluster_walktrap(g, weights = E(g)$weight)
      }
      
      if(!is.null(comm)) {
        mod <- modularity(g, membership(comm), weights = E(g)$weight)
        results[[method]] <- list(
          comm = comm,
          modularity = mod
        )
        message(paste("    -", method, "algorithm: found", max(membership(comm)), 
                    "communities, modularity =", round(mod, 3)))
      }
    }, error = function(e) {
      message(paste("    -", method, "algorithm failed:", e$message))
    })
  }
  
  # If all methods failed, create singleton communities
  if(length(results) == 0) {
    message("    - All community detection algorithms failed, creating singleton communities...")
    membership_vec <- 1:vcount(g)
    comm <- make_clusters(g, membership = membership_vec)
    results[["singleton"]] <- list(
      comm = comm,
      modularity = 0
    )
  }
  
  # Select best result based on modularity
  best_method <- names(results)[which.max(sapply(results, function(x) x$modularity))]
  best_comm <- results[[best_method]]$comm
  message(paste("    - Selected", best_method, "algorithm as best community detection method"))
  
  # Extract communities as lists of markers
  communities <- list()
  
  for(i in 1:max(membership(best_comm))) {
    # Get nodes in this community
    community_nodes <- which(membership(best_comm) == i)
    
    # Skip if empty
    if(length(community_nodes) == 0) {
      next
    }
    
    # Get node names
    node_names <- V(g)$name[community_nodes]
    
    # Extract species and chromosome from node names
    species_chr <- strsplit(node_names, ":")
    species <- sapply(species_chr, function(x) x[1])
    chromosomes <- sapply(species_chr, function(x) paste(x[-1], collapse=":"))
    
    # Get groups
    groups <- V(g)$group[community_nodes]
    
    # Create community data frame
    community_df <- data.frame(
      species = species,
      chromosome = chromosomes,
      group = groups,
      community = i,
      stringsAsFactors = FALSE
    )
    
    communities[[i]] <- community_df
  }
  
  # Filter small communities
  small_communities <- which(sapply(communities, nrow) < config$small_community_threshold)
  if(length(small_communities) > 0) {
    communities <- communities[-small_communities]
    message(paste("    - Filtered out", length(small_communities), "small communities"))
  }
  
  message(paste("    - Found", length(communities), "communities"))
  
  return(list(
    communities = communities,
    modularity = results[[best_method]]$modularity,
    method = best_method
  ))
}

get_community_confidence <- function(comm_markers) {
  # Calculate confidence score for a community
  # Based on number of species and consistency of groups
  
  # Number of unique species
  n_species <- length(unique(comm_markers$species))
  
  # Consistency of groups (proportion of markers from the most common group)
  group_counts <- table(comm_markers$group)
  max_group_count <- max(group_counts)
  group_consistency <- max_group_count / nrow(comm_markers)
  
  # Combine metrics using configurable weights
  confidence <- (config$conservation_weight * group_consistency) + 
              (config$species_coverage_weight * min(1, n_species / 5))
  
  return(confidence)
}

get_most_common_pattern <- function(group_markers) {
  # Get the most common chromosome pattern for a group
  
  # Count occurrences of each chromosome
  chr_counts <- table(group_markers$chromosome)
  
  # Sort by frequency
  sorted_chrs <- names(sort(chr_counts, decreasing = TRUE))
  
  # Return the most common chromosomes (up to 3)
  return(sorted_chrs[1:min(3, length(sorted_chrs))])
}

# ==============================
# Event Reversal Function (New)
# ==============================

reconstruct_by_event_reversal <- function(context, conserved_groups, event_results) {
  # Reconstruct ancestral states by reversing evolutionary events
  message("    - Reconstructing using event reversal logic...")
  
  # Initialize result
  result <- list()
  
  # Get all events relevant to this node and its descendants
  relevant_events <- context$node_events
  
  # Skip if no events
  if(nrow(relevant_events) == 0) {
    message("      - No relevant events found for event reversal")
    return(result)
  }
  
  # Process fusion events (indicating that ancestral state had more chromosomes)
  fusion_events <- relevant_events[relevant_events$event_type == "fusion", ]
  
  if(nrow(fusion_events) > 0) {
    message(paste("      - Processing", nrow(fusion_events), "fusion events for reversal"))
    
    for(i in 1:nrow(fusion_events)) {
      event <- fusion_events[i, ]
      
      # If event involves the current node as parent, it means two chromosomes in this node
      # were fused in the child node
      if(event$parent_node == context$node_id) {
        # Create ancestral groups based on the fusion event
        group_name <- paste0("reverse_fusion_", i)
        
        # Extract group if available
        group <- if(!is.na(event$group)) event$group else paste0("Group_", i)
        
        # Create two ancestral chromosomes from this fusion event
        for(j in 1:2) {
          result[[paste0(group_name, "_", j)]] <- list(
            group = group,
            conservation_score = 0.7,  # Moderate confidence based on event reversal
            coverage = NA,  # Coverage not applicable for event-based reconstruction
            confidence = 0.7,
            method = "event_reversal",
            chromosomes = paste0("ancestral_", group, "_", j),
            event_type = "fusion_reversal",
            event_id = i
          )
        }
      }
    }
  }
  
  # Process fission events (indicating that ancestral state had fewer chromosomes)
  fission_events <- relevant_events[relevant_events$event_type == "fission", ]
  
  if(nrow(fission_events) > 0) {
    message(paste("      - Processing", nrow(fission_events), "fission events for reversal"))
    
    for(i in 1:nrow(fission_events)) {
      event <- fission_events[i, ]
      
      # If event involves the current node as parent, it means one chromosome in this node
      # was split in the child node
      if(event$parent_node == context$node_id) {
        # Create ancestral group based on the fission event
        group_name <- paste0("reverse_fission_", i)
        
        # Extract group if available
        group <- if(!is.na(event$group)) event$group else paste0("Group_", i)
        
        result[[group_name]] <- list(
          group = group,
          conservation_score = 0.65,  # Slightly lower confidence than fusion reversal
          coverage = NA,
          confidence = 0.65,
          method = "event_reversal",
          chromosomes = paste0("ancestral_", group),
          event_type = "fission_reversal",
          event_id = i
        )
      }
    }
  }
  
  # 修复后的全基因组复制(WGD)事件检测代码
  if(config$allow_wgd_events) {
    # 获取所有子节点
    child_nodes <- context$descendants
    
    for(child_node in child_nodes) {
      # 安全获取子节点的染色体数量
      child_count <- NA
      
      # 先检查event_results的结构
      if(!is.null(event_results$ancestral_counts) && 
         !is.null(event_results$ancestral_counts$node_counts)) {
        
        # 尝试多种可能的数据结构
        if(is.list(event_results$ancestral_counts$node_counts)) {
          # 如果是列表结构
          child_node_str <- as.character(child_node)
          for(i in seq_along(event_results$ancestral_counts$node_counts)) {
            if(any(names(event_results$ancestral_counts$node_counts[[i]]) == child_node_str)) {
              child_count <- event_results$ancestral_counts$node_counts[[i]][child_node_str]
              break
            }
          }
        } else if(is.vector(event_results$ancestral_counts$node_counts) && 
                 !is.null(names(event_results$ancestral_counts$node_counts))) {
          # 如果是命名向量
          child_node_str <- as.character(child_node)
          if(child_node_str %in% names(event_results$ancestral_counts$node_counts)) {
            child_count <- event_results$ancestral_counts$node_counts[child_node_str]
          }
        } else if(is.data.frame(event_results$ancestral_counts$node_counts)) {
          # 如果是数据框
          child_idx <- which(event_results$ancestral_counts$node_counts$node_id == child_node)
          if(length(child_idx) > 0) {
            child_count <- event_results$ancestral_counts$node_counts$count[child_idx]
          }
        }
      }
      
      # 检测可能的全基因组复制事件
      if(!is.na(child_count) && !is.na(context$chr_count) && 
         (child_count / context$chr_count) >= config$max_chromosome_ratio) {
        message(paste("      - Detected potential WGD event from node", context$node_id, "to", child_node))
        
        # 使用已有的结果作为基础创建可能的前WGD状态
        group_names <- unique(sapply(result, function(x) x$group))
        
        if(length(group_names) > 0) {
          for(group in group_names) {
            result[[paste0("pre_wgd_", group)]] <- list(
              group = group,
              conservation_score = 0.6,  # 较低的信度用于WGD推断
              coverage = NA,
              confidence = 0.6,
              method = "event_reversal",
              chromosomes = paste0("ancestral_pre_wgd_", group),
              event_type = "wgd_reversal",
              event_id = paste0(context$node_id, "_to_", child_node)
            )
          }
        }
      }
    }
  }
  
  message(paste("      - Identified", length(result), "potential ancestral groups through event reversal"))
  
  return(result)
}

# ==============================
# Bottom-Up Reconstruction Functions
# ==============================

reconstruct_bottom_up <- function(node_contexts, linkage_markers, tree, conserved_groups) {
  # Perform bottom-up reconstruction along the phylogenetic tree
  message("Performing bottom-up reconstruction along the phylogenetic tree...")
  
  # Initialize result
  result <- list()
  
  # Process nodes in bottom-up order
  for(idx in node_contexts$bottom_up_order) {
    context <- node_contexts$contexts[[idx]]
    node_id <- context$node_id
    
    message(paste("  - Processing node", node_id, "(", context$node_label, ")..."))
    
    # Filter markers for species under this node
    node_species <- context$tip_species
    node_markers <- linkage_markers[linkage_markers$species %in% node_species, ]
    
    # Use parsimony to reconstruct ancestral linkage
    parsimony_linkage <- reconstruct_by_parsimony(node_markers, context, conserved_groups)
    
    # Use network analysis to reconstruct ancestral linkage
    network_linkage <- reconstruct_by_network(node_markers, context, conserved_groups)
    
    # Use event reversal to reconstruct ancestral linkage (NEW)
    event_linkage <- reconstruct_by_event_reversal(context, conserved_groups, node_contexts)
    
    # Integrate with child nodes that have already been processed
    child_integration <- integrate_child_nodes(result, context, conserved_groups)
    
    # Combine results from different methods
    expected_count <- context$chr_count
    combined_linkage <- combine_reconstruction_methods(
      parsimony_linkage, network_linkage, event_linkage, child_integration, expected_count
    )
    
    # Add to result
    result[[as.character(node_id)]] <- list(
      node_id = node_id,
      node_label = context$node_label,
      linkage_groups = combined_linkage$groups,
      confidence = combined_linkage$confidence,
      method = "bottom_up",
      taxon = context$taxon
    )
    
    # Apply taxon-specific adjustment if available
    rate_factor <- config$taxon_specific_rates[[result[[as.character(node_id)]]$taxon]]
    if(!is.null(rate_factor)) {
      # Adjust confidence based on taxon-specific evolutionary rate
      for(group_name in names(result[[as.character(node_id)]]$linkage_groups)) {
        result[[as.character(node_id)]]$linkage_groups[[group_name]]$confidence <- 
          result[[as.character(node_id)]]$linkage_groups[[group_name]]$confidence * 
          (1 + (1 - rate_factor) * 0.2)  # Small adjustment based on rate
      }
      result[[as.character(node_id)]]$confidence <- 
        sapply(result[[as.character(node_id)]]$linkage_groups, function(g) g$confidence)
    }
  }
  
  return(result)
}

reconstruct_by_parsimony <- function(markers, context, conserved_groups) {
  # Reconstruct ancestral linkage using the parsimony principle
  message("    - Reconstructing using parsimony principle...")
  
  # Initialize result
  result <- list()
  
  # Get all conserved groups
  if(length(conserved_groups$conserved_groups) == 0) {
    message("      - No conserved groups available for parsimony reconstruction")
    return(result)
  }
  
  # Get species under this node
  tip_species <- context$tip_species
  
  # Process each conserved group
  for(i in seq_along(conserved_groups$conserved_groups)) {
    group <- conserved_groups$conserved_groups[[i]]
    group_name <- names(conserved_groups$conserved_groups)[i]
    
    # Count how many species under this node have this group
    species_with_group <- names(group$species_chromosomes)
    coverage <- sum(tip_species %in% species_with_group) / length(tip_species)
    
    # Apply parsimony threshold
    if(coverage >= config$parsimony_coverage_threshold) {
      # This group likely existed in the ancestor
      result[[paste0("parsimony_", group_name)]] <- list(
        group = group_name,
        conservation_score = group$conservation_score,
        coverage = coverage,
        confidence = group$conservation_score * coverage,
        method = "parsimony",
        chromosomes = paste0("ancestral_", group_name)
      )
    }
  }
  
  message(paste("      - Found", length(result), "ancestral linkage groups using parsimony"))
  
  return(result)
}

reconstruct_by_network <- function(markers, context, conserved_groups) {
  # Reconstruct ancestral linkage using network analysis
  message("    - Reconstructing using network analysis...")
  
  # Initialize result
  result <- list()
  
  # Create marker network
  network <- create_marker_network(markers)
  
  # Skip if empty network
  if(vcount(network) == 0) {
    message("      - Empty network, no groups identified")
    return(result)
  }
  
  # Detect communities
  communities <- detect_communities(network)
  
  # Skip if no communities found
  if(length(communities$communities) == 0) {
    message("      - No communities detected in network")
    return(result)
  }
  
  # Create linkage groups from communities
  for(i in seq_along(communities$communities)) {
    comm <- communities$communities[[i]]
    
    # Get most common group in this community
    group_counts <- table(comm$group)
    if(length(group_counts) > 0) {
      main_group <- names(sort(group_counts, decreasing = TRUE))[1]
    } else {
      main_group <- paste0("Community", i)
    }
    
    # Calculate community confidence
    confidence <- get_community_confidence(comm)
    
    # Adjust confidence based on modularity
    adjusted_confidence <- confidence * (0.5 + (0.5 * communities$modularity))
    
    # Create linkage group
    result[[paste0("network_", i)]] <- list(
      group = main_group,
      conservation_score = ifelse(main_group %in% names(conserved_groups$conserved_groups),
                               conserved_groups$conserved_groups[[main_group]]$conservation_score,
                               0.5),
      coverage = length(unique(comm$species)) / length(context$tip_species),
      confidence = adjusted_confidence,
      method = paste0("network_", communities$method),
      chromosomes = paste0("ancestral_network_", i)
    )
  }
  
  message(paste("      - Found", length(result), "ancestral linkage groups using network analysis"))
  
  return(result)
}

integrate_child_nodes <- function(current_results, context, conserved_groups) {
  # Integrate already processed child nodes
  message("    - Integrating information from child nodes...")
  
  # Initialize result
  result <- list()
  
  # Get child nodes
  child_nodes <- context$descendants[context$descendants > length(conserved_groups$species)]
  
  # Skip if no child nodes or no child with results
  if(length(child_nodes) == 0) {
    message("      - No internal child nodes to integrate")
    return(result)
  }
  
  # Filter child nodes that have already been processed
  processed_children <- child_nodes[sapply(child_nodes, function(node) as.character(node) %in% names(current_results))]
  
  if(length(processed_children) == 0) {
    message("      - No processed child nodes to integrate")
    return(result)
  }
  
  message(paste("      - Integrating", length(processed_children), "child nodes"))
  
  # Collect all linkage groups from children
  child_groups <- list()
  
  for(child_id in processed_children) {
    child_result <- current_results[[as.character(child_id)]]
    
    # Skip if no linkage groups
    if(length(child_result$linkage_groups) == 0) {
      next
    }
    
    for(group_name in names(child_result$linkage_groups)) {
      group <- child_result$linkage_groups[[group_name]]
      
      # Use original group identifier
      orig_group <- group$group
      
      # Add to collection with reduced confidence
      if(orig_group %in% names(child_groups)) {
        # Update if better confidence
        if(group$confidence > child_groups[[orig_group]]$confidence) {
          child_groups[[orig_group]] <- list(
            group = orig_group,
            conservation_score = group$conservation_score,
            coverage = group$coverage,
            confidence = group$confidence * 0.95,  # Small reduction in confidence
            method = "child_integration",
            chromosomes = paste0("ancestral_", orig_group)
          )
        }
      } else {
        # Add new group
        child_groups[[orig_group]] <- list(
          group = orig_group,
          conservation_score = group$conservation_score,
          coverage = group$coverage,
          confidence = group$confidence * 0.95,  # Small reduction in confidence
          method = "child_integration",
          chromosomes = paste0("ancestral_", orig_group)
        )
      }
    }
  }
  
  # Apply events between child nodes and current node
  for(child_id in processed_children) {
    # Get events between this child and parent
    events <- context$node_events[
      (context$node_events$parent_node == context$node_id & context$node_events$child_node == child_id) |
      (context$node_events$parent_node == child_id & context$node_events$child_node == context$node_id),
    ]
    
    # Apply event logic to adjust groups
    if(nrow(events) > 0) {
      message(paste("      - Found", nrow(events), "events between node", context$node_id, "and child", child_id))
      
      # Process each event
      for(i in 1:nrow(events)) {
        event <- events[i, ]
        
        if(event$event_type == "fusion" && event$parent_node == context$node_id) {
          # Fusion in child means two groups in parent
          if(!is.na(event$group) && event$group %in% names(child_groups)) {
            # Split this group into two for the parent
            orig_group <- child_groups[[event$group]]
            
            for(j in 1:2) {
              result[[paste0("fusion_parent_", event$group, "_", j)]] <- list(
                group = event$group,
                conservation_score = orig_group$conservation_score * 0.9,
                coverage = orig_group$coverage,
                confidence = orig_group$confidence * config$fission_confidence_penalty,
                method = "fusion_reversal_integration",
                chromosomes = paste0("ancestral_", event$group, "_", j)
              )
            }
            
            # Remove original group since we've split it
            child_groups[[event$group]] <- NULL
          }
        } else if(event$event_type == "fission" && event$parent_node == context$node_id) {
          # Fission in child means combined group in parent
          if(!is.na(event$group)) {
            # Try to find the split groups
            split_groups <- grep(paste0("^", event$group, "_"), names(child_groups), value = TRUE)
            
            if(length(split_groups) >= 2) {
              # Merge these groups for the parent
              merged_confidence <- mean(sapply(split_groups, function(g) child_groups[[g]]$confidence))
              
              result[[paste0("fission_parent_", event$group)]] <- list(
                group = event$group,
                conservation_score = mean(sapply(split_groups, function(g) child_groups[[g]]$conservation_score)),
                coverage = mean(sapply(split_groups, function(g) child_groups[[g]]$coverage)),
                confidence = merged_confidence * config$fusion_confidence_penalty,
                method = "fission_reversal_integration",
                chromosomes = paste0("ancestral_", event$group)
              )
              
              # Remove the split groups
              child_groups[split_groups] <- NULL
            }
          }
        }
      }
    }
  }
  
  # Add remaining child groups to result
  for(group_name in names(child_groups)) {
    result[[paste0("integrated_", group_name)]] <- child_groups[[group_name]]
  }
  
  message(paste("      - Integrated", length(result), "linkage groups from child nodes"))
  
  return(result)
}

combine_reconstruction_methods <- function(parsimony_result, network_result, event_result, child_integration, expected_count) {
  # Combine different reconstruction methods for bottom-up approach
  
  # Collect all groups
  all_groups <- c(parsimony_result, network_result, event_result, child_integration)
  
  # If no groups, return empty list
  if(length(all_groups) == 0) {
    return(list(
      groups = list(),
      confidence = list()
    ))
  }
  
  # Get group names
  group_names <- sapply(all_groups, function(g) g$group)
  
  # Find duplicates
  duplicated_groups <- group_names[duplicated(group_names)]
  
  # Initialize combined results
  combined_groups <- list()
  confidence_scores <- list()
  
  # Process unique groups
  unique_groups <- unique(group_names)
  
  for(group in unique_groups) {
    # Find all instances of this group
    instances <- which(group_names == group)
    
    if(length(instances) == 1) {
      # Single instance - add as is
      combined_groups[[paste0("combined_", group)]] <- all_groups[[instances]]
      confidence_scores[[paste0("combined_", group)]] <- all_groups[[instances]]$confidence
    } else {
      # Multiple instances - combine
      
      # Get all methods that found this group
      methods <- sapply(all_groups[instances], function(g) g$method)
      
      # Calculate combined confidence (boost if found by multiple methods)
      combined_confidence <- mean(sapply(all_groups[instances], function(g) g$confidence)) * 
                            min(1.5, 1 + (length(instances) * 0.1))
      
      # Get chromosomes from all instances
      all_chrs <- unlist(lapply(all_groups[instances], function(g) g$chromosomes))
      
      combined_groups[[paste0("combined_", group)]] <- list(
        group = group,
        conservation_score = max(sapply(all_groups[instances], function(g) g$conservation_score)),
        coverage = max(sapply(all_groups[instances], function(g) g$coverage)),
        confidence = combined_confidence,
        method = paste(methods, collapse="+"),
        chromosomes = unique(all_chrs)
      )
      
      confidence_scores[[paste0("combined_", group)]] <- combined_confidence
    }
  }
  
  # Adjust number of groups to match expected chromosome count
  if(!is.na(expected_count) && length(combined_groups) > expected_count) {
    # Sort by confidence
    sorted_idx <- order(unlist(confidence_scores), decreasing = TRUE)
    
    # Keep top groups
    keep_idx <- sorted_idx[1:min(length(sorted_idx), expected_count)]
    combined_groups <- combined_groups[keep_idx]
    confidence_scores <- confidence_scores[keep_idx]
  }
  
  return(list(
    groups = combined_groups,
    confidence = confidence_scores
  ))
}

# ==============================
# Top-Down Reconstruction Functions
# ==============================

reconstruct_top_down <- function(node_contexts, linkage_markers, tree, conserved_groups) {
  # Perform top-down reconstruction along the phylogenetic tree
  message("Performing top-down reconstruction along the phylogenetic tree...")
  
  # Initialize result
  result <- list()
  
  # Get root node
  root_node <- node_contexts$contexts[[node_contexts$top_down_order[1]]]$node_id
  root_context <- node_contexts$contexts[[node_contexts$top_down_order[1]]]
  
  # Initialize root node first
  root_linkage <- initialize_root_linkage(root_context, conserved_groups)
  
  result[[as.character(root_node)]] <- list(
    node_id = root_node,
    node_label = root_context$node_label,
    linkage_groups = root_linkage$groups,
    confidence = root_linkage$confidence,
    method = "top_down",
    taxon = root_context$taxon
  )
  
  # Process other nodes in top-down order
  for(idx in node_contexts$top_down_order[-1]) {
    context <- node_contexts$contexts[[idx]]
    node_id <- context$node_id
    
    message(paste("  - Processing node", node_id, "(", context$node_label, ")..."))
    
    # Get parent node results
    parent_id <- context$parent_id
    
    if(!is.na(parent_id) && as.character(parent_id) %in% names(result)) {
      parent_linkage <- result[[as.character(parent_id)]]$linkage_groups
      
      # Apply events to create linkage groups for this node
      node_linkage <- apply_events_to_linkage(parent_linkage, context, conserved_groups)
      
      # Add to result
      result[[as.character(node_id)]] <- list(
        node_id = node_id,
        node_label = context$node_label,
        linkage_groups = node_linkage$groups,
        confidence = node_linkage$confidence,
        method = "top_down",
        taxon = context$taxon
      )
      
      # Apply taxon-specific adjustment if available
      rate_factor <- config$taxon_specific_rates[[result[[as.character(node_id)]]$taxon]]
      if(!is.null(rate_factor)) {
        # Adjust confidence based on taxon-specific evolutionary rate
        for(group_name in names(result[[as.character(node_id)]]$linkage_groups)) {
          result[[as.character(node_id)]]$linkage_groups[[group_name]]$confidence <- 
            result[[as.character(node_id)]]$linkage_groups[[group_name]]$confidence * 
            (1 + (1 - rate_factor) * 0.2)  # Small adjustment based on rate
        }
        result[[as.character(node_id)]]$confidence <- 
          sapply(result[[as.character(node_id)]]$linkage_groups, function(g) g$confidence)
      }
    } else {
      # No parent results, skip this node
      message("    - No parent node results available, skipping...")
    }
  }
  
  return(result)
}

initialize_root_linkage <- function(root_context, conserved_groups) {
  # Initialize linkage groups for the root node
  message("  - Initializing root node linkage groups...")
  
  # Initialize results
  result <- list()
  confidence <- list()
  
  # Sort conserved groups by conservation score
  sorted_indices <- order(sapply(conserved_groups$conserved_groups, function(g) g$conservation_score), decreasing = TRUE)
  
  # Determine how many groups to include
  target_count <- ifelse(!is.na(root_context$chr_count), 
                       root_context$chr_count, 
                       min(10, length(conserved_groups$conserved_groups)))
  
  # Select top groups
  top_groups <- min(length(sorted_indices), target_count)
  
  for(i in 1:top_groups) {
    idx <- sorted_indices[i]
    group <- conserved_groups$conserved_groups[[idx]]
    group_name <- names(conserved_groups$conserved_groups)[idx]
    
    result[[paste0("root_", i)]] <- list(
      group = group_name,
      conservation_score = group$conservation_score,
      coverage = length(unique(names(group$species_chromosomes))) / length(root_context$tip_species),
      confidence = group$conservation_score,
      method = "root_initialization",
      chromosomes = paste0("ancestral_", i)
    )
    
    confidence[[paste0("root_", i)]] <- group$conservation_score
  }
  
  message(paste("      - Initialized", length(result), "linkage groups for root node"))
  
  return(list(
    groups = result,
    confidence = confidence
  ))
}

apply_events_to_linkage <- function(parent_linkage, context, conserved_groups) {
  # Apply evolutionary events to parent linkage groups
  message("    - Applying evolutionary events to parent linkage groups...")
  
  # Start with parent linkage
  result <- parent_linkage
  confidence <- sapply(parent_linkage, function(g) g$confidence)
  
  # Get chromosome events on this branch
  node_events <- context$node_events
  
  # If no events, inherit parent linkage with slightly reduced confidence
  if(nrow(node_events) == 0) {
    message("      - No specific events on this branch, inheriting parent linkage")
    
    # Reduce confidence slightly for inherited groups
    for(i in seq_along(result)) {
      result[[i]]$confidence <- result[[i]]$confidence * 0.95  # Small reduction in confidence
      result[[i]]$method <- "inherited"
      confidence[i] <- result[[i]]$confidence
    }
    
    return(list(
      groups = result,
      confidence = confidence
    ))
  }
  
  # Process fusion events
  fusion_events <- node_events[node_events$event_type == "fusion", ]
  if(nrow(fusion_events) > 0) {
    message(paste("      - Processing", nrow(fusion_events), "fusion events"))
    
    # For each fusion event, try to find the affected groups and merge them
    for(i in 1:nrow(fusion_events)) {
      event <- fusion_events[i, ]
      
      # If group information is available, try to find affected linkage groups
      if(!is.na(event$group)) {
        # Find linkage groups matching this group
        affected_idx <- which(sapply(result, function(g) g$group == event$group))
        
        if(length(affected_idx) >= 2) {
          # Merge the first two groups
          merged_name <- paste0("merged_", names(result)[affected_idx[1]], "_", names(result)[affected_idx[2]])
          result[[merged_name]] <- list(
            group = event$group,
            conservation_score = mean(sapply(result[affected_idx[1:2]], function(g) g$conservation_score)),
            coverage = mean(sapply(result[affected_idx[1:2]], function(g) g$coverage)),
            confidence = mean(sapply(result[affected_idx[1:2]], function(g) g$confidence)) * config$fusion_confidence_penalty,
            method = "fusion_event",
            chromosomes = c(result[[affected_idx[1]]]$chromosomes, result[[affected_idx[2]]]$chromosomes)
          )
          
          # Update confidence
          confidence[[merged_name]] <- result[[merged_name]]$confidence
          
          # Remove the original groups
          result <- result[-affected_idx[1:2]]
          confidence <- confidence[-affected_idx[1:2]]
        }
      }
    }
  }
  
  # Process fission events
  fission_events <- node_events[node_events$event_type == "fission", ]
  if(nrow(fission_events) > 0) {
    message(paste("      - Processing", nrow(fission_events), "fission events"))
    
    # For each fission event, try to find the affected group and split it
    for(i in 1:nrow(fission_events)) {
      event <- fission_events[i, ]
      
      # If group information is available, try to find affected linkage group
      if(!is.na(event$group)) {
        # Find linkage group matching this group
        affected_idx <- which(sapply(result, function(g) g$group == event$group))
        
        if(length(affected_idx) > 0) {
          # Split the group into two
          for(j in 1:2) {
            split_name <- paste0("split", j, "_", names(result)[affected_idx[1]])
            result[[split_name]] <- list(
              group = event$group,
              conservation_score = result[[affected_idx[1]]]$conservation_score * 0.9,
              coverage = result[[affected_idx[1]]]$coverage,
              confidence = result[[affected_idx[1]]]$confidence * config$fission_confidence_penalty,
              method = "fission_event",
              chromosomes = paste0(result[[affected_idx[1]]]$chromosomes, "_split", j)
            )
            
            # Update confidence
            confidence[[split_name]] <- result[[split_name]]$confidence
          }
          
          # Remove the original group
          result <- result[-affected_idx[1]]
          confidence <- confidence[-affected_idx[1]]
        }
      }
    }
  }
  
  # Check for WGD events if enabled
  if(config$allow_wgd_events) {
    wgd_events <- node_events[node_events$event_type == "wgd", ]
    if(nrow(wgd_events) > 0) {
      message(paste("      - Processing", nrow(wgd_events), "whole genome duplication events"))
      
      # Create copy of all linkage groups
      original_result <- result
      original_confidence <- confidence
      
      for(group_name in names(original_result)) {
        wgd_name <- paste0("wgd_", group_name)
        result[[wgd_name]] <- original_result[[group_name]]
        result[[wgd_name]]$method <- "wgd_event"
        result[[wgd_name]]$confidence <- original_result[[group_name]]$confidence * 0.85
        confidence[[wgd_name]] <- result[[wgd_name]]$confidence
      }
    }
  }
  
  message(paste("      - Final linkage group count after events:", length(result)))
  
  return(list(
    groups = result,
    confidence = confidence
  ))
}

# ==============================
# Integration Functions
# ==============================

integrate_reconstruction_approaches <- function(bottom_up_linkage, top_down_linkage, node_contexts) {
  # Integrate bottom-up and top-down approaches
  message("Integrating bottom-up and top-down reconstruction approaches...")
  
  # Initialize integrated results
  integrated_linkage <- list()
  
  # Process each node
  for(context in node_contexts$contexts) {
    node_id <- context$node_id
    node_label <- context$node_label
    
    message(paste("  - Integrating results for node", node_id, "(", node_label, ")..."))
    
    # Get bottom-up and top-down results for this node
    bu_result <- bottom_up_linkage[[as.character(node_id)]]
    td_result <- top_down_linkage[[as.character(node_id)]]
    
    # If both approaches have results, integrate them
    if(!is.null(bu_result) && !is.null(td_result)) {
      integrated_groups <- integrate_node_approaches(bu_result$linkage_groups, td_result$linkage_groups, context)
      
      integrated_linkage[[as.character(node_id)]] <- list(
        node_id = node_id,
        node_label = node_label,
        linkage_groups = integrated_groups$groups,
        confidence = integrated_groups$confidence,
        method = "integrated",
        bottom_up = bu_result,
        top_down = td_result,
        agreement = integrated_groups$agreement
      )
    } else if(!is.null(bu_result)) {
      # Only bottom-up results available
      message("    - Only bottom-up results available for this node")
      integrated_linkage[[as.character(node_id)]] <- list(
        node_id = node_id,
        node_label = node_label,
        linkage_groups = bu_result$linkage_groups,
        confidence = sapply(bu_result$linkage_groups, function(g) g$confidence),
        method = "bottom_up_only",
        bottom_up = bu_result,
        top_down = NULL,
        agreement = NA
      )
    } else if(!is.null(td_result)) {
      # Only top-down results available
      message("    - Only top-down results available for this node")
      integrated_linkage[[as.character(node_id)]] <- list(
        node_id = node_id,
        node_label = node_label,
        linkage_groups = td_result$linkage_groups,
        confidence = sapply(td_result$linkage_groups, function(g) g$confidence),
        method = "top_down_only",
        bottom_up = NULL,
        top_down = td_result,
        agreement = NA
      )
    } else {
      # No results available
      message("    - No reconstruction results available for this node")
      integrated_linkage[[as.character(node_id)]] <- list(
        node_id = node_id,
        node_label = node_label,
        linkage_groups = list(),
        confidence = list(),
        method = "none",
        bottom_up = NULL,
        top_down = NULL,
        agreement = NA
      )
    }
    
    # Adjust number of groups to match expected chromosome count if necessary
    if(!is.na(context$chr_count) && 
       length(integrated_linkage[[as.character(node_id)]]$linkage_groups) != context$chr_count) {
      message(paste("    - Adjusting group count from", 
                  length(integrated_linkage[[as.character(node_id)]]$linkage_groups), 
                  "to match expected chromosome count of", context$chr_count))
      
      integrated_linkage[[as.character(node_id)]] <- adjust_group_count(
        integrated_linkage[[as.character(node_id)]],
        context$chr_count
      )
    }
  }
  
  message(paste("  - Completed integration for", length(integrated_linkage), "nodes"))
  
  return(integrated_linkage)
}

integrate_node_approaches <- function(bottom_up_groups, top_down_groups, context) {
  # Integrate bottom-up and top-down approaches for a single node
  
  # Initialize result
  result <- list()
  confidence <- list()
  
  # Get unique groups from both approaches
  bu_groups <- sapply(bottom_up_groups, function(g) g$group)
  td_groups <- sapply(top_down_groups, function(g) g$group)
  
  # Count common groups
  common_groups <- intersect(bu_groups, td_groups)
  common_count <- length(common_groups)
  
  # Calculate agreement between approaches
  total_groups <- length(unique(c(bu_groups, td_groups)))
  agreement <- ifelse(total_groups > 0, common_count / total_groups, 0)
  
  message(paste("    - Agreement between approaches:", round(agreement * 100, 1), "%"))
  
  # Process common groups first (highest confidence)
  for(group in common_groups) {
    # Find corresponding groups in bottom-up and top-down
    bu_idx <- which(bu_groups == group)[1]
    td_idx <- which(td_groups == group)[1]
    
    # Create integrated group with boosted confidence
    bu_conf <- bottom_up_groups[[bu_idx]]$confidence
    td_conf <- top_down_groups[[td_idx]]$confidence
    int_conf <- (bu_conf + td_conf) / 2 * config$common_group_boost
    
    # Create integrated group
    name <- paste0("integrated_", group)
    result[[name]] <- list(
      group = group,
      conservation_score = max(bottom_up_groups[[bu_idx]]$conservation_score, 
                           top_down_groups[[td_idx]]$conservation_score),
      coverage = max(bottom_up_groups[[bu_idx]]$coverage, 
                  ifelse(is.null(top_down_groups[[td_idx]]$coverage), 0, top_down_groups[[td_idx]]$coverage)),
      confidence = int_conf,
      method = "integrated_common",
      chromosomes = unique(c(bottom_up_groups[[bu_idx]]$chromosomes, 
                           top_down_groups[[td_idx]]$chromosomes))
    )
    
    confidence[[name]] <- int_conf
  }
  
  # Process unique bottom-up groups
  unique_bu_groups <- setdiff(bu_groups, common_groups)
  
  for(group in unique_bu_groups) {
    # Find corresponding group
    idx <- which(bu_groups == group)[1]
    
    # Add with slightly reduced confidence
    name <- paste0("bu_", group)
    result[[name]] <- bottom_up_groups[[idx]]
    result[[name]]$confidence <- result[[name]]$confidence * config$unique_group_penalty
    result[[name]]$method <- "integrated_bu_only"
    
    confidence[[name]] <- result[[name]]$confidence
  }
  
  # Process unique top-down groups
  unique_td_groups <- setdiff(td_groups, common_groups)
  
  for(group in unique_td_groups) {
    # Find corresponding group
    idx <- which(td_groups == group)[1]
    
    # Add with slightly reduced confidence
    name <- paste0("td_", group)
    result[[name]] <- top_down_groups[[idx]]
    result[[name]]$confidence <- result[[name]]$confidence * config$unique_group_penalty
    result[[name]]$method <- "integrated_td_only"
    
    confidence[[name]] <- result[[name]]$confidence
  }
  
  return(list(
    groups = result,
    confidence = confidence,
    agreement = agreement
  ))
}

adjust_group_count <- function(node_result, target_count) {
  # Adjust number of linkage groups to match expected chromosome count
  current_count <- length(node_result$linkage_groups)
  
  # If already matching, return unchanged
  if(current_count == target_count) {
    return(node_result)
  }
  
  # Sort groups by confidence
  conf_values <- sapply(node_result$linkage_groups, function(g) g$confidence)
  sorted_idx <- order(conf_values, decreasing = TRUE)
  
  if(current_count > target_count) {
    # Too many groups, remove lowest confidence ones
    keep_idx <- sorted_idx[1:target_count]
    node_result$linkage_groups <- node_result$linkage_groups[keep_idx]
    node_result$confidence <- node_result$confidence[keep_idx]
  } else {
    # Too few groups, add inferred groups
    for(i in 1:(target_count - current_count)) {
      name <- paste0("inferred_", i)
      node_result$linkage_groups[[name]] <- list(
        group = paste0("Inferred", i),
        conservation_score = 0.4,
        coverage = NA,
        confidence = 0.4,
        method = "count_adjustment",
        chromosomes = paste0("ancestral_inferred_", i)
      )
      
      node_result$confidence[[name]] <- 0.4
    }
  }
  
  return(node_result)
}

# ==============================
# Validation Functions
# ==============================

validate_ancestral_linkage <- function(ancestral_linkage, node_contexts, conserved_groups) {
  # Validate the ancestral linkage reconstruction
  message("Validating ancestral linkage reconstruction...")
  
  # Initialize validation dataframe
  validation_df <- data.frame(
    node_id = integer(),
    node_label = character(),
    expected_chr_count = integer(),
    reconstructed_chr_count = integer(),
    count_match = logical(),
    agreement_score = numeric(),
    avg_confidence = numeric(),
    avg_conservation = numeric(),
    event_consistency = logical(),
    stringsAsFactors = FALSE
  )
  
  # Process each node
  for(node_id in names(ancestral_linkage)) {
    node_result <- ancestral_linkage[[node_id]]
    
    # Find context for this node
    context_idx <- which(sapply(node_contexts$contexts, function(c) as.character(c$node_id) == node_id))
    
    if(length(context_idx) > 0) {
      context <- node_contexts$contexts[[context_idx]]
      
      # Count reconstructed groups
      reconstructed_count <- length(node_result$linkage_groups)
      
      # Check if count matches expected
      count_match <- ifelse(!is.na(context$chr_count), 
                          reconstructed_count == context$chr_count, 
                          NA)
      
      # Calculate average confidence
      avg_confidence <- mean(unlist(sapply(node_result$linkage_groups, function(g) g$confidence)))
      
      # Calculate average conservation score
      avg_conservation <- mean(unlist(sapply(node_result$linkage_groups, function(g) g$conservation_score)))
      
      # Check event consistency
      event_consistency <- check_event_consistency(node_result, context, ancestral_linkage)
      
      # Add to validation dataframe
      validation_df <- rbind(validation_df, data.frame(
        node_id = as.integer(node_id),
        node_label = node_result$node_label,
        expected_chr_count = context$chr_count,
        reconstructed_chr_count = reconstructed_count,
        count_match = count_match,
        agreement_score = ifelse(is.null(node_result$agreement), NA, node_result$agreement),
        avg_confidence = avg_confidence,
        avg_conservation = avg_conservation,
        event_consistency = event_consistency,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  message(paste("  - Validation completed for", nrow(validation_df), "nodes"))
  
  return(validation_df)
}

check_event_consistency <- function(node_result, context, ancestral_linkage) {
  # Check if reconstruction is consistent with detected events
  
  # If no events or insufficient data, skip consistency check
  if(is.null(context$node_events) || nrow(context$node_events) == 0) {
    return(NA)
  }
  
  # Get parent node
  parent_id <- context$parent_id
  
  if(is.na(parent_id) || !as.character(parent_id) %in% names(ancestral_linkage)) {
    return(NA)
  }
  
  parent_result <- ancestral_linkage[[as.character(parent_id)]]
  
  # Count groups in current and parent node
  current_count <- length(node_result$linkage_groups)
  parent_count <- length(parent_result$linkage_groups)
  
  # Get events between parent and current node
  events <- context$node_events[context$node_events$parent_node == parent_id & 
                               context$node_events$child_node == context$node_id, ]
  
  if(nrow(events) == 0) {
    events <- context$node_events[context$node_events$parent_node == context$node_id & 
                                 context$node_events$child_node == parent_id, ]
  }
  
  if(nrow(events) == 0) {
    return(NA)
  }
  
  # Count fusion and fission events
  fusion_count <- sum(events$event_type == "fusion")
  fission_count <- sum(events$event_type == "fission")
  
  # Calculate expected count change
  expected_change <- fission_count - fusion_count
  
  # Check if actual change matches expected change
  actual_change <- current_count - parent_count
  
  return(actual_change == expected_change)
}

# ==============================
# Visualization Functions
# ==============================

# 修复后的 visualize_tree_with_linkage_groups 函数
visualize_tree_with_linkage_groups <- function(tree, ancestral_linkage, output_file) {
  # Visualize linkage groups on the phylogenetic tree
  message("Generating phylogenetic tree visualization with linkage groups...")
  
  # 修复：确保我们有正确数量的节点标签
  tree_nnodes <- length(tree$tip.label) + tree$Nnode
  
  # 为每个节点创建唯一的标签
  all_labels <- character(tree_nnodes)
  all_labels[1:length(tree$tip.label)] <- tree$tip.label
  
  # 处理内部节点标签
  for(i in (length(tree$tip.label) + 1):tree_nnodes) {
    node_idx <- i - length(tree$tip.label)
    if(is.null(tree$node.label) || length(tree$node.label) < node_idx || 
       is.na(tree$node.label[node_idx]) || tree$node.label[node_idx] == "") {
      all_labels[i] <- paste0("Node", i)
    } else {
      all_labels[i] <- tree$node.label[node_idx]
    }
  }
  
  # 创建节点数据框
  node_data <- data.frame(
    node = 1:tree_nnodes,
    label = all_labels,
    stringsAsFactors = FALSE
  )
  
  # 添加连锁群信息
  node_data$linkage_count <- 0
  node_data$avg_confidence <- NA
  
  for(node_id in names(ancestral_linkage)) {
    node_idx <- as.integer(node_id)
    if(node_idx > 0 && node_idx <= nrow(node_data)) {
      linkage_groups <- ancestral_linkage[[node_id]]$linkage_groups
      node_data$linkage_count[node_idx] <- length(linkage_groups)
      
      confidence_values <- sapply(linkage_groups, function(g) g$confidence)
      if(length(confidence_values) > 0) {
        node_data$avg_confidence[node_idx] <- mean(confidence_values)
      }
    }
  }
  
  # 创建带注释的树
  p <- ggtree(tree) %<+% node_data
  
  # 添加节点标签
  p <- p + geom_nodelab(aes(label=label)) +
    geom_tiplab(aes(label=label))
  
  # 为内部节点添加连锁群计数
  p <- p + geom_text(aes(label=ifelse(node > length(tree$tip.label) & linkage_count > 0, 
                                paste0(linkage_count, " groups"), ""), 
                    x=x+0.1), size=3, color="blue")
  
  # 使用置信度作为内部节点的颜色
  p <- p + geom_nodepoint(aes(fill=avg_confidence), shape=21, size=5) +
    scale_fill_gradient(low="yellow", high="green", na.value="white", 
                      name="Confidence")
  
  # 保存图形
  pdf(output_file, width=10, height=8)
  print(p)
  dev.off()
  
  message(paste("  - Tree visualization saved to", output_file))
  
  return(p)
}

visualize_node_linkage_groups <- function(node_id, ancestral_linkage, output_dir) {
  # Visualize linkage groups for a specific node
  
  # Create node-specific output directory
  node_dir <- file.path(output_dir, "node_linkage_maps")
  if(!dir.exists(node_dir)) {
    dir.create(node_dir, recursive = TRUE)
  }
  
  # Get linkage groups for this node
  node_result <- ancestral_linkage[[node_id]]
  
  if(is.null(node_result) || length(node_result$linkage_groups) == 0) {
    return(NULL)
  }
  
  # Create output file
  output_file <- file.path(node_dir, paste0("node_", node_id, "_linkage_map.pdf"))
  
  # Prepare data for visualization
  linkage_data <- data.frame(
    group_id = character(),
    group_name = character(),
    confidence = numeric(),
    method = character(),
    stringsAsFactors = FALSE
  )
  
  for(i in seq_along(node_result$linkage_groups)) {
    group <- node_result$linkage_groups[[i]]
    linkage_data <- rbind(linkage_data, data.frame(
      group_id = names(node_result$linkage_groups)[i],
      group_name = group$group,
      confidence = group$confidence,
      method = group$method,
      stringsAsFactors = FALSE
    ))
  }
  
  # Create basic plot
  p <- ggplot(linkage_data, aes(x=group_name, y=confidence, fill=method)) +
    geom_col() +
    coord_flip() +
    labs(title = paste0("Linkage Groups for ", node_result$node_label), 
         x = "Group", 
         y = "Confidence") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2")
  
  # Save plot
  pdf(output_file, width=8, height=6)
  print(p)
  dev.off()
  
  return(output_file)
}

generate_linkage_heatmap <- function(ancestral_linkage, integrated_data, conserved_groups, output_file) {
  # Generate heatmap showing conservation of linkage groups across nodes
  message("Generating linkage group conservation heatmap...")
  
  # Extract all unique groups
  all_groups <- unique(unlist(lapply(ancestral_linkage, function(node) {
    sapply(node$linkage_groups, function(g) g$group)
  })))
  
  # Create matrix for heatmap
  nodes <- names(ancestral_linkage)
  node_labels <- sapply(ancestral_linkage, function(n) n$node_label)
  
  # Create empty matrix
  heatmap_matrix <- matrix(0, nrow=length(nodes), ncol=length(all_groups))
  rownames(heatmap_matrix) <- node_labels
  colnames(heatmap_matrix) <- all_groups
  
  # Fill matrix with confidence values
  for(i in seq_along(nodes)) {
    node_id <- nodes[i]
    node_result <- ancestral_linkage[[node_id]]
    
    for(j in seq_along(all_groups)) {
      group <- all_groups[j]
      
      # Find this group in node results
      group_idx <- which(sapply(node_result$linkage_groups, function(g) g$group) == group)
      
      if(length(group_idx) > 0) {
        # Use highest confidence if multiple matches
        heatmap_matrix[i, j] <- max(sapply(node_result$linkage_groups[group_idx], function(g) g$confidence))
      }
    }
  }
  
  # Generate heatmap
  pdf(output_file, width=10, height=8)
  pheatmap::pheatmap(heatmap_matrix, 
                    cluster_rows = TRUE, 
                    cluster_cols = TRUE, 
                    main = "Linkage Group Conservation Across Nodes",
                    color = colorRampPalette(c("white", "blue", "darkblue"))(100))
  dev.off()
  
  message(paste("  - Linkage heatmap saved to", output_file))
  
  return(heatmap_matrix)
}

# ==============================
# Results Saving Functions
# ==============================

save_ancestral_linkage <- function(ancestral_linkage, validation_df, integrated_data, output_dir) {
  # Save ancestral linkage reconstruction results
  message("Saving ancestral linkage reconstruction results...")
  
  # Create output directory if it doesn't exist
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message(paste("Created output directory:", output_dir))
  }
  
  # Save RData file with complete results
  rdata_file <- file.path(output_dir, "ancestral_linkage_groups.RData")
  saveRDS(ancestral_linkage, rdata_file)
  
  # Create summary table
  summary_data <- data.frame(
    node_id = integer(),
    node_label = character(),
    linkage_group_count = integer(),
    avg_confidence = numeric(),
    method = character(),
    agreement = numeric(),
    stringsAsFactors = FALSE
  )
  
  for(node_id in names(ancestral_linkage)) {
    node_result <- ancestral_linkage[[node_id]]
    
    summary_data <- rbind(summary_data, data.frame(
      node_id = as.integer(node_id),
      node_label = node_result$node_label,
      linkage_group_count = length(node_result$linkage_groups),
      avg_confidence = mean(unlist(sapply(node_result$linkage_groups, function(g) g$confidence))),
      method = node_result$method,
      agreement = ifelse(is.null(node_result$agreement), NA, node_result$agreement),
      stringsAsFactors = FALSE
    ))
  }
  
  # Save summary table
  summary_file <- file.path(output_dir, "ancestral_linkage_summary.csv")
  write.csv(summary_data, summary_file, row.names = FALSE)
  
  # Save validation results
  validation_file <- file.path(output_dir, "ancestral_linkage_validation.csv")
  write.csv(validation_df, validation_file, row.names = FALSE)
  
  # Generate report
  report_file <- file.path(output_dir, "ancestral_linkage_report.txt")
  report_result <- generate_linkage_report(
    ancestral_linkage, 
    summary_data, 
    validation_df, 
    integrated_data,
    report_file
  )
  
  return(list(
    rdata_file = rdata_file,
    summary_file = summary_file,
    validation_file = validation_file,
    report_file = report_file
  ))
}

generate_linkage_report <- function(ancestral_linkage, summary_data, validation_df, integrated_data, report_file) {
  # Generate detailed report of ancestral linkage reconstruction
  
  sink(report_file)
  
  cat("Ancestral Linkage Group Reconstruction Report\n")
  cat("===========================================\n\n")
  
  # Summary statistics
  cat("1. Overall Summary\n")
  cat("----------------\n")
  cat(paste("Total nodes analyzed:", length(ancestral_linkage), "\n"))
  cat(paste("Average linkage groups per node:", round(mean(summary_data$linkage_group_count), 1), "\n"))
  cat(paste("Overall average confidence:", round(mean(summary_data$avg_confidence), 3), "\n"))
  
  # Method distribution
  method_counts <- table(summary_data$method)
  cat("\nMethod distribution:\n")
  for(method in names(method_counts)) {
    cat(paste(" -", method, ":", method_counts[method], "nodes\n"))
  }
  
  # Validation summary
  cat("\n2. Validation Results\n")
  cat("-----------------\n")
  
  if(nrow(validation_df) > 0) {
    # Chromosome count match rate
    count_match_rate <- sum(validation_df$count_match, na.rm = TRUE) / sum(!is.na(validation_df$count_match))
    cat(paste("Chromosome count match rate:", round(count_match_rate * 100, 1), "%\n"))
    
    # Event consistency rate
    event_consistency_rate <- sum(validation_df$event_consistency, na.rm = TRUE) / sum(!is.na(validation_df$event_consistency))
    cat(paste("Event consistency rate:", round(event_consistency_rate * 100, 1), "%\n"))
    
    # Average agreement
    avg_agreement <- mean(validation_df$agreement_score, na.rm = TRUE)
    cat(paste("Average agreement between methods:", round(avg_agreement * 100, 1), "%\n"))
  } else {
    cat("No validation data available.\n")
  }
  
  # Node by node detail
  cat("\n3. Node Details\n")
  cat("-------------\n")
  
  for(node_id in names(ancestral_linkage)) {
    node_data <- ancestral_linkage[[node_id]]
    val_data <- validation_df[validation_df$node_id == as.integer(node_id), ]
    
    cat(paste("\nNode", node_id, "(", node_data$node_label, "):\n"))
    cat(paste(" - Linkage groups:", length(node_data$linkage_groups), "\n"))
    cat(paste(" - Reconstruction method:", node_data$method, "\n"))
    cat(paste(" - Average confidence:", round(mean(unlist(sapply(node_data$linkage_groups, function(g) g$confidence))), 3), "\n"))
    
    if(length(node_data$linkage_groups) > 0) {
      cat(" - Top groups by confidence:\n")
      
      # Sort by confidence
      conf_values <- sapply(node_data$linkage_groups, function(g) g$confidence)
      sorted_idx <- order(conf_values, decreasing = TRUE)
      top_idx <- sorted_idx[1:min(length(sorted_idx), 5)]
      
      for(idx in top_idx) {
        group <- node_data$linkage_groups[[idx]]
        cat(paste("   * Group", group$group, "(", names(node_data$linkage_groups)[idx], "):",
                round(group$confidence, 3), "-", group$method, "\n"))
      }
    }
    
    # Validation info if available
    if(nrow(val_data) > 0) {
      cat(" - Validation:\n")
      cat(paste("   * Expected chromosome count:", 
                ifelse(is.na(val_data$expected_chr_count), "Unknown", val_data$expected_chr_count), "\n"))
      cat(paste("   * Count match:", 
                ifelse(is.na(val_data$count_match), "Unknown", val_data$count_match), "\n"))
      cat(paste("   * Average conservation score:", round(val_data$avg_conservation, 3), "\n"))
      cat(paste("   * Event consistency:", 
                ifelse(is.na(val_data$event_consistency), "Unknown", val_data$event_consistency), "\n"))
    }
    
    # List linkage groups
    if(length(node_data$linkage_groups) > 0) {
      cat(" - Linkage groups:\n")
      
      # Sort by confidence
      group_confidence <- sapply(node_data$linkage_groups, function(g) g$confidence)
      sorted_groups <- names(node_data$linkage_groups)[order(unlist(group_confidence), decreasing = TRUE)]
      
      for(group_name in sorted_groups) {
        group <- node_data$linkage_groups[[group_name]]
        
        cat(paste0("   * ", group_name, " (", group$group, "): confidence = ", 
                 round(group$confidence, 2), "\n"))
        cat(paste0("     Method: ", group$method, "\n"))
        cat(paste0("     Chromosomes: ", paste(group$chromosomes, collapse=", "), "\n"))
      }
    }
    
    cat("\n")
  }
  
  cat("4. Evolutionary Insights\n")
  cat("---------------------\n")
  # Identify trends and patterns across the tree
  
  # Calculate average linkage group count by tree depth
  if(nrow(summary_data) > 0) {
    # Add depth information
    tree_depths <- get_node_depths(integrated_data$tree)
    summary_data$depth <- sapply(summary_data$node_id, function(id) tree_depths[as.character(id)])
    
    # Group by depth
    depth_summary <- aggregate(linkage_group_count ~ depth, summary_data, mean)
    
    if(nrow(depth_summary) > 0) {
      cat("Average linkage group count by tree depth:\n")
      for(i in 1:nrow(depth_summary)) {
        cat(paste(" - Depth", depth_summary$depth[i], ":", round(depth_summary$linkage_group_count[i], 1), "\n"))
      }
    }
  }
  
  cat("\n5. Conclusions\n")
  cat("-------------\n")
  if(length(ancestral_linkage) > 0) {
    # Calculate overall confidence
    overall_confidence <- mean(summary_data$avg_confidence)
    
    cat(paste("Overall reconstruction confidence:", round(overall_confidence, 3), "\n"))
    
    # Calculate agreement between approaches
    if(sum(!is.na(summary_data$agreement)) > 0) {
      overall_agreement <- mean(summary_data$agreement, na.rm = TRUE)
      cat(paste("Average agreement between bottom-up and top-down approaches:", 
                round(overall_agreement * 100, 1), "%\n"))
    }
    
    # Most reliable reconstructions
    high_confidence <- summary_data[summary_data$avg_confidence >= 0.8,]
    if(nrow(high_confidence) > 0) {
      cat(paste("Nodes with high confidence reconstruction (>= 0.8):", nrow(high_confidence), "\n"))
      for(i in 1:min(5, nrow(high_confidence))) {
        cat(paste(" -", high_confidence$node_label[i], "(confidence:", 
                round(high_confidence$avg_confidence[i], 2), ")\n"))
      }
      if(nrow(high_confidence) > 5) {
        cat("   ... and others\n")
      }
    }
  }
  
  sink()
  
  message(paste("  - Ancestral linkage report saved to", report_file))
  
  # 删除对未定义变量的引用
  # 不返回 rdata_file 等不存在的变量
  return(report_file)
}

# ==============================
# Main Function
# ==============================

reconstruct_ancestral_linkage_groups <- function(phase1_file, phase2_file, phase3_file, output_dir, custom_config = NULL) {
  # Main function for ancestral linkage group reconstruction
  message("Starting ancestral linkage group reconstruction...")
  
  # Apply custom configuration if provided
  if(!is.null(custom_config)) {
    # Update configuration with custom parameters
    for(param_name in names(custom_config)) {
      config[[param_name]] <- custom_config[[param_name]]
    }
    message("Applied custom configuration parameters")
  }
  
  # Load data from previous phases
  prev_data <- load_previous_phases(phase1_file, phase2_file, phase3_file)
  
  # Prepare node reconstruction contexts
  node_contexts <- prepare_node_reconstruction(
    prev_data$integrated_data$tree, 
    prev_data$event_results
  )
  
  # Create linkage markers
  linkage_markers <- create_linkage_markers(
    prev_data$conserved_groups, 
    prev_data$integrated_data
  )
  
  # Perform bottom-up reconstruction
  bottom_up_linkage <- reconstruct_bottom_up(
    node_contexts, 
    linkage_markers, 
    prev_data$integrated_data$tree, 
    prev_data$conserved_groups
  )
  
  # Perform top-down reconstruction
  top_down_linkage <- reconstruct_top_down(
    node_contexts, 
    linkage_markers, 
    prev_data$integrated_data$tree, 
    prev_data$conserved_groups
  )
  
  # Integrate both approaches
  integrated_linkage <- integrate_reconstruction_approaches(
    bottom_up_linkage, 
    top_down_linkage, 
    node_contexts
  )
  
  # Validate reconstruction
  validation_df <- validate_ancestral_linkage(
    integrated_linkage, 
    node_contexts, 
    prev_data$conserved_groups
  )
  
  # Save results
  saved_files <- save_ancestral_linkage(
    integrated_linkage, 
    validation_df, 
    prev_data$integrated_data, 
    output_dir
  )
  
  # Create visualizations
  tree_plot <- tryCatch({
    visualize_tree_with_linkage_groups(
      prev_data$integrated_data$tree, 
      integrated_linkage, 
      file.path(output_dir, "ancestral_linkage_tree.pdf")
    )
  }, error = function(e) {
    message(paste("Warning: Could not create tree visualization:", e$message))
    message("This is likely due to a mismatch between tree node labels and linkage data.")
    message("The analysis results are still valid, but the visualization could not be generated.")
    NULL
  })
  
  heatmap_matrix <- generate_linkage_heatmap(
    integrated_linkage, 
    prev_data$integrated_data, 
    prev_data$conserved_groups, 
    file.path(output_dir, "ancestral_linkage_heatmap.pdf")
  )
  
  # Generate node-specific visualizations
  for(node_id in names(integrated_linkage)) {
    visualize_node_linkage_groups(
      node_id, 
      integrated_linkage, 
      output_dir
    )
  }
  
  return(list(
    ancestral_linkage = integrated_linkage,
    validation = validation_df,
    bottom_up = bottom_up_linkage,
    top_down = top_down_linkage,
    saved_files = saved_files,
    tree_plot = tree_plot,
    heatmap_matrix = heatmap_matrix
  ))
}

# ==============================
# Command-line Interface
# ==============================

if(!interactive()) {
  # Parse command-line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  if(length(args) < 3) {
    cat("Usage: Rscript ancestral_chromosomes_phase4_linkage.R <integrated_data_file> <conserved_groups_file> <chromosome_events_file> [output_dir] [config_file]\n")
    cat("Example: Rscript ancestral_chromosomes_phase4_linkage.R results/integrated_data.RData results/phase2/conserved_chromosome_groups.RData results/phase3/chromosome_events.RData results/phase4\n")
    quit(status = 1)
  }
  
  # Parse arguments
  phase1_file <- args[1]
  phase2_file <- args[2]
  phase3_file <- args[3]
  
  # Set default output directory if not provided
  output_dir <- ifelse(length(args) >= 4, args[4], "results/phase4")
  
  # Load custom config if provided
  custom_config <- NULL
  if(length(args) >= 5) {
    config_file <- args[5]
    if(file.exists(config_file)) {
      custom_config <- tryCatch({
        readRDS(config_file)
      }, error = function(e) {
        message(paste("Warning: Could not load config file:", e$message))
        NULL
      })
    }
  }
  
  # Run reconstruction
  result <- reconstruct_ancestral_linkage_groups(
    phase1_file, 
    phase2_file, 
    phase3_file, 
    output_dir,
    custom_config
  )
}
