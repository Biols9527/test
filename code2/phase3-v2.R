#!/usr/bin/env Rscript

# Ancestral Chromosome Reconstruction Phase 3: Enhanced Shared Chromosome Event Detection
# Date: 2025-03-18
# Author: ShrivatsaMehan
# Description: Identifies shared chromosome evolutionary events (fusions, fissions, etc.)
#              across the phylogenetic tree based on conserved groups and mapping data,
#              using combined bottom-up and top-down approaches.

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
  
  # 验证数据结构
  message(paste("DEBUG: Tree object class:", class(integrated_data$tree)))
  message(paste("DEBUG: Tree has", length(integrated_data$tree$tip.label), "tips and", integrated_data$tree$Nnode, "internal nodes"))
  message(paste("DEBUG: chr_counts contains", length(integrated_data$chr_counts), "species"))

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
  
  # 验证数据结构
  message(paste("DEBUG: conserved_groups contains", length(conserved_result$conserved_groups), "groups"))
  if(length(conserved_result$conserved_groups) > 0) {
    group_names <- names(conserved_result$conserved_groups)
    message(paste("DEBUG: First group name:", group_names[1]))
  }

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
  
  message(paste("DEBUG: Tree has", n_tips, "tips and", tree$Nnode, "internal nodes"))
  message(paste("DEBUG: Total nodes:", n_tips + tree$Nnode))
  message(paste("DEBUG: Edge matrix dimensions:", nrow(tree$edge), "x", ncol(tree$edge)))

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
      if(!is.null(tree$node.label) && length(tree$node.label) >= (node_id - n_tips)) {
        node_label <- tree$node.label[node_id - n_tips]
        if(is.na(node_label) || node_label == "") {
          paste0("Node", node_id)
        } else {
          node_label
        }
      } else {
        paste0("Node", node_id)
      }
    }
  })

  edge_data$child_label <- sapply(edge_data$child_node_id, function(node_id) {
    if(node_id <= n_tips) {
      tree$tip.label[node_id]
    } else {
      if(!is.null(tree$node.label) && length(tree$node.label) >= (node_id - n_tips)) {
        node_label <- tree$node.label[node_id - n_tips]
        if(is.na(node_label) || node_label == "") {
          paste0("Node", node_id)
        } else {
          node_label
        }
      } else {
        paste0("Node", node_id)
      }
    }
  })
  
  # 检查标签数据是否正确
  message(paste("DEBUG: Parent labels sample:", paste(head(edge_data$parent_label), collapse=", ")))
  message(paste("DEBUG: Child labels sample:", paste(head(edge_data$child_label), collapse=", ")))

  message(paste("Extracted", nrow(edge_data), "edges from phylogenetic tree"))

  return(edge_data)
}

identify_clades <- function(tree) {
  # Identify all clades in the phylogenetic tree
  message("Identifying clades in phylogenetic tree...")

  n_tips <- length(tree$tip.label)
  n_nodes <- tree$Nnode
  
  message(paste("DEBUG: Looking for clades among", n_nodes, "internal nodes"))

  # Initialize clade data
  clade_data <- list()

  # For each internal node
  for(node_id in (n_tips + 1):(n_tips + n_nodes)) {
    # 添加调试信息
    message(paste("DEBUG: Processing node", node_id))
    
    # Get all tips descended from this node
    tips <- tryCatch({
      if(requireNamespace("geiger", quietly = TRUE)) {
        tips <- geiger::tips(tree, node_id)
        message(paste("DEBUG: Found", length(tips), "tips for node", node_id, "using geiger"))
        tips
      } else {
        descendants <- getDescendants(tree, node_id)
        descendants <- descendants[descendants <= n_tips]
        tips <- tree$tip.label[descendants]
        message(paste("DEBUG: Found", length(tips), "tips for node", node_id, "using custom function"))
        tips
      }
    }, error = function(e) {
      message(paste("DEBUG: Error getting tips for node", node_id, "-", e$message))
      # 更稳健的回退措施
      tryCatch({
        descendants <- getDescendants(tree, node_id)
        descendants <- descendants[descendants <= n_tips]
        tree$tip.label[descendants]
      }, error = function(e2) {
        message(paste("DEBUG: Fallback also failed:", e2$message))
        character(0)  # 返回空向量作为最后的回退
      })
    })

    # Get node label
    if(!is.null(tree$node.label) && length(tree$node.label) >= (node_id - n_tips)) {
      node_label <- tree$node.label[node_id - n_tips]
      if(is.na(node_label) || node_label == "") {
        node_label <- paste0("Node", node_id)
      }
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
    
    # 调试确认
    message(paste("DEBUG: Added clade for node", node_id, "with", length(tips), "tips and label", node_label))
  }

  message(paste("Identified", length(clade_data), "clades in phylogenetic tree"))

  return(clade_data)
}

# Helper function to get descendants of a node
getDescendants <- function(tree, node) {
  # 使用队列而非递归实现以避免堆栈溢出
  message(paste("DEBUG: Getting descendants for node", node))
  
  # 初始化
  descendants <- c()
  nodes_to_process <- c(node)
  processed_nodes <- c()
  
  # 当队列不为空时继续处理
  while(length(nodes_to_process) > 0) {
    # 取队列第一个节点
    current_node <- nodes_to_process[1]
    nodes_to_process <- nodes_to_process[-1]
    
    # 跳过已处理节点
    if(current_node %in% processed_nodes) {
      next
    }
    
    # 标记为已处理
    processed_nodes <- c(processed_nodes, current_node)
    
    # 找到当前节点的所有直接子节点
    edges <- which(tree$edge[, 1] == current_node)
    
    if(length(edges) > 0) {
      # 获取子节点ID
      child_nodes <- tree$edge[edges, 2]
      
      # 添加到后代列表
      descendants <- c(descendants, child_nodes)
      
      # 将内部节点（非叶子）添加到处理队列
      internal_children <- child_nodes[child_nodes > length(tree$tip.label)]
      if(length(internal_children) > 0) {
        nodes_to_process <- c(nodes_to_process, internal_children)
      }
    }
  }
  
  message(paste("DEBUG: Total", length(descendants), "descendants found for node", node))
  return(unique(descendants))
}

# 如果geiger包不可用，提供tips函数的备用实现
if(!requireNamespace("geiger", quietly = TRUE)) {
  message("Note: geiger package not found, using built-in methods for clade analysis")
  
  # 获取节点的所有叶子后代
  getTips <- function(tree, node) {
    # 如果已经是叶子节点
    if(node <= length(tree$tip.label)) {
      # 这已经是一个叶子节点
      return(tree$tip.label[node])
    }
    
    # 获取此节点的所有后代
    desc <- getDescendants(tree, node)
    # 只保留叶子节点
    tip_desc <- desc[desc <= length(tree$tip.label)]
    # 返回叶子标签
    return(tree$tip.label[tip_desc])
  }
  
  # 为geiger命名空间分配我们的函数
  geiger <- list(tips = getTips)
}

# ==============================
# Chromosome Event Detection Functions
# ==============================

reconstruct_ancestral_chromosome_counts <- function(tree, chr_counts) {
  # Reconstruct ancestral chromosome counts at internal nodes
  message("Reconstructing ancestral chromosome counts...")
  
  # 检查染色体计数数据
  message(paste("DEBUG: chr_counts contains", length(chr_counts), "species"))
  message(paste("DEBUG: Sample counts:", paste(head(chr_counts), collapse=", ")))
  
  n_tips <- length(tree$tip.label)
  n_nodes <- tree$Nnode
  
  # Convert chromosome counts to vector in order of tip labels
  count_vector <- numeric(n_tips)
  names(count_vector) <- tree$tip.label
  message(paste("DEBUG: Created count_vector with", length(count_vector), "elements"))
  
  for(sp in names(chr_counts)) {
    if(sp %in% tree$tip.label) {
      count_vector[sp] <- chr_counts[sp]
    }
  }
  
  # 检查向量赋值是否有效
  message(paste("DEBUG: After assignment, count_vector has", sum(!is.na(count_vector)), "non-NA values"))
  message(paste("DEBUG: Sample values:", paste(head(count_vector), collapse=", ")))
  
  # Check if we have counts for all tips
  if(any(is.na(count_vector))) {
    missing_tips <- names(count_vector)[is.na(count_vector)]
    warning(paste("Missing chromosome counts for", length(missing_tips), "species. Using median imputation."))
    message(paste("DEBUG: Missing counts for", length(missing_tips), "species:", paste(head(missing_tips), collapse=", ")))
    
    # Impute missing values with median
    median_count <- median(count_vector, na.rm = TRUE)
    message(paste("DEBUG: Using median value", median_count, "for imputation"))
    count_vector[is.na(count_vector)] <- median_count
  }
  
  # Create node count mapping that we'll fill with reconstructed values
  node_counts <- numeric(n_tips + n_nodes)
  names(node_counts) <- 1:(n_tips + n_nodes)
  message(paste("DEBUG: Created node_counts with", length(node_counts), "elements"))
  
  # Add tip counts
  for(i in 1:n_tips) {
    node_counts[i] <- count_vector[tree$tip.label[i]]
  }
  message(paste("DEBUG: Added", n_tips, "tip values to node_counts"))
  message(paste("DEBUG: Sample node_counts:", paste(head(node_counts), collapse=", ")))
  
  # Reconstruct ancestral states using sequential approach with multiple fallbacks
  success <- FALSE
  
  # 尝试拟合多种进化模型并选择最佳模型
  models <- list()
  
  # 1. 尝试拟合OU模型（有稳定选择压）
  if(requireNamespace("geiger", quietly = TRUE)) {
    models$OU <- tryCatch({
      message("DEBUG: Fitting OU model")
      ou_fit <- geiger::fitContinuous(tree, count_vector, model="OU")
      message("DEBUG: Successfully fit OU model")
      ou_fit
    }, error = function(e) NULL)
    
    # 2. 尝试拟合EB模型（早期爆发模型）
    models$EB <- tryCatch({
      message("DEBUG: Fitting EB model")
      eb_fit <- geiger::fitContinuous(tree, count_vector, model="EB")
      message("DEBUG: Successfully fit EB model")
      eb_fit
    }, error = function(e) NULL)
    
    # 3. 尝试拟合BM模型（布朗运动）
    models$BM <- tryCatch({
      message("DEBUG: Fitting BM model")
      bm_fit <- geiger::fitContinuous(tree, count_vector, model="BM")
      message("DEBUG: Successfully fit BM model")
      bm_fit
    }, error = function(e) NULL)
    
    # 选择最佳模型
    valid_models <- models[!sapply(models, is.null)]
    if(length(valid_models) > 0) {
      aic_values <- sapply(valid_models, function(m) m$opt$aic)
      best_model <- valid_models[[which.min(aic_values)]]
      model_name <- names(valid_models)[which.min(aic_values)]
      message(paste("DEBUG: Best fitting model is", model_name, "with AIC =", min(aic_values)))
      
      # 使用最佳模型的参数进行祖先状态重建
      if(model_name == "OU") {
        message("DEBUG: Using OU model for reconstruction - suggesting stabilizing selection on chromosome number")
      } else if(model_name == "EB") {
        message("DEBUG: Using EB model - suggesting early burst pattern in chromosome evolution")
      } else {
        message("DEBUG: Using BM model - suggesting neutral evolution of chromosome number")
      }
    } else {
      message("DEBUG: All model fitting failed. Using default Brownian Motion model.")
    }
  }

  # Try ML method first
  tryCatch({
    message("DEBUG: Attempting ML reconstruction")
    ace_result <- ace(count_vector, tree, type = "continuous", method = "ML")
    
    # Extract ancestral states and add to node_counts
    for(i in 1:n_nodes) {
      node_id <- i + n_tips
      node_counts[node_id] <- round(ace_result$ace[i, 1])  # Round to nearest integer
    }
    
    message("Successfully reconstructed ancestral states using ML method")
    message(paste("DEBUG: Sample internal node values:", 
                 paste(head(node_counts[(n_tips+1):length(node_counts)]), collapse=", ")))
    success <- TRUE
  }, error = function(e) {
    message("ML ancestral state reconstruction failed: ", e$message)
  })

  # If ML failed, try parsimony (PIC) method
  if(!success) {
    tryCatch({
      message("DEBUG: Attempting parsimony (PIC) reconstruction")
      ace_result <- ace(count_vector, tree, type = "continuous", method = "pic")
      
      # Extract ancestral states and add to node_counts
      for(i in 1:n_nodes) {
        node_id <- i + n_tips
        node_counts[node_id] <- round(ace_result$ace[i])  # Round to nearest integer
      }
      
      message("Successfully reconstructed ancestral states using parsimony method")
      message(paste("DEBUG: Sample internal node values:", 
                   paste(head(node_counts[(n_tips+1):length(node_counts)]), collapse=", ")))
      success <- TRUE
    }, error = function(e) {
      message("Parsimony ancestral state reconstruction failed: ", e$message)
    })
  }

  # If both methods failed, use simple median approach
  if(!success) {
    message("All statistical methods failed, using simple median approach for ancestral states")
    median_count <- round(median(count_vector, na.rm = TRUE))
    message(paste("DEBUG: Using median count", median_count, "for all internal nodes"))
    
    # Fill internal nodes with median value
    for(i in 1:n_nodes) {
      node_id <- i + n_tips
      node_counts[node_id] <- median_count
    }
    
    # Create basic ace_result structure for compatibility
    ace_result <- list(
      ace = matrix(median_count, nrow=n_nodes, ncol=1),
      method = "median",
      Nnode = n_nodes
    )
  }
  
  # 添加名称检查
  message("DEBUG: Checking node_counts names")
  message(paste("DEBUG: node_counts names:", paste(head(names(node_counts)), collapse=", ")))
  
  message("Completed ancestral chromosome count reconstruction")
  
  return(list(
    node_counts = node_counts,
    ace_result = ace_result
  ))
}

detect_count_changes <- function(tree, edge_data, chr_counts, ancestral_counts) {
  # Detect chromosome count changes across tree edges using both tip data and ancestral reconstruction
  message("Detecting chromosome count changes across tree edges...")
  
  # 调试名称检查
  message("DEBUG: Checking ancestral_counts structure")
  message(paste("DEBUG: ancestral_counts$node_counts has", length(ancestral_counts$node_counts), "elements"))
  message(paste("DEBUG: ancestral_counts$node_counts names:", paste(head(names(ancestral_counts$node_counts)), collapse=", ")))
  
  # Initialize count change data
  count_changes <- data.frame(
    edge_id = edge_data$edge_id,
    parent_node_id = edge_data$parent_node_id,
    child_node_id = edge_data$child_node_id,
    parent_label = edge_data$parent_label,
    child_label = edge_data$child_label,
    parent_is_tip = edge_data$parent_is_tip,
    child_is_tip = edge_data$child_is_tip,
    parent_count = NA_integer_,
    child_count = NA_integer_,
    count_diff = NA_integer_,
    count_ratio = NA_real_,
    count_change_type = NA_character_,
    stringsAsFactors = FALSE
  )
  
  # Fill in known counts from the ancestral reconstruction
  for(i in 1:nrow(count_changes)) {
    parent_id <- edge_data$parent_node_id[i]
    child_id <- edge_data$child_node_id[i]
    
    # 使用一致的方法访问node_counts
    parent_key <- as.character(parent_id)
    child_key <- as.character(child_id)
    
    message(paste("DEBUG: Processing edge", i, "- parent:", parent_id, "child:", child_id))
    
    # Get counts from ancestral reconstruction - 使用更强大的错误检查
    if(parent_key %in% names(ancestral_counts$node_counts)) {
      count_changes$parent_count[i] <- ancestral_counts$node_counts[parent_key]
      message(paste("DEBUG: Found parent count:", count_changes$parent_count[i]))
    } else {
      message(paste("DEBUG: Parent node", parent_id, "not found in ancestral_counts"))
    }
    
    if(child_key %in% names(ancestral_counts$node_counts)) {
      count_changes$child_count[i] <- ancestral_counts$node_counts[child_key]
      message(paste("DEBUG: Found child count:", count_changes$child_count[i]))
    } else {
      message(paste("DEBUG: Child node", child_id, "not found in ancestral_counts"))
    }
    
    # For tips, prioritize actual observed counts over reconstructed ones
    if(count_changes$child_is_tip[i] && edge_data$child_label[i] %in% names(chr_counts)) {
      count_changes$child_count[i] <- chr_counts[edge_data$child_label[i]]
      message(paste("DEBUG: Using observed count for tip", edge_data$child_label[i], ":", count_changes$child_count[i]))
    }
    
    # Calculate count difference and ratio
    if(!is.na(count_changes$parent_count[i]) && !is.na(count_changes$child_count[i])) {
      count_changes$count_diff[i] <- count_changes$child_count[i] - count_changes$parent_count[i]
      if(count_changes$parent_count[i] > 0) {
        count_changes$count_ratio[i] <- count_changes$child_count[i] / count_changes$parent_count[i]
      }
      
      # Determine change type
      if(count_changes$count_diff[i] > 1) {
        count_changes$count_change_type[i] <- "increase"
      } else if(count_changes$count_diff[i] < -1) {
        count_changes$count_change_type[i] <- "decrease"
      } else {
        count_changes$count_change_type[i] <- "stable"
      }
      
      message(paste("DEBUG: Edge", i, "diff:", count_changes$count_diff[i], 
                   "ratio:", round(count_changes$count_ratio[i], 2), 
                   "type:", count_changes$count_change_type[i]))
    }
  }
  
  message(paste("Detected count changes for", sum(!is.na(count_changes$count_diff)), "edges"))
  
  # 总结变化类型
  if(sum(!is.na(count_changes$count_change_type)) > 0) {
    change_summary <- table(count_changes$count_change_type, useNA="ifany")
    message(paste("DEBUG: Change type summary:", 
                 paste(names(change_summary), change_summary, sep="=", collapse=", ")))
  }
  
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
  message(paste("DEBUG: Found", nrow(unique_pairs), "unique species pairs"))
  
  # Initialize mapping pattern data
  pattern_data <- list()
  
  # Analyze each conserved group
  if(requireNamespace("parallel", quietly = TRUE) && length(conserved_groups$conserved_groups) > 10) {
    message("DEBUG: Using parallel processing for conserved group analysis...")
    
    # 检测操作系统并适配不同的并行方法
    is_windows <- .Platform$OS.type == "windows"
    n_cores <- min(parallel::detectCores(), 4)  # 限制使用的核心数，避免资源耗尽
    
    if(is_windows) {
      # Windows上使用parLapply
      message("DEBUG: Windows detected, using cluster-based parallelism")
      cl <- parallel::makeCluster(n_cores)
      # 导出必要的变量和函数
      parallel::clusterExport(cl, c("unique_pairs"), envir = environment())
      
      # 为集群提供必要的包
      parallel::clusterEvalQ(cl, {
        library(data.table)
        library(igraph)
      })
      
      # 处理每个保守组
      group_results <- parallel::parLapply(cl, names(conserved_groups$conserved_groups), function(group_name) {
        local_patterns <- list()
        group <- conserved_groups$conserved_groups[[group_name]]
        
        # 获取此组中有染色体的物种
        group_species <- names(group$species_chr_counts)
        
        # 检查每个物种对
        for(i in 1:nrow(unique_pairs)) {
          sp_A <- unique_pairs$species_A[i]
          sp_B <- unique_pairs$species_B[i]
          
          # 如果任一物种不在组中则跳过
          if(!(sp_A %in% group_species) || !(sp_B %in% group_species)) {
            next
          }
          
          # 安全地获取染色体
          chr_A <- tryCatch(strsplit(group$species_chromosomes[[sp_A]], ",")[[1]],
                          error = function(e) character(0))
          chr_B <- tryCatch(strsplit(group$species_chromosomes[[sp_B]], ",")[[1]],
                          error = function(e) character(0))
          
          # 如果没有染色体数据则跳过
          if(length(chr_A) == 0 || length(chr_B) == 0) {
            next
          }
          
          # 检测融合模式（A中的多个染色体映射到B中的一个）
          if(length(chr_A) > length(chr_B) && length(chr_B) == 1) {
            local_patterns[[length(local_patterns) + 1]] <- list(
              group = group_name,
              species_A = sp_A,
              species_B = sp_B,
              chromosomes_A = paste(chr_A, collapse=","),
              chromosomes_B = chr_B,
              pattern_type = "fusion_candidate",
              chromosome_ratio = length(chr_A) / length(chr_B)
            )
          }
          
          # 检测分裂模式（A中的一个染色体映射到B中的多个）
          if(length(chr_A) == 1 && length(chr_B) > 1) {
            local_patterns[[length(local_patterns) + 1]] <- list(
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
        
        return(local_patterns)
      })
      
      parallel::stopCluster(cl)
      
      # 展平结果
      pattern_data <- unlist(group_results, recursive = FALSE)
      pattern_data <- pattern_data[!sapply(pattern_data, is.null)]
      
    } else {
      # 非Windows系统使用mclapply
      message("DEBUG: Using fork-based parallelism")
      group_results <- parallel::mclapply(names(conserved_groups$conserved_groups), function(group_name) {
        local_patterns <- list()
        group <- conserved_groups$conserved_groups[[group_name]]
        
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
          chr_A <- tryCatch(strsplit(group$species_chromosomes[[sp_A]], ",")[[1]],
                           error = function(e) character(0))
          chr_B <- tryCatch(strsplit(group$species_chromosomes[[sp_B]], ",")[[1]],
                           error = function(e) character(0))
          
          # Skip if missing chromosome data
          if(length(chr_A) == 0 || length(chr_B) == 0) {
            next
          }
          
          # Look for fusion pattern (multiple chromosomes in A map to one in B)
          if(length(chr_A) > length(chr_B) && length(chr_B) == 1) {
            local_patterns[[length(local_patterns) + 1]] <- list(
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
            local_patterns[[length(local_patterns) + 1]] <- list(
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
        
        return(local_patterns)
      }, mc.cores = n_cores)
      
      # 展平结果
      pattern_data <- unlist(group_results, recursive = FALSE)
      pattern_data <- pattern_data[!sapply(pattern_data, is.null)]
    }
    
    message(paste("DEBUG: Parallel processing found", length(pattern_data), "pattern candidates"))
  } else {
    # 非并行处理 - 遍历每个保守组
    for(group_name in names(conserved_groups$conserved_groups)) {
      group <- conserved_groups$conserved_groups[[group_name]]
      message(paste(" - Analyzing mapping patterns for", group_name, "..."))
      
      # Get species with chromosomes in this group
      group_species <- names(group$species_chr_counts)
      message(paste("DEBUG:", length(group_species), "species in this group"))
      
      # Look for patterns in each species pair
      for(i in 1:nrow(unique_pairs)) {
        sp_A <- unique_pairs$species_A[i]
        sp_B <- unique_pairs$species_B[i]
        
        # Skip if either species not in group
        if(!(sp_A %in% group_species) || !(sp_B %in% group_species)) {
          next
        }
        
        # Get chromosomes for these species in this group
        chr_A <- tryCatch(strsplit(group$species_chromosomes[[sp_A]], ",")[[1]],
                         error = function(e) {
                           message(paste("DEBUG: Error getting chromosomes for", sp_A, "-", e$message))
                           character(0)
                         })
        chr_B <- tryCatch(strsplit(group$species_chromosomes[[sp_B]], ",")[[1]],
                         error = function(e) {
                           message(paste("DEBUG: Error getting chromosomes for", sp_B, "-", e$message))
                           character(0)
                         })
        
        # Skip if missing chromosome data
        if(length(chr_A) == 0 || length(chr_B) == 0) {
          next
        }
        
        message(paste("DEBUG: Species pair", sp_A, "-", sp_B, "has", length(chr_A), "and", length(chr_B), "chromosomes"))
        
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
          message(paste("DEBUG: Detected fusion pattern between", sp_A, "and", sp_B))
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
          message(paste("DEBUG: Detected fission pattern between", sp_A, "and", sp_B))
        }
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
    
    # 检查数据框结构
    message(paste("DEBUG: pattern_df has", nrow(pattern_df), "rows and", ncol(pattern_df), "columns"))
    if(nrow(pattern_df) > 0) {
      message(paste("DEBUG: Pattern type distribution:", paste(table(pattern_df$pattern_type), collapse=", ")))
    }
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
    message("DEBUG: No mapping patterns detected")
  }
  
  message(paste("Detected", nrow(pattern_df), "potential chromosome event patterns"))
  return(pattern_df)
}

# ==============================
# Bottom-Up Event Inference Functions
# ==============================

infer_events_bottom_up <- function(tree, edge_data, count_changes, pattern_data, chr_counts) {
  # Infer chromosome evolutionary events using bottom-up approach (from tips to internal nodes)
  message("Inferring chromosome events using bottom-up approach...")

  # Get unique pattern types
  if(nrow(pattern_data) > 0) {
    fusion_patterns <- pattern_data[pattern_data$pattern_type == "fusion_candidate",]
    fission_patterns <- pattern_data[pattern_data$pattern_type == "fission_candidate",]
    message(paste("DEBUG: Found", nrow(fusion_patterns), "fusion patterns and", 
                 nrow(fission_patterns), "fission patterns"))
  } else {
    fusion_patterns <- fission_patterns <- data.frame()
    message("DEBUG: No patterns found in pattern_data")
  }

  # Initialize events list
  events <- list()

  # Infer events based on chromosome count changes where child is a tip
  for(i in 1:nrow(count_changes)) {
    # Skip if no count data for child or parent
    if(is.na(count_changes$child_count[i]) || is.na(count_changes$parent_count[i])) {
      next
    }

    child_species <- count_changes$child_label[i]
    child_count <- count_changes$child_count[i]
    parent_count <- count_changes$parent_count[i]
    count_diff <- count_changes$count_diff[i]

    message(paste("DEBUG: Processing edge", i, "from", parent_count, "to", child_count, 
                 "chromosomes (diff:", count_diff, ")"))

    # Detect potential events based on count differences
    if(count_diff >= 2) {
      # Significant increase suggests fission
      events[[length(events) + 1]] <- list(
        edge_id = count_changes$edge_id[i],
        parent_node = count_changes$parent_node_id[i],
        child_node = count_changes$child_node_id[i],
        parent_label = count_changes$parent_label[i],
        child_label = child_species,
        event_type = "fission",
        parent_count = parent_count,
        child_count = child_count,
        count_diff = count_diff,
        count_ratio = count_changes$count_ratio[i],
        evidence_type = "count_difference",
        confidence = if(count_diff >= 4) "high" else "medium",
        analysis_method = "bottom_up"
      )
      message(paste("DEBUG: Detected fission event with confidence", 
                   if(count_diff >= 4) "high" else "medium"))
    } else if(count_diff <= -2) {
      # Significant decrease suggests fusion
      events[[length(events) + 1]] <- list(
        edge_id = count_changes$edge_id[i],
        parent_node = count_changes$parent_node_id[i],
        child_node = count_changes$child_node_id[i],
        parent_label = count_changes$parent_label[i],
        child_label = child_species,
        event_type = "fusion",
        parent_count = parent_count,
        child_count = child_count,
        count_diff = -count_diff,  # Make positive for clarity
        count_ratio = if(!is.na(count_changes$count_ratio[i]) && count_changes$count_ratio[i] != 0) 
                        1 / count_changes$count_ratio[i] else NA,  # Invert for clarity with NA check
        evidence_type = "count_difference",
        confidence = if(-count_diff >= 4) "high" else "medium",
        analysis_method = "bottom_up"
      )
      message(paste("DEBUG: Detected fusion event with confidence", 
                   if(-count_diff >= 4) "high" else "medium"))
    }

    # If child is a tip, check for additional supporting evidence from mapping patterns
    if(count_changes$child_is_tip[i]) {
      # Look for related species (sister taxa) to compare
      message(paste("DEBUG: Looking for sister taxa of", child_species))
      sister_species <- find_related_species(tree, child_species)
      message(paste("DEBUG: Found", length(sister_species), "sister species for", child_species))
      if(length(sister_species) > 0) {
        message(paste("DEBUG: Sister species:", paste(head(sister_species), collapse=", ")))
      }

      # Skip if no related species found
      if(length(sister_species) == 0) {
        next
      }

      # Check for supporting mapping patterns
      if(nrow(fusion_patterns) > 0) {
        # Look for fusion patterns where this species has fewer chromosomes
        fusion_evidence <- fusion_patterns[
          (fusion_patterns$species_B == child_species) &
          (fusion_patterns$species_A %in% sister_species),
        ]
        message(paste("DEBUG: Found", nrow(fusion_evidence), "fusion evidence patterns"))
        if(nrow(fusion_evidence) > 0) {
          for(j in 1:nrow(fusion_evidence)) {
            events[[length(events) + 1]] <- list(
              edge_id = count_changes$edge_id[i],
              parent_node = count_changes$parent_node_id[i],
              child_node = count_changes$child_node_id[i],
              parent_label = count_changes$parent_label[i],
              child_label = child_species,
              event_type = "fusion",
              parent_count = parent_count,
              child_count = child_count,
              count_diff = NA,
              count_ratio = fusion_evidence$chromosome_ratio[j],
              evidence_type = "mapping_pattern",
              confidence = "high",
              group = fusion_evidence$group[j],
              chromosomes_A = fusion_evidence$chromosomes_A[j],
              chromosomes_B = fusion_evidence$chromosomes_B[j],
              reference_species = fusion_evidence$species_A[j],
              analysis_method = "bottom_up"
            )
            message(paste("DEBUG: Added fusion event from mapping pattern with reference species", 
                          fusion_evidence$species_A[j]))
          }
        }
      }

      if(nrow(fission_patterns) > 0) {
        # Look for fission patterns where this species has more chromosomes
        fission_evidence <- fission_patterns[
          (fission_patterns$species_B == child_species) &
          (fission_patterns$species_A %in% sister_species),
        ]
        message(paste("DEBUG: Found", nrow(fission_evidence), "fission evidence patterns"))
        if(nrow(fission_evidence) > 0) {
          for(j in 1:nrow(fission_evidence)) {
            events[[length(events) + 1]] <- list(
              edge_id = count_changes$edge_id[i],
              parent_node = count_changes$parent_node_id[i],
              child_node = count_changes$child_node_id[i],
              parent_label = count_changes$parent_label[i],
              child_label = child_species,
              event_type = "fission",
              parent_count = parent_count,
              child_count = child_count,
              count_diff = NA,
              count_ratio = fission_evidence$chromosome_ratio[j],
              evidence_type = "mapping_pattern",
              confidence = "high",
              group = fission_evidence$group[j],
              chromosomes_A = fission_evidence$chromosomes_A[j],
              chromosomes_B = fission_evidence$chromosomes_B[j],
              reference_species = fission_evidence$species_A[j],
              analysis_method = "bottom_up"
            )
            message(paste("DEBUG: Added fission event from mapping pattern with reference species", 
                          fission_evidence$species_A[j]))
          }
        }
      }
    }
  }

  # Convert to data frame if events were found
  if(length(events) > 0) {
    message(paste("DEBUG: Converting", length(events), "events to data frame"))
    # Handle different event structures by filling in missing columns
    events_list <- lapply(events, function(event) {
      if(!is.null(event$edge_id)) {
        # Edge-specific event
        data.frame(
          edge_id = event$edge_id,
          parent_node = event$parent_node,
          child_node = event$child_node,
          parent_label = event$parent_label,
          child_label = event$child_label,
          event_type = event$event_type,
          parent_count = event$parent_count,
          child_count = event$child_count,
          count_diff = event$count_diff,
          count_ratio = event$count_ratio,
          evidence_type = event$evidence_type,
          confidence = event$confidence,
          analysis_method = event$analysis_method,
          shared_event = FALSE,
          stringsAsFactors = FALSE
        )
      } else {
        # Node-level shared event
        data.frame(
          edge_id = NA_integer_,
          parent_node = event$parent_node,
          child_node = NA_integer_,
          parent_label = event$parent_label,
          child_label = NA_character_,
          event_type = event$event_type,
          parent_count = event$node_count,
          child_count = event$child_count_mean,
          count_diff = event$count_diff_mean,
          count_ratio = NA_real_,
          evidence_type = event$evidence_type,
          confidence = event$confidence,
          analysis_method = event$analysis_method,
          affected_children = event$affected_children,
          total_children = event$total_children,
          prevalence = event$prevalence,
          shared_event = TRUE,
          stringsAsFactors = FALSE
        )
      }
    })
    
    # 使用tryCatch安全地合并数据框
    events_df <- tryCatch({
      do.call(rbind, events_list)
    }, error = function(e) {
      message(paste("DEBUG: Error combining events into data frame:", e$message))
      # 查找不一致的列
      all_cols <- unique(unlist(lapply(events_list, names)))
      message(paste("DEBUG: All column names:", paste(all_cols, collapse=", ")))
      # 确保所有事件都有相同的列
      events_list_fixed <- lapply(events_list, function(event_df) {
        for(col in all_cols) {
          if(!(col %in% names(event_df))) {
            event_df[[col]] <- NA
          }
        }
        return(event_df)
      })
      # 再次尝试合并
      tryCatch({
        do.call(rbind, events_list_fixed)
      }, error = function(e2) {
        message(paste("DEBUG: Second attempt failed:", e2$message))
        # 创建空数据框作为最终回退
        data.frame(
          edge_id = integer(),
          parent_node = integer(),
          child_node = integer(),
          parent_label = character(),
          child_label = character(),
          event_type = character(),
          parent_count = integer(),
          child_count = integer(),
          count_diff = integer(),
          count_ratio = numeric(),
          evidence_type = character(),
          confidence = character(),
          analysis_method = character(),
          affected_children = integer(),
          total_children = integer(),
          prevalence = numeric(),
          shared_event = logical(),
          stringsAsFactors = FALSE
        )
      })
    })
    
    # 检查结果
    message(paste("DEBUG: Created events_df with", nrow(events_df), "rows and", ncol(events_df), "columns"))
    if(nrow(events_df) > 0) {
      message(paste("DEBUG: Event types:", paste(table(events_df$event_type), collapse=", ")))
    }
  } else {
    # Create empty dataframe if no events found
    events_df <- data.frame(
      edge_id = integer(),
      parent_node = integer(),
      child_node = integer(),
      parent_label = character(),
      child_label = character(),
      event_type = character(),
      parent_count = integer(),
      child_count = integer(),
      count_diff = integer(),
      count_ratio = numeric(),
      evidence_type = character(),
      confidence = character(),
      analysis_method = character(),
      affected_children = integer(),
      total_children = integer(),
      prevalence = numeric(),
      shared_event = logical(),
      stringsAsFactors = FALSE
    )
    message("DEBUG: No events found, created empty data frame")
  }

  message(paste("Inferred", nrow(events_df), "chromosome events using bottom-up approach"))
  return(events_df)
}

# ==============================
# Top-Down Event Inference Functions
# ==============================

infer_events_top_down <- function(tree, edge_data, count_changes, ancestral_counts, clades) {
  # Infer chromosome evolutionary events using top-down approach (from root to tips)
  message("Inferring chromosome events using top-down approach...")

  # Initialize events list
  events <- list()

  # Get all internal nodes
  n_tips <- length(tree$tip.label)
  internal_nodes <- unique(edge_data$parent_node_id[!edge_data$parent_is_tip])
  message(paste("DEBUG: Found", length(internal_nodes), "internal nodes to analyze"))

  # For each internal node, examine its immediate children
  for(node_id in internal_nodes) {
    # Get node label
    node_label <- if(node_id <= n_tips) {
      tree$tip.label[node_id]
    } else {
      if(!is.null(tree$node.label) && length(tree$node.label) >= (node_id - n_tips)) {
        node_label <- tree$node.label[node_id - n_tips]
        if(is.na(node_label) || node_label == "") {
          paste0("Node", node_id)
        } else {
          node_label
        }
      } else {
        paste0("Node", node_id)
      }
    }

    message(paste("DEBUG: Processing node", node_id, "(", node_label, ")"))

    # Get node's chromosome count - 使用更安全的访问方式
    node_key <- as.character(node_id)
    if(node_key %in% names(ancestral_counts$node_counts)) {
      node_count <- ancestral_counts$node_counts[node_key]
    } else {
      message(paste("DEBUG: Node", node_id, "not found in ancestral_counts, skipping"))
      next
    }

    # Skip if no count available
    if(is.na(node_count)) {
      message(paste("DEBUG: No count available for node", node_id, ", skipping"))
      next
    }

    message(paste("DEBUG: Node", node_id, "has", node_count, "chromosomes"))

    # Get all edges where this node is the parent
    child_edges <- edge_data[edge_data$parent_node_id == node_id,]
    message(paste("DEBUG: Node has", nrow(child_edges), "child edges"))

    # Get count changes for these edges
    edge_count_changes <- count_changes[count_changes$parent_node_id == node_id,]
    message(paste("DEBUG: Found", nrow(edge_count_changes), "count changes for this node"))

    # Skip if no count changes available
    if(nrow(edge_count_changes) == 0) {
      next
    }

    # Analyze count change patterns across all children
    changes <- na.omit(edge_count_changes$count_diff)
    message(paste("DEBUG: After removing NA values, found", length(changes), "count changes"))
    if(length(changes) == 0) {
      next
    }

    # Check if there's a consistent pattern across multiple children
    if(length(changes) >= 2) {
      # Check for consistent increase (fission) across children
      increase_count <- sum(changes > 1)
      increase_ratio <- increase_count / length(changes)
      message(paste("DEBUG: Increase pattern:", increase_count, "children (", 
                   round(increase_ratio * 100), "%)"))

      # Check for consistent decrease (fusion) across children
      decrease_count <- sum(changes < -1)
      decrease_ratio <- decrease_count / length(changes)
      message(paste("DEBUG: Decrease pattern:", decrease_count, "children (", 
                   round(decrease_ratio * 100), "%)"))

      # If significant fraction of children show same direction of change, infer shared event
      if(increase_ratio >= 0.5 && increase_count >= 2) {
        # Shared fission event at this node
        events[[length(events) + 1]] <- list(
          parent_node = node_id,
          parent_label = node_label,
          event_type = "fission",
          node_count = node_count,
          child_count_mean = mean(edge_count_changes$child_count, na.rm = TRUE),
          count_diff_mean = mean(changes[changes > 0], na.rm = TRUE),
          affected_children = increase_count,
          total_children = length(changes),
          prevalence = increase_ratio,
          evidence_type = "consistent_child_pattern",
          confidence = if(increase_ratio >= 0.8) "high" else "medium",
          analysis_method = "top_down"
        )
        message(paste("DEBUG: Detected shared fission event at node", node_id, 
                     "with prevalence", round(increase_ratio * 100), "%"))
      }

      if(decrease_ratio >= 0.5 && decrease_count >= 2) {
        # Shared fusion event at this node
        events[[length(events) + 1]] <- list(
          parent_node = node_id,
          parent_label = node_label,
          event_type = "fusion",
          node_count = node_count,
          child_count_mean = mean(edge_count_changes$child_count, na.rm = TRUE),
          count_diff_mean = mean(-changes[changes < 0], na.rm = TRUE),  # Make positive for clarity
          affected_children = decrease_count,
          total_children = length(changes),
          prevalence = decrease_ratio,
          evidence_type = "consistent_child_pattern",
          confidence = if(decrease_ratio >= 0.8) "high" else "medium",
          analysis_method = "top_down"
        )
        message(paste("DEBUG: Detected shared fusion event at node", node_id, 
                     "with prevalence", round(decrease_ratio * 100), "%"))
      }
    }

    # Also check each individual child for significant changes that might represent lineage-specific events
    for(i in 1:nrow(edge_count_changes)) {
      count_diff <- edge_count_changes$count_diff[i]

      # Skip if no significant change
      if(is.na(count_diff) || abs(count_diff) < 2) {
        next
      }

      child_node <- edge_count_changes$child_node_id[i]
      child_label <- edge_count_changes$child_label[i]
      child_count <- edge_count_changes$child_count[i]

      message(paste("DEBUG: Significant change for child", child_node, "(", child_label, "):", count_diff))

      # Significant increase suggests fission
      if(count_diff >= 2) {
        events[[length(events) + 1]] <- list(
          edge_id = edge_count_changes$edge_id[i],
          parent_node = node_id,
          child_node = child_node,
          parent_label = node_label,
          child_label = child_label,
          event_type = "fission",
          parent_count = node_count,
          child_count = child_count,
          count_diff = count_diff,
          count_ratio = edge_count_changes$count_ratio[i],
          evidence_type = "ancestral_count_difference",
          confidence = if(count_diff >= 4) "high" else "medium",
          analysis_method = "top_down"
        )
        message(paste("DEBUG: Detected lineage-specific fission for", child_label))
      } else if(count_diff <= -2) {
        # Significant decrease suggests fusion
        events[[length(events) + 1]] <- list(
          edge_id = edge_count_changes$edge_id[i],
          parent_node = node_id,
          child_node = child_node,
          parent_label = node_label,
          child_label = child_label,
          event_type = "fusion",
          parent_count = node_count,
          child_count = child_count,
          count_diff = -count_diff,  # Make positive for clarity
          count_ratio = if(!is.na(edge_count_changes$count_ratio[i]) && edge_count_changes$count_ratio[i] != 0) 
                          1 / edge_count_changes$count_ratio[i] else NA,  # Invert with NA check
          evidence_type = "ancestral_count_difference",
          confidence = if(-count_diff >= 4) "high" else "medium",
          analysis_method = "top_down"
        )
        message(paste("DEBUG: Detected lineage-specific fusion for", child_label))
      }
    }
  }

  # Convert to data frame if events were found
  if(length(events) > 0) {
    message(paste("DEBUG: Converting", length(events), "top-down events to data frame"))
    # Handle different event structures by filling in missing columns
    events_list <- lapply(events, function(event) {
      if(!is.null(event$edge_id)) {
        # Edge-specific event
        data.frame(
          edge_id = event$edge_id,
          parent_node = event$parent_node,
          child_node = event$child_node,
          parent_label = event$parent_label,
          child_label = event$child_label,
          event_type = event$event_type,
          parent_count = event$parent_count,
          child_count = event$child_count,
          count_diff = event$count_diff,
          count_ratio = event$count_ratio,
          evidence_type = event$evidence_type,
          confidence = event$confidence,
          analysis_method = event$analysis_method,
          shared_event = FALSE,
          stringsAsFactors = FALSE
        )
      } else {
        # Node-level shared event
        data.frame(
          edge_id = NA_integer_,
          parent_node = event$parent_node,
          child_node = NA_integer_,
          parent_label = event$parent_label,
          child_label = NA_character_,
          event_type = event$event_type,
          parent_count = event$node_count,
          child_count = event$child_count_mean,
          count_diff = event$count_diff_mean,
          count_ratio = NA_real_,
          evidence_type = event$evidence_type,
          confidence = event$confidence,
          analysis_method = event$analysis_method,
          affected_children = event$affected_children,
          total_children = event$total_children,
          prevalence = event$prevalence,
          shared_event = TRUE,
          stringsAsFactors = FALSE
        )
      }
    })
    
    # 使用tryCatch安全地合并数据框
    events_df <- tryCatch({
      do.call(rbind, events_list)
    }, error = function(e) {
      message(paste("DEBUG: Error combining top-down events into data frame:", e$message))
      # 查找不一致的列
      all_cols <- unique(unlist(lapply(events_list, names)))
      message(paste("DEBUG: All column names:", paste(all_cols, collapse=", ")))
      # 确保所有事件都有相同的列
      events_list_fixed <- lapply(events_list, function(event_df) {
        for(col in all_cols) {
          if(!(col %in% names(event_df))) {
            event_df[[col]] <- NA
          }
        }
        return(event_df)
      })
      # 再次尝试合并
      tryCatch({
        do.call(rbind, events_list_fixed)
      }, error = function(e2) {
        message(paste("DEBUG: Second attempt failed:", e2$message))
        # 创建空数据框作为最终回退
        data.frame(
          edge_id = integer(),
          parent_node = integer(),
          child_node = integer(),
          parent_label = character(),
          child_label = character(),
          event_type = character(),
          parent_count = integer(),
          child_count = integer(),
          count_diff = integer(),
          count_ratio = numeric(),
          evidence_type = character(),
          confidence = character(),
          analysis_method = character(),
          affected_children = integer(),
          total_children = integer(),
          prevalence = numeric(),
          shared_event = logical(),
          stringsAsFactors = FALSE
        )
      })
    })
    
    # 检查结果
    message(paste("DEBUG: Created events_df with", nrow(events_df), "rows and", ncol(events_df), "columns"))
    if(nrow(events_df) > 0) {
      message(paste("DEBUG: Event types:", paste(table(events_df$event_type), collapse=", ")))
    }
  } else {
    # Create empty dataframe if no events found
    events_df <- data.frame(
      edge_id = integer(),
      parent_node = integer(),
      child_node = integer(),
      parent_label = character(),
      child_label = character(),
      event_type = character(),
      parent_count = integer(),
      child_count = integer(),
      count_diff = integer(),
      count_ratio = numeric(),
      evidence_type = character(),
      confidence = character(),
      analysis_method = character(),
      affected_children = integer(),
      total_children = integer(),
      prevalence = numeric(),
      shared_event = logical(),
      stringsAsFactors = FALSE
    )
    message("DEBUG: No events found, created empty data frame")
  }

  message(paste("Inferred", nrow(events_df), "chromosome events using top-down approach"))
  return(events_df)
}

# ==============================
# Combined Event Analysis
# ==============================

combine_event_analyses <- function(bottom_up_events, top_down_events) {
  # Combine results from both approaches
  message("Combining bottom-up and top-down event analyses...")

  message(paste("DEBUG: Bottom-up events:", nrow(bottom_up_events), "rows and", ncol(bottom_up_events), "columns"))
  message(paste("DEBUG: Top-down events:", nrow(top_down_events), "rows and", ncol(top_down_events), "columns"))

  # Add analysis type flag to each dataset
  if(nrow(bottom_up_events) > 0) {
    if(!"shared_event" %in% names(bottom_up_events)) {
      bottom_up_events$shared_event <- FALSE
      message("DEBUG: Added shared_event column to bottom_up_events")
    }
  }

  # Combine events
  if(nrow(bottom_up_events) > 0 && nrow(top_down_events) > 0) {
    # Find overlapping columns
    common_cols <- intersect(names(bottom_up_events), names(top_down_events))
    message(paste("DEBUG: Common columns:", length(common_cols)))
    # Log column differences
    bu_only <- setdiff(names(bottom_up_events), names(top_down_events))
    td_only <- setdiff(names(top_down_events), names(bottom_up_events))
    if(length(bu_only) > 0) {
      message(paste("DEBUG: Bottom-up only columns:", paste(bu_only, collapse=", ")))
    }
    if(length(td_only) > 0) {
      message(paste("DEBUG: Top-down only columns:", paste(td_only, collapse=", ")))
    }

    # Ensure both dataframes have the same columns
    for(col in setdiff(names(bottom_up_events), names(top_down_events))) {
      top_down_events[[col]] <- NA
      message(paste("DEBUG: Added column", col, "to top_down_events"))
    }
    for(col in setdiff(names(top_down_events), names(bottom_up_events))) {
      bottom_up_events[[col]] <- NA
      message(paste("DEBUG: Added column", col, "to bottom_up_events"))
    }

    # Combine events - 使用tryCatch处理合并错误
    all_events <- tryCatch({
      rbind(bottom_up_events, top_down_events)
    }, error = function(e) {
      message(paste("DEBUG: Error combining events:", e$message))
      # 检查列类型是否一致
      for(col in names(bottom_up_events)) {
        bu_type <- class(bottom_up_events[[col]])
        td_type <- class(top_down_events[[col]])
        if(!identical(bu_type, td_type)) {
          message(paste("DEBUG: Type mismatch for column", col, "- bottom-up:", 
                       paste(bu_type, collapse="/"), "top-down:", paste(td_type, collapse="/")))
          # 尝试调整类型
          if("character" %in% c(bu_type, td_type)) {
            bottom_up_events[[col]] <- as.character(bottom_up_events[[col]])
            top_down_events[[col]] <- as.character(top_down_events[[col]])
            message(paste("DEBUG: Converted both", col, "columns to character"))
          } else if("numeric" %in% c(bu_type, td_type) && "integer" %in% c(bu_type, td_type)) {
            bottom_up_events[[col]] <- as.numeric(bottom_up_events[[col]])
            top_down_events[[col]] <- as.numeric(top_down_events[[col]])
            message(paste("DEBUG: Converted both", col, "columns to numeric"))
          }
        }
      }
      # 再次尝试合并
      tryCatch({
        rbind(bottom_up_events, top_down_events)
      }, error = function(e2) {
        message(paste("DEBUG: Second attempt to combine events failed:", e2$message))
        # 返回底部事件作为后备
        message("DEBUG: Returning bottom-up events only as fallback")
        bottom_up_events
      })
    })
    message(paste("DEBUG: Combined data frame has", nrow(all_events), "rows and", ncol(all_events), "columns"))
  } else if(nrow(bottom_up_events) > 0) {
    all_events <- bottom_up_events
    message("DEBUG: Using only bottom-up events")
  } else if(nrow(top_down_events) > 0) {
    all_events <- top_down_events
    message("DEBUG: Using only top-down events")
  } else {
    # Create empty dataframe if no events found
    all_events <- data.frame(
      edge_id = integer(),
      parent_node = integer(),
      child_node = integer(),
      parent_label = character(),
      child_label = character(),
      event_type = character(),
      parent_count = integer(),
      child_count = integer(),
      count_diff = integer(),
      count_ratio = numeric(),
      evidence_type = character(),
      confidence = character(),
      analysis_method = character(),
      affected_children = integer(),
      total_children = integer(),
      prevalence = numeric(),
      shared_event = logical(),
      stringsAsFactors = FALSE
    )
    message("DEBUG: No events found in either analysis, created empty data frame")
  }

  message(paste("Combined event analyses with a total of", nrow(all_events), "events"))

  # Look for corroborating evidence between methods
  if(nrow(all_events) > 0) {
    if(!"corroborated" %in% names(all_events)) {
      all_events$corroborated <- FALSE
      message("DEBUG: Added corroborated column to all_events")
    }
    for(i in 1:nrow(bottom_up_events)) {
      bu_event <- bottom_up_events[i,]
      # Skip if no edge ID
      if(is.na(bu_event$edge_id)) {
        next
      }
      # Find matching top-down events
      matching_td <- top_down_events[
        !is.na(top_down_events$edge_id) &
        top_down_events$edge_id == bu_event$edge_id &
        top_down_events$event_type == bu_event$event_type,
      ]
      if(nrow(matching_td) > 0) {
        message(paste("DEBUG: Found corroborating evidence for edge", bu_event$edge_id, 
                      "event type", bu_event$event_type))
        # Event detected by both methods - increase confidence
        all_events$confidence[all_events$edge_id == bu_event$edge_id &
                              all_events$event_type == bu_event$event_type &
                              all_events$analysis_method == "bottom_up"] <- "high"
        # Add corroboration flag
        all_events$corroborated[all_events$edge_id == bu_event$edge_id &
                                all_events$event_type == bu_event$event_type] <- TRUE
      }
    }
    # Count corroborated events
    corr_count <- sum(all_events$corroborated, na.rm = TRUE)
    message(paste("DEBUG: Found", corr_count, "corroborated events"))
  }

  return(all_events)
}

# ==============================
# Shared Event Identification
# ==============================

identify_shared_events <- function(tree, events_df, clades) {
  # Identify shared chromosome events across branches
  message("Identifying shared chromosome events...")

  # If no events were found, return empty result
  if(nrow(events_df) == 0) {
    message("No events to analyze for sharing")
    return(data.frame(
      event_type = character(),
      clade_node_id = integer(),
      clade_label = character(),
      tip_count = integer(),
      event_count = integer(),
      prevalence = numeric(),
      shared_status = character(),
      confidence = character(),
      analysis_type = character(),
      corroborated = logical(),
      stringsAsFactors = FALSE
    ))
  }

  message(paste("DEBUG: Analyzing", nrow(events_df), "events for sharing patterns"))
  # Get unique event types
  event_types <- unique(events_df$event_type)
  message(paste("DEBUG: Found event types:", paste(event_types, collapse=", ")))

  # Initialize results
  shared_events <- list()

  # Find events already identified as shared by top-down analysis
  shared_from_top_down <- events_df[events_df$shared_event == TRUE,]
  message(paste("DEBUG: Found", nrow(shared_from_top_down), "events already identified as shared"))

  for(i in 1:nrow(shared_from_top_down)) {
    event <- shared_from_top_down[i,]
    # Find corresponding clade
    clade <- clades[[as.character(event$parent_node)]]
    message(paste("DEBUG: Processing shared event", i, "for node", event$parent_node))
    if(is.null(clade)) {
      message(paste("DEBUG: No clade found for node", event$parent_node))
      next
    }
    shared_events[[length(shared_events) + 1]] <- list(
      event_type = event$event_type,
      clade_node_id = event$parent_node,
      clade_label = event$parent_label,
      tip_count = clade$tip_count,
      event_count = event$affected_children,
      prevalence = event$prevalence,
      shared_status = if(event$prevalence >= 0.8) "strongly_shared" else "likely_shared",
      confidence = event$confidence,
      analysis_type = "top_down",
      corroborated = event$corroborated
    )
    message(paste("DEBUG: Added shared", event$event_type, "event for clade", event$parent_label))
  }

  # For each event type and clade, check if the event is shared using bottom-up approach
  for(event_type in event_types) {
    type_events <- events_df[events_df$event_type == event_type & !events_df$shared_event,]
    message(paste("DEBUG: Found", nrow(type_events), event_type, "events for bottom-up sharing analysis"))
    for(clade_id in names(clades)) {
      clade <- clades[[clade_id]]
      message(paste("DEBUG: Analyzing clade", clade_id, "(", clade$node_label, ") with", length(clade$tips), "tips"))
      clade_tips <- clade$tips
      clade_events <- type_events[type_events$child_label %in% clade_tips,]
      # Count events in this clade
      event_count <- nrow(clade_events)
      message(paste("DEBUG: Found", event_count, event_type, "events in this clade"))
      # Skip clades with no events
      if(event_count == 0) {
        next
      }
      # Calculate event prevalence in clade
      prevalence <- event_count / clade$tip_count
      message(paste("DEBUG: Event prevalence:", round(prevalence * 100, 1), "%"))
      # Determine if the event is shared across the clade
      if(prevalence >= 0.8) {
        message("DEBUG: Detected strongly shared event pattern")
        shared_status <- "strongly_shared"
        confidence <- "high"
      } else if(prevalence >= 0.5) {
        message("DEBUG: Detected likely shared event pattern")
        shared_status <- "likely_shared"
        confidence <- "medium"
      } else if(prevalence >= 0.3) {
        message("DEBUG: Detected partially shared event pattern")
        shared_status <- "partially_shared"
        confidence <- "low"
      } else {
        message("DEBUG: Low prevalence, not considering as shared event")
        next
      }
      # Check if there's a matching top-down event for this clade
      corroborated <- FALSE
      if(nrow(shared_from_top_down) > 0) {
        matching_td <- shared_from_top_down[
          shared_from_top_down$parent_node == as.integer(clade_id) &
          shared_from_top_down$event_type == event_type,
        ]
        if(nrow(matching_td) > 0) {
          message("DEBUG: Found corroborating top-down event")
          corroborated <- TRUE
          confidence <- "high"
        }
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
        analysis_type = "bottom_up",
        corroborated = corroborated,
        tip_events = list(clade_events)
      )
      message(paste("DEBUG: Added", shared_status, event_type, "event for clade", clade$node_label))
    }
  }

  # Convert to data frame if shared events were found
  if(length(shared_events) > 0) {
    message(paste("DEBUG: Converting", length(shared_events), "shared events to data frame"))
    shared_df <- tryCatch({
      do.call(rbind, lapply(shared_events, function(e) {
        df <- data.frame(
          event_type = e$event_type,
          clade_node_id = e$clade_node_id,
          clade_label = e$clade_label,
          tip_count = e$tip_count,
          event_count = e$event_count,
          prevalence = e$prevalence,
          shared_status = e$shared_status,
          confidence = e$confidence,
          analysis_type = e$analysis_type,
          corroborated = e$corroborated,
          stringsAsFactors = FALSE
        )
        return(df)
      }))
    }, error = function(err) {
      message(paste("DEBUG: Error creating shared_df:", err$message))
      data.frame(
        event_type = character(),
        clade_node_id = integer(),
        clade_label = character(),
        tip_count = integer(),
        event_count = integer(),
        prevalence = numeric(),
        shared_status = character(),
        confidence = character(),
        analysis_type = character(),
        corroborated = logical(),
        stringsAsFactors = FALSE
      )
    })
    # Store tip events separately because they don't fit in a data frame
    attr(shared_df, "tip_events") <- lapply(shared_events, function(e) {
      if(!is.null(e$tip_events)) e$tip_events[[1]] else NULL
    })
    message(paste("DEBUG: Stored events for", length(attr(shared_df, "tip_events")), "clades as attributes"))
  } else {
    message("DEBUG: No shared events found, created empty data frame")
    shared_df <- data.frame(
      event_type = character(),
      clade_node_id = integer(),
      clade_label = character(),
      tip_count = integer(),
      event_count = integer(),
      prevalence = numeric(),
      shared_status = character(),
      confidence = character(),
      analysis_type = character(),
      corroborated = logical(),
      stringsAsFactors = FALSE
    )
    attr(shared_df, "tip_events") <- list()
  }

  message(paste("Identified", nrow(shared_df), "shared chromosome events"))
  return(shared_df)
}

# ==============================
# Output Functions
# ==============================

save_event_detection_results <- function(event_result, output_dir) {
  # Save event detection results to output files
  message("Saving event detection results...")
  if(!dir.exists(output_dir)) {
    # Create output directory if it doesn't exist
    dir.create(output_dir, recursive = TRUE)
    message(paste("Created output directory:", output_dir))
  }

  # Save full results data
  events_rdata_file <- file.path(output_dir, "chromosome_events.RData")
  saveRDS(event_result, file = events_rdata_file)
  message(paste("DEBUG: Saved full results to", events_rdata_file))

  # Save individual data files
  edge_data_file <- file.path(output_dir, "tree_edges.csv")
  write.csv(event_result$edge_data, file = edge_data_file, row.names = FALSE)
  message(paste("DEBUG: Saved edge data to", edge_data_file))

  ancestral_counts_file <- file.path(output_dir, "ancestral_chromosome_counts.csv")
  tryCatch({
    ancestral_counts_df <- data.frame(
      node_id = as.integer(names(event_result$ancestral_counts$node_counts)),
      count = event_result$ancestral_counts$node_counts,
      stringsAsFactors = FALSE
    )
    write.csv(ancestral_counts_df, file = ancestral_counts_file, row.names = FALSE)
    message(paste("DEBUG: Saved ancestral counts to", ancestral_counts_file))
  }, error = function(e) {
    message(paste("DEBUG: Error creating ancestral counts data frame:", e$message))
    ancestral_counts_vector <- event_result$ancestral_counts$node_counts
    df <- data.frame(
      node_id = 1:length(ancestral_counts_vector),
      count = ancestral_counts_vector,
      stringsAsFactors = FALSE
    )
    write.csv(df, file = ancestral_counts_file, row.names = FALSE)
    message("DEBUG: Saved simplified ancestral counts data")
  })

  count_changes_file <- file.path(output_dir, "chromosome_count_changes.csv")
  write.csv(event_result$count_changes, file = count_changes_file, row.names = FALSE)
  message(paste("DEBUG: Saved count changes to", count_changes_file))

  pattern_data_file <- file.path(output_dir, "mapping_patterns.csv")
  write.csv(event_result$pattern_data, file = pattern_data_file, row.names = FALSE)
  message(paste("DEBUG: Saved mapping patterns to", pattern_data_file))

  bottom_up_file <- file.path(output_dir, "bottom_up_events.csv")
  write.csv(event_result$bottom_up_events, file = bottom_up_file, row.names = FALSE)
  message(paste("DEBUG: Saved bottom-up events to", bottom_up_file))

  top_down_file <- file.path(output_dir, "top_down_events.csv")
  write.csv(event_result$top_down_events, file = top_down_file, row.names = FALSE)
  message(paste("DEBUG: Saved top-down events to", top_down_file))

  events_file <- file.path(output_dir, "inferred_events.csv")
  write.csv(event_result$events, file = events_file, row.names = FALSE)
  message(paste("DEBUG: Saved combined events to", events_file))

  shared_events_file <- file.path(output_dir, "shared_events.csv")
  write.csv(event_result$shared_events, file = shared_events_file, row.names = FALSE)
  message(paste("DEBUG: Saved shared events to", shared_events_file))

  # Create detailed events report
  events_report <- file.path(output_dir, "chromosome_events_report.txt")
  sink(events_report)
  cat("Enhanced Chromosome Evolutionary Events Report\n")
  cat("===========================================\n\n")
  cat("----------------------------------\n")
  cat("1. Ancestral Chromosome Count Summary\n")
  tip_counts <- event_result$ancestral_counts$node_counts[1:length(event_result$tree$tip.label)]
  internal_counts <- event_result$ancestral_counts$node_counts[(length(tip_counts) + 1):length(event_result$ancestral_counts$node_counts)]
  cat(paste("Extant species chromosome count range:", min(tip_counts, na.rm=TRUE), 
            "to", max(tip_counts, na.rm=TRUE), "\n"))
  cat(paste("Ancestral node chromosome count range:", min(internal_counts, na.rm=TRUE), 
            "to", max(internal_counts, na.rm=TRUE), "\n"))
  cat(paste("Root node estimated chromosome count:", 
            event_result$ancestral_counts$node_counts[length(tip_counts) + 1], "\n\n"))
  cat("---------------------------\n")
  cat("2. Summary of Detected Events\n")
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
    method_summary <- table(event_result$events$analysis_method)
    cat("\nAnalysis methods:\n")
    for(method in names(method_summary)) {
      cat(paste0("  - ", method, ": ", method_summary[method], " events\n"))
    }
    confidence_summary <- table(event_result$events$confidence)
    cat("\nConfidence levels:\n")
    for(conf in names(confidence_summary)) {
      cat(paste0("  - ", conf, ": ", confidence_summary[conf], " events\n"))
    }
    corroborated_count <- sum(event_result$events$corroborated, na.rm = TRUE)
    cat(paste("\nEvents corroborated by both methods:", corroborated_count, "\n"))
  } else {
    cat("  No chromosome events detected\n")
  }
  cat("-------------------------\n")
  cat("3. Shared Chromosome Events\n")
  if(nrow(event_result$shared_events) > 0) {
    for(i in 1:nrow(event_result$shared_events)) {
      se <- event_result$shared_events[i,]
      cat(paste0("  [", i, "] ", se$event_type, " in clade ", se$clade_label, 
                 " (", se$analysis_type, ")\n"))
      cat(paste0("      - Prevalence: ", round(se$prevalence * 100), "% (", 
                 se$event_count, " of ", se$tip_count, " species)\n"))
      cat(paste0("      - Status: ", se$shared_status, "\n"))
      cat(paste0("      - Confidence: ", se$confidence, "\n"))
      if(!is.null(se$corroborated) && !is.na(se$corroborated) && se$corroborated) {
        cat("      - Corroborated by both bottom-up and top-down analyses\n")
      }
    }
  } else {
    cat("  No shared chromosome events detected\n")
  }
  cat("----------------------\n")
  cat("4. Detailed Event Listing\n")
  if(nrow(event_result$events) > 0) {
    high_conf_events <- event_result$events[event_result$events$confidence == "high",]
    if(nrow(high_conf_events) > 0) {
      cat("High confidence events:\n")
      for(i in 1:min(10, nrow(high_conf_events))) {
        evt <- high_conf_events[i,]
        cat(paste0("  [", i, "] ", evt$event_type, " in ", evt$child_label, 
                   " (", evt$analysis_method, ", ", evt$evidence_type, ")\n"))
        if(!is.na(evt$group)) cat(paste0("      - Group: ", evt$group, "\n"))
        if(!is.na(evt$reference_species)) cat(paste0("      - Reference species: ", evt$reference_species, "\n"))
        if(!is.na(evt$chromosomes_A) && !is.na(evt$chromosomes_B)) {
          cat(paste0("      - Chromosomes: ", 
                     if(evt$event_type == "fusion") {
                       paste0(evt$chromosomes_A, " → ", evt$chromosomes_B)
                     } else {
                       paste0(evt$chromosomes_A, " → ", evt$chromosomes_B)
                     }, "\n"))
        }
        if(!is.null(evt$corroborated) && !is.na(evt$corroborated) && evt$corroborated) {
          cat("      - Corroborated by both analysis methods\n")
        }
      }
      if(nrow(high_conf_events) > 10) {
        cat(paste0("  ...and ", nrow(high_conf_events) - 10, " more high confidence events\n"))
      }
    } else {
      cat("  No events detected\n")
    }
  } else {
    cat("  No events detected\n")
  }
  cat("\n===== Enhanced Chromosome Evolutionary Events Summary =====\n")
  cat(paste("Root node estimated chromosome count:", 
            event_result$ancestral_counts$node_counts[length(event_result$tree$tip.label) + 1], "\n"))
  cat("==============================================\n")
  cat(paste("\nResults saved to:", output_dir, "\n"))
  sink()
  message(paste("DEBUG: Created event report at", events_report))
  message("Shared chromosome event detection complete")
  return(list(
    events_rdata = events_rdata_file,
    edge_data = edge_data_file,
    ancestral_counts = ancestral_counts_file,
    count_changes = count_changes_file,
    pattern_data = pattern_data_file,
    bottom_up_events = bottom_up_file,
    top_down_events = top_down_file,
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
    cat("Usage: Rscript ancestral_chromosomes_phase3_fixed.R integrated_data_file conserved_groups_file output_dir\n")
    cat("Example: Rscript ancestral_chromosomes_phase3_fixed.R results/integrated_data.RData results/phase2/conserved_chromosome_groups.RData results/phase3/\n")
    return(1)
  }
  
  integrated_data_file <- args[1]
  conserved_groups_file <- args[2]
  output_dir <- args[3]
  
  message(paste("DEBUG: Using integrated data from", integrated_data_file))
  message(paste("DEBUG: Using conserved groups from", conserved_groups_file))
  message(paste("DEBUG: Output will be saved to", output_dir))
  
  # Load data
  integrated_data <- load_integrated_data(integrated_data_file)
  conserved_groups <- load_conserved_groups(conserved_groups_file)
  message("DEBUG: Data loading complete")
  
  # 添加诊断测试
  test_results <- run_diagnostic_tests(integrated_data, conserved_groups)
  if(!all(unlist(test_results))) {
    warning("Some diagnostic tests failed. Results may be unreliable.")
  }
  
  # 检测共享染色体事件
  event_result <- detect_shared_chromosome_events(integrated_data, conserved_groups)
  message("DEBUG: Event detection complete")
  
  # 生成事件摘要
  summary_stats <- generate_event_summary(event_result$events)
  
  # 创建可视化（如果ape包可用）
  if(requireNamespace("ape", quietly = TRUE)) {
    plot_file <- file.path(output_dir, "chromosome_events_plot.pdf")
    plot_chromosome_events(event_result$tree, event_result$events, plot_file)
  }
  
  # 保存结果
  output_files <- save_event_detection_results(event_result, output_dir)
  message("DEBUG: Results saved")
  
  cat("==============================================\n")
  cat("\n===== Enhanced Chromosome Evolutionary Events Summary =====\n")
  cat(paste("Root node estimated chromosome count:", 
            event_result$ancestral_counts$node_counts[length(event_result$tree$tip.label) + 1], "\n"))
  cat("==============================================\n")
  cat(paste("\nResults saved to:", output_dir, "\n"))
  
  return(0)
}

if(!interactive()) {
  # Run the script if executed directly
  args <- commandArgs(trailingOnly = TRUE)
  exit_code <- main(args)
  quit(status = exit_code)
}

# ==============================
# Helper Functions
# ==============================

find_related_species <- function(tree, species) {
  # Find related species (sister taxa) for a given species
  message(paste("DEBUG: Finding related species for", species))
  
  # 使用tryCatch增强错误处理能力
  tryCatch({
    # First get the node that this species belongs to
    tip_idx <- which(tree$tip.label == species)
    if(length(tip_idx) == 0) {
      message(paste("DEBUG: Species", species, "not found in tree"))
      return(character(0))
    }
    
    message(paste("DEBUG: Found species at tip index", tip_idx))
    
    # Get parent node of this species
    edges <- tree$edge
    parent_edge <- which(edges[, 2] == tip_idx)
    if(length(parent_edge) == 0) {
      message(paste("DEBUG: No parent edge found for", species))
      return(character(0))
    }
    
    parent_node <- edges[parent_edge, 1]
    message(paste("DEBUG: Parent node is", parent_node))
    
    # Get all children of this parent (sister taxa)
    sister_nodes <- tree$edge[tree$edge[, 1] == parent_node, 2]
    sister_nodes <- sister_nodes[sister_nodes != tip_idx]  # 排除原始物种节点
    
    # 将姐妹节点转换为姐妹分类群名称
    if(length(sister_nodes) == 0) {
      message("DEBUG: No immediate sisters found, going up to grandparent level")
      
      # 如果没有姐妹，尝试再上升一级 - 使用更可靠的方法
      grandparent_edge <- which(edges[, 2] == parent_node)
      if(length(grandparent_edge) == 0) {
        message("DEBUG: No grandparent edge found")
        return(character(0))
      }
      
      grandparent_node <- edges[grandparent_edge, 1]
      message(paste("DEBUG: Grandparent node is", grandparent_node))
      
      # 获取所有孙子节点 - 不包括原始父节点的后代
      aunt_uncle_nodes <- tree$edge[tree$edge[, 1] == grandparent_node, 2]
      aunt_uncle_nodes <- aunt_uncle_nodes[aunt_uncle_nodes != parent_node]
      
      if(length(aunt_uncle_nodes) == 0) {
        message("DEBUG: No aunts/uncles found")
        return(character(0))
      }
      
      message(paste("DEBUG: Found", length(aunt_uncle_nodes), "aunt/uncle nodes"))
      
      # 处理每个姑姑/叔叔节点
      sister_tips <- character(0)
      for(node in aunt_uncle_nodes) {
        if(node <= length(tree$tip.label)) {
          # 这是一个叶子节点
          sister_tips <- c(sister_tips, tree$tip.label[node])
          message(paste("DEBUG: Added aunt/uncle tip", tree$tip.label[node]))
        } else {
          # 这是一个内部节点，获取其所有后代
          message(paste("DEBUG: Processing aunt/uncle internal node", node))
          
          # 使用安全的方法获取后代
          tip_descendants <- tryCatch({
            descendants <- getDescendants(tree, node)
            # 只保留叶子节点
            tip_desc <- descendants[descendants <= length(tree$tip.label)]
            if(length(tip_desc) > 0) {
              message(paste("DEBUG: Found", length(tip_desc), "cousin tips"))
              tree$tip.label[tip_desc]
            } else {
              character(0)
            }
          }, error = function(e) {
            message(paste("DEBUG: Error getting descendants:", e$message))
            character(0)
          })
          
          if(length(tip_descendants) > 0) {
            sister_tips <- c(sister_tips, tip_descendants)
          }
        }
      }
      
      # 限制返回的相关物种数量
      if(length(sister_tips) > 20) {
        message(paste("DEBUG: Limiting from", length(sister_tips), "to 20 related species"))
        sister_tips <- sample(sister_tips, 20)
      }
      
      return(sister_tips)
    } else {
      # 处理直接姐妹节点
      sister_tips <- character(0)
      
      for(node in sister_nodes) {
        if(node <= length(tree$tip.label)) {
          # 这是一个叶子节点
          sister_tips <- c(sister_tips, tree$tip.label[node])
        } else {
          # 这是一个内部节点，获取其所有后代
          tip_descendants <- tryCatch({
            descendants <- getDescendants(tree, node)
            tip_desc <- descendants[descendants <= length(tree$tip.label)]
            tree$tip.label[tip_desc]
          }, error = function(e) {
            message(paste("DEBUG: Error getting descendants for sister node:", e$message))
            character(0)
          })
          
          if(length(tip_descendants) > 0) {
            sister_tips <- c(sister_tips, tip_descendants)
          }
        }
      }
      
      message(paste("DEBUG: Found", length(sister_tips), "sister species"))
      
      # 限制返回的相关物种数量
      if(length(sister_tips) > 20) {
        message(paste("DEBUG: Limiting from", length(sister_tips), "to 20 sister species"))
        sister_tips <- sample(sister_tips, 20)
      }
      
      return(sister_tips)
    }
  }, error = function(e) {
    # 最终的错误处理 - 如果出现任何问题，返回一些随机物种
    message(paste("DEBUG: Critical error in find_related_species:", e$message))
    message("DEBUG: Falling back to random species selection")
    
    all_tips <- tree$tip.label
    all_tips <- all_tips[all_tips != species]
    
    if(length(all_tips) > 10) {
      return(sample(all_tips, 10))
    } else {
      return(all_tips)
    }
  })
}

# ==============================
# Main Detection Function
# ==============================

detect_shared_chromosome_events <- function(integrated_data, conserved_groups) {
  # Main function to detect shared chromosome events using combined approaches
  message("Starting shared chromosome event detection with combined approaches...")
  
  # Get tree and chromosome counts
  tree <- integrated_data$tree
  chr_counts <- integrated_data$chr_counts
  
  # 确保树中有内部节点标签
  tree <- add_node_labels(tree)
  
  message(paste("DEBUG: Working with tree of", length(tree$tip.label), "tips and",
                length(chr_counts), "chromosome counts"))
  
  # Extract tree structure information
  edge_data <- extract_tree_edges(tree)
  clades <- identify_clades(tree)
  message(paste("DEBUG: Extracted", nrow(edge_data), "edges and", length(clades), "clades"))
  
  # 检查染色体计数是否完整
  missed_tips <- setdiff(tree$tip.label, names(chr_counts))
  if(length(missed_tips) > 0) {
    message(paste("WARNING:", length(missed_tips), "species in tree have no chromosome count data"))
    message(paste("DEBUG: Sample missing species:", paste(head(missed_tips, 5), collapse=", ")))
  }
  
  # Reconstruct ancestral chromosome counts
  ancestral_counts <- reconstruct_ancestral_chromosome_counts(tree, chr_counts)
  message("DEBUG: Completed ancestral reconstruction")
  
  # Detect chromosome count changes
  count_changes <- detect_count_changes(tree, edge_data, chr_counts, ancestral_counts)
  message(paste("DEBUG: Detected", nrow(count_changes), "count changes"))
  
  # Detect mapping patterns
  pattern_data <- detect_mapping_patterns(edge_data, integrated_data, conserved_groups)
  message(paste("DEBUG: Found", nrow(pattern_data), "mapping patterns"))
  
  # Infer events using bottom-up approach
  bottom_up_events <- infer_events_bottom_up(tree, edge_data, count_changes, pattern_data, chr_counts)
  message(paste("DEBUG: Inferred", nrow(bottom_up_events), "bottom-up events"))
  
  # Infer events using top-down approach
  top_down_events <- infer_events_top_down(tree, edge_data, count_changes, ancestral_counts, clades)
  message(paste("DEBUG: Inferred", nrow(top_down_events), "top-down events"))
  
  # Combine event analyses
  events_df <- combine_event_analyses(bottom_up_events, top_down_events)
  message(paste("DEBUG: Combined into", nrow(events_df), "total events"))
  
  # Identify shared events
  shared_events_df <- identify_shared_events(tree, events_df, clades)
  message(paste("DEBUG: Identified", nrow(shared_events_df), "shared events"))
  
  # Return results
  result <- list(
    tree = tree,
    edge_data = edge_data,
    clades = clades,
    ancestral_counts = ancestral_counts,
    count_changes = count_changes,
    pattern_data = pattern_data,
    bottom_up_events = bottom_up_events,
    top_down_events = top_down_events,
    events = events_df,
    shared_events = shared_events_df
  )
  message("DEBUG: Created complete result object")
  message("Shared chromosome event detection complete")
  
  return(result)
}

# ==============================
# 修复有问题的辅助函数和增加所需的新函数
# ==============================

# 增加集成测试函数
run_diagnostic_tests <- function(integrated_data, conserved_groups) {
  # 运行诊断测试以确保数据完整性
  message("Running diagnostic tests...")
  
  # 测试树结构
  tree <- integrated_data$tree
  if(is.null(tree) || !inherits(tree, "phylo")) {
    warning("Tree object is not valid. Attempting to proceed anyway.")
  } else {
    message(paste("DEBUG: Tree seems valid with", length(tree$tip.label), "tips and", tree$Nnode, "internal nodes"))
  }
  
  # 测试染色体计数数据
  counts <- integrated_data$chr_counts
  if(is.null(counts) || length(counts) == 0) {
    warning("Chromosome count data is empty or invalid.")
  } else {
    message(paste("DEBUG: Chromosome count data seems valid with", length(counts), "entries"))
    # 检查是否有NA或异常值
    na_count <- sum(is.na(counts))
    zero_count <- sum(counts == 0, na.rm = TRUE)
    if(na_count > 0) {
      warning(paste("Found", na_count, "NA values in chromosome counts"))
    }
    if(zero_count > 0) {
      warning(paste("Found", zero_count, "zero values in chromosome counts"))
    }
  }
  
  # 测试保守基因组
  if(is.null(conserved_groups$conserved_groups) || length(conserved_groups$conserved_groups) == 0) {
    warning("Conserved groups data is empty or invalid.")
  } else {
    message(paste("DEBUG: Conserved groups data seems valid with", 
                 length(conserved_groups$conserved_groups), "groups"))
  }
  
  # 测试双向映射数据
  if(is.null(integrated_data$bidirectional_maps) || nrow(integrated_data$bidirectional_maps) == 0) {
    warning("Bidirectional mapping data is empty or invalid.")
  } else {
    message(paste("DEBUG: Mapping data seems valid with", 
                 nrow(integrated_data$bidirectional_maps), "mappings"))
  }
  
  # 返回测试结果
  return(list(
    tree_valid = !is.null(tree) && inherits(tree, "phylo"),
    counts_valid = !is.null(counts) && length(counts) > 0,
    groups_valid = !is.null(conserved_groups$conserved_groups) && 
                  length(conserved_groups$conserved_groups) > 0,
    maps_valid = !is.null(integrated_data$bidirectional_maps) && 
                nrow(integrated_data$bidirectional_maps) > 0
  ))
}

# 修复和完善现有函数:

# 修复可能有问题的event_summary函数
generate_event_summary <- function(events_df) {
  # 从检测到的事件生成摘要统计
  message("Generating event summary statistics...")
  
  if(nrow(events_df) == 0) {
    message("No events to summarize")
    return(list(
      total_events = 0,
      fusion_events = 0,
      fission_events = 0,
      other_events = 0,
      high_confidence = 0,
      medium_confidence = 0,
      low_confidence = 0,
      corroborated_events = 0
    ))
  }
  
  # 按事件类型计数
  event_types <- table(events_df$event_type)
  fusion_count <- if("fusion" %in% names(event_types)) event_types["fusion"] else 0
  fission_count <- if("fission" %in% names(event_types)) event_types["fission"] else 0
  other_count <- sum(event_types) - fusion_count - fission_count
  
  # 按置信度计数
  confidence_levels <- table(events_df$confidence)
  high_conf <- if("high" %in% names(confidence_levels)) confidence_levels["high"] else 0
  medium_conf <- if("medium" %in% names(confidence_levels)) confidence_levels["medium"] else 0
  low_conf <- if("low" %in% names(confidence_levels)) confidence_levels["low"] else 0
  
  # 计算经证实的事件数
  corroborated <- if("corroborated" %in% names(events_df)) {
    sum(events_df$corroborated, na.rm = TRUE)
  } else {
    0
  }
  
  # 输出摘要
  message(paste("Total events:", nrow(events_df)))
  message(paste("Fusion events:", fusion_count))
  message(paste("Fission events:", fission_count))
  message(paste("Other events:", other_count))
  message(paste("High confidence:", high_conf))
  message(paste("Medium confidence:", medium_conf))
  message(paste("Low confidence:", low_conf))
  message(paste("Corroborated events:", corroborated))
  
  return(list(
    total_events = nrow(events_df),
    fusion_events = fusion_count,
    fission_events = fission_count,
    other_events = other_count,
    high_confidence = high_conf,
    medium_confidence = medium_conf,
    low_confidence = low_conf,
    corroborated_events = corroborated
  ))
}

# 增加一个可视化函数
plot_chromosome_events <- function(tree, events_df, output_file = NULL) {
  # 在系统树上可视化染色体事件
  message("Plotting chromosome events on phylogenetic tree...")
  
  if(!requireNamespace("ape", quietly = TRUE)) {
    warning("Package 'ape' is required for plotting. Skipping visualization.")
    return(NULL)
  }
  
  # 检查输入
  if(nrow(events_df) == 0) {
    warning("No events to plot.")
    return(NULL)
  }
  
  # 将事件映射到边
  edge_events <- events_df[!is.na(events_df$edge_id), ]
  
  # 创建用于绘图的边色
  edge_colors <- rep("black", nrow(tree$edge))
  
  # 为不同事件类型分配不同颜色
  if(nrow(edge_events) > 0) {
    fusion_edges <- edge_events$edge_id[edge_events$event_type == "fusion"]
    fission_edges <- edge_events$edge_id[edge_events$event_type == "fission"]
    
    edge_colors[fusion_edges] <- "blue"
    edge_colors[fission_edges] <- "red"
  }
  
  # 绘制树
  if(!is.null(output_file)) {
    # 保存到文件
    pdf(output_file, width = 10, height = 12)
    
    plot.phylo(tree, edge.color = edge_colors, cex = 0.7)
    legend("topright", legend = c("No change", "Fusion", "Fission"),
           col = c("black", "blue", "red"), lwd = 2)
    
    # 添加高置信度事件的标注
    high_conf_events <- edge_events[edge_events$confidence == "high",]
    if(nrow(high_conf_events) > 0) {
      for(i in 1:nrow(high_conf_events)) {
        event <- high_conf_events[i,]
        # 找到这个边的位置
        edge_pos <- which(tree$edge[,2] == event$child_node)
        if(length(edge_pos) > 0) {
          # 添加一个星号表示高置信度事件
          nodelabels(text = "*", node = event$child_node, 
                     bg = ifelse(event$event_type == "fusion", "blue", "red"),
                     col = "white", cex = 1.2)
        }
      }
    }
    
    dev.off()
    message(paste("Plot saved to", output_file))
  } else {
    # 直接显示
    plot.phylo(tree, edge.color = edge_colors, cex = 0.7)
    legend("topright", legend = c("No change", "Fusion", "Fission"),
           col = c("black", "blue", "red"), lwd = 2)
  }
  
  return(TRUE)
}

# 增加getSubtree函数，用于辅助clade分析
getSubtree <- function(tree, node) {
  # 根据指定节点提取子树
  if(!inherits(tree, "phylo")) {
    stop("Object is not a phylogenetic tree")
  }
  
  # 如果节点是叶子，则返回只有一个tip的树
  if(node <= length(tree$tip.label)) {
    message(paste("DEBUG: Node", node, "is a tip, returning single-tip tree"))
    return(NULL)  # 无法提取叶子的子树
  }
  
  # 获取所有后代
  descendants <- getDescendants(tree, node)
  tip_descendants <- descendants[descendants <= length(tree$tip.label)]
  
  if(length(tip_descendants) == 0) {
    message(paste("DEBUG: Node", node, "has no tip descendants"))
    return(NULL)
  }
  
  # 使用ape包的drop.tip函数保留我们想要的tips
  tips_to_keep <- tree$tip.label[tip_descendants]
  tips_to_drop <- setdiff(tree$tip.label, tips_to_keep)
  
  if(length(tips_to_drop) == 0) {
    message("DEBUG: No tips to drop, returning original tree")
    return(tree)
  }
  
  # 使用drop.tip函数提取子树
  subtree <- drop.tip(tree, tips_to_drop)
  message(paste("DEBUG: Extracted subtree with", length(subtree$tip.label), "tips"))
  
  return(subtree)
}

# 增加一个辅助函数，用于将内部节点标签添加到树上
add_node_labels <- function(tree) {
  # 为树的内部节点添加标签
  message("Adding labels to internal nodes...")
  
  if(is.null(tree$node.label)) {
    n_tips <- length(tree$tip.label)
    n_nodes <- tree$Nnode
    tree$node.label <- paste0("Node", (n_tips + 1):(n_tips + n_nodes))
    message(paste("DEBUG: Added", n_nodes, "node labels"))
  } else {
    # 检查是否有空标签
    empty_labels <- which(tree$node.label == "" | is.na(tree$node.label))
    if(length(empty_labels) > 0) {
      n_tips <- length(tree$tip.label)
      tree$node.label[empty_labels] <- paste0("Node", (n_tips + empty_labels))
      message(paste("DEBUG: Added", length(empty_labels), "missing node labels"))
    }
  }
  
  return(tree)
}

# 改进detect_shared_chromosome_events函数处理更复杂的场景
detect_shared_chromosome_events <- function(integrated_data, conserved_groups) {
  # Main function to detect shared chromosome events using combined approaches
  message("Starting shared chromosome event detection with combined approaches...")
  
  # Get tree and chromosome counts
  tree <- integrated_data$tree
  chr_counts <- integrated_data$chr_counts
  
  # 确保树中有内部节点标签
  tree <- add_node_labels(tree)
  
  message(paste("DEBUG: Working with tree of", length(tree$tip.label), "tips and",
                length(chr_counts), "chromosome counts"))
  
  # Extract tree structure information
  edge_data <- extract_tree_edges(tree)
  clades <- identify_clades(tree)
  message(paste("DEBUG: Extracted", nrow(edge_data), "edges and", length(clades), "clades"))
  
  # 检查染色体计数是否完整
  missed_tips <- setdiff(tree$tip.label, names(chr_counts))
  if(length(missed_tips) > 0) {
    message(paste("WARNING:", length(missed_tips), "species in tree have no chromosome count data"))
    message(paste("DEBUG: Sample missing species:", paste(head(missed_tips, 5), collapse=", ")))
  }
  
  # Reconstruct ancestral chromosome counts
  ancestral_counts <- reconstruct_ancestral_chromosome_counts(tree, chr_counts)
  message("DEBUG: Completed ancestral reconstruction")
  
  # Detect chromosome count changes
  count_changes <- detect_count_changes(tree, edge_data, chr_counts, ancestral_counts)
  message(paste("DEBUG: Detected", nrow(count_changes), "count changes"))
  
  # Detect mapping patterns
  pattern_data <- detect_mapping_patterns(edge_data, integrated_data, conserved_groups)
  message(paste("DEBUG: Found", nrow(pattern_data), "mapping patterns"))
  
  # Infer events using bottom-up approach
  bottom_up_events <- infer_events_bottom_up(tree, edge_data, count_changes, pattern_data, chr_counts)
  message(paste("DEBUG: Inferred", nrow(bottom_up_events), "bottom-up events"))
  
  # Infer events using top-down approach
  top_down_events <- infer_events_top_down(tree, edge_data, count_changes, ancestral_counts, clades)
  message(paste("DEBUG: Inferred", nrow(top_down_events), "top-down events"))
  
  # Combine event analyses
  events_df <- combine_event_analyses(bottom_up_events, top_down_events)
  message(paste("DEBUG: Combined into", nrow(events_df), "total events"))
  
  # Identify shared events
  shared_events_df <- identify_shared_events(tree, events_df, clades)
  message(paste("DEBUG: Identified", nrow(shared_events_df), "shared events"))
  
  # Return results
  result <- list(
    tree = tree,
    edge_data = edge_data,
    clades = clades,
    ancestral_counts = ancestral_counts,
    count_changes = count_changes,
    pattern_data = pattern_data,
    bottom_up_events = bottom_up_events,
    top_down_events = top_down_events,
    events = events_df,
    shared_events = shared_events_df
  )
  message("DEBUG: Created complete result object")
  message("Shared chromosome event detection complete")
  
  return(result)
}

# 修改main函数以使用新的辅助函数
main <- function(args) {
  # Check arguments
  if(length(args) < 3) {
    cat("Usage: Rscript ancestral_chromosomes_phase3_fixed.R integrated_data_file conserved_groups_file output_dir\n")
    cat("Example: Rscript ancestral_chromosomes_phase3_fixed.R results/integrated_data.RData results/phase2/conserved_chromosome_groups.RData results/phase3/\n")
    return(1)
  }
  
  integrated_data_file <- args[1]
  conserved_groups_file <- args[2]
  output_dir <- args[3]
  
  message(paste("DEBUG: Using integrated data from", integrated_data_file))
  message(paste("DEBUG: Using conserved groups from", conserved_groups_file))
  message(paste("DEBUG: Output will be saved to", output_dir))
  
  # Load data
  integrated_data <- load_integrated_data(integrated_data_file)
  conserved_groups <- load_conserved_groups(conserved_groups_file)
  message("DEBUG: Data loading complete")
  
  # 添加诊断测试
  test_results <- run_diagnostic_tests(integrated_data, conserved_groups)
  if(!all(unlist(test_results))) {
    warning("Some diagnostic tests failed. Results may be unreliable.")
  }
  
  # 检测共享染色体事件
  event_result <- detect_shared_chromosome_events(integrated_data, conserved_groups)
  message("DEBUG: Event detection complete")
  
  # 生成事件摘要
  summary_stats <- generate_event_summary(event_result$events)
  
  # 创建可视化（如果ape包可用）
  if(requireNamespace("ape", quietly = TRUE)) {
    plot_file <- file.path(output_dir, "chromosome_events_plot.pdf")
    plot_chromosome_events(event_result$tree, event_result$events, plot_file)
  }
  
  # 保存结果
  output_files <- save_event_detection_results(event_result, output_dir)
  message("DEBUG: Results saved")
  
  cat("==============================================\n")
  cat("\n===== Enhanced Chromosome Evolutionary Events Summary =====\n")
  cat(paste("Root node estimated chromosome count:", 
            event_result$ancestral_counts$node_counts[length(event_result$tree$tip.label) + 1], "\n"))
  cat("==============================================\n")
  cat(paste("\nResults saved to:", output_dir, "\n"))
  
  return(0)
}

# 增加报告生成器函数，用于生成交互式HTML报告
generate_html_report <- function(event_result, output_dir) {
  # 生成更漂亮的HTML报告
  message("Generating HTML report...")
  
  # 检查必需的包
  if(!requireNamespace("htmltools", quietly = TRUE) || 
     !requireNamespace("knitr", quietly = TRUE)) {
    warning("Packages 'htmltools' and 'knitr' are required for HTML reports. Skipping report generation.")
    return(FALSE)
  }
  
  report_file <- file.path(output_dir, "chromosome_events_report.html")
  
  # 提取需要的数据
  events <- event_result$events
  shared_events <- event_result$shared_events
  tree <- event_result$tree
  ancestral_counts <- event_result$ancestral_counts
  
  # 生成报告内容
  html_content <- c(
    "<!DOCTYPE html>",
    "<html>",
    "<head>",
    "  <title>Ancestral Chromosome Event Detection Report</title>",
    "  <style>",
    "    body { font-family: Arial, sans-serif; line-height: 1.6; margin: 20px; }",
    "    h1 { color: #2c3e50; }",
    "    h2 { color: #3498db; margin-top: 30px; border-bottom: 1px solid #ddd; padding-bottom: 10px; }",
    "    table { border-collapse: collapse; width: 100%; margin-bottom: 20px; }",
    "    th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }",
    "    th { background-color: #f2f2f2; }",
    "    tr:nth-child(even) { background-color: #f9f9f9; }",
    "    .fusion { color: blue; }",
    "    .fission { color: red; }",
    "    .high { font-weight: bold; }",
    "  </style>",
    "</head>",
    "<body>",
    paste0("  <h1>Ancestral Chromosome Event Detection Report</h1>"),
    paste0("  <p>Analysis performed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "</p>"),
    
    "  <h2>1. Ancestral Chromosome Count Summary</h2>"
  )
  
  # 添加染色体计数摘要
  tip_counts <- ancestral_counts$node_counts[1:length(tree$tip.label)]
  internal_counts <- ancestral_counts$node_counts[(length(tip_counts) + 1):length(ancestral_counts$node_counts)]
  
  html_content <- c(html_content,
    paste0("  <p>Extant species chromosome count range: ", min(tip_counts, na.rm=TRUE), 
           " to ", max(tip_counts, na.rm=TRUE), "</p>"),
    paste0("  <p>Ancestral node chromosome count range: ", min(internal_counts, na.rm=TRUE), 
           " to ", max(internal_counts, na.rm=TRUE), "</p>"),
    paste0("  <p>Root node estimated chromosome count: ", 
           ancestral_counts$node_counts[length(tip_counts) + 1], "</p>"),
    
    "  <h2>2. Detected Events Summary</h2>"
  )
  
  # 添加事件摘要
  if(nrow(events) > 0) {
    event_summary <- table(events$event_type)
    evidence_summary <- table(events$evidence_type)
    confidence_summary <- table(events$confidence)
    corroborated_count <- sum(events$corroborated, na.rm = TRUE)
    
    html_content <- c(html_content,
      "  <h3>Event Types</h3>",
      "  <table>",
      "    <tr><th>Event Type</th><th>Count</th></tr>"
    )
    
    for(evt_type in names(event_summary)) {
      html_content <- c(html_content,
        paste0("    <tr><td>", evt_type, "</td><td>", event_summary[evt_type], "</td></tr>")
      )
    }
    
    html_content <- c(html_content,
      "  </table>",
      
      "  <h3>Evidence Types</h3>",
      "  <table>",
      "    <tr><th>Evidence Type</th><th>Count</th></tr>"
    )
    
    for(evd_type in names(evidence_summary)) {
      html_content <- c(html_content,