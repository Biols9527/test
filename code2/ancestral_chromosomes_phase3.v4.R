#!/usr/bin/env Rscript

# Ancestral Chromosome Reconstruction Phase 3: Enhanced Shared Chromosome Event Detection
# Date: 2025-03-21
# Description: Identifies shared chromosome evolutionary events (fusions, fissions, etc.)
#              across the phylogenetic tree based on conserved groups and mapping data,
#              using combined bottom-up and top-down approaches with statistical validation.

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

load_genomic_data <- function(genomic_data_file) {
  # Load genomic evidence data if available
  if(!is.null(genomic_data_file) && file.exists(genomic_data_file)) {
    message("Loading genomic evidence data...")
    genomic_data <- readRDS(genomic_data_file)
    message(paste("Loaded genomic data for", length(unique(genomic_data$species)), "species"))
    return(genomic_data)
  } else {
    message("Genomic evidence file not found, proceeding without genomic integration")
    return(NULL)
  }
}

# ==============================
# Tree Analysis Functions
# ==============================

prepare_phylogenetic_tree <- function(tree) {
  # 预处理系统发育树，确保其适合祖先状态重建分析
  message("Preparing phylogenetic tree...")
  
  # 检查树是否有效
  if(is.null(tree) || !inherits(tree, "phylo")) {
    stop("Invalid phylogenetic tree object")
  }
  
  # 创建树的副本避免修改原始对象
  tree_copy <- tree
  has_issues <- FALSE
  
  # 检查并修复根节点问题
  if(!is.rooted(tree_copy)) {
    message("Tree is not rooted. Attempting to root at midpoint...")
    has_issues <- TRUE
    
    if(requireNamespace("phytools", quietly = TRUE)) {
      # 尝试使用中点重建(更好的方法)
      tree_copy <- tryCatch({
        phytools::midpoint.root(tree_copy)
      }, error = function(e) {
        message("Midpoint rooting failed, rooting at first tip...")
        root(tree_copy, outgroup = tree_copy$tip.label[1], resolve.root = TRUE)
      })
    } else {
      # 如果phytools不可用，使用简单方法
      tree_copy <- root(tree_copy, outgroup = tree_copy$tip.label[1], resolve.root = TRUE)
    }
    message("Tree has been rooted")
  }
  
  # 检查并解决多分叉问题
  if(!is.binary(tree_copy)) {
    message("Tree contains multifurcations. Resolving to binary tree...")
    has_issues <- TRUE
    tree_copy <- multi2di(tree_copy)
    message("Multifurcations in the tree have been resolved")
  }
  
  # 确保所有分支长度都是有效的正数
  if(is.null(tree_copy$edge.length)) {
    message("Tree has no branch lengths. Adding uniform branch lengths...")
    has_issues <- TRUE
    tree_copy$edge.length <- rep(1, nrow(tree_copy$edge))
  } else {
    # 检查零或近零分支长度
    zero_branches <- which(tree_copy$edge.length <= 1e-8)
    if(length(zero_branches) > 0) {
      message(paste("Found", length(zero_branches), "branches with zero or near-zero lengths. Fixing..."))
      has_issues <- TRUE
      
      # 将0或近0的分支长度替换为小的正值(0.0001)
      tree_copy$edge.length[zero_branches] <- 0.0001
      
      # 计算替换后的总长度，以确认修复成功
      total_length <- sum(tree_copy$edge.length)
      message(paste("After fixing, total tree length:", round(total_length, 4)))
    }
  }
  
  # 确保节点编号连续且有效
  if(max(tree_copy$edge) != (length(tree_copy$tip.label) + tree_copy$Nnode)) {
    message("Tree has non-sequential node numbering. Renumbering...")
    has_issues <- TRUE
    tree_copy <- ape::reorder.phylo(tree_copy, "postorder")
  }
  
  # 检查分支长度的分布情况，发现异常分支
  branch_stats <- summary(tree_copy$edge.length)
  if(branch_stats["Max."] > 1000 * branch_stats["Median"]) {
    message("WARNING: Tree has extremely long branches relative to median length")
    message(paste("Branch length summary - Min:", branch_stats["Min."], 
                  "Median:", branch_stats["Median"], 
                  "Max:", branch_stats["Max."]))
  }
  
  if(has_issues) {
    message("Tree has been corrected and is now suitable for analysis")
  } else {
    message("Tree is already in suitable format")
  }
  
  return(tree_copy)
}

extract_tree_edges <- function(tree) {
  # Extract edge information from the phylogenetic tree
  message("Extracting tree edges and nodes...")
  
  # Get number of species (tips)
  n_tips <- length(tree$tip.label)
  
  # Use data.table for better performance if available
  if(requireNamespace("data.table", quietly = TRUE)) {
    # Create edge data
    edge_data <- data.table::data.table(
      edge_id = 1:nrow(tree$edge),
      parent_node_id = tree$edge[, 1],
      child_node_id = tree$edge[, 2],
      edge_length = tree$edge.length
    )
    
    # Use efficient vectorized operations
    edge_data[, parent_is_tip := parent_node_id <= n_tips]
    edge_data[, child_is_tip := child_node_id <= n_tips]
    
    # Create node label mapping
    node_labels <- character(max(tree$edge))
    node_labels[1:n_tips] <- tree$tip.label
    
    if(!is.null(tree$node.label)) {
      node_labels[(n_tips+1):length(node_labels)] <- tree$node.label
    } else {
      node_labels[(n_tips+1):length(node_labels)] <- paste0("Node", (n_tips+1):length(node_labels))
    }
    
    # Vectorized label assignment
    edge_data[, parent_label := node_labels[parent_node_id]]
    edge_data[, child_label := node_labels[child_node_id]]
    
    # Convert back to data.frame for interface consistency
    edge_data <- as.data.frame(edge_data)
  } else {
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
  }
  
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
perform_ml_reconstruction <- function(tree, count_vector) {
  # 专门处理ML方法的祖先状态重建，增强健壮性
  message("Attempting ML ancestral state reconstruction...")
  
  # 验证输入
  if(length(count_vector) != length(tree$tip.label)) {
    message(paste("Warning: count_vector length", length(count_vector), 
                "doesn't match number of tips", length(tree$tip.label)))
    # 确保count_vector与树的tip.label匹配
    valid_counts <- count_vector[names(count_vector) %in% tree$tip.label]
    if(length(valid_counts) == 0) {
      stop("No valid chromosome counts match tree tip labels")
    }
    
    # 重建完整的count_vector
    full_vector <- numeric(length(tree$tip.label))
    names(full_vector) <- tree$tip.label
    for(tip in tree$tip.label) {
      if(tip %in% names(count_vector)) {
        full_vector[tip] <- count_vector[tip]
      } else {
        full_vector[tip] <- NA
      }
    }
    count_vector <- full_vector
  }
  
  # 处理缺失值
  na_count <- sum(is.na(count_vector))
  if(na_count > 0) {
    message(paste("Note:", na_count, "missing values in count_vector"))
    # 这里不做插补，让ace函数自己处理
  }
  
  # 尝试ML方法
  tryCatch({
    # 使用ace函数的标准调用，但显式设置所有参数
    ml_result <- ace(count_vector, tree, type = "continuous", method = "ML",
                    model = "BM", scaled = FALSE, kappa = 1, corStruct = NULL)
    
    # 验证结果
    if(is.null(ml_result$ace) || length(ml_result$ace) == 0) {
      message("ML reconstruction returned empty results")
      return(NULL)
    }
    
    message("ML ancestral state reconstruction successful")
    return(ml_result)
  }, error = function(e) {
    message(paste("ML reconstruction error:", e$message))
    
    # 尝试备选方法
    tryCatch({
      message("Trying ML with restricted maximum likelihood...")
      ml_result <- ace(count_vector, tree, type = "continuous", method = "REML")
      
      if(is.null(ml_result$ace) || length(ml_result$ace) == 0) {
        message("REML reconstruction returned empty results")
        return(NULL)
      }
      
      message("REML ancestral state reconstruction successful")
      return(ml_result)
    }, error = function(e2) {
      message(paste("REML reconstruction error:", e2$message))
      return(NULL)
    })
  })
}

perform_bayesian_reconstruction <- function(tree, count_vector) {
  # 专门处理贝叶斯方法的祖先状态重建，增强健壮性
  message("Attempting Bayesian ancestral state reconstruction...")
  
  if(!requireNamespace("phytools", quietly = TRUE)) {
    message("Package 'phytools' required for Bayesian reconstruction but not available")
    return(NULL)
  }
  
  # 验证树是否适合贝叶斯分析
  if(any(tree$edge.length <= 1e-8)) {
    message("Tree contains zero-length branches, not suitable for Bayesian analysis")
    return(NULL)
  }
  
  # 确保计数向量格式正确
  if(!all(names(count_vector) %in% tree$tip.label)) {
    message("Warning: Some count_vector names don't match tree tip labels")
    # 确保计数向量与树匹配
    fixed_counts <- count_vector[names(count_vector) %in% tree$tip.label]
    if(length(fixed_counts) < 5) { # 至少需要5个有效数据点
      message("Too few valid tip values for Bayesian analysis")
      return(NULL)
    }
    count_vector <- fixed_counts
  }
  
  # 简化MCMC参数，提高成功率
  ngen <- 30000        # 减少代数以避免长时间运行
  samplefreq <- 100    # 采样频率
  burnin <- 0.2        # burn-in比例
  control_params <- list(
    sig2 = var(count_vector, na.rm = TRUE), # 使用数据方差作为初始值
    quiet = TRUE                           # 静默运行
  )
  
  # 尝试贝叶斯MCMC
  tryCatch({
    message(paste("Running Bayesian MCMC with", ngen, "generations..."))
    
    bayes_result <- phytools::anc.Bayes(
      tree = tree,
      x = count_vector,
      ngen = ngen,
      samplefreq = samplefreq,
      control = control_params
    )
    
    # 验证结果
    if(is.null(bayes_result$mcmc) || ncol(bayes_result$mcmc) < 2) {
      message("Bayesian analysis returned invalid results")
      return(NULL)
    }
    
    message("Bayesian ancestral state reconstruction successful")
    return(bayes_result)
  }, error = function(e) {
    message(paste("Bayesian reconstruction error:", e$message))
    return(NULL)
  })
}

reconstruct_ancestral_chromosome_counts <- function(tree, chr_counts, method = "ML") {
  # 祖先染色体计数重建，整合各种方法并增强健壮性
  message(paste("Reconstructing ancestral chromosome counts using", method, "method..."))
  
  # 预处理进化树，确保其符合分析要求
  tree <- prepare_phylogenetic_tree(tree)
  
  n_tips <- length(tree$tip.label)
  n_nodes <- tree$Nnode
  
  # 构建计数向量
  count_vector <- numeric(n_tips)
  names(count_vector) <- tree$tip.label
  
  for(sp in names(chr_counts)) {
    if(sp %in% tree$tip.label) {
      count_vector[sp] <- chr_counts[sp]
    }
  }
  
  # 测试系统发育信号
  phylo_signal <- FALSE
  if(requireNamespace("phytools", quietly = TRUE)) {
    tryCatch({
      lambda_test <- phytools::phylosig(tree, count_vector, method="lambda", test=TRUE)
      message(paste("Phylogenetic signal test (Pagel's lambda):", 
                    round(lambda_test$lambda, 3), 
                    "p-value:", round(lambda_test$P, 4)))
      if(lambda_test$P < 0.05) {
        phylo_signal <- TRUE 
        message("Significant phylogenetic signal detected, using phylogenetic methods")
      } else {
        message("WARNING: Low phylogenetic signal, ancestral reconstruction may be less reliable")
      }
    }, error = function(e) {
      message("Phylogenetic signal test failed, proceeding with standard methods")
    })
  }
  
  # 系统发育插补缺失值
  if(any(is.na(count_vector))) {
    missing_tips <- names(count_vector)[is.na(count_vector)]
    message(paste("Missing chromosome counts for", length(missing_tips), "species. Using phylogenetic imputation."))
    
    for(sp in missing_tips) {
      # 获取近缘物种
      related_species <- find_related_species(tree, sp)
      related_counts <- count_vector[names(count_vector) %in% related_species]
      
      if(length(na.omit(related_counts)) > 0) {
        # 使用近缘物种的中位数
        count_vector[sp] <- median(related_counts, na.rm=TRUE)
        message(paste(" - Imputed", sp, "with value", count_vector[sp], 
                    "based on", length(na.omit(related_counts)), "related species"))
      } else {
        # 如果没有近缘物种数据，使用全局中位数
        count_vector[sp] <- median(count_vector, na.rm=TRUE)
        message(paste(" - Imputed", sp, "with global median", count_vector[sp]))
      }
    }
  }
  
  # 创建节点计数映射
  node_counts <- numeric(n_tips + n_nodes)
  names(node_counts) <- 1:(n_tips + n_nodes)
  
  # 添加叶节点计数
  for(i in 1:n_tips) {
    node_counts[i] <- count_vector[tree$tip.label[i]]
  }
  
  # 祖先状态重建方法选择
  chosen_method <- method
  ace_result <- NULL
  success <- FALSE
  
  # 贝叶斯方法
  if(method == "bayesian") {
    bayes_result <- perform_bayesian_reconstruction(tree, count_vector)
    
    if(!is.null(bayes_result)) {
      # 处理MCMC样本
      mcmc_samples <- bayes_result$mcmc
      
      # 应用burn-in
      burnin_samples <- floor(nrow(mcmc_samples) * 0.2)
      if(burnin_samples > 0) {
        mcmc_samples <- mcmc_samples[-(1:burnin_samples),]
      }
      
      # 创建后验摘要
      node_posteriors <- list()
      
      # 处理内部节点
      for(i in 1:n_nodes) {
        node_id <- i + n_tips
        node_col <- paste0("node", i)
        
        if(node_col %in% colnames(mcmc_samples)) {
          # 获取后验分布
          posterior <- mcmc_samples[, node_col]
          
          # 存储点估计(中位数)
          node_counts[node_id] <- round(median(posterior, na.rm = TRUE))
          
          # 存储后验信息
          node_posteriors[[as.character(node_id)]] <- list(
            mean = mean(posterior, na.rm = TRUE),
            median = median(posterior, na.rm = TRUE),
            ci95 = quantile(posterior, c(0.025, 0.975), na.rm = TRUE)
          )
        }
      }
      
      message("Successfully reconstructed ancestral states using Bayesian method")
      success <- TRUE
      result <- list(
        node_counts = node_counts,
        posteriors = node_posteriors,
        method = "bayesian",
        phylo_signal = phylo_signal,
        mcmc_result = bayes_result
      )
      return(result)
    } else {
      message("Bayesian method failed, falling back to ML")
      chosen_method <- "ML"
    }
  }
  
  # ML方法
  if(chosen_method == "ML" || chosen_method == "REML") {
    ml_result <- perform_ml_reconstruction(tree, count_vector)
    
    if(!is.null(ml_result)) {
      # 提取祖先状态
      for(i in 1:n_nodes) {
        node_id <- i + n_tips
        if("ace" %in% names(ml_result)) {
          # 标准的ace结果格式 (矩阵型)
          if(is.matrix(ml_result$ace) && nrow(ml_result$ace) >= i) {
            node_counts[node_id] <- round(ml_result$ace[i, 1])
          } 
          # 向量格式的ace结果
          else if(length(ml_result$ace) >= i) {
            node_counts[node_id] <- round(ml_result$ace[i])
          }
        }
      }
      
      message("Successfully reconstructed ancestral states using", ml_result$method, "method")
      success <- TRUE
      result <- list(
        node_counts = node_counts,
        ace_result = ml_result,
        method = ml_result$method,
        phylo_signal = phylo_signal
      )
      return(result)
    } else {
      message("ML methods failed, falling back to parsimony")
      chosen_method <- "parsimony"
    }
  }
  
  # 简约法
  if(chosen_method == "parsimony" || !success) {
    tryCatch({
      ace_result <- ace(count_vector, tree, type = "continuous", method = "pic")
      
      # 提取祖先状态
      for(i in 1:n_nodes) {
        node_id <- i + n_tips
        if(length(ace_result$ace) >= i) {
          node_counts[node_id] <- round(ace_result$ace[i])
        }
      }
      
      message("Successfully reconstructed ancestral states using parsimony method")
      success <- TRUE
      result <- list(
        node_counts = node_counts,
        ace_result = ace_result,
        method = "parsimony",
        phylo_signal = phylo_signal
      )
      return(result)
    }, error = function(e) {
      message("Parsimony ancestral state reconstruction failed:", e$message)
    })
  }
  
  # 如果所有统计方法都失败，使用基于拓扑的插值
  if(!success) {
    message("All statistical methods failed, using topology-based interpolation")
    
    # 从叶节点开始向上填充
    for(node_id in (n_tips+1):(n_tips+n_nodes)) {
      # 获取该节点的所有直接子节点
      child_nodes <- tree$edge[tree$edge[,1] == node_id, 2]
      
      # 获取子节点的染色体数
      child_counts <- node_counts[as.character(child_nodes)]
      
      # 使用子节点的加权中位数（按分支长度加权）
      child_weights <- rep(1, length(child_counts))
      
      # 如果有分支长度信息，使用其作为权重的倒数
      if(!is.null(tree$edge.length)) {
        branch_lengths <- tree$edge.length[tree$edge[,1] == node_id]
        if(length(branch_lengths) > 0 && all(!is.na(branch_lengths))) {
          child_weights <- 1/branch_lengths
        }
      }
      
      node_counts[as.character(node_id)] <- round(weighted.median(child_counts, child_weights, na.rm=TRUE))
    }
    
    # 创建基本的ace_result以保持兼容性
    ace_result <- list(
      ace = matrix(node_counts[(n_tips+1):(n_tips+n_nodes)], nrow=n_nodes, ncol=1),
      method = "topology_based",
      Nnode = n_nodes
    )
    
    result <- list(
      node_counts = node_counts,
      ace_result = ace_result,
      method = "topology_based",
      phylo_signal = phylo_signal
    )
  }
  
  # 验证结果的生物学合理性
  min_anc_count <- min(node_counts[(n_tips+1):(n_tips+n_nodes)], na.rm=TRUE)
  max_anc_count <- max(node_counts[(n_tips+1):(n_tips+n_nodes)], na.rm=TRUE)
  
  if(min_anc_count < 3) {
    message("WARNING: Reconstructed ancestral counts include very low values (<3), which may be biologically implausible")
  }
  
  if(max_anc_count > 100) {
    message("WARNING: Reconstructed ancestral counts include very high values (>100), which may indicate reconstruction artifacts")
  }
  
  message("Completed ancestral chromosome count reconstruction")
  
  return(result)
}

# Helper function: weighted median
weighted.median <- function(x, w, na.rm = FALSE) {
  if(na.rm) {
    valid <- !is.na(x) & !is.na(w)
    x <- x[valid]
    w <- w[valid]
  }
  
  if(length(x) == 0) return(NA)
  if(length(x) == 1) return(x)
  
  # Sort values
  order_idx <- order(x)
  x <- x[order_idx]
  w <- w[order_idx]
  
  # Calculate cumulative weights
  cum_w <- cumsum(w)
  total_w <- sum(w)
  
  # Find median
  median_idx <- which(cum_w >= total_w/2)[1]
  return(x[median_idx])
}

detect_count_changes <- function(tree, edge_data, chr_counts, ancestral_counts) {
  # Detect chromosome count changes with adaptive thresholds
  message("Detecting chromosome count changes with adaptive thresholds...")
  
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
    branch_length = edge_data$edge_length,
    standardized_change_rate = NA_real_,
    stringsAsFactors = FALSE
  )
  
  # Fill in known counts
  for(i in 1:nrow(count_changes)) {
    parent_id <- edge_data$parent_node_id[i]
    child_id <- edge_data$child_node_id[i]
    
    # Get counts from ancestral reconstruction
    if(as.character(parent_id) %in% names(ancestral_counts$node_counts)) {
      count_changes$parent_count[i] <- ancestral_counts$node_counts[as.character(parent_id)]
    }
    
    if(as.character(child_id) %in% names(ancestral_counts$node_counts)) {
      count_changes$child_count[i] <- ancestral_counts$node_counts[as.character(child_id)]
    }
    
    # For tips, prioritize actual observed counts
    if(count_changes$child_is_tip[i] && edge_data$child_label[i] %in% names(chr_counts)) {
      count_changes$child_count[i] <- chr_counts[edge_data$child_label[i]]
    }
    
    # Calculate count difference and ratio
    if(!is.na(count_changes$parent_count[i]) && !is.na(count_changes$child_count[i])) {
      count_changes$count_diff[i] <- count_changes$child_count[i] - count_changes$parent_count[i]
      
      if(count_changes$parent_count[i] > 0) {
        count_changes$count_ratio[i] <- count_changes$child_count[i] / count_changes$parent_count[i]
      }
      
      # Calculate standardized change rate
      if(!is.na(count_changes$branch_length[i]) && count_changes$branch_length[i] > 0) {
        count_changes$standardized_change_rate[i] <- abs(count_changes$count_diff[i]) / count_changes$branch_length[i]
      }
    }
  }
  
  # Calculate data-driven thresholds
  valid_diffs <- count_changes$count_diff[!is.na(count_changes$count_diff)]
  
  if(length(valid_diffs) > 4) {
    # Minor threshold - keep ability to detect single chromosome changes
    minor_threshold <- 1
    
    # Medium threshold - use 75% quantile or standard deviation
    sd_threshold <- sd(abs(valid_diffs), na.rm=TRUE)
    quant75_threshold <- quantile(abs(valid_diffs), 0.75, na.rm=TRUE)
    medium_threshold <- max(minor_threshold + 1, min(sd_threshold, quant75_threshold))
    
    # Major threshold - use 90% quantile
    major_threshold <- max(medium_threshold + 1, quantile(abs(valid_diffs), 0.90, na.rm=TRUE))
    
    message(paste("Adaptive thresholds - Minor:", minor_threshold, 
                "Medium:", round(medium_threshold, 1), 
                "Major:", round(major_threshold, 1)))
  } else {
    # Not enough data, use reasonable defaults
    minor_threshold <- 1
    medium_threshold <- 2
    major_threshold <- 4
    message("Insufficient data for adaptive thresholds, using defaults (1, 2, 4)")
  }
  
  # Apply finer change classification
  for(i in 1:nrow(count_changes)) {
    diff <- count_changes$count_diff[i]
    if(is.na(diff)) next
    
    if(diff > major_threshold) {
      count_changes$count_change_type[i] <- "major_increase"
    } else if(diff > medium_threshold) {
      count_changes$count_change_type[i] <- "medium_increase"
    } else if(diff > minor_threshold) {
      count_changes$count_change_type[i] <- "minor_increase"
    } else if(diff < -major_threshold) {
      count_changes$count_change_type[i] <- "major_decrease"
    } else if(diff < -medium_threshold) {
      count_changes$count_change_type[i] <- "medium_decrease"
    } else if(diff < -minor_threshold) {
      count_changes$count_change_type[i] <- "minor_decrease"
    } else {
      count_changes$count_change_type[i] <- "stable"
    }
  }
  
  # Summarize change types
  change_summary <- table(count_changes$count_change_type)
  if(length(change_summary) > 0) {
    message("Change distribution:")
    for(change_type in names(change_summary)) {
      message(paste(" -", change_type, ":", change_summary[change_type]))
    }
  }
  
  message(paste("Detected count changes for", sum(!is.na(count_changes$count_diff)), "edges"))
  
  return(list(
    count_changes = count_changes,
    thresholds = list(
      minor = minor_threshold,
      medium = medium_threshold, 
      major = major_threshold
    )
  ))
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
# Bottom-Up Event Inference Functions
# ==============================

infer_events_bottom_up <- function(tree, edge_data, count_changes, pattern_data, chr_counts) {
  # Infer chromosome evolutionary events using bottom-up approach (from tips to internal nodes)
  message("Inferring chromosome events using bottom-up approach...")
  
  # Get the actual count_changes data
  count_changes_df <- count_changes$count_changes
  thresholds <- count_changes$thresholds
  
  # Get unique pattern types
  if(nrow(pattern_data) > 0) {
    fusion_patterns <- pattern_data[pattern_data$pattern_type == "fusion_candidate",]
    fission_patterns <- pattern_data[pattern_data$pattern_type == "fission_candidate",]
  } else {
    fusion_patterns <- fission_patterns <- data.frame()
  }
  
  # Initialize events list
  events <- list()
  
  # Map change types to event types
  change_type_map <- list(
    major_increase = "fission",
    medium_increase = "fission",
    minor_increase = "minor_fission",
    major_decrease = "fusion",
    medium_decrease = "fusion",
    minor_decrease = "minor_fusion",
    stable = "stable"
  )
  
  # Infer events based on chromosome count changes
  for(i in 1:nrow(count_changes_df)) {
    # Skip if no count data
    if(is.na(count_changes_df$child_count[i]) || is.na(count_changes_df$parent_count[i])) {
      next
    }
    
    child_species <- count_changes_df$child_label[i]
    child_count <- count_changes_df$child_count[i]
    parent_count <- count_changes_df$parent_count[i]
    count_diff <- count_changes_df$count_diff[i]
    change_type <- count_changes_df$count_change_type[i]
    
    # Skip stable changes
    if(change_type == "stable") {
      next
    }
    
    # Determine event type and confidence from change type
    event_type <- change_type_map[[change_type]]
    if(is.null(event_type)) event_type <- "unknown" # Fallback
    
    confidence <- if(grepl("major", change_type)) "high" else if(grepl("medium", change_type)) "medium" else "low"
    
    # Create event record
    events[[length(events) + 1]] <- list(
      edge_id = count_changes_df$edge_id[i],
      parent_node = count_changes_df$parent_node_id[i],
      child_node = count_changes_df$child_node_id[i],
      parent_label = count_changes_df$parent_label[i],
      child_label = child_species,
      event_type = event_type,
      parent_count = parent_count,
      child_count = child_count,
      count_diff = abs(count_diff),
      count_ratio = if(count_diff > 0) count_changes_df$count_ratio[i] else 1/count_changes_df$count_ratio[i],
      evidence_type = "count_difference",
      confidence = confidence,
      analysis_method = "bottom_up",
      group = NA_character_,
      chromosomes_A = NA_character_,
      chromosomes_B = NA_character_,
      reference_species = NA_character_
    )
    
    # If child is a tip, check for supporting mapping patterns
    if(count_changes_df$child_is_tip[i]) {
      # Find related species for comparison
      sister_species <- find_related_species(tree, child_species)
      
      # Skip if no related species
      if(length(sister_species) == 0) {
        next
      }
      
      # Check for fusion evidence
      if(nrow(fusion_patterns) > 0 && grepl("decrease", change_type)) {
        fusion_evidence <- fusion_patterns[
          (fusion_patterns$species_B == child_species) & 
            (fusion_patterns$species_A %in% sister_species),
        ]
        
        if(nrow(fusion_evidence) > 0) {
          for(j in 1:nrow(fusion_evidence)) {
            events[[length(events) + 1]] <- list(
              edge_id = count_changes_df$edge_id[i],
              parent_node = count_changes_df$parent_node_id[i],
              child_node = count_changes_df$child_node_id[i],
              parent_label = count_changes_df$parent_label[i],
              child_label = child_species,
              event_type = "fusion",
              parent_count = parent_count,
              child_count = child_count,
              count_diff = NA_integer_,
              count_ratio = fusion_evidence$chromosome_ratio[j],
              evidence_type = "mapping_pattern",
              confidence = "high",
              analysis_method = "bottom_up",
              group = fusion_evidence$group[j],
              chromosomes_A = fusion_evidence$chromosomes_A[j],
              chromosomes_B = fusion_evidence$chromosomes_B[j],
              reference_species = fusion_evidence$species_A[j]
            )
          }
        }
      }
      
      # Check for fission evidence
      if(nrow(fission_patterns) > 0 && grepl("increase", change_type)) {
        fission_evidence <- fission_patterns[
          (fission_patterns$species_B == child_species) & 
            (fission_patterns$species_A %in% sister_species),
        ]
        
        if(nrow(fission_evidence) > 0) {
          for(j in 1:nrow(fission_evidence)) {
            events[[length(events) + 1]] <- list(
              edge_id = count_changes_df$edge_id[i],
              parent_node = count_changes_df$parent_node_id[i],
              child_node = count_changes_df$child_node_id[i],
              parent_label = count_changes_df$parent_label[i],
              child_label = child_species,
              event_type = "fission",
              parent_count = parent_count,
              child_count = child_count,
              count_diff = NA_integer_,
              count_ratio = fission_evidence$chromosome_ratio[j],
              evidence_type = "mapping_pattern",
              confidence = "high",
              analysis_method = "bottom_up",
              group = fission_evidence$group[j],
              chromosomes_A = fission_evidence$chromosomes_A[j],
              chromosomes_B = fission_evidence$chromosomes_B[j],
              reference_species = fission_evidence$species_A[j]
            )
          }
        }
      }
    }
  }
  
  # Convert to data frame
  if(length(events) > 0) {
    # Create a data frame directly from the list of events
    events_df <- do.call(rbind, lapply(events, function(event) {
      as.data.frame(event, stringsAsFactors = FALSE)
    }))
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
      group = character(),
      chromosomes_A = character(),
      chromosomes_B = character(),
      reference_species = character(),
      stringsAsFactors = FALSE
    )
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
  
  # Get the actual count_changes data
  count_changes_df <- count_changes$count_changes
  
  # Initialize events list
  events <- list()
  
  # Get all internal nodes
  n_tips <- length(tree$tip.label)
  internal_nodes <- unique(edge_data$parent_node_id[!edge_data$parent_is_tip])
  
  # For each internal node
  for(node_id in internal_nodes) {
    # Get node label
    node_label <- if(node_id <= n_tips) {
      tree$tip.label[node_id]
    } else {
      if(!is.null(tree$node.label)) {
        tree$node.label[node_id - n_tips]
      } else {
        paste0("Node", node_id)
      }
    }
    
    # Get node's chromosome count
    node_count <- ancestral_counts$node_counts[as.character(node_id)]
    
    # Skip if no count available
    if(is.na(node_count)) {
      next
    }
    
    # Get all edges where this node is the parent
    child_edges <- edge_data[edge_data$parent_node_id == node_id,]
    
    # Get count changes for these edges
    edge_count_changes <- count_changes_df[count_changes_df$parent_node_id == node_id,]
    
    # Skip if no count changes available
    if(nrow(edge_count_changes) == 0) {
      next
    }
    
    # Analyze change type patterns across children
    change_types <- na.omit(edge_count_changes$count_change_type)
    
    if(length(change_types) == 0) {
      next
    }
    
    # Check for consistent pattern across multiple children
    if(length(change_types) >= 2) {
      # Count increases (fissions)
      increase_types <- c("major_increase", "medium_increase")
      increase_count <- sum(change_types %in% increase_types)
      increase_ratio <- increase_count / length(change_types)
      
      # Count decreases (fusions)
      decrease_types <- c("major_decrease", "medium_decrease")
      decrease_count <- sum(change_types %in% decrease_types)
      decrease_ratio <- decrease_count / length(change_types)
      
      # Detect shared fission
      if(increase_ratio >= 0.4 && increase_count >= 2) {
        # Get the mean increase
        increase_diffs <- edge_count_changes$count_diff[edge_count_changes$count_change_type %in% increase_types]
        mean_increase <- mean(increase_diffs, na.rm = TRUE)
        
        # Create shared fission event
        events[[length(events) + 1]] <- list(
          edge_id = NA_integer_,
          parent_node = node_id,
          child_node = NA_integer_,
          parent_label = node_label,
          child_label = NA_character_,
          event_type = "fission",
          parent_count = node_count,
          child_count = mean(edge_count_changes$child_count[edge_count_changes$count_change_type %in% increase_types], na.rm = TRUE),
          count_diff = mean_increase,
          count_ratio = NA_real_,
          evidence_type = "consistent_child_pattern",
          confidence = if(increase_ratio >= 0.7) "high" else "medium",
          analysis_method = "top_down",
          affected_children = increase_count,
          total_children = length(change_types),
          prevalence = increase_ratio,
          shared_event = TRUE,
          group = NA_character_,
          chromosomes_A = NA_character_,
          chromosomes_B = NA_character_,
          reference_species = NA_character_
        )
      }
      
      # Detect shared fusion
      if(decrease_ratio >= 0.4 && decrease_count >= 2) {
        # Get the mean decrease
        decrease_diffs <- edge_count_changes$count_diff[edge_count_changes$count_change_type %in% decrease_types]
        mean_decrease <- mean(decrease_diffs, na.rm = TRUE)
        
        # Create shared fusion event
        events[[length(events) + 1]] <- list(
          edge_id = NA_integer_,
          parent_node = node_id,
          child_node = NA_integer_,
          parent_label = node_label,
          child_label = NA_character_,
          event_type = "fusion",
          parent_count = node_count,
          child_count = mean(edge_count_changes$child_count[edge_count_changes$count_change_type %in% decrease_types], na.rm = TRUE),
          count_diff = -mean_decrease,  # Make positive for clarity
          count_ratio = NA_real_,
          evidence_type = "consistent_child_pattern",
          confidence = if(decrease_ratio >= 0.7) "high" else "medium",
          analysis_method = "top_down",
          affected_children = decrease_count,
          total_children = length(change_types),
          prevalence = decrease_ratio,
          shared_event = TRUE,
          group = NA_character_,
          chromosomes_A = NA_character_,
          chromosomes_B = NA_character_,
          reference_species = NA_character_
        )
      }
    }
    
    # Also detect lineage-specific events
    for(i in 1:nrow(edge_count_changes)) {
      change_type <- edge_count_changes$count_change_type[i]
      
      # Skip stable changes
      if(is.na(change_type) || change_type == "stable") {
        next
      }
      
      child_node <- edge_count_changes$child_node_id[i]
      child_label <- edge_count_changes$child_label[i]
      child_count <- edge_count_changes$child_count[i]
      count_diff <- edge_count_changes$count_diff[i]
      
      # Determine event type from change type
      if(grepl("increase", change_type)) {
        event_type <- if(grepl("minor", change_type)) "minor_fission" else "fission"
        events[[length(events) + 1]] <- list(
          edge_id = edge_count_changes$edge_id[i],
          parent_node = node_id,
          child_node = child_node,
          parent_label = node_label,
          child_label = child_label,
          event_type = event_type,
          parent_count = node_count,
          child_count = child_count,
          count_diff = count_diff,
          count_ratio = edge_count_changes$count_ratio[i],
          evidence_type = "ancestral_count_difference",
          confidence = if(grepl("major", change_type)) "high" else if(grepl("medium", change_type)) "medium" else "low",
          analysis_method = "top_down",
          affected_children = NA_integer_,
          total_children = NA_integer_,
          prevalence = NA_real_,
          shared_event = FALSE,
          group = NA_character_,
          chromosomes_A = NA_character_,
          chromosomes_B = NA_character_,
          reference_species = NA_character_
        )
      } else if(grepl("decrease", change_type)) {
        event_type <- if(grepl("minor", change_type)) "minor_fusion" else "fusion"
        events[[length(events) + 1]] <- list(
          edge_id = edge_count_changes$edge_id[i],
          parent_node = node_id,
          child_node = child_node,
          parent_label = node_label,
          child_label = child_label,
          event_type = event_type,
          parent_count = node_count,
          child_count = child_count,
          count_diff = -count_diff,  # Make positive for clarity
          count_ratio = 1 / edge_count_changes$count_ratio[i],  # Invert for clarity
          evidence_type = "ancestral_count_difference",
          confidence = if(grepl("major", change_type)) "high" else if(grepl("medium", change_type)) "medium" else "low",
          analysis_method = "top_down",
          affected_children = NA_integer_,
          total_children = NA_integer_,
          prevalence = NA_real_,
          shared_event = FALSE,
          group = NA_character_,
          chromosomes_A = NA_character_,
          chromosomes_B = NA_character_,
          reference_species = NA_character_
        )
      }
    }
  }
  
  # Convert to data frame
  if(length(events) > 0) {
    # Create a data frame directly from the list of events
    events_df <- do.call(rbind, lapply(events, function(event) {
      as.data.frame(event, stringsAsFactors = FALSE)
    }))
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
      group = character(),
      chromosomes_A = character(),
      chromosomes_B = character(),
      reference_species = character(),
      stringsAsFactors = FALSE
    )
  }
  
  message(paste("Inferred", nrow(events_df), "chromosome events using top-down approach"))
  
  return(events_df)
}

# ==============================
# Combined Event Analysis
# ==============================

combine_event_analyses <- function(bottom_up_events, top_down_events) {
  # 合并自下而上和自上而下的事件分析结果，确保列结构兼容
  message("Combining bottom-up and top-down event analyses with enhanced integration...")
  
  # 首先检查我们是否有有效的数据框
  have_bottom_up <- !is.null(bottom_up_events) && inherits(bottom_up_events, "data.frame") && nrow(bottom_up_events) > 0
  have_top_down <- !is.null(top_down_events) && inherits(top_down_events, "data.frame") && nrow(top_down_events) > 0
  
  # 如果两者都为空，返回空数据框
  if(!have_bottom_up && !have_top_down) {
    message("No events to combine")
    return(data.frame(
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
      analysis_orig = character(),
      shared_event = logical(),
      affected_children = integer(),
      total_children = integer(),
      prevalence = numeric(),
      group = character(),
      chromosomes_A = character(),
      chromosomes_B = character(),
      reference_species = character(),
      stringsAsFactors = FALSE
    ))
  }
  
  # 创建一个函数来标准化数据框
  standardize_df <- function(df, ref_df = NULL) {
    # 标准化单个数据框的列
    if(nrow(df) == 0) return(df)
    
    # 如果没有提供参考数据框，添加一些基本列
    if(is.null(ref_df)) {
      if(!"shared_event" %in% names(df)) df$shared_event <- FALSE
      if(!"analysis_orig" %in% names(df)) df$analysis_orig <- NA_character_
      
      return(df)
    }
    
    # 所有可能的列（合并两个数据框的所有列）
    all_columns <- unique(c(names(df), names(ref_df)))
    
    # 对齐列：添加缺失的列并应用正确的类型
    for(col in all_columns) {
      if(!(col %in% names(df))) {
        # 如果列在df中缺失，添加具有适当类型的NA值列
        if(col %in% names(ref_df)) {
          # 根据ref_df中的类型确定
          col_class <- class(ref_df[[col]])
          
          if("integer" %in% col_class) {
            df[[col]] <- NA_integer_
          } else if("numeric" %in% col_class) {
            df[[col]] <- NA_real_
          } else if("logical" %in% col_class) {
            df[[col]] <- NA
          } else {
            df[[col]] <- NA_character_
          }
        } else {
          # 如果在ref_df中也缺失，根据列名判断
          if(grepl("id$|node$|count$|children$", col)) {
            df[[col]] <- NA_integer_
          } else if(grepl("ratio$|prevalence$", col)) {
            df[[col]] <- NA_real_
          } else if(grepl("^is_|event$", col)) {
            df[[col]] <- NA
          } else {
            df[[col]] <- NA_character_
          }
        }
      }
    }
    
    # 确保列的顺序匹配
    df <- df[, all_columns]
    
    return(df)
  }
  
  # 处理单个数据框情况
  if(have_bottom_up && !have_top_down) {
    message("Only bottom-up events available")
    bottom_up_events <- standardize_df(bottom_up_events)
    bottom_up_events$analysis_orig <- "bottom_up"
    return(bottom_up_events)
  }
  
  if(!have_bottom_up && have_top_down) {
    message("Only top-down events available")
    top_down_events <- standardize_df(top_down_events)
    top_down_events$analysis_orig <- "top_down"
    return(top_down_events)
  }
  
  # 如果两个数据框都存在，先给每个添加起源标记
  bottom_up_events$analysis_orig <- "bottom_up"
  top_down_events$analysis_orig <- "top_down"
  
  # 标准化双方的数据框
  message("Standardizing data frames before combining...")
  
  # 通过相互参考标准化
  bu_standardized <- standardize_df(bottom_up_events, top_down_events)
  td_standardized <- standardize_df(top_down_events, bu_standardized)
  
  # 重新检查一次以确保完全对齐
  all_columns <- unique(c(names(bu_standardized), names(td_standardized)))
  bu_standardized <- bu_standardized[, all_columns]
  td_standardized <- td_standardized[, all_columns]
  
  # 输出列信息以进行调试
  message(paste("Bottom-up events data frame has", ncol(bu_standardized), "columns"))
  message(paste("Top-down events data frame has", ncol(td_standardized), "columns"))
  
  # 使用rbind合并标准化的数据框
  if(ncol(bu_standardized) != ncol(td_standardized)) {
    # 最后的安全措施 - 使用更可靠的data.table方法
    if(requireNamespace("data.table", quietly = TRUE)) {
      message("Using data.table::rbindlist for safer merging")
      all_events <- data.table::rbindlist(list(bu_standardized, td_standardized), fill = TRUE)
      all_events <- as.data.frame(all_events)
    } else {
      stop(paste("Cannot combine data frames with different column counts:",
                "Bottom-up has", ncol(bu_standardized), "columns,",
                "Top-down has", ncol(td_standardized), "columns"))
    }
  } else {
    # 标准rbind应该可以安全工作
    all_events <- rbind(bu_standardized, td_standardized)
  }
  
  # 初始化佐证和趋同列
  all_events$corroborated <- FALSE
  if(!"convergent_evolution" %in% names(all_events)) all_events$convergent_evolution <- FALSE
  if(!"convergent_group" %in% names(all_events)) all_events$convergent_group <- NA_character_
  
  # 寻找两种方法都检测到的事件
  for(i in 1:nrow(bu_standardized)) {
    bu_event <- bu_standardized[i,]
    
    # 跳过无边界ID的事件
    if(is.na(bu_event$edge_id)) {
      next
    }
    
    # 查找匹配的自上而下事件
    matching_td <- td_standardized[
      !is.na(td_standardized$edge_id) &
      td_standardized$edge_id == bu_event$edge_id &
      td_standardized$event_type == bu_event$event_type,
    ]
    
    if(nrow(matching_td) > 0) {
      # 事件被两种方法都检测到 - 增加置信度
      all_events$confidence[all_events$edge_id == bu_event$edge_id & 
                        all_events$event_type == bu_event$event_type] <- "high"
      
      # 添加佐证标志
      all_events$corroborated[all_events$edge_id == bu_event$edge_id & 
                          all_events$event_type == bu_event$event_type] <- TRUE
    }
  }
  
  # 检测趋同进化模式
  if(requireNamespace("stats", quietly = TRUE) && nrow(all_events) >= 5) {
    tryCatch({
      # 分析每种主要事件类型
      for(evt_type in c("fusion", "fission")) {
        # 过滤这种类型的事件
        type_events <- all_events[all_events$event_type == evt_type & !is.na(all_events$count_diff),]
        
        if(nrow(type_events) >= 3) {
          # 按计数差异聚类事件
          # 准备数据用于聚类
          event_diffs <- type_events$count_diff
          names(event_diffs) <- 1:nrow(type_events)
          
          # 层次聚类
          clust_result <- stats::hclust(stats::dist(event_diffs), method="complete")
          
          # 剪裁树得到聚类
          k <- min(3, nrow(type_events) - 1)
          clusters <- stats::cutree(clust_result, k = k)
          
          # 检查每个聚类
          for(cluster_id in unique(clusters)) {
            # 获取此聚类中的事件
            cluster_events <- type_events[as.integer(names(clusters[clusters == cluster_id])),]
            
            # 如果多个事件在不同物种中有相似变化
            if(nrow(cluster_events) >= 2) {
              species_set <- unique(cluster_events$child_label)
              
              if(length(species_set) >= 2) {
                # 标记为潜在的趋同进化
                group_id <- paste0("convergent_", evt_type, "_", cluster_id)
                
                # 更新convergent_evolution和convergent_group列
                for(i in 1:nrow(all_events)) {
                  if(all_events$child_label[i] %in% species_set && 
                     all_events$event_type[i] == evt_type) {
                    all_events$convergent_evolution[i] <- TRUE
                    all_events$convergent_group[i] <- group_id
                  }
                }
              }
            }
          }
        }
      }
    }, error = function(e) {
      message("Convergent evolution detection failed: ", e$message)
    })
  }
  
  # 汇总趋同进化发现
  convergent_count <- sum(all_events$convergent_evolution, na.rm=TRUE)
  if(convergent_count > 0) {
    message(paste("Detected", convergent_count, "events potentially involved in convergent evolution"))
    
    # 显示每个趋同组
    convergent_groups <- unique(all_events$convergent_group[!is.na(all_events$convergent_group)])
    for(group in convergent_groups) {
      group_events <- all_events[all_events$convergent_group == group,]
      species_involved <- unique(group_events$child_label)
      message(paste(" -", group, ":", length(species_involved), "species:", 
                  paste(species_involved[1:min(3, length(species_involved))], collapse=", "),
                  ifelse(length(species_involved) > 3, "...", "")))
    }
  }
  
  message(paste("Combined event analyses with a total of", nrow(all_events), "events"))
  
  return(all_events)
}

identify_shared_events <- function(tree, events_df, clades) {
  # Identify shared chromosome events across branches with dynamic thresholds
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
  
  # Find events already identified as shared by top-down analysis
  shared_from_top_down <- events_df[events_df$shared_event == TRUE,]
  
  if(nrow(shared_from_top_down) > 0) {
    for(i in 1:nrow(shared_from_top_down)) {
      event <- shared_from_top_down[i,]
      
      # Find corresponding clade
      clade <- clades[[as.character(event$parent_node)]]
      
      if(is.null(clade)) {
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
    }
  }
  
  # Get unique event types
  event_types <- unique(events_df$event_type)
  
  # For each event type and clade, check using dynamic thresholds
  for(event_type in event_types) {
    # Filter events of this type (exclude minor events for sharing analysis)
    if(event_type %in% c("minor_fusion", "minor_fission", "stable")) {
      next
    }
    
    type_events <- events_df[events_df$event_type == event_type & !events_df$shared_event,]
    
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
      
      # Calculate event prevalence
      prevalence <- event_count / clade$tip_count
      
      # Dynamic threshold based on clade size
      threshold <- adjust_threshold(clade$tip_count)
      
      # Determine if event is shared across the clade
      if(prevalence >= threshold) {
        # Calculate sharing status and confidence
        if(prevalence >= 0.8) {
          shared_status <- "strongly_shared"
          confidence <- "high"
        } else if(prevalence >= 0.5) {
          shared_status <- "likely_shared"
          confidence <- "medium"
        } else {
          shared_status <- "weakly_shared"
          confidence <- "low"
        }
        
        # Statistical validation using binomial test
        # Assume background rate of this event type is 10%
        if(requireNamespace("stats", quietly = TRUE)) {
          # Count total events of this type
          total_type_events <- nrow(type_events)
          total_tips <- length(tree$tip.label)
          background_rate <- max(0.1, total_type_events / total_tips)
          
          # Binomial test
          p_value <- tryCatch({
            binom.test(event_count, clade$tip_count, p=background_rate, alternative="greater")$p.value
          }, error = function(e) NA)
          
          # Adjust confidence based on p-value
          if(!is.na(p_value)) {
            if(p_value < 0.01) {
              confidence <- "high"
            } else if(p_value < 0.05) {
              confidence <- "medium"
            } else if(p_value >= 0.1) {
              confidence <- "low"
            }
          }
        }
        
        # Check if there's a matching top-down event for this clade
        corroborated <- FALSE
        if(nrow(shared_from_top_down) > 0) {
          matching_td <- shared_from_top_down[
            shared_from_top_down$parent_node == as.integer(clade_id) &
            shared_from_top_down$event_type == event_type,
          ]
          
          if(nrow(matching_td) > 0) {
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
      }
    }
  }
  
  # Convert to data frame
  if(length(shared_events) > 0) {
    shared_df <- do.call(rbind, lapply(shared_events, function(e) {
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
    
    # Store tip events separately because they don't fit in a data frame
    attr(shared_df, "tip_events") <- lapply(shared_events, function(e) {
      if(!is.null(e$tip_events)) e$tip_events[[1]] else NULL
    })
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
      analysis_type = character(),
      corroborated = logical(),
      stringsAsFactors = FALSE
    )
    attr(shared_df, "tip_events") <- list()
  }
  
  message(paste("Identified", nrow(shared_df), "shared chromosome events"))
  
  return(shared_df)
}

# Helper function for dynamic threshold based on clade size
adjust_threshold <- function(clade_size) {
  # Small clades need higher threshold, large clades can use lower threshold
  if(clade_size < 5) return(0.9)  # Very small clades require near-universal sharing
  if(clade_size < 10) return(0.8) # Small clades
  if(clade_size < 30) return(0.6) # Medium clades
  if(clade_size < 100) return(0.5) # Large clades
  return(0.4) # Very large clades
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

# ==============================
# Genomic Evidence Integration
# ==============================

integrate_genomic_evidence <- function(events_df, genomic_data) {
  # Integrate genomic evidence with chromosome events if available
  if(is.null(genomic_data)) {
    message("No genomic data provided, skipping genomic integration")
    return(events_df)
  }
  
  message("Integrating genomic evidence with chromosome events...")
  
  # Ensure we have events to integrate
  if(nrow(events_df) == 0) {
    message("No events to integrate with genomic data")
    return(events_df)
  }
  
  # Add genomic support columns
  events_df$genomic_support <- "none"
  events_df$support_type <- NA_character_
  events_df$support_details <- NA_character_
  
  # Map event types to genomic evidence types
  event_to_evidence_map <- list(
    "fusion" = c("fusion", "robertsonian", "centromere_shift"),
    "fission" = c("fission", "chromosome_break", "neo_centromere"),
    "minor_fusion" = c("fusion", "robertsonian"),
    "minor_fission" = c("fission", "chromosome_break")
  )
  
  # For each event involving tip species, check for genomic evidence
  for(i in 1:nrow(events_df)) {
    if(!is.na(events_df$child_label[i]) && events_df$child_is_tip[i]) {
      species <- events_df$child_label[i]
      event_type <- events_df$event_type[i]
      
      # Skip if species not in genomic data
      if(!(species %in% genomic_data$species)) {
        next
      }
      
      # Get evidence types to look for
      evidence_types <- event_to_evidence_map[[event_type]]
      if(is.null(evidence_types)) {
        next
      }
      
      # Extract relevant genomic evidence
      sp_evidence <- genomic_data[genomic_data$species == species & 
                                  genomic_data$evidence_type %in% evidence_types,]
      
      if(nrow(sp_evidence) > 0) {
        # Found supporting genomic evidence
        events_df$genomic_support[i] <- "confirmed"
        events_df$support_type[i] <- paste(unique(sp_evidence$evidence_type), collapse=",")
        events_df$support_details[i] <- paste(sp_evidence$details[1:min(3, nrow(sp_evidence))], collapse="; ")
        
        # Upgrade confidence based on genomic evidence
        events_df$confidence[i] <- "high"
      }
    }
  }
  
  # Summarize genomic integration
  support_count <- sum(events_df$genomic_support == "confirmed", na.rm=TRUE)
  message(paste("Added genomic support for", support_count, "events"))
  
  return(events_df)
}

# ==============================
# Main Event Detection Function
# ==============================

detect_shared_chromosome_events <- function(integrated_data, conserved_groups, anc_method = "ML", genomic_data = NULL) {
  # Main function for chromosome event detection using combined approaches
  message("Starting shared chromosome event detection with combined approaches...")
  
  # Get tree and chromosome counts
  tree <- integrated_data$tree
  chr_counts <- integrated_data$chr_counts
  
  # Extract tree structure information
  edge_data <- extract_tree_edges(tree)
  clades <- identify_clades(tree)
  
  # Reconstruct ancestral chromosome counts using specified method
  ancestral_counts <- reconstruct_ancestral_chromosome_counts(tree, chr_counts, method = anc_method)
  
  # Detect chromosome count changes
  count_changes <- detect_count_changes(tree, edge_data, chr_counts, ancestral_counts)
  
  # Detect mapping patterns
  pattern_data <- detect_mapping_patterns(edge_data, integrated_data, conserved_groups)
  
  # Infer events using bottom-up approach (from tips to internal nodes)
  bottom_up_events <- infer_events_bottom_up(tree, edge_data, count_changes, pattern_data, chr_counts)
  
  # Infer events using top-down approach (from root to tips)
  top_down_events <- infer_events_top_down(tree, edge_data, count_changes, ancestral_counts, clades)
  
  # Combine event analyses
  events_df <- combine_event_analyses(bottom_up_events, top_down_events)
  
  # Integrate genomic evidence if available
  if(!is.null(genomic_data)) {
    events_df <- integrate_genomic_evidence(events_df, genomic_data)
  }
  
  # Identify shared events
  shared_events_df <- identify_shared_events(tree, events_df, clades)
  
  # Return results
  result <- list(
    tree = tree,
    edge_data = edge_data,
    clades = clades,
    ancestral_counts = ancestral_counts,
    count_changes = count_changes$count_changes,
    pattern_data = pattern_data,
    bottom_up_events = bottom_up_events,
    top_down_events = top_down_events,
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
  
  # Ancestral counts
  ancestral_counts_file <- file.path(output_dir, "ancestral_chromosome_counts.csv")
  ancestral_counts_df <- data.frame(
    node_id = as.integer(names(event_result$ancestral_counts$node_counts)),
    count = event_result$ancestral_counts$node_counts,
    stringsAsFactors = FALSE
  )
  write.csv(ancestral_counts_df, file = ancestral_counts_file, row.names = FALSE)
  
  # Count changes
  count_changes_file <- file.path(output_dir, "chromosome_count_changes.csv")
  write.csv(event_result$count_changes, file = count_changes_file, row.names = FALSE)
  
  # Mapping patterns
  pattern_data_file <- file.path(output_dir, "mapping_patterns.csv")
  write.csv(event_result$pattern_data, file = pattern_data_file, row.names = FALSE)
  
  # Bottom-up events
  bottom_up_file <- file.path(output_dir, "bottom_up_events.csv")
  write.csv(event_result$bottom_up_events, file = bottom_up_file, row.names = FALSE)
  
  # Top-down events
  top_down_file <- file.path(output_dir, "top_down_events.csv")
  write.csv(event_result$top_down_events, file = top_down_file, row.names = FALSE)
  
  # Combined events
  events_file <- file.path(output_dir, "inferred_events.csv")
  write.csv(event_result$events, file = events_file, row.names = FALSE)
  
  # Shared events
  shared_events_file <- file.path(output_dir, "shared_events.csv")
  write.csv(event_result$shared_events, file = shared_events_file, row.names = FALSE)
  
  # Create JSON files for visualization if package available
  if(requireNamespace("jsonlite", quietly = TRUE)) {
    # Events JSON
    events_json_file <- file.path(output_dir, "events.json")
    jsonlite::write_json(event_result$events, events_json_file, pretty=TRUE)
    
    # Shared events JSON
    shared_json_file <- file.path(output_dir, "shared_events.json")
    jsonlite::write_json(event_result$shared_events, shared_json_file, pretty=TRUE)
  }
  
  # Create detailed events report
  events_report <- file.path(output_dir, "chromosome_events_report.txt")
  
  sink(events_report)
  cat("Enhanced Chromosome Evolutionary Events Report\n")
  cat("===========================================\n\n")
  
  cat("1. Ancestral Chromosome Count Summary\n")
  cat("----------------------------------\n")
  tip_counts <- event_result$ancestral_counts$node_counts[1:length(event_result$tree$tip.label)]
  internal_counts <- event_result$ancestral_counts$node_counts[(length(tip_counts) + 1):length(event_result$ancestral_counts$node_counts)]
  
  cat(paste("Extant species chromosome count range:", min(tip_counts, na.rm=TRUE), 
            "to", max(tip_counts, na.rm=TRUE), "\n"))
  cat(paste("Ancestral node chromosome count range:", min(internal_counts, na.rm=TRUE), 
            "to", max(internal_counts, na.rm=TRUE), "\n"))
  cat(paste("Root node estimated chromosome count:", 
            event_result$ancestral_counts$node_counts[length(tip_counts) + 1], "\n\n"))
  cat(paste("Ancestral state reconstruction method:", 
            event_result$ancestral_counts$method, "\n\n"))
  
  cat("2. Summary of Detected Events\n")
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
    
    # Report corroborated events
    corroborated_count <- sum(event_result$events$corroborated, na.rm = TRUE)
    cat(paste("\nEvents corroborated by both methods:", corroborated_count, "\n"))
    
    # Report convergent evolution if detected
    convergent_count <- sum(event_result$events$convergent_evolution, na.rm = TRUE)
    if(convergent_count > 0) {
      cat(paste("\nEvents involved in potential convergent evolution:", convergent_count, "\n"))
      
      convergent_groups <- unique(event_result$events$convergent_group[!is.na(event_result$events$convergent_group)])
      cat(paste("Found", length(convergent_groups), "convergent evolution patterns\n"))
    }
    
    # Report genomic evidence if available
    if("genomic_support" %in% names(event_result$events)) {
      genomic_support_count <- sum(event_result$events$genomic_support == "confirmed", na.rm = TRUE)
      cat(paste("\nEvents with genomic evidence support:", genomic_support_count, "\n"))
    }
  } else {
    cat("  No chromosome events detected\n")
  }
  cat("\n")
  
  cat("3. Shared Chromosome Events\n")
  cat("-------------------------\n")
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
      cat("\n")
    }
  } else {
    cat("  No shared chromosome events detected\n")
  }
  cat("\n")
  
  cat("4. Detailed Event Listing\n")
  cat("----------------------\n")
  if(nrow(event_result$events) > 0) {
    high_conf_events <- event_result$events[event_result$events$confidence == "high",]
    
    if(nrow(high_conf_events) > 0) {
      cat("High confidence events:\n")
      for(i in 1:min(10, nrow(high_conf_events))) {
        evt <- high_conf_events[i,]
        cat(paste0("  [", i, "] ", evt$event_type, " in ", evt$child_label, 
                   " (", evt$analysis_method, ", ", evt$evidence_type, ")\n"))
        
        if(!is.na(evt$evidence_type) && evt$evidence_type == "mapping_pattern") {
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
        }
        
        # Report corroboration
        if(!is.null(evt$corroborated) && !is.na(evt$corroborated) && evt$corroborated) {
          cat("      - Corroborated by both analysis methods\n")
        }
        
        # Report convergent evolution
        if(!is.null(evt$convergent_evolution) && !is.na(evt$convergent_evolution) && evt$convergent_evolution) {
          cat(paste0("      - Part of convergent evolution group: ", evt$convergent_group, "\n"))
        }
        
        # Report genomic support if available
        if(!is.null(evt$genomic_support) && !is.na(evt$genomic_support) && evt$genomic_support == "confirmed") {
          cat(paste0("      - Genomic evidence: ", evt$support_type, "\n"))
          if(!is.na(evt$support_details)) cat(paste0("      - Details: ", evt$support_details, "\n"))
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
    cat("Usage: Rscript ancestral_chromosomes_phase3_optimized.R integrated_data_file conserved_groups_file output_dir [anc_method] [genomic_data_file]\n")
    cat("Example: Rscript ancestral_chromosomes_phase3_optimized.R results/integrated_data.RData results/phase2/conserved_chromosome_groups.RData results/phase3/ ML\n")
    cat("Options:\n")
    cat("  integrated_data_file: RData file with integrated mapping data from phase 1\n")
    cat("  conserved_groups_file: RData file with conserved chromosome groups from phase 2\n")
    cat("  output_dir: Directory for saving results\n")
    cat("  anc_method: Method for ancestral state reconstruction (ML, REML, bayesian) [default: ML]\n")
    cat("  genomic_data_file: Optional RData file with genomic evidence for event validation\n")
    return(1)
  }
  
  # Parse arguments
  integrated_data_file <- args[1]
  conserved_groups_file <- args[2]
  output_dir <- args[3]
  
  # Optional arguments
  anc_method <- if(length(args) >= 4) args[4] else "ML"
  genomic_data_file <- if(length(args) >= 5) args[5] else NULL
  
  # Load data
  integrated_data <- load_integrated_data(integrated_data_file)
  conserved_groups <- load_conserved_groups(conserved_groups_file)
  
  # Load genomic data if available
  genomic_data <- load_genomic_data(genomic_data_file)
  
  # Detect shared chromosome events
  event_result <- detect_shared_chromosome_events(integrated_data, conserved_groups, anc_method, genomic_data)
  
  # Save results
  output_files <- save_event_detection_results(event_result, output_dir)
  
  # Print basic summary
  cat("\n===== Enhanced Chromosome Evolutionary Events Summary =====\n")
  
  # Summarize ancestral reconstruction
  cat("Ancestral chromosome counts:\n")
  cat(paste("  Root node count:", round(event_result$ancestral_counts$node_counts[length(event_result$tree$tip.label) + 1]), "\n"))
  cat(paste("  Method used:", event_result$ancestral_counts$method, "\n"))
  
  # Summarize events
  if(nrow(event_result$events) > 0) {
    # Summarize events by type
    event_types <- table(event_result$events$event_type)
    cat("\nDetected events by type:\n")
    for(evt_type in names(event_types)) {
      cat(paste0("  - ", evt_type, ": ", event_types[evt_type], "\n"))
    }
    
    # Summarize by analysis method
    method_types <- table(event_result$events$analysis_method)
    cat("\nEvents by analysis method:\n")
    for(method in names(method_types)) {
      cat(paste0("  - ", method, ": ", method_types[method], "\n"))
    }
    
    # Corroborated events
    corroborated_count <- sum(event_result$events$corroborated, na.rm = TRUE)
    cat(paste("\nEvents corroborated by both methods:", corroborated_count, "\n"))
    
    # Convergent evolution
    convergent_count <- sum(event_result$events$convergent_evolution, na.rm = TRUE)
    if(convergent_count > 0) {
      cat(paste("\nEvents in convergent evolution patterns:", convergent_count, "\n"))
    }
    
    # Genomic support
    if("genomic_support" %in% names(event_result$events)) {
      genomic_support_count <- sum(event_result$events$genomic_support == "confirmed", na.rm = TRUE)
      if(genomic_support_count > 0) {
        cat(paste("\nEvents with genomic evidence support:", genomic_support_count, "\n"))
      }
    }
  } else {
    cat("No chromosome events detected\n")
  }
  
  # Summarize shared events
  if(nrow(event_result$shared_events) > 0) {
    cat("\nShared events by type:\n")
    shared_by_type <- table(event_result$shared_events$event_type)
    for(evt_type in names(shared_by_type)) {
      cat(paste0("  - ", evt_type, ": ", shared_by_type[evt_type], "\n"))
    }
    
    # By analysis method
    shared_by_method <- table(event_result$shared_events$analysis_type)
    cat("\nShared events by analysis method:\n")
    for(method in names(shared_by_method)) {
      cat(paste0("  - ", method, ": ", shared_by_method[method], "\n"))
    }
    
    # Corroborated shared events
    corroborated_shared <- sum(event_result$shared_events$corroborated, na.rm = TRUE)
    cat(paste("\nShared events corroborated by both methods:", corroborated_shared, "\n"))
    
    # List top shared events
    cat("\nTop shared events by prevalence:\n")
    top_shared <- event_result$shared_events[order(-event_result$shared_events$prevalence),]
    for(i in 1:min(3, nrow(top_shared))) {
      se <- top_shared[i,]
      cat(paste0("  [", i, "] ", se$event_type, " in clade ", se$clade_label,
                 " (", round(se$prevalence * 100), "% prevalence, ", 
                 ifelse(!is.null(se$corroborated) && !is.na(se$corroborated) && se$corroborated, "corroborated", "single method"), ")\n"))
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
