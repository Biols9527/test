#!/usr/bin/env Rscript

# ==============================
# Ancestral Chromosome Reconstruction Phase 5: Synteny Data Processing for Invertebrates
# Date: 2025-03-18
# Author: Biols9527
# Description: Processes pairwise synteny data with specific optimizations for invertebrate genomes
# ==============================

# Required packages
suppressPackageStartupMessages({
  library(data.table)  # For efficient data handling
  library(ape)         # For tree manipulation
  library(igraph)      # For network analysis
  library(ggplot2)     # For visualization
  library(reshape2)    # For data reshaping
})

# ==============================
# Invertebrate-Optimized Configuration
# ==============================

load_config <- function(config_file) {
  # Default configuration optimized for invertebrate genomes
  default_config <- list(
    # 无脊椎动物优化的共线性处理参数
    min_block_size = 2,                 # 更低的阈值，捕获更小的保守区块 (脊椎动物通常为3-5)
    min_block_significance = 0.1,       # 放宽显著性阈值
    max_gap = 15,                       # 允许更大的间隙，适应重排频繁的基因组
    normalize_coordinates = TRUE,       # 坐标标准化
    
    # 针对频繁重排的过滤参数
    min_species_coverage = 0.2,         # 降低物种覆盖要求
    max_distance_cutoff = 1.0,          # 增大进化距离阈值，因为即使较远的物种也可能保留微共线性
    
    # 事件检测参数调整
    fusion_detection_threshold = 0.5,   # 降低检测阈值，适应更复杂的融合模式
    fission_detection_threshold = 0.5,  # 降低检测阈值，适应更复杂的裂变模式
    translocation_threshold = 0.3,      # 降低转位事件阈值，更敏感地检测重排
    inversion_threshold = 0.4,          # 检测倒位的阈值
    
    # 无脊椎动物特定参数
    detect_micro_synteny = TRUE,        # 专门检测微共线性区块
    micro_synteny_min_size = 2,         # 微共线性最小基因数
    micro_synteny_max_gap = 5,          # 微共线性允许的最大间隙
    
    # 重排分析参数
    detect_rearrangement_hotspots = TRUE,  # 检测重排热点
    min_rearrangements_for_hotspot = 2,    # 判定为热点的最小重排数
    
    # 多层次分析
    hierarchical_synteny = TRUE,        # 分层次分析共线性(宏共线性和微共线性)
    
    # 可视化选项
    generate_visualizations = TRUE,
    max_dotplots = 30,                  # 增加点图数量以覆盖更多物种对
    
    # 日志记录选项
    verbose = TRUE,                     # 是否输出详细日志
    log_level = "INFO",                 # 日志等级: DEBUG, INFO, WARNING, ERROR
    
    # 物种名映射
    species_name_mapping = list()
  )
  
  # 加载自定义配置
  if(!is.null(config_file) && file.exists(config_file)) {
    custom_config <- tryCatch({
      readRDS(config_file)
    }, error = function(e) {
      log_message(paste("加载配置文件时出错:", e$message), "ERROR")
      return(NULL)
    })
    
    if(!is.null(custom_config)) {
      for(param in names(custom_config)) {
        default_config[[param]] <- custom_config[[param]]
      }
      log_message("已加载自定义配置", "INFO")
    }
  }
  
  return(default_config)
}

# ==============================
# 日志记录函数
# ==============================

log_message <- function(message, level = "INFO", config = NULL) {
  # 如果提供了config并且verbose为FALSE，则只记录ERROR级别的消息
  if(!is.null(config) && !is.null(config$verbose) && !config$verbose && level != "ERROR") {
    return(invisible(NULL))
  }
  
  # 获取config中的log_level，如果没有则默认为"INFO"
  log_level <- if(!is.null(config) && !is.null(config$log_level)) config$log_level else "INFO"
  
  # 日志级别优先级: DEBUG < INFO < WARNING < ERROR
  level_priority <- list("DEBUG" = 1, "INFO" = 2, "WARNING" = 3, "ERROR" = 4)
  
  # 仅当消息级别高于或等于配置的级别时才输出
  if(is.null(level_priority[[level]]) || is.null(level_priority[[log_level]]) || 
     level_priority[[level]] < level_priority[[log_level]]) {
    return(invisible(NULL))
  }
  
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  formatted_message <- paste0("[", level, " ", timestamp, "] ", message)
  
  if(level == "ERROR") {
    cat(formatted_message, "\n", file = stderr())
  } else {
    cat(formatted_message, "\n")
  }
  
  return(invisible(NULL))
}

# ==============================
# 数据加载函数
# ==============================

load_synteny_files <- function(synteny_dir, tree_file, config) {
  log_message("加载共线性数据文件...", "INFO", config)
  
  # 读取系统发育树
  tree <- tryCatch({
    read.tree(tree_file)
  }, error = function(e) {
    log_message(paste("无法加载系统发育树:", tree_file, "-", e$message), "ERROR", config)
    return(NULL)
  })
  
  if(is.null(tree)) {
    stop("无法加载系统发育树，处理终止")
  }
  
  # 获取树中的所有物种名
  valid_species <- tree$tip.label
  log_message(paste("系统发育树中包含", length(valid_species), "个物种"), "INFO", config)
  
  # 检查目录是否存在
  if(!dir.exists(synteny_dir)) {
    log_message(paste("指定的共线性数据目录不存在:", synteny_dir), "ERROR", config)
    stop("共线性数据目录不存在")
  }
  
  # 列出所有符合模式的共线性文件
  synteny_files <- list.files(synteny_dir, pattern = "\\.tsv$", full.names = TRUE)
  
  if(length(synteny_files) == 0) {
    log_message(paste("在目录中未找到共线性数据文件:", synteny_dir), "ERROR", config)
    stop("未找到共线性数据文件")
  }
  
  log_message(paste("找到", length(synteny_files), "个潜在的共线性数据文件"), "INFO", config)
  
  # 处理每个文件
  synteny_list <- list()
  processed_count <- 0
  skipped_count <- 0
  
  for(file in synteny_files) {
    # 从文件名提取物种对
    filename <- basename(file)
    log_message(paste("处理文件:", filename), "DEBUG", config)
    
    species_pair <- sub("\\.tsv$", "", filename)
    
    # 检查是否包含"-"分隔符
    if(!grepl("-", species_pair)) {
      log_message(paste("跳过文件: 名称格式不符合'物种1-物种2'模式:", filename), "WARNING", config)
      skipped_count <- skipped_count + 1
      next
    }
    
    # 分割物种名
    species <- strsplit(species_pair, "-")[[1]]
    
    if(length(species) != 2) {
      log_message(paste("跳过文件: 无法从名称中提取两个物种:", filename), "WARNING", config)
      skipped_count <- skipped_count + 1
      next
    }
    
    species1 <- species[1]
    species2 <- species[2]
    
    # 应用物种名映射（如果存在）
    if(length(config$species_name_mapping) > 0) {
      if(species1 %in% names(config$species_name_mapping)) {
        log_message(paste("映射物种名:", species1, "->", config$species_name_mapping[[species1]]), "DEBUG", config)
        species1 <- config$species_name_mapping[[species1]]
      }
      if(species2 %in% names(config$species_name_mapping)) {
        log_message(paste("映射物种名:", species2, "->", config$species_name_mapping[[species2]]), "DEBUG", config)
        species2 <- config$species_name_mapping[[species2]]
      }
    }
    
    # 检查物种名是否在进化树中
    if(!(species1 %in% valid_species)) {
      log_message(paste("跳过文件: 物种名不在进化树中:", species1), "WARNING", config)
      skipped_count <- skipped_count + 1
      next
    }
    if(!(species2 %in% valid_species)) {
      log_message(paste("跳过文件: 物种名不在进化树中:", species2), "WARNING", config)
      skipped_count <- skipped_count + 1
      next
    }
    
    # 加载数据
    data <- tryCatch({
      df <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      log_message(paste("文件", filename, "包含", nrow(df), "行数据"), "DEBUG", config)
      df
    }, error = function(e) {
      log_message(paste("加载文件时出错:", file, "-", e$message), "ERROR", config)
      return(NULL)
    })
    
    if(is.null(data)) {
      skipped_count <- skipped_count + 1
      next
    }
    
    if(nrow(data) == 0) {
      log_message(paste("跳过空文件:", filename), "WARNING", config)
      skipped_count <- skipped_count + 1
      next
    }
    
    # 检查所需列
    required_cols <- c("x", "y", "chr1", "chr2")
    if(!all(required_cols %in% colnames(data))) {
      missing_cols <- required_cols[!required_cols %in% colnames(data)]
      log_message(paste("文件缺少所需列:", paste(missing_cols, collapse=", "), "-", filename), "WARNING", config)
      skipped_count <- skipped_count + 1
      next
    }
    
    # 对列进行初步检查
    if(any(is.na(data$x)) || any(is.na(data$y)) || any(is.na(data$chr1)) || any(is.na(data$chr2))) {
      log_message(paste("文件包含NA值:", filename, 
                       "- NA值计数: x=", sum(is.na(data$x)), 
                       "y=", sum(is.na(data$y)), 
                       "chr1=", sum(is.na(data$chr1)), 
                       "chr2=", sum(is.na(data$chr2))), "WARNING", config)
      
      # 移除包含NA值的行
      data <- data[!is.na(data$x) & !is.na(data$y) & !is.na(data$chr1) & !is.na(data$chr2), ]
      
      if(nrow(data) == 0) {
        log_message(paste("删除NA值后文件为空:", filename), "WARNING", config)
        skipped_count <- skipped_count + 1
        next
      }
      
      log_message(paste("移除NA值后剩余", nrow(data), "行数据"), "DEBUG", config)
    }
    
    # 应用配置中的过滤
    if(config$min_block_size > 1) {
      # 按染色体对分组
      chr_pairs <- paste0(data$chr1, "_", data$chr2)
      blocks <- split(data, chr_pairs)
      
      log_message(paste("文件包含", length(blocks), "个潜在染色体对"), "DEBUG", config)
      
      # 仅保留足够大小的块
      before_count <- length(blocks)
      blocks <- blocks[sapply(blocks, nrow) >= config$min_block_size]
      after_count <- length(blocks)
      
      if(after_count < before_count) {
        log_message(paste("过滤掉", before_count - after_count, "个小于min_block_size的染色体对"), "DEBUG", config)
      }
      
      if(length(blocks) == 0) {
        log_message(paste("文件中没有满足最小大小要求的共线性块:", filename), "WARNING", config)
        skipped_count <- skipped_count + 1
        next
      }
      
      # 重组块
      data <- do.call(rbind, blocks)
    }
    
    # 计算进化距离
    distance <- NA
    tryCatch({
      distance <- cophenetic.phylo(tree)[species1, species2]
      log_message(paste("物种", species1, "和", species2, "之间的进化距离:", distance), "DEBUG", config)
    }, error = function(e) {
      log_message(paste("无法计算物种间的进化距离:", species1, "和", species2, "-", e$message), "WARNING", config)
    })
    
    # 添加到数据集
    synteny_list[[species_pair]] <- list(
      species1 = species1,
      species2 = species2,
      data = data,
      file = file,
      distance = distance,
      
      # 计算基本统计
      gene_count = nrow(data),
      chr1_count = length(unique(data$chr1)),
      chr2_count = length(unique(data$chr2)),
      chr_pairs = length(unique(paste0(data$chr1, "_", data$chr2)))
    )
    
    processed_count <- processed_count + 1
    
    log_message(paste("处理文件:", filename, "-", 
                    nrow(data), "个共线性基因匹配,",
                    length(unique(data$chr1)), "个染色体(物种1),",
                    length(unique(data$chr2)), "个染色体(物种2)"), "INFO", config)
  }
  
  log_message(paste("成功加载了", processed_count, "个共线性数据集,", 
                  "跳过了", skipped_count, "个文件"), "INFO", config)
  
  # 如果没有加载任何文件，发出警告
  if(processed_count == 0) {
    log_message("未能成功加载任何共线性数据文件!", "ERROR", config)
    stop("未能成功加载任何共线性数据文件!")
  }
  
  # 添加系统树信息
  result <- list(
    synteny_data = synteny_list,
    tree = tree,
    valid_species = valid_species,
    file_count = processed_count,
    skipped_count = skipped_count
  )
  
  return(result)
}

# ==============================
# 共线性数据验证和可视化
# ==============================

validate_synteny_quality <- function(synteny_results, config) {
  log_message("验证共线性数据质量...", "INFO", config)
  
  # 初始化结果
  quality_results <- list(
    overall_quality = "good",
    issues = character(),
    statistics = list()
  )
  
  if(length(synteny_results$synteny_data) == 0) {
    quality_results$overall_quality <- "poor"
    quality_results$issues <- c(quality_results$issues, "未加载任何共线性数据")
    log_message("数据质量检查失败:未加载任何共线性数据", "ERROR", config)
    return(quality_results)
  }
  
  # 计算总体统计
  total_genes <- sum(sapply(synteny_results$synteny_data, function(x) x$gene_count))
  total_species <- length(unique(c(
    sapply(synteny_results$synteny_data, function(x) x$species1),
    sapply(synteny_results$synteny_data, function(x) x$species2)
  )))
  
  quality_results$statistics$total_genes <- total_genes
  quality_results$statistics$total_species <- total_species
  quality_results$statistics$total_pairs <- length(synteny_results$synteny_data)
  
  # 检查覆盖率 - 需要的物种配对数量
  required_pairs <- total_species * (total_species - 1) / 2
  coverage <- length(synteny_results$synteny_data) / required_pairs
  quality_results$statistics$coverage <- coverage
  
  if(coverage < 0.3) {
    quality_results$overall_quality <- "warning"
    quality_results$issues <- c(quality_results$issues, "物种配对覆盖率低")
    log_message(paste("警告: 物种配对覆盖率低 (", round(coverage*100, 1), "%)"), "WARNING", config)
  }
  
  # 检查共线性基因数量
  avg_genes <- total_genes / length(synteny_results$synteny_data)
  quality_results$statistics$avg_genes_per_pair <- avg_genes
  
  if(avg_genes < 10) {
    quality_results$overall_quality <- "warning"
    quality_results$issues <- c(quality_results$issues, "平均共线性基因数量较低")
    log_message(paste("警告: 平均共线性基因数量较低 (", round(avg_genes, 1), ")"), "WARNING", config)
  }
  
  # 创建各数据集的质量评分
  dataset_quality <- lapply(synteny_results$synteny_data, function(dataset) {
    # 基于共线性基因数量的质量评分
    gene_score <- min(1, dataset$gene_count / 100)
    
    # 基于覆盖染色体的质量评分
    chr_coverage <- dataset$chr_pairs / (dataset$chr1_count * dataset$chr2_count)
    chr_score <- min(1, chr_coverage * 2)  # 加倍以处理预期的稀疏性
    
    # 综合评分
    overall_score <- (gene_score + chr_score) / 2
    
    list(
      gene_count = dataset$gene_count,
      chr_coverage = chr_coverage,
      quality_score = overall_score
    )
  })
  
  quality_results$dataset_quality <- dataset_quality
  
  # 输出结果
  log_message(paste("共线性数据质量评估:", quality_results$overall_quality), "INFO", config)
  if(length(quality_results$issues) > 0) {
    log_message("发现的问题:", "INFO", config)
    for(issue in quality_results$issues) {
      log_message(paste(" -", issue), "INFO", config)
    }
  }
  log_message(paste("总体统计: 共", total_genes, "个共线性基因,", 
                  total_species, "个物种,", 
                  length(synteny_results$synteny_data), "个物种对"), "INFO", config)
  
  return(quality_results)
}

visualize_synteny_overview <- function(synteny_results, output_file = NULL, config) {
  # 创建数据概览
  if(length(synteny_results$synteny_data) == 0) {
    log_message("无数据可视化", "WARNING", config)
    return(NULL)
  }
  
  log_message("生成共线性数据概览可视化...", "INFO", config)
  
  # 提取每个物种对的信息
  overview_data <- tryCatch({
    data.frame(
      species_pair = names(synteny_results$synteny_data),
      species1 = sapply(synteny_results$synteny_data, function(x) x$species1),
      species2 = sapply(synteny_results$synteny_data, function(x) x$species2),
      gene_count = sapply(synteny_results$synteny_data, function(x) x$gene_count),
      chr1_count = sapply(synteny_results$synteny_data, function(x) x$chr1_count),
      chr2_count = sapply(synteny_results$synteny_data, function(x) x$chr2_count),
      chr_pairs = sapply(synteny_results$synteny_data, function(x) x$chr_pairs),
      distance = sapply(synteny_results$synteny_data, function(x) x$distance),
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    log_message(paste("创建概览数据时出错:", e$message), "ERROR", config)
    return(NULL)
  })
  
  if(is.null(overview_data) || nrow(overview_data) == 0) {
    log_message("无法创建概览数据", "ERROR", config)
    return(NULL)
  }
  
  # 创建图形
  tryCatch({
    if(!is.null(output_file)) {
      pdf(output_file, width=12, height=8)
      log_message(paste("创建PDF文件:", output_file), "DEBUG", config)
    }
    
    # 绘制基因数量条形图
    par(mar=c(10, 4, 4, 2) + 0.1) # 增加底部边距以容纳长标签
    barplot(overview_data$gene_count, 
            names.arg = overview_data$species_pair,
            las = 2, # 垂直标签
            cex.names = 0.7,
            col = "lightblue",
            main = "物种对共线性基因数量",
            ylab = "共线性基因数量")
    
    # 如果输出到文件，关闭设备
    if(!is.null(output_file)) {
      dev.off()
      log_message(paste("共线性数据概览已保存到:", output_file), "INFO", config)
    }
  }, error = function(e) {
    log_message(paste("创建概览图形时出错:", e$message), "ERROR", config)
    if(!is.null(output_file) && dev.cur() > 1) dev.off()
  })
  
  return(overview_data)
}

# ==============================
# 共线性块检测功能
# ==============================

extract_synteny_blocks <- function(synteny_data, config) {
  log_message("提取共线性块 (无脊椎动物优化)...", "INFO", config)
  
  # 初始化块数据库
  blocks <- list()
  block_id <- 1
  micro_blocks <- list()
  micro_block_id <- 1
  
  for(pair_name in names(synteny_data)) {
    pair_data <- synteny_data[[pair_name]]
    species1 <- pair_data$species1
    species2 <- pair_data$species2
    data <- pair_data$data
    
    log_message(paste("处理物种对:", species1, "-", species2, 
                     "（", nrow(data), "个共线性基因）"), "DEBUG", config)
    
    # 分组数据按染色体对
    chr_pairs <- paste0(data$chr1, "_", data$chr2)
    unique_pairs <- unique(chr_pairs)
    
    log_message(paste("发现", length(unique_pairs), "个染色体对"), "DEBUG", config)
    
    for(pair in unique_pairs) {
      # 获取此染色体对的数据
      pair_indices <- which(chr_pairs == pair)
      pair_data <- data[pair_indices, ]
      
      # 提取染色体
      chrs <- strsplit(pair, "_")[[1]]
      chr1 <- chrs[1]
      chr2 <- chrs[2]
      
      if(nrow(pair_data) < 2) {
        log_message(paste("跳过仅有", nrow(pair_data), "个基因的染色体对:", pair), "DEBUG", config)
        next
      }
      
      log_message(paste("处理染色体对:", pair, "（", nrow(pair_data), "个基因）"), "DEBUG", config)
      
      # 提取坐标
      x_coords <- pair_data$x
      y_coords <- pair_data$y
      
      # 应用改进的聚类算法，更适合检测重排基因组中的共线性
      clusters <- improved_synteny_clustering(x_coords, y_coords, config)
      
      # 检查聚类结果
      if(length(unique(clusters)) == 0 || all(clusters == 0)) {
        log_message(paste("对染色体对", pair, "的聚类未产生任何块"), "DEBUG", config)
        next
      }
      
      log_message(paste("染色体对", pair, "被分为", length(unique(clusters[clusters != 0])), "个块"), "DEBUG", config)
      
      # 处理每个聚类作为共线性块
      for(cluster_id in unique(clusters)) {
        # 跳过分配为0的点（未聚类）
        if(cluster_id == 0) next
        
        cluster_indices <- which(clusters == cluster_id)
        cluster_data <- pair_data[cluster_indices, ]
        
        # 日志记录块大小
        log_message(paste("块", cluster_id, "包含", nrow(cluster_data), "个基因"), "DEBUG", config)
        
        # 检测是否为宏共线性块
        if(nrow(cluster_data) >= config$min_block_size) {
          # 计算块内基因的顺序一致性 (检测倒位)
          order_consistency <- calc_order_consistency(cluster_data$x, cluster_data$y)
          is_inverted <- order_consistency < 0
          
          # 日志块信息
          log_message(paste("添加宏共线性块:", species1, ":", chr1, "-", species2, ":", chr2, 
                           "(", nrow(cluster_data), "个基因, 顺序一致性=", round(order_consistency, 2), 
                           ifelse(is_inverted, ", 倒位", ""), ")"), "DEBUG", config)
          
          # 创建块条目
          blocks[[paste0("block_", block_id)]] <- list(
            id = block_id,
            species1 = species1,
            species2 = species2,
            chr1 = chr1,
            chr2 = chr2,
            genes = nrow(cluster_data),
            x_min = min(cluster_data$x),
            x_max = max(cluster_data$x),
            y_min = min(cluster_data$y),
            y_max = max(cluster_data$y),
            order_consistency = order_consistency,
            is_inverted = is_inverted,
            block_type = "macro",
            data = cluster_data
          )
          
          block_id <- block_id + 1
        }
        # 检测微共线性块
        else if(config$detect_micro_synteny && 
                nrow(cluster_data) >= config$micro_synteny_min_size) {
          
          # 计算顺序一致性
          order_consistency <- calc_order_consistency(cluster_data$x, cluster_data$y)
          
          log_message(paste("添加微共线性块:", species1, ":", chr1, "-", species2, ":", chr2, 
                           "(", nrow(cluster_data), "个基因, 顺序一致性=", round(order_consistency, 2), ")"), 
                     "DEBUG", config)
          
          micro_blocks[[paste0("micro_", micro_block_id)]] <- list(
            id = micro_block_id,
            species1 = species1,
            species2 = species2,
            chr1 = chr1,
            chr2 = chr2,
            genes = nrow(cluster_data),
            x_min = min(cluster_data$x),
            x_max = max(cluster_data$x),
            y_min = min(cluster_data$y),
            y_max = max(cluster_data$y),
            order_consistency = order_consistency,
            is_inverted = order_consistency < 0,
            block_type = "micro",
            data = cluster_data
          )
          
          micro_block_id <- micro_block_id + 1
        }
      }
    }
  }
  
  log_message(paste("提取了", length(blocks), "个宏共线性块和", 
                   length(micro_blocks), "个微共线性块"), "INFO", config)
  
  return(list(
    macro_blocks = blocks,
    micro_blocks = micro_blocks
  ))
}

# 改进的共线性聚类算法，针对无脊椎动物的重排特点优化
improved_synteny_clustering <- function(x_coords, y_coords, config) {
  # 初始化聚类
  n_points <- length(x_coords)
  clusters <- rep(0, n_points)
  current_cluster <- 1
  
  # 特殊情况处理
  if(n_points == 0) return(clusters)
  if(n_points == 1) {
    clusters[1] <- current_cluster
    return(clusters)
  }
  
  log_message(paste("对", n_points, "个点执行改进共线性聚类 (max_gap=", config$max_gap, ")"), "DEBUG", config)
  
  # 计算所有点对之间的距离矩阵
  dist_matrix <- matrix(NA, nrow=n_points, ncol=n_points)
  for(i in 1:n_points) {
    for(j in 1:n_points) {
      # 使用曼哈顿距离可以更好地检测沿轴分布的共线性
      dist_matrix[i,j] <- abs(x_coords[i] - x_coords[j]) + 
                           abs(y_coords[i] - y_coords[j])
    }
  }
  
  # 使用密度聚类方法，更适合捕获无规则形状的共线性块
  # 这里使用简化版本的DBSCAN
  visited <- rep(FALSE, n_points)
  
  for(i in 1:n_points) {
    if(visited[i]) next
    
    # 标记当前点为已访问
    visited[i] <- TRUE
    
    # 查找邻居 (距离小于阈值的点)
    neighbors <- which(dist_matrix[i,] <= config$max_gap)
    
    if(length(neighbors) >= config$micro_synteny_min_size) {
      # 创建新簇
      clusters[i] <- current_cluster
      
      # 扩展簇
      j <- 1
      while(j <= length(neighbors)) {
        next_pt <- neighbors[j]
        
        if(!visited[next_pt]) {
          visited[next_pt] <- TRUE
          new_neighbors <- which(dist_matrix[next_pt,] <= config$max_gap)
          
          if(length(new_neighbors) >= config$micro_synteny_min_size) {
            # 合并新邻居
            new_neighbors <- setdiff(new_neighbors, neighbors)
            neighbors <- c(neighbors, new_neighbors)
          }
        }
        
        # 如果点还没被分配到簇，将其添加到当前簇
        if(clusters[next_pt] == 0) {
          clusters[next_pt] <- current_cluster
        }
        
        j <- j + 1
      }
      
      log_message(paste("创建簇", current_cluster, "包含", sum(clusters == current_cluster), "个点"), "DEBUG", config)
      current_cluster <- current_cluster + 1
    } 
    # 对孤立点的特殊处理：如果附近有其他点，形成微共线性块
    else if(length(neighbors) >= 2) {
      clusters[i] <- current_cluster
      for(n in neighbors) {
        if(clusters[n] == 0) {
          clusters[n] <- current_cluster
        }
      }
      log_message(paste("创建小簇", current_cluster, "包含", sum(clusters == current_cluster), "个点"), "DEBUG", config)
      current_cluster <- current_cluster + 1
    }
  }
  
  # 对未分配的点进行后续处理 - 分配给最近的簇或创建新簇
  unassigned <- which(clusters == 0)
  if(length(unassigned) > 0) {
    log_message(paste("处理", length(unassigned), "个未分配点"), "DEBUG", config)
    for(i in unassigned) {
      # 寻找距离最近的已分配点
      min_dist <- Inf
      closest_cluster <- 0
      
      for(j in setdiff(1:n_points, unassigned)) {
        if(dist_matrix[i,j] < min_dist) {
          min_dist <- dist_matrix[i,j]
          closest_cluster <- clusters[j]
        }
      }
      
      # 如果最近点足够近，分配到同一簇
      if(min_dist <= config$max_gap * 1.5) {
        clusters[i] <- closest_cluster
        log_message(paste("将点", i, "分配到最近的簇", closest_cluster, "(距离=", 
                         round(min_dist, 2), ")"), "DEBUG", config)
      } else {
        # 否则创建新簇
        current_cluster <- current_cluster + 1
        clusters[i] <- current_cluster
        log_message(paste("为孤立点", i, "创建新簇", current_cluster), "DEBUG", config)
      }
    }
  }
  
  log_message(paste("聚类完成，共", length(unique(clusters[clusters > 0])), "个簇,", 
                   sum(clusters == 0), "个未分配点"), "DEBUG", config)
  
  return(clusters)
}

# 计算序列顺序一致性，用于检测倒位
calc_order_consistency <- function(x_coords, y_coords) {
  if(length(x_coords) < 2) return(0)
  
  # 计算斯皮尔曼相关系数
  correlation <- tryCatch({
    cor(x_coords, y_coords, method="spearman")
  }, error = function(e) {
    # 如果计算相关性失败（例如，所有值相同），返回0
    return(0)
  })
  
  return(correlation)
}

# ==============================
# 事件检测函数
# ==============================

detect_synteny_events <- function(synteny_data, chr_networks, tree, config) {
  log_message("检测共线性事件(无脊椎动物优化)...", "INFO", config)
  
  # 初始化事件列表
  events <- list(
    fusions = data.frame(),
    fissions = data.frame(),
    translocations = data.frame(),
    inversions = data.frame(),   # 倒位事件
    complex_rearrangements = data.frame(),  # 复杂重排事件
    conserved = data.frame(),
    hotspots = data.frame()      # 重排热点
  )
  
  # 检测融合事件 (多个染色体映射到单个染色体)
  fusion_rows <- list()
  fusion_count <- 0
  
  for(pair_name in names(synteny_data)) {
    pair_data <- synteny_data[[pair_name]]
    species1 <- pair_data$species1
    species2 <- pair_data$species2
    data <- pair_data$data
    
    log_message(paste("分析物种对", species1, "-", species2, "的融合事件"), "DEBUG", config)
    
    # 计算每个染色体的共线性基因数
    chr1_counts <- table(data$chr1)
    chr2_counts <- table(data$chr2)
    
    # 寻找物种2中潜在融合目标 (一个染色体有多个配对伙伴)
    for(chr2 in names(chr2_counts)) {
      partners <- data[data$chr2 == chr2, "chr1"]
      partner_table <- table(partners)
      
      if(length(partner_table) > 1) {
        # 物种2中的潜在融合
        # 使用更灵活的阈值，适应无脊椎动物的小块共线性
        major_partners <- names(partner_table)[partner_table >= config$min_block_size - 1]
        
        if(length(major_partners) > 1) {
          # 计算融合置信度 - 无脊椎动物修正版
          total_genes <- sum(partner_table[major_partners])
          avg_per_partner <- total_genes / length(major_partners)
          
          # 使用灵活的基因数阈值
          min_genes_threshold <- max(3, min(5, config$min_block_size * 1.5))
          
          # 包括更多潜在的融合事件
          if(avg_per_partner >= min_genes_threshold || 
             (total_genes >= min_genes_threshold * 2 && length(major_partners) >= 3)) {
            
            # 应用逆向距离加权，优先考虑系统发育上更近的物种
            distance_factor <- 1.0
            tryCatch({
              if(!is.na(pair_data$distance)) {
                distance <- pair_data$distance
                distance_factor <- max(0.5, 1 - (distance / 2))
                log_message(paste("物种间距离因子:", round(distance_factor, 3)), "DEBUG", config)
              }
            }, error = function(e) {})
            
            # 计算事件的置信度
            confidence_score <- min(1.0, (avg_per_partner / 8) * distance_factor)
            
            # 日志记录潜在融合事件
            log_message(paste("检测到潜在融合事件:", paste(major_partners, collapse=","), 
                             "->", chr2, "(", species1, "->", species2, ")",
                             "基因数:", total_genes, "置信度:", round(confidence_score, 3)), 
                       "DEBUG", config)
            
            fusion_count <- fusion_count + 1
            fusion_rows[[fusion_count]] <- data.frame(
              species_ancestor = species1,
              species_derived = species2,
              chr_ancestor = paste(major_partners, collapse=","),
              chr_derived = chr2,
              genes_involved = total_genes,
              partner_count = length(major_partners),
              event_type = "fusion",
              confidence = confidence_score,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
    
    # 寻找物种1中的潜在融合 (类似逻辑应用于物种1)
    for(chr1 in names(chr1_counts)) {
      partners <- data[data$chr1 == chr1, "chr2"]
      partner_table <- table(partners)
      
      if(length(partner_table) > 1) {
        # 物种1中的潜在融合
        major_partners <- names(partner_table)[partner_table >= config$min_block_size - 1]
        
        if(length(major_partners) > 1) {
          # 计算融合置信度
          total_genes <- sum(partner_table[major_partners])
          avg_per_partner <- total_genes / length(major_partners)
          
          # 使用灵活的阈值
          min_genes_threshold <- max(3, min(5, config$min_block_size * 1.5))
          
          if(avg_per_partner >= min_genes_threshold || 
             (total_genes >= min_genes_threshold * 2 && length(major_partners) >= 3)) {
            
            # 应用系统发育距离加权
            distance_factor <- 1.0
            tryCatch({
              if(!is.na(pair_data$distance)) {
                distance <- pair_data$distance
                distance_factor <- max(0.5, 1 - (distance / 2))
              }
            }, error = function(e) {})
            
            # 计算置信度
            confidence_score <- min(1.0, (avg_per_partner / 8) * distance_factor)
            
            log_message(paste("检测到潜在融合事件:", paste(major_partners, collapse=","), 
                             "->", chr1, "(", species2, "->", species1, ")",
                             "基因数:", total_genes, "置信度:", round(confidence_score, 3)), 
                       "DEBUG", config)
            
            fusion_count <- fusion_count + 1
            fusion_rows[[fusion_count]] <- data.frame(
              species_ancestor = species2,
              species_derived = species1,
              chr_ancestor = paste(major_partners, collapse=","),
              chr_derived = chr1,
              genes_involved = total_genes,
              partner_count = length(major_partners),
              event_type = "fusion",
              confidence = confidence_score,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }
  
  if(length(fusion_rows) > 0) {
    events$fusions <- do.call(rbind, fusion_rows)
    log_message(paste("检测到", nrow(events$fusions), "个潜在融合事件"), "INFO", config)
  } else {
    log_message("未检测到融合事件", "INFO", config)
  }
  
  # 检测裂变事件 (基于融合事件的转换)
  fission_rows <- list()
  fission_count <- 0
  
  if(nrow(events$fusions) > 0) {
    log_message("从融合事件转换裂变事件...", "DEBUG", config)
    # 将融合事件转换为裂变事件 (它们是互逆的)
    for(i in 1:nrow(events$fusions)) {
      fusion <- events$fusions[i,]
      
      fission_count <- fission_count + 1
      fission_rows[[fission_count]] <- data.frame(
        species_ancestor = fusion$species_derived,
        species_derived = fusion$species_ancestor,
        chr_ancestor = fusion$chr_derived,
        chr_derived = fusion$chr_ancestor,
        genes_involved = fusion$genes_involved,
        partner_count = fusion$partner_count,
        event_type = "fission",
        confidence = fusion$confidence * 0.9,  # 稍微降低置信度
        stringsAsFactors = FALSE
      )
      
      log_message(paste("转换为裂变事件:", fusion$chr_derived, "->", 
                       fusion$chr_ancestor, "(", fusion$species_derived, "->", 
                       fusion$species_ancestor, ")"), "DEBUG", config)
    }
    
    events$fissions <- do.call(rbind, fission_rows)
    log_message(paste("转换为", nrow(events$fissions), "个潜在裂变事件"), "INFO", config)
  } else {
    log_message("无融合事件可转换为裂变事件", "DEBUG", config)
  }
  
  # 检测倒位事件
  inversion_rows <- list()
  inversion_id <- 1
  
  # 通过检查相邻基因对的方向一致性来检测倒位
  log_message("检测倒位事件...", "DEBUG", config)
  
  for(pair_name in names(synteny_data)) {
    pair_data <- synteny_data[[pair_name]]
    species1 <- pair_data$species1
    species2 <- pair_data$species2
    data <- pair_data$data
    
    # 分组数据按染色体对
    chr_pairs <- paste0(data$chr1, "_", data$chr2)
    unique_pairs <- unique(chr_pairs)
    
    for(pair in unique_pairs) {
      # 获取此染色体对的数据
      pair_indices <- which(chr_pairs == pair)
      pair_data <- data[pair_indices, ]
      
      # 至少需要多个基因才能检测倒位
      if(nrow(pair_data) >= 4) {
        # 提取染色体
        chrs <- strsplit(pair, "_")[[1]]
        chr1 <- chrs[1]
        chr2 <- chrs[2]
        
        # 计算顺序相关性
        correlation <- cor(pair_data$x, pair_data$y, method="spearman")
        
        # 检测强烈的负相关，表明倒位
        if(correlation < -config$inversion_threshold) {
          log_message(paste("检测到倒位事件:", species1, ":", chr1, "<->", 
                           species2, ":", chr2, "相关性:", round(correlation, 3)), 
                     "DEBUG", config)
          
          inversion_rows[[inversion_id]] <- data.frame(
            species1 = species1,
            species2 = species2,
            chr1 = chr1,
            chr2 = chr2,
            gene_count = nrow(pair_data),
            correlation = correlation,
            event_type = "inversion",
            confidence = min(1.0, abs(correlation)),
            stringsAsFactors = FALSE
          )
          inversion_id <- inversion_id + 1
        }
      }
    }
  }
  
  if(length(inversion_rows) > 0) {
    events$inversions <- do.call(rbind, inversion_rows)
    log_message(paste("检测到", nrow(events$inversions), "个潜在倒位事件"), "INFO", config)
  } else {
    log_message("未检测到倒位事件", "INFO", config)
  }
  
  # 检测保守共线性关系
  log_message("检测保守的共线性关系...", "DEBUG", config)
  conserved_rows <- list()
  cons_id <- 1
  
  for(pair_name in names(synteny_data)) {
    pair_data <- synteny_data[[pair_name]]
    species1 <- pair_data$species1
    species2 <- pair_data$species2
    data <- pair_data$data
    
    # 计算每个染色体对的基因数
    chr_pairs <- paste0(data$chr1, "_", data$chr2)
    pair_counts <- table(chr_pairs)
    
    for(pair in names(pair_counts)) {
      count <- pair_counts[[pair]]
      chrs <- strsplit(pair, "_")[[1]]
      
      # 跳过小型关联
      if(count < config$min_block_size) next
      
      # 检查这是否是两染色体的主要关联
      chr1 <- chrs[1]
      chr2 <- chrs[2]
      
      chr1_partners <- data[data$chr1 == chr1, "chr2"]
      chr1_partner_counts <- table(chr1_partners)
      
      chr2_partners <- data[data$chr2 == chr2, "chr1"]
      chr2_partner_counts <- table(chr2_partners)
      
      # 检查这些染色体是否是互相的主要伙伴
      is_primary_partner <- FALSE
      if(length(chr1_partner_counts) > 0 && length(chr2_partner_counts) > 0) {
        if(names(chr1_partner_counts)[which.max(chr1_partner_counts)] == chr2 && 
           names(chr2_partner_counts)[which.max(chr2_partner_counts)] == chr1) {
          is_primary_partner <- TRUE
        }
      }
      
      # 使用无脊椎动物优化阈值
      threshold <- config$min_block_size * 1.5
      
      if(is_primary_partner && count >= threshold) {
        log_message(paste("检测到保守的共线性关系:", species1, ":", chr1, "<->", 
                         species2, ":", chr2, "基因数:", count), "DEBUG", config)
        
        conserved_rows[[cons_id]] <- data.frame(
          id = paste0("consv_", cons_id),
          species1 = species1,
          species2 = species2,
          chr1 = chr1,
          chr2 = chr2,
          gene_count = count,
          confidence = min(1.0, count / 15),  # 无脊椎动物优化的置信度计算
          stringsAsFactors = FALSE
        )
        cons_id <- cons_id + 1
      }
    }
  }
  
  if(length(conserved_rows) > 0) {
    events$conserved <- do.call(rbind, conserved_rows)
    log_message(paste("检测到", nrow(events$conserved), "个保守共线性关系"), "INFO", config)
  } else {
    log_message("未检测到保守共线性关系", "INFO", config)
  }
  
  # 检测重排热点
  if(config$detect_rearrangement_hotspots) {
    log_message("检测染色体重排热点...", "DEBUG", config)
    hotspot_rows <- detect_rearrangement_hotspots(
      events$fusions, events$fissions, events$inversions, config
    )
    
    if(length(hotspot_rows) > 0) {
      events$hotspots <- do.call(rbind, hotspot_rows)
      log_message(paste("检测到", nrow(events$hotspots), "个重排热点"), "INFO", config)
    } else {
      log_message("未检测到重排热点", "INFO", config)
    }
  }
  
  log_message("事件检测完成 (无脊椎动物优化)", "INFO", config)
  return(events)
}

# 检测重排热点的函数
detect_rearrangement_hotspots <- function(fusions, fissions, inversions, config) {
  # 初始化热点列表
  hotspot_rows <- list()
  hotspot_id <- 1
  
  # 如果没有足够的事件数据，返回空列表
  total_events <- nrow(fusions) + nrow(fissions) + nrow(inversions)
  if(total_events < config$min_rearrangements_for_hotspot) {
    log_message(paste("事件总数(", total_events, ")低于热点检测阈值(", 
                     config$min_rearrangements_for_hotspot, ")"), "DEBUG", config)
    return(hotspot_rows)
  }
  
  log_message(paste("从", total_events, "个事件中搜索重排热点..."), "DEBUG", config)
  
  # 收集所有染色体和涉及的重排事件
  all_chrs <- data.frame()
  
  # 处理融合事件
  if(nrow(fusions) > 0) {
    fusion_chrs <- rbind(
      data.frame(
        species = fusions$species_ancestor,
        chromosome = fusions$chr_ancestor,
        event_type = "fusion_source",
        stringsAsFactors = FALSE
      ),
      data.frame(
        species = fusions$species_derived,
        chromosome = fusions$chr_derived,
        event_type = "fusion_target",
        stringsAsFactors = FALSE
      )
    )
    all_chrs <- rbind(all_chrs, fusion_chrs)
    log_message(paste("添加", nrow(fusion_chrs), "条染色体-融合事件记录"), "DEBUG", config)
  }
  
  # 处理裂变事件
  if(nrow(fissions) > 0) {
    fission_chrs <- rbind(
      data.frame(
        species = fissions$species_ancestor,
        chromosome = fissions$chr_ancestor,
        event_type = "fission_source",
        stringsAsFactors = FALSE
      ),
      data.frame(
        species = fissions$species_derived,
        chromosome = fissions$chr_derived,
        event_type = "fission_target",
        stringsAsFactors = FALSE
      )
    )
    all_chrs <- rbind(all_chrs, fission_chrs)
    log_message(paste("添加", nrow(fission_chrs), "条染色体-裂变事件记录"), "DEBUG", config)
  }
  
  # 处理倒位事件
  if(nrow(inversions) > 0) {
    inversion_chrs <- rbind(
      data.frame(
        species = inversions$species1,
        chromosome = inversions$chr1,
        event_type = "inversion",
        stringsAsFactors = FALSE
      ),
      data.frame(
        species = inversions$species2,
        chromosome = inversions$chr2,
        event_type = "inversion",
        stringsAsFactors = FALSE
      )
    )
    all_chrs <- rbind(all_chrs, inversion_chrs)
    log_message(paste("添加", nrow(inversion_chrs), "条染色体-倒位事件记录"), "DEBUG", config)
  }
  
  # 计算每个染色体涉及的重排事件数
  if(nrow(all_chrs) > 0) {
    chr_ids <- paste(all_chrs$species, all_chrs$chromosome, sep=":")
    chr_counts <- table(chr_ids)
    
    log_message(paste("分析", length(chr_counts), "个染色体的事件频率"), "DEBUG", config)
    
    # 识别高于阈值的热点
    for(chr in names(chr_counts)) {
      if(chr_counts[chr] >= config$min_rearrangements_for_hotspot) {
        # 解析物种和染色体
        parts <- strsplit(chr, ":")[[1]]
        species <- parts[1]
        chromosome <- paste(parts[-1], collapse=":")  # 处理染色体ID中可能包含冒号的情况
        
        # 获取此染色体涉及的所有事件类型
        events_involved <- all_chrs[all_chrs$species == species & 
                                     all_chrs$chromosome == chromosome, "event_type"]
        event_types <- paste(unique(events_involved), collapse=",")
        
        log_message(paste("检测到重排热点:", species, ":", chromosome, 
                         "涉及", chr_counts[chr], "个事件:", event_types), "DEBUG", config)
        
        # 添加到热点列表
        hotspot_rows[[hotspot_id]] <- data.frame(
          species = species,
          chromosome = chromosome,
          event_count = chr_counts[chr],
          event_types = event_types,
          stringsAsFactors = FALSE
        )
        hotspot_id <- hotspot_id + 1
      }
    }
  }
  
  log_message(paste("检测到", length(hotspot_rows), "个重排热点"), "INFO", config)
  return(hotspot_rows)
}

# ==============================
# 事件到节点映射函数
# ==============================

map_events_to_nodes <- function(synteny_events, tree, phase4_results = NULL, config) {
  log_message("将共线性事件映射到祖先节点(无脊椎动物优化)...", "INFO", config)
  
  # 获取树中的所有末端物种
  tip_species <- tree$tip.label
  
  # 获取所有内部节点
  n_tips <- length(tip_species)
  internal_nodes <- (n_tips + 1):(n_tips + tree$Nnode)
  
  log_message(paste("进化树包含", n_tips, "个物种和", tree$Nnode, "个内部节点"), "DEBUG", config)
  
  # 初始化节点映射
  node_mapping <- list()
  
  for(node_id in internal_nodes) {
    # 获取节点的后代
    descendants <- get_descendants(tree, node_id)
    
    # 获取此节点下的末端物种
    tip_descendants <- get_tip_descendants(tree, node_id)
    node_species <- tree$tip.label[tip_descendants]
    
    log_message(paste("处理节点", node_id, "下的", length(node_species), "个物种"), "DEBUG", config)
    
    # 找出与此节点相关的事件
    node_events <- list(
      fusions = data.frame(),
      fissions = data.frame(),
      inversions = data.frame(),
      conserved = list()
    )
    
    # 检查融合事件
    if(nrow(synteny_events$fusions) > 0) {
      # 寻找涉及此节点下物种的事件
      relevant_fusions <- synteny_events$fusions[
        synteny_events$fusions$species_ancestor %in% node_species &
          synteny_events$fusions$species_derived %in% node_species,
      ]
      
      if(nrow(relevant_fusions) > 0) {
        node_events$fusions <- relevant_fusions
        log_message(paste("节点", node_id, "相关融合事件:", nrow(relevant_fusions)), "DEBUG", config)
      }
    }
    
    # 检查裂变事件
    if(nrow(synteny_events$fissions) > 0) {
      relevant_fissions <- synteny_events$fissions[
        synteny_events$fissions$species_ancestor %in% node_species &
          synteny_events$fissions$species_derived %in% node_species,
      ]
      
      if(nrow(relevant_fissions) > 0) {
        node_events$fissions <- relevant_fissions
        log_message(paste("节点", node_id, "相关裂变事件:", nrow(relevant_fissions)), "DEBUG", config)
      }
    }
    
    # 检查倒位事件
    if(nrow(synteny_events$inversions) > 0) {
      relevant_inversions <- synteny_events$inversions[
        (synteny_events$inversions$species1 %in% node_species &
          synteny_events$inversions$species2 %in% node_species),
      ]
      
      if(nrow(relevant_inversions) > 0) {
        node_events$inversions <- relevant_inversions
        log_message(paste("节点", node_id, "相关倒位事件:", nrow(relevant_inversions)), "DEBUG", config)
      }
    }
    
    # 检查保守共线性
    if(nrow(synteny_events$conserved) > 0) {
      relevant_conserved <- synteny_events$conserved[
        synteny_events$conserved$species1 %in% node_species &
          synteny_events$conserved$species2 %in% node_species,
      ]
      
      if(nrow(relevant_conserved) > 0) {
        # 按染色体对分组
        cons_groups <- split(relevant_conserved, 
                             paste0(relevant_conserved$species1, "_", relevant_conserved$chr1))
        
        node_events$conserved <- cons_groups
        log_message(paste("节点", node_id, "相关保守共线性关系:", 
                         length(cons_groups), "组"), "DEBUG", config)
      }
    }
    
    # 添加到节点映射（如果找到任何相关事件）
    if(nrow(node_events$fusions) > 0 || 
       nrow(node_events$fissions) > 0 || 
       nrow(node_events$inversions) > 0 || 
       length(node_events$conserved) > 0) {
      
      # 为节点命名（如果可能）
      node_label <- if(!is.null(tree$node.label) && length(tree$node.label) >= node_id - n_tips) {
        tree$node.label[node_id - n_tips]
      } else {
        paste("Node", node_id)
      }
      
      node_mapping[[as.character(node_id)]] <- list(
        node_id = node_id,
        node_label = node_label,
        species = node_species,
        events = node_events
      )
      
      log_message(paste("向节点", node_id, "(", node_label, ")添加事件映射"), "DEBUG", config)
    }
  }
  
  log_message(paste("将共线性事件映射到", length(node_mapping), "个祖先节点"), "INFO", config)
  
  return(node_mapping)
}

# ==============================
# 数据处理和报告生成函数
# ==============================

process_synteny_data <- function(synteny_dir, tree_file, output_dir, config_file = NULL) {
  log_message("开始处理共线性数据...", "INFO")
  
  # 创建输出目录
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    log_message(paste("创建输出目录:", output_dir), "INFO")
  }
  
  # 加载配置
  config <- load_config(config_file)
  
  # 读取共线性数据
  synteny_results <- tryCatch({
    load_synteny_files(synteny_dir, tree_file, config)
  }, error = function(e) {
    log_message(paste("加载共线性数据时出错:", e$message), "ERROR", config)
    stop(paste("加载共线性数据时出错:", e$message))
  })
  
  if(is.null(synteny_results) || length(synteny_results$synteny_data) == 0) {
    log_message("未能成功加载任何共线性数据，处理终止", "ERROR", config)
    stop("未能成功加载任何共线性数据")
  }
  
  # 验证数据质量
  quality_results <- validate_synteny_quality(synteny_results, config)
  
  # 创建数据概览可视化
  if(config$generate_visualizations) {
    overview_plot <- visualize_synteny_overview(
      synteny_results, 
      file.path(output_dir, "synteny_overview.pdf"),
      config
    )
  }
  
  # 保存处理过的数据
  synteny_processed <- list(
    data = synteny_results,
    quality = quality_results,
    config = config,
    processing_date = Sys.time()
  )
  
  saveRDS(synteny_processed, file.path(output_dir, "synteny_processed.RData"))
  log_message(paste("处理后的共线性数据已保存到:", file.path(output_dir, "synteny_processed.RData")), "INFO", config)
  
  return(synteny_processed)
}

# 生成分析报告
generate_synteny_report <- function(analysis_results, output_dir, config = NULL) {
  if(is.null(config)) {
    config <- analysis_results$synteny_data$config
    if(is.null(config)) config <- list(verbose = TRUE)
  }
  
  log_message("生成共线性分析报告...", "INFO", config)
  
  # 创建报告文件
  report_file <- file.path(output_dir, "synteny_analysis_report.txt")
  
  # 打开文件写入
  sink(report_file)
  
  cat("=======================================================\n")
  cat("共线性分析报告 (无脊椎动物优化版)\n")
  cat(paste0("日期: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n"))
  cat(paste0("用户: ", Sys.info()["user"], "\n"))
  cat("=======================================================\n\n")
  
  # 摘要统计
  cat("1. 共线性数据摘要\n")
  cat("-------------------------\n")
  
  synteny_results <- analysis_results$synteny_data
  
  if(length(synteny_results$synteny_data) > 0) {
    cat(paste0("共线性数据集总数: ", length(synteny_results$synteny_data), "\n"))
    
    total_genes <- sum(sapply(synteny_results$synteny_data, function(x) x$gene_count))
    cat(paste0("共线性基因总数: ", total_genes, "\n"))
    
    all_species <- unique(unlist(lapply(synteny_results$synteny_data, function(x) c(x$species1, x$species2))))
    cat(paste0("物种数量: ", length(all_species), "\n"))
    cat(paste0("包含物种: ", paste(all_species, collapse=", "), "\n\n"))
    
    # 打印按共线性基因数排序的顶部物种对
    if(length(synteny_results$synteny_data) > 0) {
      gene_counts <- sapply(synteny_results$synteny_data, function(x) x$gene_count)
      sorted_indices <- order(gene_counts, decreasing = TRUE)
      top_pairs <- sorted_indices[1:min(5, length(sorted_indices))]
      
      cat("共线性基因数量最多的物种对:\n")
      for(i in top_pairs) {
        pair_name <- names(synteny_results$synteny_data)[i]
        pair_data <- synteny_results$synteny_data[[pair_name]]
        cat(paste0("  - ", pair_data$species1, " <-> ", pair_data$species2,
                 ": ", pair_data$gene_count, " 个基因, ",
                 pair_data$chr1_count, " 个染色体(物种1), ",
                 pair_data$chr2_count, " 个染色体(物种2)\n"))
      }
    }
  }
  
  cat("\n")
  
  # 共线性块统计
  if(!is.null(analysis_results$blocks)) {
    blocks <- analysis_results$blocks
    
    cat("2. 共线性块摘要\n")
    cat("-------------------------\n")
    
    # 宏共线性块
    macro_count <- length(blocks$macro_blocks)
    cat(paste0("宏共线性块数量: ", macro_count, "\n"))
    
    if(macro_count > 0) {
      block_sizes <- sapply(blocks$macro_blocks, function(x) x$genes)
      cat(paste0("宏共线性块大小分布:\n"))
      cat(paste0("  - 最小: ", min(block_sizes), " 基因\n"))
      cat(paste0("  - 最大: ", max(block_sizes), " 基因\n"))
      cat(paste0("  - 平均: ", round(mean(block_sizes), 1), " 基因\n"))
      cat(paste0("  - 中位数: ", median(block_sizes), " 基因\n"))
      
      # 最大的宏共线性块
      cat("\n最大的宏共线性块:\n")
      sorted_indices <- order(block_sizes, decreasing = TRUE)
      top_indices <- sorted_indices[1:min(5, length(sorted_indices))]
      
      for(i in top_indices) {
        block <- blocks$macro_blocks[[i]]
        cat(paste0("  - ", block$species1, ":", block$chr1, 
                 " <-> ", block$species2, ":", block$chr2,
                 " (", block$genes, " 基因, 顺序一致性: ", 
                 round(block$order_consistency, 2), ")\n"))
      }
    }
    
    # 微共线性块
    micro_count <- length(blocks$micro_blocks)
    cat(paste0("\n微共线性块数量: ", micro_count, "\n"))
    
    if(micro_count > 0) {
      micro_sizes <- sapply(blocks$micro_blocks, function(x) x$genes)
      cat(paste0("微共线性块大小分布:\n"))
      cat(paste0("  - 最小: ", min(micro_sizes), " 基因\n"))
      cat(paste0("  - 最大: ", max(micro_sizes), " 基因\n"))
      cat(paste0("  - 平均: ", round(mean(micro_sizes), 1), " 基因\n"))
    }
  }
  
  cat("\n")
  
  # 事件检测摘要
  if(!is.null(analysis_results$events)) {
    events <- analysis_results$events
    
    cat("3. 事件检测摘要\n")
    cat("-------------------------\n")
    
    cat(paste0("融合事件数量: ", nrow(events$fusions), "\n"))
    if(nrow(events$fusions) > 0) {
      cat("置信度最高的融合事件:\n")
      sorted_fusions <- events$fusions[order(events$fusions$confidence, decreasing = TRUE),]
      top_fusions <- head(sorted_fusions, min(5, nrow(sorted_fusions)))
      
      for(i in 1:nrow(top_fusions)) {
        cat(paste0("  - ", top_fusions$chr_ancestor[i], " in ", top_fusions$species_ancestor[i],
                 " -> ", top_fusions$chr_derived[i], " in ", top_fusions$species_derived[i],
                 " (置信度: ", round(top_fusions$confidence[i], 2), 
                 ", 涉及", top_fusions$genes_involved[i], "个基因)\n"))
      }
    }
    
    cat(paste0("\n裂变事件数量: ", nrow(events$fissions), "\n"))
    if(nrow(events$fissions) > 0) {
      cat("置信度最高的裂变事件:\n")
      sorted_fissions <- events$fissions[order(events$fissions$confidence, decreasing = TRUE),]
      top_fissions <- head(sorted_fissions, min(5, nrow(sorted_fissions)))
      
      for(i in 1:nrow(top_fissions)) {
        cat(paste0("  - ", top_fissions$chr_ancestor[i], " in ", top_fissions$species_ancestor[i],
                 " -> ", top_fissions$chr_derived[i], " in ", top_fissions$species_derived[i],
                 " (置信度: ", round(top_fissions$confidence[i], 2), 
                 ", 涉及", top_fissions$genes_involved[i], "个基因)\n"))
      }
    }
    
    cat(paste0("\n倒位事件数量: ", nrow(events$inversions), "\n"))
    if(nrow(events$inversions) > 0) {
      cat("最显著的倒位事件:\n")
      sorted_inversions <- events$inversions[order(abs(events$inversions$correlation), decreasing = TRUE),]
      top_inversions <- head(sorted_inversions, min(5, nrow(sorted_inversions)))
      
      for(i in 1:nrow(top_inversions)) {
        cat(paste0("  - ", top_inversions$species1[i], ":", top_inversions$chr1[i],
                 " <-> ", top_inversions$species2[i], ":", top_inversions$chr2[i],
                 " (相关性: ", round(top_inversions$correlation[i], 2), 
                 ", 涉及", top_inversions$gene_count[i], "个基因)\n"))
      }
    }
    
    cat(paste0("\n保守共线性关系数量: ", nrow(events$conserved), "\n"))
    if(nrow(events$conserved) > 0) {
      cat("基因数量最多的保守共线性关系:\n")
      sorted_conserved <- events$conserved[order(events$conserved$gene_count, decreasing = TRUE),]
      top_conserved <- head(sorted_conserved, min(5, nrow(sorted_conserved)))
      
      for(i in 1:nrow(top_conserved)) {
        cat(paste0("  - ", top_conserved$species1[i], ":", top_conserved$chr1[i],
                 " <-> ", top_conserved$species2[i], ":", top_conserved$chr2[i],
                 " (", top_conserved$gene_count[i], "个基因, 置信度: ", 
                 round(top_conserved$confidence[i], 2), ")\n"))
      }
    }
    
    cat(paste0("\n重排热点数量: ", nrow(events$hotspots), "\n"))
    if(nrow(events$hotspots) > 0) {
      cat("事件最多的重排热点:\n")
      sorted_hotspots <- events$hotspots[order(events$hotspots$event_count, decreasing = TRUE),]
      top_hotspots <- head(sorted_hotspots, min(5, nrow(sorted_hotspots)))
      
      for(i in 1:nrow(top_hotspots)) {
        cat(paste0("  - ", top_hotspots$species[i], ":", top_hotspots$chromosome[i],
                 " (涉及", top_hotspots$event_count[i], "个事件: ", 
                 top_hotspots$event_types[i], ")\n"))
      }
    }
  }
  
  cat("\n")
  
  # 节点映射摘要
  if(!is.null(analysis_results$node_mapping)) {
    node_mapping <- analysis_results$node_mapping
    
    cat("4. 系统发育树节点事件映射\n")
    cat("-------------------------\n")
    
    cat(paste0("共线性事件已映射到", length(node_mapping), "个祖先节点\n\n"))
    
    for(node_id in names(node_mapping)) {
      node_data <- node_mapping[[node_id]]
      
      cat(paste0("节点 ", node_id, " (", node_data$node_label, "):\n"))
      cat(paste0("  物种数量: ", length(node_data$species), "\n"))
      if(length(node_data$species) <= 10) {
        cat(paste0("  物种: ", paste(node_data$species, collapse=", "), "\n"))
      } else {
        cat(paste0("  物种: ", paste(node_data$species[1:5], collapse=", "), "... (等", 
                   length(node_data$species), "个物种)\n"))
      }
      
      cat(paste0("  融合事件: ", nrow(node_data$events$fusions), "\n"))
      cat(paste0("  裂变事件: ", nrow(node_data$events$fissions), "\n"))
      cat(paste0("  倒位事件: ", nrow(node_data$events$inversions), "\n"))
      cat(paste0("  保守共线性组: ", length(node_data$events$conserved), "\n"))
      cat("\n")
    }
  }
  
  # 添加执行信息
  cat("5. 执行信息\n")
  cat("-------------------------\n")
  cat(paste0("报告生成时间: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n"))
  cat(paste0("用户: ", Sys.info()["user"], "\n"))
  cat(paste0("R版本: ", R.version$version.string, "\n"))
  cat(paste0("操作系统: ", Sys.info()["sysname"], " ", Sys.info()["release"], "\n"))
  
  if(!is.null(analysis_results$config)) {
    cat("\n配置参数:\n")
    for(param_name in names(analysis_results$config)) {
      param_value <- analysis_results$config[[param_name]]
      if(!is.list(param_value) && !is.function(param_value)) {
        cat(paste0("  - ", param_name, ": ", param_value, "\n"))
      }
    }
  }
  
  # 关闭报告文件
  sink()
  
  log_message(paste("共线性分析报告已保存到:", report_file), "INFO", config)
}

# ==============================
# 辅助函数
# ==============================

get_tip_descendants <- function(tree, node) {
  # 获取内部节点的所有末端后代
  n_tips <- length(tree$tip.label)
  
  # 处理节点已经是末端的情况
  if(node <= n_tips) {
    return(node)
  }
  
  # 获取所有后代
  descendants <- c()
  to_process <- node
  
  while(length(to_process) > 0) {
    current <- to_process[1]
    to_process <- to_process[-1]
    
    # 查找子节点
    children_idx <- which(tree$edge[, 1] == current)
    if(length(children_idx) == 0) {
      next  # 没有子节点
    }
    
    children <- tree$edge[children_idx, 2]
    
    # 添加到队列或结果
    for(child in children) {
      if(child <= n_tips) {
        # 末端节点
        descendants <- c(descendants, child)
      } else {
        # 内部节点
        to_process <- c(to_process, child)
      }
    }
  }
  
  return(descendants)
}

get_descendants <- function(tree, node) {
  # 获取节点的直接后代
  children_idx <- which(tree$edge[, 1] == node)
  
  if(length(children_idx) == 0) {
    return(integer(0))  # 返回空的整数向量
  }
  
  children <- tree$edge[children_idx, 2]
  return(children)
}

# ==============================
# 主函数
# ==============================

main <- function(args) {
  # 解析命令行参数
  if(length(args) < 3) {
    cat("用法: Rscript ancestral_chromosomes_phase5_invertebrate.R <synteny_dir> <tree_file> <output_dir> [config_file]\n")
    cat("示例: Rscript ancestral_chromosomes_phase5_invertebrate.R synteny_files/ tree.nwk results/phase5\n")
    return(1)
  }
  
  synteny_dir <- args[1]
  tree_file <- args[2]
  output_dir <- args[3]
  config_file <- if(length(args) >= 4) args[4] else NULL
  
  # 加载配置
  config <- load_config(config_file)
  
  # 输出基本信息
  log_message(paste("开始无脊椎动物共线性数据处理"), "INFO", config)
  log_message(paste("共线性目录:", synteny_dir), "INFO", config)
  log_message(paste("进化树文件:", tree_file), "INFO", config)
  log_message(paste("输出目录:", output_dir), "INFO", config)
  if(!is.null(config_file)) log_message(paste("配置文件:", config_file), "INFO", config)
  
  # 处理共线性数据
  tryCatch({
    synteny_processed <- process_synteny_data(synteny_dir, tree_file, output_dir, config_file)
    
    # 检查是否有共线性数据
    if(length(synteny_processed$data$synteny_data) == 0) {
      log_message("未找到有效的共线性数据，处理终止", "ERROR", config)
      return(1)
    }
    
    # 提取共线性块
    log_message("开始提取共线性块...", "INFO", config)
    blocks <- extract_synteny_blocks(synteny_processed$data$synteny_data, config)
    
    # 检测事件
    log_message("开始检测共线性事件...", "INFO", config)
    events <- detect_synteny_events(
      synteny_processed$data$synteny_data, 
      chr_networks = NULL,  # 这里可能需要先构建网络
      tree = synteny_processed$data$tree, 
      config = config
    )
    
    # 映射事件到节点
    log_message("开始映射事件到祖先节点...", "INFO", config)
    node_mapping <- map_events_to_nodes(events, synteny_processed$data$tree, NULL, config)
    
    # 保存结果
    full_results <- list(
      synteny_data = synteny_processed$data,
      quality = synteny_processed$quality,
      blocks = blocks,
      events = events,
      node_mapping = node_mapping,
      config = config,
      processing_date = Sys.time()
    )
    
    results_file <- file.path(output_dir, "synteny_analysis_results.RData")
    saveRDS(full_results, results_file)
    log_message(paste("完整的共线性分析结果已保存到:", results_file), "INFO", config)
    
    # 生成报告
    generate_synteny_report(full_results, output_dir, config)
    
    # 生成可视化
    if(config$generate_visualizations) {
      log_message("开始生成可视化...", "INFO", config)
      # TODO: 添加更多可视化功能
      # 如 visualize_blocks, visualize_events, visualize_on_tree 等
    }
    
    log_message("共线性分析完成", "INFO", config)
    return(0)
    
  }, error = function(e) {
    log_message(paste("处理过程中发生错误:", e$message), "ERROR", config)
    return(1)
  })
}

# 如果脚本直接运行，执行main函数
if(!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  status <- main(args)
  quit(status = status)
}
