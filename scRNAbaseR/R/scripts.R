#' Load Single File from featureCounts output
#' 
#' This function loads a single scRNA-seq / RNA-seq data file, removes duplicate rows, and renames the 7th column to 'counts'.
#' 
#' @param path A string specifying the path to the input file.
#' 
#' @return A data frame containing the processed data.
#' 
#' @examples
#' load_single_file("data/matrix.txt")
#' @export
load_single_file <- function(path) {
  set.seed(123)
  
  data <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  data <- dplyr::distinct(data)
  colnames(data)[7] <- 'counts'
  
  return(data)
}



#' Load scRNA-seq / RNA-seq Data
#'
#' This function processes multiple scRNA-seq / RNA-seq files, filters based on annotation types, and concatenates data as specified.
#'
#' @param path Path to the directory containing scRNA-seq / RNA-seq matrices.
#' @param concatenation_type Method to concatenate data: 'sum', 'avr', or 'top_diff'.
#' @param annotation_sides A vector of annotation sides to include (e.g., 'three_prime_UTR', 'exon', 'five_prime_UTR').
#'
#' @return A list of data frames, each corresponding to an annotation type.
#'
#' @examples
#' load_sc_rnaseq("data/matrices/", concatenation_type = 'sum')
#' @export
load_rnaseq <- function(path, concatenation_type = 'sum', annotation_sides = c('three_prime_UTR', 'exon', 'five_prime_UTR')) {
  set.seed(123)
  `%>%` <- dplyr::`%>%`
  files <- list.files(path = "results/matrices/", pattern = "*_genes_count_matrix.txt", full.names = TRUE)
  files <- files[!grepl(pattern = 'summary', files)]
  files <- files[grepl(pattern =  paste(annotation_sides, collapse = "|"), files)]
  
  basenames <- basename(files)
  
  annotation_types <- unique(gsub("^(([^_]*_){2}).*", "\\1", basenames))
  annotation_types <- sub("_$", "", annotation_types)
  
  
  rest_names <- unique(gsub("^[^_]*_[^_]*_", "", basenames))
  
  set_names <- unique(gsub(paste(annotation_sides, collapse = "|"), "", rest_names))
  set_names <- unique(gsub('_genes_count_matrix.txt', "", set_names))
  set_names <- sub("_", "", set_names)
  
  annotation_list <- list()
  
  for (at in annotation_types) {
    
    print(at)
    
    full_df <- list()
    
    for (s in set_names) {
      print(s)
      
      annp_list <- data.frame()
      
      for (as in annotation_sides) {
        
        regex <- paste0("(?=.*", paste(c(s,as,at), collapse = ")(?=.*"), ")")
        
        file_name <- files[grepl(pattern = regex, files, perl = TRUE)]
        
        file <- load_single_file(file_name)
        
        colnames(file)[1] <- 'name'
        
        file$annotation_sides <- as
        
        
        annp_list <- rbind(annp_list, file)
        
        
        
      }
      
      
      if (concatenation_type == 'sum') {
        
        res_df <- annp_list %>%
          dplyr::group_by(name) %>%
          dplyr::summarize(counts = sum(counts, na.rm = TRUE), .groups = 'drop') %>%
          dplyr::select(name, counts)
        
        res_df$counts <- as.numeric(res_df$counts)
        
        
        
      } else if (concatenation_type == 'avg') {
        
        res_df <- annp_list %>%
          dplyr::group_by(name) %>%
          dplyr::summarize(counts = mean(counts, na.rm = TRUE), .groups = 'drop') %>%
          dplyr::select(name, counts) %>%
          dplyr::distinct()
        
        res_df$counts <- as.numeric(res_df$counts)
        
        
      } else if (concatenation_type == 'top_diff') {
        
        res_df <- annp_list %>%
          dplyr::group_by(name) %>%
          dplyr::arrange(Length, .by_group = TRUE) %>%
          dplyr::slice_max(Length, n = 1) %>%
          dplyr::summarize(counts = mean(counts, na.rm = TRUE), .groups = 'drop') %>%
          dplyr::select(name, counts)
        
        res_df$counts <- as.numeric(res_df$counts)
        
      }
      
      
      full_df[[s]] <- res_df
      
    }
    
    merged_df <- NULL
    
    for (i in names(full_df)) {
      df <- as.data.frame(full_df[[i]])
      
      colnames(df)[2] <- i
      
      rownames(df) <- df$name
      df <- df[, -which(names(df) == "name"), drop = FALSE]
      
      if (is.null(merged_df)) {
        merged_df <- df
      } else {
        merged_df <- cbind(merged_df, df)
      }
    }
    
    
    merged_df[is.na(merged_df)] <- 0
    
    annotation_list[[at]] <- merged_df
    
  }
  
  return(annotation_list)
  
}



#' Count Genes
#'
#' This function calculates the number of expressed genes per sample / cell.
#'
#' @param df A data frame or matrix of gene expression counts, where rows represent genes and columns represent samples or cells.
#' @param count A character string specifying the counting method. Acceptable values are:
#'   \describe{
#'     \item{'counts'}{Calculate the total counts for each sample or cell (default).}
#'     \item{'genes'}{Count the number of expressed genes (non-zero counts) for each sample or cell.}
#'   }
#' @return A data frame with two columns:
#'   \describe{
#'     \item{set_name}{The names of the samples or cells.}
#'     \item{n}{The total counts or the number of expressed genes, depending on the `count` parameter.}
#'   }
#'   
#' @return A data frame with the count of expressed genes per sample.
#'
#' @examples
#' count_genes(df, count = 'counts')
#' @export
count_genes <- function(df, count = 'counts') {
  set.seed(123)
  
  
  if (count == 'counts') {
    
    sum_col <- colSums(df)
    result_df <- data.frame(
      set_name = names(sum_col),
      n = sum_col              
    )
    
    return(result_df)
    
    
  } else  if (count == 'genes') {
    
    df[df > 0] = 1
    
    sum_col <- colSums(df)
    result_df <- data.frame(
      set_name = names(sum_col),
      n = sum_col              
    )
    
    return(result_df)
    
  } else {
    
    
    print("Invalid `counts` parameter. The `counts` parameter should be included in either `counts` or `genes`.")
    
  } 
  
  
}


#' Cells Threshold
#'
#' This function calculates thresholds for the number of expressed genes per cell.
#'
#' @param data A data frame of gene counts.
#'
#' @return A list of thresholds at various percentiles.
#'
#' @examples
#' cells_threshold(data_frame)
#' @export
cells_threshold <- function(data) {
  set.seed(123)
  thresh_list <- list()
  
  thresh_list[['perc.10']] <- quantile(data$n, 0.10)
  thresh_list[['perc.25']] <- quantile(data$n, 0.25)
  thresh_list[['perc.33']] <- quantile(data$n, 0.33)
  thresh_list[['perc.50']] <- quantile(data$n, 0.50)
  thresh_list[['perc.66']] <- quantile(data$n, 0.66)
  thresh_list[['perc.75']] <- quantile(data$n, 0.75)
  thresh_list[['perc.90']] <- quantile(data$n, 0.90)
  
  return((thresh_list))
  
  
}



#' Genes per Cell / Set Plot
#'
#' This function generates a plot of the number of expressed genes per cell.
#'
#' @param data A data frame of gene counts.
#' @param cut_point Optional cutoff point for visualizing the threshold. Default is NaN. If NaN, the threshold is set on median.
#'
#' @return A ggplot object.
#'
#' @examples
#' genes_per_set_plot(data_frame)
#' @export
genes_per_cell_plot <- function(data, cut_point = NaN) {
  set.seed(123)
  library(ggplot2)
  theme_set(theme_bw())
  
  if (is.na(cut_point)) {
    cut_point <- quantile(genes$n, 0.5)
  }
  
  plot = ggplot2::ggplot(genes, aes(x = set_name, y = n)) +
    ggplot2::geom_point(size = 3) + 
    ggplot2::geom_segment(aes(x = set_name, 
                     xend = set_name, 
                     y = 0, 
                     yend = n)) + 
    ggplot2::geom_hline(yintercept = cut_point, linetype = "dashed", color = "red") + 
    ggplot2::labs(x = "Name", 
         y = "Genes") + 
    ggplot2::coord_flip() +  
    ggplot2::theme(axis.text.x = element_text(vjust = 0.6))
  
  return(plot)
}




#' Percent Features
#'
#' This function calculates the percentage of features / feature counts compared to the remaining features / feature counts in a set or cell.

#'
#' @param df A data frame of gene counts obtained from load_rnaseq()
#'
#' @return A data frame with the percent features / features counts.
#'
#' @examples
#' features_percentage(data, features_list)
#' @export
features_percentage <- function(data, features_list) {
  set.seed(123)
  
  sum_counts <- colSums(data)
  sum_features_counts <- colSums(data[toupper(rownames(data)) %in% features_list, ])
  
  data[data > 0] = 1
  
  sum_genes <- colSums(data)
  sum_genes_counts <- colSums(data[toupper(rownames(data)) %in% features_list, ])
  
  
  result_df <- data.frame(
    set_name = names(sum_counts),
    n_features_counts = sum_features_counts,   
    total_features_counts = sum_counts,
    n_features = sum_genes_counts,   
    total_features = sum_genes 
  )
  
  result_df$perc_features_counts <- round((result_df$n_features_counts / result_df$total_features_counts) * 100, 2)
  result_df$perc_features <- round((result_df$n_features / result_df$total_features) * 100, 2)
  
  
  return(result_df)
  
}






#' Features percent per Cell / Set Plot
#'
#' This function generates a plot of the percent of features / features counts per cell /set.
#'
#' @param data A data frame with percentage of features / features counts..
#' @param cut_point Optional cutoff point for visualizing the threshold. Default is 0.
#'
#' @return A ggplot object.
#'
#' @examples
#' features_perc_plot(data, count = 'counts', cut_point = 0)
#' @export
features_perc_plot <- function(data, count = 'counts', cut_point = 0) {
  set.seed(123)
  theme_set(theme_bw())
  
  if (count == 'counts') {
    
    plot = ggplot2::ggplot(data, aes(x = set_name, y = perc_features_counts)) +
      ggplot2::geom_point(size = 3, color = 'orange') + 
      ggplot2::geom_segment(aes(x = set_name, 
                       xend = set_name, 
                       y = 0, 
                       yend = perc_features_counts), color = 'orange') + 
      ggplot2::geom_hline(yintercept = cut_point, linetype = "dashed", color = "red") + 
      ggplot2::labs(x = "Name", 
           y = "Features counts [%]") + 
      ggplot2::coord_flip() +  
      ggplot2::theme(axis.text.x = element_text(vjust = 0.6))
    
  } else if (count == 'genes') {
    
    plot = ggplot2::ggplot(data, aes(x = set_name, y = perc_features)) +
      ggplot2::geom_point(size = 3, color = 'orange') + 
      ggplot2::geom_segment(aes(x = set_name, 
                       xend = set_name, 
                       y = 0, 
                       yend = perc_features), color = 'orange') + 
      ggplot2::geom_hline(yintercept = cut_point, linetype = "dashed", color = "red") + 
      ggplot2::labs(x = "Name", 
           y = "Features [%]") + 
      ggplot2::coord_flip() +  
      ggplot2::theme(axis.text.x = element_text(vjust = 0.6))
    
  }
  
  
  return(plot)
}








#' Select cells based on a gene count threshold
#' 
#' Filters cells in the dataset based on the minimum number of genes they express.
#' 
#' @param data A matrix or data frame of gene expression data, where rows represent genes and columns represent cells.
#' @param threshold Numeric value specifying the minimum number of genes required for a cell to be retained. Default is NaN.
#' @return A filtered dataset retaining only the cells meeting the threshold. If no threshold is provided, a message is printed, and no operation is performed.
#' @examples
#' filtered_data <- select_cells(data, threshold = 5000)
#' @export
select_cells <- function(data, threshold = NaN) {
  set.seed(123)
  if (!is.na(threshold)) {
    
    counts <- count_genes(data)
    counts <- counts[counts$n >= threshold, ]
    
    
    data = data[, counts$set_name]
    return(data)
    
    
  } else {
    
    print("Provide 'threshold' value!")
    
  }
  
  
}




#' Select Cells Based on Feature Percentage Data
#'
#' This function filters cells in a dataset based on the percentage occurrence of features 
#' or feature counts calculated using `features_percentage()`. It retains only the cells 
#' that meet a specified threshold.
#'
#' @param data A matrix or data frame of gene expression data, where rows represent genes 
#'   and columns represent cells.
#' @param features_perc A data frame containing feature percentage data, including a column 
#'   `set_name` for cell identifiers and `perc_features_counts` or `perc_features` for percentage values.
#' @param count A character string specifying the type of data to filter by. Acceptable values 
#'   are 'counts' (for filtering based on feature counts) or 'genes' (for filtering based on 
#'   gene occurrence percentages). Default is 'counts'.
#' @param threshold A numeric value specifying the maximum percentage threshold for a cell 
#'   to be retained. Only cells with a percentage below or equal to this threshold will be 
#'   included. If no threshold is provided (default is NaN), the function will print a message 
#'   and perform no operation.
#' @return A filtered dataset (matrix or data frame) containing only the cells that meet 
#'   the specified threshold criteria. If an invalid `counts` parameter is provided or 
#'   no threshold is specified, a message will be printed instead.
#' @examples
#' filtered_data <- select_on_features(data = reduced_data, features_perc, counts = 'counts', threshold = 20)
#' @export
select_on_features <- function(data, features_perc, count = 'counts', threshold = NaN) {
  set.seed(123)
  if (!is.na(threshold)) {
    
    if (count == 'counts') {
      
      data <- data[, toupper(colnames(data)) %in% toupper(features_perc$set_name[features_perc$perc_features_counts <= threshold])]
      return(data)
      
    } else  if (count == 'genes') {
      
      data <- data[, toupper(colnames(data)) %in% toupper(features_perc$set_name[features_perc$perc_features <= threshold])]
      return(data)
      
    } else {
      
      
      print("Invalid `count` parameter. The `count` parameter should be included in either `counts` or `genes`.")
      
    }
    
    
    
  } else {
    
    print("Provide 'threshold' value!")
    
  }
  
  
}



#' Normalize gene expression data
#' 
#' Applies log2 normalization to gene expression data based on gene counts and a scaling factor.
#' 
#' @param data A matrix or data frame of raw gene expression data.
#' @param type A type of normalization based on number of positive expressed genes or total counts per cell / set. Default is 'counts'
#' @param factor A numeric scaling factor for normalization. Default is 1,000,000.
#' @return A normalized dataset with log2-transformed values.
#' @examples
#' normalized_data <- normalize_data(data)
#' @export
normalize_data <- function(data, type = 'counts',  factor = 1000000) {
  set.seed(123)
  
  if (type %in% c('counts', 'genes')) {
    
    counts <- count_genes (df = data, count = type) 
    
    for  (r in row.names(counts)) {
      
      
      data[,r] =  log2(((data[,r]/as.integer(counts[r,]$n)) * factor) + 1) 
      
    }
    
  } else {
    
    print("Invalid `type` provided. The `type` parameter should be either 'counts' or 'genes'")
    
  }
  
  
  data[is.na(data)] <- 0
  
  return(data)
}



#' Calculate variance and occurrence of genes
#' 
#' Computes variance, mean, and percentage occurrence of each gene across samples.
#' 
#' @param data A matrix or data frame of gene expression data.
#' @param min Minimum threshold for occurrence percentage calculation. Default is 0.5.
#' @return A data frame containing variance, mean, and occurrence percentage for each gene.
#' @examples
#' gene_variance_data <- genes_variance(data)
#' @export
genes_variance <- function(data, min = 0.5) {
  set.seed(123)
  `%>%` <- dplyr::`%>%`
  gene_variance <- apply(data, 1, var)
  gene_mean <- apply(data, 1, mean)
  pct_occurrence <- apply(data, 1, function(row) mean(row > min) * 100)
  
  
  variance_df <- data.frame(gene = rownames(data), variance = gene_variance, avg = gene_mean, pct_occurrence = pct_occurrence)
  
  variance_df <- variance_df %>%
    dplyr::arrange(desc(variance))
  
  
  return(variance_df)
  
  
}





#' Plot variable genes
#' 
#' Generates a scatter plot of gene variance and mean, highlighting the most variable genes.
#' 
#' @param var_data A data frame of variance and mean values for genes.
#' @param side Character vector specifying the type of genes to plot ('equal', 'variable'). Default includes both.
#' @param n_top Integer specifying the number of top genes to label. Default is 20.
#' @return A ggplot2 scatter plot object.
#' @examples
#' var_gene_plot <- var_plot(var_data)
#' @export
var_plot <- function(var_data, side = c('equal', 'variable'), n_top = 20) {
  set.seed(123)
  `%>%` <- dplyr::`%>%`
  var_genes2 <- var_genes[order(var_genes$variance, var_genes$avg, decreasing = c(TRUE, FALSE)), ]
  
  var_genes2$color = 'bisque'
  
  subset1 <- var_genes2 %>%
    dplyr::filter(variance > quantile(variance, 0.9)) %>%
    dplyr::filter(avg > quantile(avg, 0.10))
  
  
  subset1 <- subset1[order(subset1$variance, subset1$avg, decreasing = TRUE), ]
  subset1$color <- 'black'
  
  var_genes2$color[rownames(var_genes2) %in% rownames(subset1)] <- 'coral'
  
  
  subset2 <- var_genes2 %>%
    dplyr::filter(avg > quantile(avg, 0.9)) %>%
    dplyr::filter(variance < quantile(variance, 0.10))
  
  subset2 <- subset2[order(subset2$avg, subset2$variance , decreasing = c(TRUE, FALSE)), ]
  subset2$color <- 'brown'
  
  var_genes2$color[rownames(var_genes2) %in% rownames(subset2)] <- 'gold'
  
  
  if ('equal' %in% side & 'variable' %in% side) {
    
    top_genes <- rbind(subset1[1:as.integer(n_top/2), ], subset2[1:as.integer(n_top/2), ])
    
  } else if ('equal' %in% side) {
    
    top_genes <- subset2[1:as.integer(n_top), ]
    
  } else if ('variable' %in% side) {
    
    top_genes <- subset1[1:as.integer(n_top), ]
    
  }
  
  
  
  plot <- ggplot2::ggplot(var_genes2, aes(x = variance, y = avg)) +
    ggplot2::geom_point(size = 3, color = var_genes2$color) + 
    ggrepel::geom_text_repel(data = top_genes, aes(label = gene), color = top_genes$color, max.overlaps = 50) +  
    ggplot2::labs(x = "var", y = "avg") +
    ggplot2::theme_minimal()
  
  
  return(plot)
  
}





#' Extract variable or equally expressed genes
#' 
#' Selects genes based on their variance and mean expression values.
#' 
#' @param var_data A data frame of gene variance and mean values.
#' @param side A string specifying which genes to extract ('equal' or 'variable'). Default is 'equal'.
#' @return A data frame of selected genes based on the specified criterion.
#' @examples
#' variable_genes <- get_var_genes(var_data, side = 'variable')
#' @export
get_var_genes <- function(var_data, side = 'equal') {
  set.seed(123)
  `%>%` <- dplyr::`%>%`
  if (!side %in% c('equal','variable')) {
    
    print('Side must be invluded in: "equal" or "variable"')
    return(NaN)
  }
  
  var_genes2 <- var_genes[order(var_genes$variance, var_genes$avg, decreasing = c(TRUE, FALSE)), ]
  
  variable <- var_genes2 %>%
    dplyr::filter(variance > quantile(variance, 0.75)) %>%
    dplyr::filter(avg > quantile(avg, 0.25))
  
  
  variable <- variable[order(variable$variance, variable$avg, decreasing = TRUE), ]
  
  
  
  equal <- var_genes2 %>%
    dplyr::filter(avg > quantile(avg, 0.75)) %>%
    dplyr::filter(variance < quantile(variance, 0.25))
  
  equal <- equal[order(equal$avg, equal$variance , decreasing = c(TRUE, FALSE)), ]
  
  if (side == 'equal') {
    
    return(equal)
    
    
  } else  if (side == 'variable') {
    
    return(variable)
    
    
  }
  
  
}





#' Perform PCA on gene expression data
#' 
#' Applies Principal Component Analysis (PCA) to identify the main components of variation in the data.
#' 
#' @param data A matrix or data frame of gene expression data.
#' @return An object of class 'prcomp' containing PCA results.
#' @examples
#' pca_result <- get_PCA(data)
#' @export
get_PCA <- function(data) {
  set.seed(123)
  
  vars <- get_var_genes(var_data = var_data, side = 'variable')
  
  data_scaled <- scale(t(reduced_data[vars$gene,]))
  
  pca_result <- prcomp(data_scaled, center = TRUE, scale. = TRUE)
  
  return(pca_result)
  
}



#' Create a knee plot for PCA results
#' 
#' Visualizes the explained variance for each principal component to identify the optimal number of components.
#' 
#' @param pca_data An object of class 'prcomp' containing PCA results.
#' @return A ggplot2 object showing the knee plot.
#' @examples
#' knee_plot <- get_knee(pca_result)
#' @export
get_knee <- function(pca_data) {
  set.seed(123)
  
  explained_variance <- (pca_data$sdev^2) / sum(pca_data$sdev^2) * 100
  cumulative_variance <- cumsum(explained_variance)
  components <- seq_along(explained_variance)
  
  # Przygotowanie ramki danych
  knee_data <- data.frame(
    Component = components,
    Variance = explained_variance,
    CumulativeVariance = cumulative_variance
  )
  
  
  
  plot <- ggplot2::ggplot(knee_data, aes(x = Component, y = Variance)) +
    ggplot2::geom_line(color = "blue", size = 1) +
    ggplot2::geom_point(size = 3, color = "red") +
    ggplot2::labs(x = "Principal Component",
         y = "Explained Variance (%)") +
    ggplot2::theme_minimal()
  
  return(plot)
}





#' Create a PCA scatter plot
#' 
#' Plots the first two principal components from PCA results.
#' 
#' @param pca_data An object of class 'prcomp' containing PCA results.
#' @return A ggplot2 scatter plot of the first two principal components.
#' @examples
#' pca_scatter <- pca_plot(pca_result)
#' @export
pca_plot <- function(pca_data) {
  set.seed(123)
  
  pca <- as.data.frame(pca_data$x[, 1:2])  
  pca$Sample <- rownames(pca_data$x) 
  colnames(pca) <- c("PC1", "PC2", "Sample")  
  
  plot <- ggplot2::ggplot(pca, aes(x = PC1, y = PC2)) +  
    ggplot2::geom_point(size = 3, color = 'brown') + 
    ggplot2::labs(x = "Principal Component 1", 
         y = "Principal Component 2") +                                
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none") 
  
  
  return(plot)
  
}





#' Perform UMAP and DBSCAN clustering
#' 
#' Combines UMAP dimensionality reduction with DBSCAN clustering for gene expression data.
#' 
#' @param pca_data An object of class 'prcomp' containing PCA results.
#' @param pc Number of principal components to use for UMAP. Default is 10.
#' @param eps The epsilon parameter for DBSCAN. Default is 0.5.
#' @param min_dist Minimum distance for UMAP. Default is 0.01.
#' @param n_neighbors Number of neighbors for UMAP. Default is 5.
#' @param minPts Minimum points for DBSCAN. Default is 5.
#' @return A ggplot2 UMAP plot with DBSCAN clusters.
#' @examples
#' umap_plot <- umap_db_scan(pca_result)
#' @export
umap_db_scan <- function(pca_data, pc = 10, eps = 0.5, min_dist = 0.01, n_neighbors = 5, minPts = 5) {
  set.seed(123)
  
  umap_config <- umap::umap.defaults
  umap_config$n_neighbors <- n_neighbors
  umap_config$min_dist <- min_dist
  umap_config$n_components <- 2
  
  
  umap_result <- umap::umap(pca_result$x[, 1:pc], config = umap_config) 
  
  umap_data <- as.data.frame(umap_result$layout)
  colnames(umap_data) <- c("UMAP1", "UMAP2")
  dbscan_result <- dbscan::dbscan(umap_data, eps = eps, minPts = n_neighbors) 
  
  umap_data$cluster <- as.factor(dbscan_result$cluster)
  
  plot <- ggplot2::ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = cluster)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(x = "UMAP Dimension 1", 
         y = "UMAP Dimension 2", 
         color = "Cluster") +
    ggplot2::theme_minimal()
  
  return(plot)
  
}





#' Retrieve clusters from UMAP and DBSCAN results
#' 
#' Extracts cluster assignments for each sample based on UMAP and DBSCAN results.
#' 
#' @param pca_data An object of class 'prcomp' containing PCA results.
#' @param pc Number of principal components to use for UMAP. Default is 10.
#' @param eps The epsilon parameter for DBSCAN. Default is 0.5.
#' @param min_dist Minimum distance for UMAP. Default is 0.01.
#' @param n_neighbors Number of neighbors for UMAP. Default is 5.
#' @param minPts Minimum points for DBSCAN. Default is 5.
#' @return A data frame mapping sample names to cluster assignments.
#' @examples
#' clusters <- get_clusters(pca_result)
#' @export
get_clusters <- function(pca_data, pc = 10, eps = 0.5, min_dist = 0.01, n_neighbors = 5, minPts = 5) {
  set.seed(123)
  
  umap_config <- umap::umap.defaults
  umap_config$n_neighbors <- n_neighbors
  umap_config$min_dist <- min_dist
  umap_config$n_components <- 2
  
  
  umap_result <- umap::umap(pca_result$x[, 1:pc], config = umap_config) 
  
  umap_data <- as.data.frame(umap_result$layout)
  colnames(umap_data) <- c("UMAP1", "UMAP2")
  dbscan_result <- dbscan::dbscan(umap_data, eps = eps, minPts = n_neighbors) 
  
  umap_data$cluster <- as.factor(dbscan_result$cluster)
  
  ret <- data.frame(name = rownames(umap_data), cluster = umap_data$cluster)
  
  return(ret)
  
}




#' Generate cluster statistics
#' 
#' Computes differential expression statistics for genes in each cluster.
#' 
#' @param data A matrix or data frame of gene expression data.
#' @param clusters A data frame of cluster assignments for each sample.
#' @param only_pos A logical value indicating whether to retain only positively differentially expressed genes (TRUE) or include both upregulated and downregulated genes (FALSE). Default is TRUE.
#' @param min_pct A numeric value specifying the minimum percentage of cells in a cluster that must express a gene for it to be considered. Default is 0.05.
#' @return A data frame of genes with log-fold change and adjusted p-values for each cluster.
#' @examples
#' cluster_stats <- get_cluster_stats(data, clusters, only_pos = TRUE, min_pct = 0.05)
#' @export
get_cluster_stats <- function(data, clusters, only_pos = TRUE, min_pct = 0.05) {
  set.seed(123)
  
  
  col_map <- setNames(clusters$cluster, clusters$name)
  
  colnames(data) <- col_map[colnames(data)]
  
  
  results <- data.frame()


  for (c in unique(clusters$cluster)) {


    cat(paste('\n\n Cluster ->  ', c, '- searching marker genes... \n\n' ))


    tmp1 = data[,colnames(data) %in% c]
    tmp_results <- data.frame(genes = rownames(tmp1))
    tmp2 = data[,!colnames(data) %in% c]
    tmp_sum <- tmp1
    tmp_sum[tmp_sum > 0] <- 1
    tmp_results$perc <- rowSums(tmp_sum)/ncol(tmp_sum)
    rm(tmp_sum)
    t1 = rowMeans(tmp1)
    t2 = rowMeans(tmp2)
    tmp_results$avg_logFC <- log2((t1 + (min(t1[t1 > 0])/2))  / (t2 + (min(t2[t2 > 0])/2)))
    tmp_results <- tmp_results[tmp_results$perc > min_pct, , drop = FALSE]



    if (only_pos == TRUE) {

      tmp_results <- tmp_results[tmp_results$avg_logFC > 0,]
      tmp_data <- data[match(tmp_results$genes, rownames(data)), ]

      tmp_results$p_val <- apply(tmp_data, 1, function(gene_expression) {
        wilcox.test(gene_expression[colnames(tmp_data) %in% c],
                    gene_expression[!colnames(tmp_data) %in% c])$p.value})

    } else {

      tmp_results$p_val <- apply(data, 1, function(gene_expression) {
        wilcox.test(gene_expression[colnames(data) %in% c],
                    gene_expression[!colnames(data) %in% c])$p.value})
    }


    tmp_results$cluster <- c
    results <- rbind(results, tmp_results)
  }


  
  return(results)
  
  
}




#' Plot heatmap for selected genes
#' 
#' Creates a heatmap for a subset of genes based on their expression values.
#' 
#' @param display_data A matrix or data frame of normalized gene expression data.
#' @param genes A vector of gene names to include in the heatmap.
#' @return A heatmap visualization of the selected genes.
#' @examples
#' heatmap_plot <- heatmap_plot(display_data, selected_genes)
#' @export
heatmap_plot <- function(display_data, features, features_metadata = NaN, min_value = 0.5) {
  set.seed(123)
  
  
  display_mean <- display_data
  
  display_mean[display_mean <= min_value] <- NaN
  
  mean <- as.data.frame(t(colMeans(display_mean, na.rm = TRUE)))
  
  rownames(mean)[1] <- 'expressed_mean'
  
  display_data <- display_data[rownames(display_data) %in% features,]
  
  
  if (!is.na(features_metadata[1])) {
    if (length(features) == length(features_metadata)) {
      
      meta = data.frame(name = features, type = features_metadata)
      meta$con <- paste0(meta$name, ' - ', meta$type)
      
      rownames(display_data) <- meta$con[match(rownames(display_data), meta$name)]
      
    } else {
      
      print('features_metadata is not equal to features; features_metadata will not be taken into consideration!')
      
    }
  }
  
  
  display_data <- rbind(display_data, mean)
  
  pheatmap::pheatmap(display_data, 
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           color = colorRampPalette(c("blue", "white", "red"))(100),  
           show_rownames = TRUE,  
           show_colnames = TRUE   
  )
  
  
}