library(ppclust)
library(factoextra)
library(data.table)
library(tools)
library(gridExtra)
library(cluster)
library(plyr)
library(tidyverse)
library(dplyr)
library(psych)
library(ClustOfVar)


# Extracting the indices of metabolite/non-metabolite columns.
# A column is considered metabolite data if it starts with mzid.
extract_metab_indices <- function(df, pattern = '^mzid'){
  metab_indices <- grep(pattern, colnames(df)) %>% sort
  non_metab_indices <- 1:ncol(df) %>% setdiff(metab_indices)
  
  return(list(metab_indices, non_metab_indices))
  
}

# This function uses three different criteria to determine the optimal number clusters using three different methods.
# This is a wrapper function around fviz_nbclust from factoextra library.
# 
# Input:
# data: A tibble/dataframe where each row represents one metabolite or variable to be clustered.
# cluster_using_correlation: if TRUE, then the correlation values between different rows are calculated 
# and used for clustering where rows with higher absolute correlation are put in the same cluster. 
# Otherwise, rows are clustered using their values.
#
# Output: The output is a list and contains one plot per optimaly_methods.
optimal_cluster_nums <- function (data, cluster_using_correlation = TRUE,
                                  clustering_method = "kmeans", 
                                  optimality_methods = c("wss", "silhouette", "gap_stat") ){
  
  dist_mat <- NULL
  if (cluster_using_correlation){
    #data <- data %>% t()
    dist_mat <- (1 - abs( data %>% t() %>% cor(method = "pearson") )) %>% as.dist()
    
  }
  
  implemented_clustering_methods <- c("kmeans", "hclust")
  if (!clustering_method %in% implemented_clustering_methods){
    err_msg <- list()
    err_msg[[1]] <- sprintf('Clustering method %s is not implemented. The available methods are :\n', clustering_method)
    for (i in 1:length(implemented_clustering_methods) )
      err_msg[[i+1]] <- paste(implemented_clustering_methods[[i]], '\n')
    stop(paste(err_msg))
  }
  
  if (clustering_method == "kmeans"){
    clustering_func <- kmeans
  }else if (clustering_method == "hclust"){
    clustering_func <- function (x, k) {
      hcut(x, k, hc_method = "ward.D", metric = "manhattan")
    }
  }
  
  titles <- list("wss" = "Total Within Sum of Squares", 
                 "silhouette" = "Average Silhouette Width", 
                 "gap_stat" = "Gap Statistics")
  plots <- list()
  for (i in 1:length(optimality_methods)){
    cur_method <- optimality_methods[[i]]
    if (cur_method == "gap_stat"){
      plots[[i]] <- fviz_nbclust(data, FUNcluster = clustering_func, method = cur_method, nboot = 10, diss = dist_mat)
    }else{
      plots[[i]] <- fviz_nbclust(data, FUNcluster = clustering_func, method = cur_method, diss = dist_mat)
    }
    cur_title <- titles[[cur_method]]
    if (cluster_using_correlation){
      cur_title <- paste(cur_title, '\nClustered Using Coefficient values')
    }
    plots[[i]] <- plots[[i]] + ggtitle(cur_title) + theme(plot.title = element_text(hjust = 0.5))
  }
  
  return(plots)
}

# Same as optimal_cluster_nums function but it returns plots as a face grid of optimality methods on x-axis and 
# clustering using coefficient on the y-axis
optimal_cluster_nums_faceted <- function (data, clustering_method = "kmeans", 
                                          optimality_methods = c("wss", "silhouette", "gap_stat"),
                                          use_correlations_as_well = T,
                                          adjusting_covariates = NULL,
                                          kmax = 10, 
                                          do_not_grob = F
                                   ){
  

  
  implemented_clustering_methods <- c("kmeans", "hclust", "fuzzyKMeans", "hclust_partial", "varclus_hclust")
  if (!clustering_method %in% implemented_clustering_methods){
    err_msg <- list()
    err_msg[[1]] <- sprintf('Clustering method %s is not implemented. The available methods are :\n', clustering_method)
    for (i in 1:length(implemented_clustering_methods) )
      err_msg[[i+1]] <- paste(implemented_clustering_methods[[i]], '\n')
    stop(paste(err_msg))
  }
  
  if (clustering_method == "kmeans"){
    clustering_func <- clustering_func <- function (x, k) {
      kmeans(x, centers = k, iter.max = 100, nstart = 10)
    }
  }else if (clustering_method == "hclust"){
    clustering_func <- function (x, k) {
      hcut(x, k, hc_method = "ward.D2", hc_metric = "manhattan")
    }
  }else if(clustering_method == "fuzzyKMeans"){
      clustering_func <- function (x, k){
        fcm(x, k, m = 1.2)
      }
  }else if(clustering_method == "hclust_partial"){
        clustering_func <- function (x, k){
          
          x %<>% t()
          adjusting_covariates %<>% t()
          cor_mat <- partial.r(data = cbind(adjusting_covariates,x),
                               x = x %>% colnames(), 
                               y = adjusting_covariates %>% colnames()) %>% unclass()
          dist_cor_mat <- (1 - cor_mat^2) %>% as.dist()
          assignments <- dist_cor_mat %>% hclust(method = "ward.D") %>% cutree(k = k)
          list(cluster = assignments) %>% return()
          
    }
    
  }else if (clustering_method == "varclus_hclust"){
    clustering_func <- function (x, k) {
      varclus_tree <- hclustvar(X.quanti = (x %>% t) )
      cutreevar(varclus_tree, k)
    }
  }
  
  # titles <- list("wss" = "Total Within Sum of Squares", 
  #                "silhouette" = "Average Silhouette Width", 
  #                "gap_stat" = "Gap Statistics")
  if (use_correlations_as_well){
    use_correlation <- c(FALSE, TRUE)
  }else{
    use_correlation <- c(F)
  }
  plots <- list()
  counter <- 1
  for (i in 1:length(optimality_methods)){
    cur_method <- optimality_methods[[i]]
    for (j in 1:length(use_correlation)){
      dist_mat <- NULL
      if (use_correlation[[j]]){
        dist_mat <- (1 - abs( data %>% t() %>% cor(method = "pearson") )) %>% as.dist()
      }
      if (cur_method == "gap_stat"){
        plots[[counter]] <- fviz_nbclust(data, FUNcluster = clustering_func, method = cur_method, nboot = 10, diss = dist_mat, k.max = kmax)
      }else{
        plots[[counter]] <- fviz_nbclust(data, FUNcluster = clustering_func, method = cur_method, diss = dist_mat, k.max = kmax)
      }
      
      
      cur_title <- NULL
      if (use_correlation[[j]]){
        cur_title <- paste(cur_title, 'Using Correlations')
      }
      plots[[counter]] <- plots[[counter]] + ggtitle(cur_title) + theme(plot.title = element_text(hjust = 0.5))
      counter <- counter+1
      }
  }
  
  if (do_not_grob){
    return(plots)
  }
  
  if (length(optimality_methods) == 3){
  
  plots <- arrangeGrob(grobs = plots, 
                       ncol = 3, 
                       layout_matrix = cbind(c(1,2),c(3,4),c(5,6)))
  
  } else {
    if (use_correlations_as_well){
      plots <- arrangeGrob(grobs = plots, 
                           ncol = 2, 
                           layout_matrix = cbind(c(1,2),c(3,4)))
    }else{
      plots <- arrangeGrob(grobs = plots, 
                           nrow = 2)
    }
  }
                
  return(plots)
}

# Converting Tab-separated-value(TSV) to CSV.
tsv_to_csv <- function(input_file, output_file, replace = F){
  
  if (!file.exists(input_file)){
    stop(sprintf( 'Input File \n%s\n Does not exist.', input_file) )
  }
  fread(input_file) %>% write_csv(paste(tools::file_path_sans_ext(output_file), '.csv', sep = "") )
  if (replace){
    file.remove(input_file)
  }
  
}

# This function orders the color vector of a ggplot using the cluster centers 
# where the clusters are re-numbered from the lowset to highets of :
# 1- First cluster centers x values 
# 2- In case of a tie, cluster centers y-values
# The centers x and y values are bucketed so that the algorithm is robust to
# small changes in centers x-value. 
reorder_cluster_ids <- function(cur_plot, bucket_num = 5){
  
  
  # Extracting the mapping variables
  x_name <- cur_plot$layers[[1]]$mapping$x[[2]] %>% toString()
  y_name <- cur_plot$layers[[1]]$mapping$y[[2]] %>% toString()
  color_name <- cur_plot$layers[[1]]$mapping$colour[[2]] %>% toString()
  
  # Getting the centers
  centers <- cur_plot$data %>% 
    group_by(.dots=color_name) %>% 
    select(color_name, x_name, y_name) %>% 
    summarise_all(median)
  
  # Bucketing and sorting
  x_buckets <- centers %>% select(x_name) %>% as.matrix() %>% cut(bucket_num, include.lowest = T)
  y_buckets <- centers %>% select(y_name) %>% as.matrix() %>% cut(bucket_num, include.lowest = T)
  
  new_order <- order(x_buckets, y_buckets)
  
  # Changing the cluster ids
  cur_colors <- cur_plot$data[color_name] %>% as.matrix() %>% as.integer() %>% unique() %>% sort()
  cur_plot$data[color_name] <- cur_plot$data[color_name] %>% as.matrix() %>%
    mapvalues(from = cur_colors[new_order], to = cur_colors)
  return(cur_plot)
}

replace_mzid_with_compoundID <- function(data, eiconsanoids_ids, mzid_column_name){
  data <- data %>% dplyr::rename(mzid = mzid_column_name) %>% inner_join(eiconsanoids_ids, by = "mzid") %>% 
     mutate(mzid = compound_id) %>% select(-compound_id) %>% dplyr::rename(compound_id = mzid)
  return(data)
}

get_metab_data <- function(data, use_test_id = F){
  
  if (use_test_id){
    data %>% 
      select(test_id, site, timepoint, plate_nums, well_nums, starts_with("mzid")) %>%
      mutate(test_id = sprintf("ID_%d", test_id)) %>%
      unite("sample_id", test_id, site, timepoint, plate_nums, well_nums, sep = "_") %>%
      gather(key = "mzid", value = "metab_value", starts_with("mzid")) %>%
      spread(key = "sample_id", value = "metab_value")
  }else{
    data %>% 
      select(plate_nums, well_nums, starts_with("mzid")) %>%
      mutate(plate_nums = sprintf('P%02d', plate_nums), well_nums = sprintf('W%02d', well_nums) ) %>%
      unite("plate_well", plate_nums, well_nums, sep = "_") %>%
      gather(key = "mzid", value = "metab_value", starts_with("mzid")) %>%
      spread(key = "plate_well", value = "metab_value")
  }
    
}

# Helper function for mlr machine learning package to define hierarchical clustering function.
mlr_define_hclust <- function(){
  
  makeRLearner.cluster.hclust <<- function() {
    makeRLearnerCluster(
      cl = "cluster.hclust",
      package = "stats",
      par.set = makeParamSet(
        makeIntegerLearnerParam(id = "k"),
        makeDiscreteLearnerParam(id = "method", default = "complete",
                                 values = c("ward.D", "ward.D2", "single", "complete", 
                                            "average", "mcquitty", "median", "centroid")
                                 ),
        makeDiscreteLearnerParam(id = "metric", default = "euclidean",
                                 values = c("euclidean", "manhattan", "maximum", 
                                            "canberra", "binary", "minkowski")
        )
      ),
      properties = c("numerics"),
      name = "Hierarchical Clustering from R",
      short.name = "hclust-R"
    )
    
  }
    trainLearner.cluster.hclust <<- function (.learner, .task, .subset, .weights = NULL, ...) 
    {
      unnamed_arguments <- list(...)
      metric <- unnamed_arguments[["metric"]]
      if (is.null(metric)){
        metric <- .learner$par.set$pars$metric$default
      }
      
      method <- unnamed_arguments[["method"]]
      if (is.null(method)){
        method <- .learner$par.set$pars$method$default
      }
      
      
      dist_mat <- dist(getTaskData(.task, .subset), method = metric)
      hclust(dist_mat, method = method)
      
    }
    
    predictLearner.cluster.hclust <<- function (.learner, .model, .newdata, ...) 
    {
      k <- getHyperPars(.learner)[["k"]]

      cutree(.model$learner.model, k = k)
    }
  
}


# This function adds guids to a heatmap according to the cluster assignments.
# It is assumed that the cluster_indices are sorted based on the cluster sizes in a descending order.

heatmap_add_guides <- function(plot_obj, cluster_indices, line_color = "white", only_rect = F){
  
  # Finding where to add the guides. Since the cluster_indices is sorted then where the values drop, marks the boundaries of the clusters.
  # Each element is compared to the next to find where the boundaries are. To do so, a helper array that is shifted is used.
  right_shift <- cluster_indices %>% shift(n = 1)
  cluster_boundaries <- (right_shift != cluster_indices) %>% which()
  
  num_boundaries <- length(cluster_boundaries)
  cluster_boundaries <- cluster_boundaries - 1
  
  coord_min <- 1
  
  for (i in 1:num_boundaries){
    
    if (only_rect){
      coord_max <- cluster_boundaries[[i]]
      plot_obj <- plot_obj + geom_rect(xmin = coord_min, ymin = coord_min, xmax = coord_max, ymax = coord_max,
                                       color = line_color, linetype = "dashed", fill = NA)
      coord_min <- coord_max
    }else{
      plot_obj <- plot_obj + 
        geom_hline(yintercept = cluster_boundaries[[i]], color = line_color, linetype = "dashed") + 
        geom_vline(xintercept = cluster_boundaries[[i]], color = line_color, linetype = "dashed")
    }
    
  }
  
  if (only_rect){
    coord_max <- length(cluster_indices)
    plot_obj <- plot_obj + geom_rect(xmin = coord_min, ymin = coord_min, xmax = coord_max, ymax = coord_max,
                                     color = line_color, linetype = "dashed", fill = NA)
  }
  
  return(plot_obj)
  
  
}


# This function extracts metabolite columns that are missing less than metabolites_percent_missing
# from a dataframe/tibble.
extract_nonMissing_metabs <- function(data, threshold, metabolites_percent_missing, remove_factor_metabolites = T){
  all_indices <- data %>% extract_metab_indices()
  metab_indices <- all_indices[[1]]
  non_metab_indices <- all_indices[[2]]
  cur_metabs <- metabolites_percent_missing %>% filter(percent_missing < MISSING_THRESHOLD) %>% dplyr::select(mzid) %>% as.matrix() %>% dplyr::intersect(data %>% colnames)
  
  if (remove_factor_metabolites){
    cur_metabs <- data %>% dplyr::select(cur_metabs) %>% dplyr::select_if(function(x) !is.factor(x)) %>% colnames
  }
  
  data %<>% dplyr::select(non_metab_indices, cur_metabs)
}


# This function is intended to merge annotations from multiple metabolites where each metabolite
# can have multiple annotations. Also, it removes suffixes _duplicate_ and \\[M-H\\] and \\[M-H\\+Acetate\\].
prune_merge_compound_ids <- function(input_data, remove_mh = F, split_annotations = F){
  # Splitting the names using ";" delimiter and merging all into one vector
  if (split_annotations){
    input_data %<>% str_split(";") %>% unlist
  }
  # Removing _duplicate_(num)/(num) suffixes
  input_data %<>% str_replace("_duplicate_\\d+/\\d+", "")
  # Removing [M-H] and [M-H+Acetate] suffixes
  if (remove_mh){
    input_data %<>% str_replace(" \\[M-H\\]", "") %>% str_replace(" \\[M-H\\+Acetate\\]", "")
  }
   
  return(input_data)
  
}