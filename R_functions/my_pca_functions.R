#library(tidyimpute)
library(tidyverse)
library(pracma)
library(psych)

impute_mean <- function(x) {
  
  # if the vector supplied is not numeric...
  # return the original vector with a warning
  if (!is.numeric(x)) {
    warning("x is not numeric, returning original")
    return(x)
  }
  
  # this value will be used to replace NA's
  vector_mean <- mean(x, na.rm = TRUE)
  
  # read as: "replace x with the vector mean where x is NA"
  x[is.na(x)] <- vector_mean
  
  x
  
}


pca_plot_function <- function(data_processed, pov,
                              color_vector_name, is_color_factor,
                              output_POV_name, output_PC1PC2_name, add_output_path = T,
                              dot_size = 6, save_plots = T,
                              plot_pov = T,
                              PC_labels = c("PC1", "PC2"), pov_range = 1:10
                              
){
  
  plots <- list()
  
  outpath <- 'address'
  if (add_output_path){
    pov_figure_path <- file.path(outpath, output_POV_name)
  }else{
    pov_figure_path <- output_POV_name
  }
  
  plot_counter <- 1
  if (plot_pov){
    plots[[plot_counter]] <- pov %>% slice(pov_range) %>% ggplot() + geom_col(mapping = aes(x = pov_range, y = pov), fill="blue") + 
      labs(y = "Proportion of Variance", x = "PC Number") + 
      theme(text = element_text(size = 20)) + 
      scale_x_continuous(labels = function (x) sprintf("%.0f", x), breaks = pov_range) + 
      scale_y_continuous(labels = function (x) sprintf("%.0f%%", 100*x))
    plot_counter <- plot_counter + 1
    if (save_plots){
      plots[[1]] + ggsave(pov_figure_path, height = 6, width = 12)
    }
  }
  # Transformed data
  if (add_output_path){
    pc_figure_path <- file.path(outpath, output_PC1PC2_name)
  }else{
    pc_figure_path <- output_PC1PC2_name
  }
  
  
  if (is_color_factor){
    data_processed$color_vector <- factor(data_processed$color_vector)
  }
  
  plots[[plot_counter]] <- ggplot(data_processed) + 
    geom_point(mapping = aes_string(x = PC_labels[[1]], y = PC_labels[[2]], color="color_vector" ), size=dot_size) + 
    labs(color=color_vector_name) + theme(text = element_text(size=20))  

  if (save_plots){
    if (plot_pov){
      plots[[2]] + ggsave(pc_figure_path, height = 6, width = 12)
    }else{
      plots[[1]] + ggsave(pc_figure_path, height = 6, width = 12)
    }
  }
  
  return(plots)
}


pca_analysis_func <- function(data, metabolite_indices,
                              output_POV_name, output_PC1PC2_name,
                              color_vector, color_vector_name, is_color_factor,
                              keep_test_ids = NA,
                              dot_size = 6,
                              save_plots = T, add_output_path = T,
                              rotate_loadings = F,
                              adjusting_covariates = NA){
  
  
  # Removing columns with zero variances

  metab_data <- data %>% select(metabolite_indices)
  variances <- apply(metab_data, 2, var, na.rm=T)
  if (sum(variances == 0, na.rm = T) > 0){
    data_processed <- data %>% select(-metabolite_indices, 
                                      metabolite_indices[variances != 0]
    )
    metabolite_indices <- metabolite_indices[variances != 0]
  } else {
    data_processed <- data
  }
    
  metab_data <- metab_data %>% select( which(variances != 0) )
  
  # Replacing the missing values with mean
  #metab_data <- impute(metab_data, na.tools::na.mean)
  metab_data <- metab_data %>% mutate_all(impute_mean)
  

  if (rotate_loadings){
    # rawLoadings          <- data_processed.pca$rotation %*% diag(data_processed.pca$sdev)
    # rotatedLoadings      <- varimax(rawLoadings)$loadings
    # invLoadings          <- t(pracma::pinv(rotatedLoadings))
    # scores               <- scale(metab_data) %*% invLoadings
    # colnames(scores)     <- data_processed.pca$x %>% colnames
    # data_processed.pca$x <- scores
    
    pca_results <- metab_data %>% psych::principal(., nfactors = ncol(.), rotate = "varimax", scores = T)
    
    data_processed.pca <- list()
    data_processed.pca$sdev <- sqrt(pca_results$values)
    data_processed.pca$rotation <- pca_results$loadings %>% unclass()
    data_processed.pca$x <- pca_results$scores
    
    colnames(data_processed.pca$rotation) <- paste("PC", 1:ncol(metab_data), sep = "")
    colnames(data_processed.pca$x) <- paste("PC", 1:ncol(metab_data), sep = "")
    
    
  }else if(!is.na(adjusting_covariates)){
    cor_mat <- partial.r(data = cbind(adjusting_covariates,metab_data),
                         x = metab_data %>% colnames(), 
                         y = adjusting_covariates %>% colnames()) %>% unclass()
    eigen_adjusted <- cor_mat %>% eigen
    
    data_processed.pca <- list()
    data_processed.pca$sdev <- sqrt(eigen_adjusted$values)
    data_processed.pca$rotation <- eigen_adjusted$vectors
    data_processed.pca$x <- metab_data * data_processed.pca$rotation
    
    colnames(data_processed.pca$rotation) <- paste("PC", 1:ncol(metab_data), sep = "")
    colnames(data_processed.pca$x) <- paste("PC", 1:ncol(metab_data), sep = "")
    
  } else{
    data_processed.pca <- metab_data  %>% 
      prcomp(center = T, scale. = T)  
  }
  
  
  # Proportion of Variances
  pov <- as.tibble( data_processed.pca$sdev^2/sum(data_processed.pca$sdev^2) )
  colnames(pov) <- "pov"
  data_processed_names <- colnames(data_processed)
  data_processed <- cbind(data_processed, color_vector, data_processed.pca$x)
  
  
  
  data_processed_filtered <- data_processed
  if (!is.na(keep_test_ids)){
    data_processed_filtered <- data_processed %>% filter(test_id %in% keep_test_ids)
  }
  
  pca_plot_function(data_processed_filtered, pov,
                    color_vector_name, is_color_factor,
                    output_POV_name, output_PC1PC2_name, add_output_path = add_output_path,
                    dot_size = dot_size, 
                    save_plots = save_plots
  )
  
  pca_loadings <- tibble(mzid = colnames(data)[metabolite_indices] )  %>% cbind(data_processed.pca$rotation) %>% as.tibble()
  return (list(data_processed %>% as_tibble(), pov, pca_loadings))
}
