library(matrixStats)
library(pracma)

my_join_columns <- function(mz, RT){
  return(sprintf("mzid_%f_%f", mz, RT))
}

my_scale_function <- function(data, scale, center){
  data %<>% as.matrix
  if (center){
    data <- data - (colMeans(data) %>% repmat(nrow(data), 1) )
  }

  if (scale){
    data <- data / (colSds(data) %>% repmat(nrow(data), 1) )
  }

  return(data)
}

# Setting the missing values to 25% of the minimum of the overall min for each metabolite
my_impute_missing_constant <- function(x, dichotomize_threshold){
  
  out <- zeros(dim(x)[1], dim(x)[2]) %>% as.data.frame()
  
  min_indices <- 1:ncol(x)
  # Dichotomizing metabolties that are missing above a certain threshold.
  if (dichotomize_threshold > 0 & dichotomize_threshold < 1){
    missing_count <- x %>% is.na() %>% colSums2()
    missing_percent <- missing_count / nrow(x)
    
    dich_indices <- missing_percent > dichotomize_threshold
    
    if ( sum(dich_indices) > 0 ){
      dich_mat <- x[,dich_indices]
      
      dich_mat[dich_mat > 0] <- 1
      dich_mat[dich_mat %>% is.na] <- 0
      
      out[,dich_indices] <- dich_mat
    }
    min_indices <- seq(ncol(x))[!dich_indices]
    
  }else{
    
    dich_indices <- c()
  }
  
  
  mins <- x[,min_indices] %>% as.matrix() %>% colMins(na.rm = T)
  mins <- mins*0.25
  list_mins <- mins %>% as.list() %>% setNames(colnames(x[,min_indices]))
  
  
  out[,min_indices] <- x[,min_indices] %>% as_tibble %>% replace_na(list_mins)
  
  return( list(out, min_indices) )
}


#imputing the missing values with samples from a normal distribution
my_impute_missing_normal <- function(x, dichotomize_threshold = .9, mean_imp_factor = .25, sd_imp_factor = .33){
  
  set.seed(1)
  out <- zeros(dim(x)[1], dim(x)[2]) %>% as.data.frame()
  
  min_indices <- 1:ncol(x)
  # Dichotomizing metabolties that are missing above a certain threshold.
  if (dichotomize_threshold > 0 & dichotomize_threshold < 1){
    missing_count <- x %>% is.na() %>% colSums2()
    missing_percent <- missing_count / nrow(x)
    
    dich_indices <- missing_percent > dichotomize_threshold
    
    if ( sum(dich_indices) > 0 ){
      dich_mat <- x[,dich_indices]
      
      dich_mat[dich_mat > 0] <- 1
      dich_mat[dich_mat %>% is.na] <- 0
      
      out[,dich_indices] <- dich_mat
    }
    min_indices <- seq(ncol(x))[!dich_indices]
    
  }else{
    
    dich_indices <- c()
  }
  
  cl <- parallel::makeCluster(6)
  doParallel::registerDoParallel(cl)
  min_out <- foreach(i = min_indices, .combine = 'bind_cols', .packages = c('tidyverse', 'magrittr')) %dopar% {
    cur_data <- x[,i]
    new_mean <- cur_data %>% mean(na.rm = T) * mean_imp_factor
    
    new_sd_1 <- ( ( cur_data %>% min(na.rm = T) ) - new_mean ) * sd_imp_factor
    # Making it highly unlikely to generate negative values
    new_sd_2 <- ( new_mean ) * sd_imp_factor
    new_sd <- min(new_sd_1, new_sd_2)
    
    cur_missing_ind <- cur_data %>% is.na %>% which
    
    cur_data[cur_missing_ind] <- rnorm(cur_missing_ind %>% length, mean = new_mean, sd = new_sd)
    tibble(!!colnames(x)[i] := cur_data)
  }
  parallel::stopCluster(cl)
  
  # Setting any (unlikely generated!) negative values to zero
  min_out[min_out < 0] <- 0
  
  out[,min_indices] <- min_out
  colnames(out)[min_indices] <- colnames(min_out)
  
  return( list(out, min_indices) )
}


# The helper function to perform 
# 1- Log transformation
# 2- Scaling/Standardizing
# 3- Norma/Constant Imputation
normalize_data <- function(df, metabolite_indices,
                           impute_missing = T, dichotomize_threshold = -1,
                           do_log_transform = T, do_scale = T, do_center = T,
                           impute_normal = F){
  
  data <- df[,metabolite_indices] %>% as.matrix
  if(do_log_transform){
    data %<>% log(base = exp(1))
  }
  
  if (impute_missing){
    if (impute_normal){
      out <- data %>% my_impute_missing_normal(dichotomize_threshold = dichotomize_threshold)
      data <- out[[1]]
      metabolite_indices <- metabolite_indices[out[[2]]]
    }else{
      out <- data %>% my_impute_missing_constant(dichotomize_threshold = dichotomize_threshold)
      data <- out[[1]]
      metabolite_indices <- metabolite_indices[out[[2]]]  
    }
    
  }
  
  if(do_scale || do_center){
    #data <- data %>% scale(scale = do_scale, center = do_center) %>% as_tibble
    data %<>% scale(scale = do_scale, center = do_center)
  }
  
  
  colnames(data) <-  colnames(df)[metabolite_indices]
  df[,metabolite_indices] <- data %>% as_tibble
  df[,metabolite_indices] %<>% map_df(as.double)
  return(df)
  
}
