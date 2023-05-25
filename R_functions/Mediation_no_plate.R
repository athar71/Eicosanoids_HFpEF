## 
library(tidyverse) # for analysis
library(survival) 
library(psych) # for skewness
library(e1071)          
library(magrittr)
library(doParallel)
library(mediation)
library(dplyr)
# Clean the environment
rm(list = ls())

stan_lgtrans_data <- read.csv("stan_lgtrans_data.csv", header =  TRUE) %>% as_tibble()
eicosanoids_ids <- read.csv("EIC_IDs_936.csv", header =  TRUE) %>% as_tibble()
stan_lgtrans_data <- stan_lgtrans_data[,-1]
eicosanoids_ids <- eicosanoids_ids[,-1]
eicosanoids_ids <- eicosanoids_ids[which(eicosanoids_ids$mzid %in% names(stan_lgtrans_data)),]

# Making them as factor 
stan_lgtrans_data$sex <- as.factor(stan_lgtrans_data$sex )
stan_lgtrans_data$diab_jeh <- as.factor(stan_lgtrans_data$diab_jeh )
stan_lgtrans_data$present_smoke <- as.factor(stan_lgtrans_data$present_smoke )
stan_lgtrans_data$mi <- as.factor(stan_lgtrans_data$mi )
stan_lgtrans_data$htn <- as.factor(stan_lgtrans_data$htn )
stan_lgtrans_data$plate_nums <- as.factor(stan_lgtrans_data$plate_nums)


# Standardizing the covariates
covariates <- c('bmi', 'age')
for (covar in covariates){
  stan_lgtrans_data[[covar]] %<>% scale(center = T, scale = T) %>% as.vector()
}

responses <- c( 'log_pap_co_slope_time',
                'log_pcw_co_slope_time' , 'peak_vo2', 'qc_mxm_hd_cavo2_peak', 
                'bp_response', 'vevco2nadir_864pull','sbp_r', 'pp_hr')

covariates <- c('bmi', 'age', 'sex', 'htn', 'mi', 'diab_jeh', 'present_smoke')
STANDARDIZE_OUTCOME <- T

if (STANDARDIZE_OUTCOME){
  
  for (i in seq_along(responses)){
    if (!(stan_lgtrans_data[[responses[[i]]]] %>% is.factor) ){
      stan_lgtrans_data[[ responses[[i]] ]] <- scale(stan_lgtrans_data[[ responses[[i]] ]], center = T, scale = T)
    }
  }
}

my_naming_func <- function(name, data){
  if (data[[name]] %>% is.factor){
    paste(name, "1", sep = "") %>% return()
  }else{
    name %>% return()
  }
}

models <- list()
total_counter <- 1

cur_traits <- c("hfpef_phys")
#controls <- c('plate_nums')


stan_lgtrans_data[['hfpef_phys']] %<>% as.integer()
#sprintf("Starting Mediation Analysis") %>% print
for (cur_trait in cur_traits){
  
  assign(cur_trait, stan_lgtrans_data[[cur_trait]])
  for (cur_covariate in covariates){
    cur_formula <- sprintf("%s ~ %s", cur_trait, cur_covariate) %>% as.formula()
    model.test <- glm(cur_formula, family = "binomial", data = stan_lgtrans_data)
    cov_adjusted_name <- my_naming_func(cur_covariate, stan_lgtrans_data)
    if (summary(model.test)$coefficients[,4][[cov_adjusted_name]] > 0.05){
      next
    }
    for (cur_eic in eicosanoids_ids$mzid){
      
      cur_formula <- sprintf("%s ~ %s", cur_eic, cur_covariate) %>% as.formula()
      model.M <- lm(cur_formula, stan_lgtrans_data)
      
      
      if (summary(model.M)$coefficients[,4][[cov_adjusted_name]] > 0.05){
        next
      }
      
      cur_formula <- sprintf("%s ~ %s + %s", cur_trait, cur_eic, cur_covariate) %>% as.formula()
      model.Y <- glm(cur_formula, family = "binomial", data = stan_lgtrans_data)
      
      if (summary(model.Y)$coefficients[,4][[cur_eic]] > 0.05){
        next
      }
      
      #sprintf("Starting X:%s M:%s Y:%s", cur_covariate, cur_eic, cur_trait) %>% print
      cur_model <- mediation::mediate(model.M, model.Y, treat=cur_covariate, mediator=cur_eic, robustSE = T, sims = 5000)
      
      
      
      
      models[[total_counter]] <- tibble("x:covariate" = cur_covariate, "m:EIC" = cur_eic, "y:trait" = cur_trait, pvalue = cur_model$n.avg.p, mediation_effect = cur_model$n.avg, 
                                        x_m_path = c(summary(model.M)$coefficients[,1][[cov_adjusted_name]], summary(model.M)$coefficients[,2][[cov_adjusted_name]], summary(model.M)$coefficients[,4][[cov_adjusted_name]]) %>% list,
                                        m_y_path = c(summary(model.Y)$coefficients[,1][[cur_eic]], summary(model.Y)$coefficients[,2][[cur_eic]],
                                                     summary(model.Y)$coefficients[,4][[cur_eic]]) %>% list,
                                        x_y_direct_path = c(summary(model.Y)$coefficients[,1][[cov_adjusted_name]], summary(model.Y)$coefficients[,2][[cov_adjusted_name]],
                                                            summary(model.Y)$coefficients[,4][[cov_adjusted_name]]) %>% list,
                                        x_y_total_path = c(summary(model.test)$coefficients[,1][[cov_adjusted_name]], summary(model.test)$coefficients[,2][[cov_adjusted_name]],
                                                           summary(model.test)$coefficients[,4][[cov_adjusted_name]]) %>% list
      )
      
      
      
      
      total_counter <- total_counter+1
      print(total_counter)
      print(cur_covariate)
      
    }
  }
}

models %<>% reduce(bind_rows) %>% inner_join(eicosanoids_ids, by = c("m:EIC" = "mzid")) %>% dplyr::select(`x:covariate`, compound_id, "m:EIC", `y:trait`, pvalue, mediation_effect,
                                                                                                          x_m_path, m_y_path, x_y_direct_path, x_y_total_path) %>% dplyr::rename(mzid = "m:EIC") %>% arrange(pvalue) %>% dplyr::rename(mediation_pvalue = pvalue)


combine_func <- function (x) {
  sprintf("%.2g +- %.2g, pval = %.2g", x[[1]], x[[2]]*1.96, x[[3]])
}

filter_func <- function(x_m_path, m_y_path, x_y_total_path){
  x_m_beta <- x_m_path %>% map_dbl(function (x) x[[1]])
  m_y_beta <- m_y_path %>% map_dbl(function (x) x[[1]])
  x_y_beta <- x_y_total_path %>% map_dbl(function (x) x[[1]])
  
  return((x_m_beta * m_y_beta * x_y_beta) >= 0)
}

# models %>% mutate(mediation_pvalue = sprintf("%.2g", mediation_pvalue), mediation_effect = sprintf("%.2f%%", mediation_effect*100)) %>% 
#   filter(filter_func(x_m_path = x_m_path, m_y_path = m_y_path, x_y_total_path = x_y_total_path)) %>% 
#   mutate(x_m_path = map_chr(x_m_path, combine_func),
#          m_y_path = map_chr(m_y_path, combine_func),
#          x_y_direct_path = map_chr(x_y_direct_path, combine_func),
#          x_y_total_path = map_chr(x_y_total_path, combine_func)
#   ) 

extract_func <- function(x, index){
  x %>% map_dbl(function(x) x[[index]])
}
l <- models %>% filter(filter_func(x_m_path = x_m_path, m_y_path = m_y_path, x_y_total_path = x_y_total_path)) %>%
  mutate(x_m_beta = extract_func(x_m_path, 1), x_m_se = extract_func(x_m_path, 2), x_m_pval = extract_func(x_m_path, 3),
         m_y_beta = extract_func(m_y_path, 1), m_y_se = extract_func(m_y_path, 2), m_y_pval = extract_func(m_y_path, 3),
         x_y_direct_beta = extract_func(x_y_direct_path, 1), x_y_direct_se = extract_func(x_y_direct_path, 2), x_y_direct_pval = extract_func(x_y_direct_path, 3),
         x_y_total_beta = extract_func(x_y_total_path, 1), x_y_total_se = extract_func(x_y_total_path, 2), x_y_total_pval = extract_func(x_y_total_path, 3)
  ) %>% dplyr::select(-x_m_path, -m_y_path, -x_y_direct_path, -x_y_total_path) %>%  write.csv(file =  'plate_mediation_analysis_hfpef1.csv' )


ll <- models  %>% mutate(x_m_beta = extract_func(x_m_path, 1), x_m_se = extract_func(x_m_path, 2), x_m_pval = extract_func(x_m_path, 3),
                         m_y_beta = extract_func(m_y_path, 1), m_y_se = extract_func(m_y_path, 2), m_y_pval = extract_func(m_y_path, 3),
                         x_y_direct_beta = extract_func(x_y_direct_path, 1), x_y_direct_se = extract_func(x_y_direct_path, 2), x_y_direct_pval = extract_func(x_y_direct_path, 3),
                         x_y_total_beta = extract_func(x_y_total_path, 1), x_y_total_se = extract_func(x_y_total_path, 2), x_y_total_pval = extract_func(x_y_total_path, 3)
) %>% dplyr::select(-x_m_path, -m_y_path, -x_y_direct_path, -x_y_total_path) #%>% write.csv(file =  'med_analysis_model_hfpef.csv' )

