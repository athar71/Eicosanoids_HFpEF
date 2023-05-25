read_CPET_data <- function(subset_metabs_name = "eicosanoids",
                           LOG_TRANSFORM_NTPROBNP = T,
                           sampling_sites = c("S"), sampling_timepoints = c(1),
                           LOG_TRANSFORM = T,
                           IMPUTE_MISSING = T,
                           IMPUTE_NORMAL = F,
                           SCALE_AND_CENTER = T,
                           NOT_TRANSFORM = F,
                           USE_NEW_DATASET = T,
                           MISSING_THRESHOLD = -1,
                           sub_subset_metabs = F,
                           subset_metab_subset = 1:50,
                           include_marathoners = F,
                           REMOVE_PHANTOM = T,
                           extra_annotations = F){
  
  library(data.table)
  library(magrittr)
  library(plyr)
  library(dplyr)
  library(tidyverse)
  library(lubridate)
  
  eicosanoids_path <- 'address'
  
  source(file.path(eicosanoids_path, 'R_scripts', 'normalize_data.R'))
  source(file.path(eicosanoids_path, 'R_scripts', 'eicosanoids_utils.R'))
  source(file.path(eicosanoids_path, 'R_scripts', 'get_CPET_filename.R'))
  
  input_file <- get_CPET_filename(eicosanoids_path = eicosanoids_path,
                                   log_transform = LOG_TRANSFORM,
                                   impute_missing = IMPUTE_MISSING,
                                   impute_normal = IMPUTE_NORMAL,
                                   scale_and_center = SCALE_AND_CENTER,
                                   use_new_dataset = USE_NEW_DATASET,
                                   extra_annotations = extra_annotations)
  
  load(input_file)
  
  if (NOT_TRANSFORM){
    CPET_data <- CPET_data_unTransformed
    patient_data <- patient_data_unTransformed
    plates_pooled_plasma <- plates_pooled_plasma_unTransformed
    meta_df <- meta_df_unTransformed
    internal_standard_data <- internal_standard_data_unTransformed
    metabolites_percent_missing <- metabolites_percent_missing_unTransformed
    patients_percent_missing <- patients_percent_missing_unTransformed
    pheno_data <- pheno_data_unTransformed
  }
  
  

  
  select <- dplyr::select
  rename <- dplyr::rename
  
  # Exctracting the Eicosanoids
  
  CPET_data <- CPET_data %>% inner_join(pheno_data, by = "test_id")
  if (include_marathoners){
    CPET_data %<>% bind_rows(patient_data %>% filter(str_detect(site, "^M")))
  }
  all_indices <- CPET_data %>% extract_metab_indices()
  metab_indices <- all_indices[[1]]
  non_metab_indices <- all_indices[[2]]
  if (subset_metabs_name != "all_metabolites"){
    if (subset_metabs_name == "all_named"){
      subset_ids <- eicosanoids_ids %>% bind_rows(other_named_ids)
    }else if (subset_metabs_name == "eicosanoids"){
      subset_ids <- eicosanoids_ids
    } else{
      subset_ids <- other_named_ids
    }
    filtered_data <- CPET_data %>% dplyr::select(non_metab_indices, subset_ids$mzid)
    
    if (sub_subset_metabs){
      cur_subset <- subset_ids$mzid[subset_metab_subset]
      
      filtered_data <- CPET_data %>% dplyr::select(non_metab_indices, cur_subset)
      subset_ids %<>% filter(mzid %in% cur_subset)
     
      if (subset_metabs_name == "all_named"){
        subset_ids <- eicosanoids_ids %>% bind_rows(other_named_ids)
      }else if (subset_metabs_name == "eicosanoids"){
        subset_ids <- eicosanoids_ids
      } else{
        subset_ids <- other_named_ids
      }
       
    }
  }else{
    filtered_data <- CPET_data %>% dplyr::select(non_metab_indices, metab_indices)
  }

  # Extracting the relevant timepoint and site
  if (include_marathoners){
    filtered_data %<>% filter(site %in% sampling_sites | str_detect(site, "^M"), timepoint %in% sampling_timepoints)
  }else{
    filtered_data %<>% filter(site %in% sampling_sites, timepoint %in% sampling_timepoints)  
  }
  
  
  # Removing the phantom samples
  if (REMOVE_PHANTOM){
    filtered_data <- filtered_data[!duplicated(filtered_data %>% select(test_id, site, timepoint)), ]
  }

  # Reading CPET dataset
  # CPET_exercise_data <- file.path(eicosanoids_path, 'Data/190313 CPET Full Dataset.csv') %>% 
  #   fread() %>% as_tibble() %>% 
  #   dplyr::rename(test_id = study_id) %>% select(test_id, brd_lvef_ex, brd_lvef_r, brd_rvef_ex, brd_rvef_r, pap_co_slope_time, pcw_co_slope_time,
  #                                                mxm_gx_vo2kg_rest, mxm_gx_vo2_peak, mxm_gx_vo2kg_peak, wt_kg, qc_mxm_hd_cavo2_peak, smoke, prev_hf_admit, mxm_gx_rer_peak, clin_hfpef, rer,
  #                                                brd_lvef_ex, brd_lvef_r,
  #                                                brd_rvef_ex , brd_rvef_r,
  #                                                mxm_gx_hr_peak, mxm_gx_hr_rest,
  #                                                hr_rest, hr_max,
  #                                                mxm_gx_sbp_peak, mxm_gx_sbp_rest,
  #                                                sbp_r, sbp_max,
  #                                                vevco2_slope
  #   )
  
  format_date_func <- Vectorize( function (x)  {ifelse(is.na(x) | x == "", NA, x %>% str_split(pattern = "/")  %>% 
      (function(y) c(paste0("20",y[[1]][[3]]), y[[1]][[1]],y[[1]][[2]]))  %>% 
      paste(collapse = "/")) }  %>% lubridate::as_date()) 
  CPET_exercise_data <- file.path(eicosanoids_path, 'Data/200210 ESL merge.csv') %>% 
    fread() %>% as_tibble() %>% 
    dplyr::rename(test_id = study_id) %>% select(test_id, brd_lvef_ex, brd_lvef_r, brd_rvef_ex, brd_rvef_r, pap_co_slope_time, pcw_co_slope_time,
            mxm_gx_vo2kg_rest, mxm_gx_vo2_peak, mxm_gx_vo2kg_peak, wt_kg, qc_mxm_hd_cavo2_peak, smoke, prev_hf_admit, mxm_gx_rer_peak, clin_hfpef, rer,
            brd_lvef_ex, brd_lvef_r,
            brd_rvef_ex , brd_rvef_r,
            mxm_gx_hr_rest, hr_max, 
            sbp_r, sbp_max,
            dbp_r, dbp_max,
            vevco2nadir_864pull,
            glucose, insulin,
            bb, ccb, acei, arb,
            date_test,
            cv_event_yes, cv_event_time,
            hf_event_yes, hf_event_time,
            vad_oht_event_yes, vad_oht_event_time,
            any_event_yes, any_event_time,
            death, death_time
            ) %>% mutate(homa_ir = (glucose*insulin)/405, ln_homa_ir = log(homa_ir, base = exp(1))) %>%
            mutate(date_test = date_test %>% format_date_func )
  
  # Removing people who did not exercise enough
  CPET_exercise_data %<>% filter(rer >= 1 | is.na(rer))
  
  # Correcting faulty datapoint
  CPET_exercise_data[CPET_exercise_data$test_id == 759, "brd_lvef_ex"] <- 66
  
  # Reading test ids which have passed exclusion criteria
  test_ids <- file.path(eicosanoids_path, 'Data/190327 PH ID List.csv') %>% fread() %>% as_tibble() %>% dplyr::rename(test_id = study_id)
  # fitting_data <- CPET_exercise_data %>% inner_join(test_ids, by = 'test_id') %>% inner_join(filtered_data, by = 'test_id') %>% 
  #   mutate(
  #     lv_systolic_reserve = brd_lvef_ex - brd_lvef_r,
  #     rv_systolic_reserve = brd_rvef_ex - brd_rvef_r,
  #     peak_vo2 = mxm_gx_vo2_peak / wt_kg,
  #     mxm_gx_vo2kg_delta = mxm_gx_vo2kg_peak - mxm_gx_vo2kg_rest,
  #     chrono_incomp = mxm_gx_hr_peak - mxm_gx_hr_rest,
  #     bp_response = mxm_gx_sbp_peak - mxm_gx_sbp_rest,
  #     pp_hr = (mxm_gx_hr_peak / (220 - age))*100
  #     
  #   )
  
  filtered_data <- CPET_exercise_data %>% inner_join(filtered_data, by = 'test_id') %>% 
    mutate(
      lv_systolic_reserve = brd_lvef_ex - brd_lvef_r,
      rv_systolic_reserve = brd_rvef_ex - brd_rvef_r,
      peak_vo2 = mxm_gx_vo2_peak / wt_kg,
      mxm_gx_vo2kg_delta = mxm_gx_vo2kg_peak - mxm_gx_vo2kg_rest,
      chrono_incomp = hr_max - mxm_gx_hr_rest,
      bp_response = sbp_max - sbp_r,
      pp_hr = (hr_max / (220 - age))*100
      
    )
  
  fitting_data <- filtered_data %>% inner_join(test_ids, by = 'test_id') 
  
  fitting_data %<>% mutate(present_smoke = case_when(
      (smoke == 0 | smoke == 1 | is.na(smoke)) ~ 0,
       smoke == 2                              ~ 1
    ),
    past_present_smoke = case_when(
      (smoke == 0 | is.na(smoke)) ~ 0,
      (smoke == 1 | smoke == 2)   ~ 1                           
    ),
    clin_hfpef         = case_when(
      (clin_hfpef == 2  | is.na(clin_hfpef) ) ~ 0,
       prev_hf_admit == 1                     ~ 1,
       (clin_hfpef == 1 | clin_hfpef == 0)    ~ clin_hfpef %>% as.double()
    ),
    htn_med = ifelse((bb == 1)|(ccb == 1)|(acei == 1)|(arb == 1), 1, 0) %>% as.integer,
    htn_new = ifelse((htn_med == 1)|(sbp_r >= 140 & !is.na(sbp_r))|(dbp_r >= 90 & !is.na(dbp_r)), 1, 0) %>% as.integer
  )
  
  fitting_data %<>% mutate(
    log_pap_co_slope_time = pap_co_slope_time %>% log(base = exp(1)),
    log_pcw_co_slope_time = pcw_co_slope_time %>% log(base = exp(1)),
    pap_co_slope_time_dict = (pap_co_slope_time >= 3) %>% as.logical() %>% as.integer() %>% as.factor(),
    pcw_co_slope_time_dict = (pcw_co_slope_time >= 2) %>% as.logical() %>% as.integer() %>% as.factor(),
    lv_systolic_reserve_dict = (lv_systolic_reserve > 0) %>% as.logical() %>% as.integer() %>% as.factor(),
    rv_systolic_reserve_dict = (rv_systolic_reserve > 0) %>% as.logical() %>% as.integer() %>% as.factor(),
    hfpef_phys = hfpef_phys %>% as.logical() %>% as.integer() %>% as.factor(),
    present_smoke = present_smoke %>% as.integer,
    age = age %>% as.double,
    obesity = (bmi > 30) %>% as.integer(),
    insulin_resist = (homa_ir > 2.0) %>% as.integer
  )
  
  # Removig the EICs that are missing too much.
  
  if (MISSING_THRESHOLD > 0){
    fitting_data %<>% extract_nonMissing_metabs(MISSING_THRESHOLD, metabolites_percent_missing)
    CPET_data %<>% extract_nonMissing_metabs(MISSING_THRESHOLD, metabolites_percent_missing)
    filtered_data %<>% extract_nonMissing_metabs(MISSING_THRESHOLD, metabolites_percent_missing)
    eicosanoids_ids %<>% filter(mzid %in% (fitting_data %>% select(starts_with("mzid")) %>% colnames()) )
  }
  
  return(list(CPET_data, filtered_data, fitting_data,
              eicosanoids_ids, other_named_ids, metabolites_percent_missing))
  
}
