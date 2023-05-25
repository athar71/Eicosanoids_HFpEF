rm(list=ls(all.names=TRUE))
library(tidyverse)
library(data.table)
library(matrixStats)
library(readxl)
library(magrittr)
library(doParallel)

# Whether just run the preprocessing for a subset of metabolites. A positive number 
# SUBSET_METABOLITES will be randomly sampled. 
SUBSET_METABOLITES <- -1

# Use the dataset with more annotations.
USE_UPDATED_DATASET <- T

# Reading data from raw datasets as opposed to loading from .bin file.
# Most of the time it should be set to T
READ_FROM_CSV <- T
# Whether to log transform metabolite data.
LOG_TRANSFORM <- F #
# 
IMPUTE_MISSING <- F #
# Whether to Scale to have mean zero and SD of one.
SCALE_AND_CENTER <- F
# Whether to remove phantom samples. phantom samples are the exact samples that have 
# been given different study_ids and put on different mass-spec plates to check data quality.
# The metabolite information from these samples are expected to be almost the same or fairly close.
REMOVE_PHANTOM_SAMPLES <- T
LOG_TRANSFORM_NTPROBNP <- T

IMPUTE_NORMAL <- F

dichotomize_threshold <- -1

ADD_EXTRA_ANNOTATIONS <- T

NOT_TRANSFORM <- F
if (NOT_TRANSFORM){
  LOG_TRANSFORM <- T
  IMPUTE_MISSING <- F
  SCALE_AND_CENTER <- F
}

all_params <- list()

all_params[['USE_UPDATED_DATASET']] <- USE_UPDATED_DATASET
all_params[['READ_FROM_CSV']] <- READ_FROM_CSV
all_params[['LOG_TRANSFORM']] <- LOG_TRANSFORM
all_params[['IMPUTE_MISSING']] <- IMPUTE_MISSING
all_params[['SCALE_AND_CENTER']] <- SCALE_AND_CENTER
all_params[['REMOVE_PHANTOM_SAMPLES']] <- REMOVE_PHANTOM_SAMPLES
all_params[['LOG_TRANSFORM_NTPROBNP']] <- LOG_TRANSFORM_NTPROBNP
all_params[['IMPUTE_NORMAL']] <- IMPUTE_NORMAL

all_params[['IMPUTE_NORMAL']] <- IMPUTE_NORMAL
all_params[['dichotomize_threshold']] <- dichotomize_threshold
all_params[['NOT_TRANSFORM']] <- NOT_TRANSFORM

eic_path <- '/PHShome/aj438/eicosanoids/2018.10 Data Sharing/'

setwd(eic_path)
data_path <- file.path(eic_path, "Data")
source(file.path(eic_path, 'R_scripts', 'normalize_data.R'))
if (READ_FROM_CSV){
  if (USE_UPDATED_DATASET){
    input_data_file <- file.path(data_path, 'CPET1', 'ProcessedDataRawDeadducted_Labelled.csv')  
  }else{
    input_data_file <- file.path(data_path, '2018.10.19_CPET_Data_Extraction_Final.csv')  
  }
	
	input_key_file <- file.path(data_path, '2018.10.19_CPET_Sample_Key_Run_Order_cleaned.csv')
	input_pooled_plasma_file <- file.path(data_path, '2018.10.19_CPET_On-Plate_PooledPlasma.csv')
	input_internal_standard_file <- file.path(data_path, '2018.10.19_CPET_InternalStandards_CV_cleaned.csv')
	input_eicosanoids_file <- file.path(data_path, '2018.10.19_CPET_Data_Extraction_Final_Separate_Categories.xlsx')
	input_pheno_file <- file.path(data_path, 'HFpEF_pheno.csv')

	input_data <- as.tibble(fread(input_data_file))
	keys <- as.tibble(fread(input_key_file))
	pooled_plasma_raw_data <- as.tibble(fread(input_pooled_plasma_file, select = c(1, 2, 7:174)) )
	internal_standard_data <- as.tibble(fread(input_internal_standard_file, select = c(1, 5:5059)) )
	eicosanoids_data <- read_excel(input_eicosanoids_file, 5)
	pheno_data <- as.tibble(fread(input_pheno_file)) %>% dplyr::rename(test_id = study_id)
} else{
	load( file.path(data_path, './Input_data.bin') )
}

#######################################################################################
########################## Patient and CPET data preparation ########################## 
#######################################################################################

# Dropping duplicated columns
input_data <- distinct(input_data)

if (ADD_EXTRA_ANNOTATIONS){
  new_data <- file.path(data_path, 'ProcessedDataNormDeadductedV2_MzRtPeakId(20210120T1225).csv') %>% 
    fread %>% as_tibble
  input_data <- new_data %>% select(MZ, RT, Identity) %>% dplyr::rename(Identity_new = Identity) %>% 
    right_join(input_data, by = c("MZ", "RT")) %>% mutate(Identity = Identity_new) %>% select(-Identity_new)
}

if (SUBSET_METABOLITES > 0){
  rand_indices <- sample.int(nrow(input_data), SUBSET_METABOLITES)
  input_data <- input_data[rand_indices,]
}

# Extracting meta-data
if (USE_UPDATED_DATASET){
  metabolite_ids <- input_data %>% transmute(metabolite_ids = my_join_columns(MZ, RT))
  col_id <- grep('RT_Error', colnames(input_data))
  meta_df <- input_data %>% select(1:col_id)
  meta_df <- add_column(meta_df, metabolite_ids$metabolite_ids, .before='MZ_lower_lim_including')
  
  # Extracting metabolite information
  col_id <- grep('a180316_LCMS_EIC_CPET_Plate_02_02_A02_Samp_227460033_647_S_01_d', colnames(input_data))
  tmp_data <- colnames(input_data)[col_id:ncol(input_data)] %>% str_match('CPET_Plate_(\\d{2})_(\\d{2})')
  colnames(input_data)[col_id:ncol(input_data)] <- paste(tmp_data[,2], tmp_data[,3], sep = "_")
  data <- input_data %>% select(col_id:ncol(input_data))
  # Transposing the data to have plate ids as rows
  data <- add_column(data, metabolite_ids$metabolite_ids, .before='02_02')
}else{
  metabolite_ids <- input_data %>% transmute(metabolite_ids = my_join_columns(`m/z`, RT))
  col_id <- grep('Compound IDs', colnames(input_data))
  meta_df <- input_data %>% select(1:col_id)
  meta_df <- add_column(meta_df, metabolite_ids$metabolite_ids, .before='Method')
  
  # Extracting metabolite information
  col_id <- grep('01_02', colnames(input_data)) 
  data <- input_data %>% select(col_id:ncol(input_data))
  # Transposing the data to have plate ids as rows
  data <- add_column(data, metabolite_ids$metabolite_ids, .before='01_02')
}



data %<>%  gather(key = "plate_well", value = "value", 2:ncol(data)) %>%
  spread(`metabolite_ids$metabolite_ids`, value) %>%
  dplyr::rename("EIC ID" = "plate_well")

#Merging the key ids
patient_data <- data %>% left_join(keys, by="EIC ID")
patient_data %<>% select(`Sample ID`, 1:(ncol(patient_data)-1))

#Separating the EIC ID to plate and well num
patient_data %<>% separate("EIC ID", 
                           into=c("plate_nums", "well_nums"),
                           convert=T
                                    )

#Seperating the test ids(patien id), site and the time point(rest=0, peak=1 and recovery=2)
#site: (S=superior vena cava, P=pulmonary artery, A=arterial)
#time point(rest=1, peak=2 and recovery=3)

patient_data %<>% separate("Sample ID", 
                           into=c("test_id", "site", "timepoint"), 
                           fill="left", convert = TRUE
                           )

# Phantom samples are defined as the sample points having the same tiple of test_id, site and timepoint.
# The phantom samples are intended for data QC and should not be included in the modelling stage.
if (REMOVE_PHANTOM_SAMPLES){
  # Removing the phantom samples
  patient_data <- patient_data[!duplicated(patient_data %>% select(test_id, site, timepoint)), ]
}

#Extracting the CPET data
#CPET data set has test id between 1-999 and the site being either S/P/A
CPET_data <-  patient_data %>% filter(
                    test_id > 0 &
                    test_id < 1000 &
                    site %in% c('A', 'S', 'P') )


# Calculating Percent missing of each metabolite
#tmp_data <- patient_data %>% filter(timepoint == 1, site == "S")
tmp_data <- patient_data 
cur_nrow <- tmp_data %>% nrow
metabolites_percent_missing <- 
  tmp_data %>% select(starts_with('mzid')) %>% 
  is.na() %>% colSums() %>% 
  ( function(x) {tibble(mzid = names(x), percent_missing = x/cur_nrow * 100)} )


patients_percent_missing <- tmp_data %>% select(starts_with('mzid')) %>% is.na() %>% rowSums2()
patients_percent_missing <- tibble(test_id = tmp_data$test_id, 
                                            percent_missing = patients_percent_missing/(tmp_data %>% select(starts_with("mzid")) %>% ncol)*100 )

####################################################################################
########################## Pooled plasma data preparation ##########################
####################################################################################

# There are 3 spots per plate for each metabolite. 
# The spots are :  1, 48 and 96

#Dropping duplicated values
pooled_plasma_raw_data <- distinct(pooled_plasma_raw_data)
#Adjusting the column names to contain only the first portion

colnames(pooled_plasma_raw_data)[3:dim(pooled_plasma_raw_data)[2]] <- 
  sapply(colnames(pooled_plasma_raw_data)[3:dim(pooled_plasma_raw_data)[2]], function(x) substr(x, 1, 5))

#Creating metabolite identifiers and adding to the data
pooled_plasma_metabolite_ids <-  pooled_plasma_raw_data %>% transmute(metabolite_ids = my_join_columns(`mz`, RT))
pooled_plasma_raw_data <- add_column(pooled_plasma_raw_data, pooled_plasma_metabolite_ids$metabolite_ids, .after = 'RT')
pooled_plasma_raw_data$mz <- NULL
pooled_plasma_raw_data$RT <- NULL
#transposing the data
plates_pooled_plasma <- pooled_plasma_raw_data %>% 
  gather(key = "plate_well", value = "value", 2:ncol(.)) %>%
  spread(`pooled_plasma_metabolite_ids$metabolite_ids`, value) %>%
  dplyr::rename("EIC ID" = "plate_well")


plates_pooled_plasma <- plates_pooled_plasma %>% separate("EIC ID", c("plate_nums", "well_nums"), convert = T)

# Removing samples where the well_nums is not in c(1, 48, 96)
plates_pooled_plasma <- plates_pooled_plasma %>% filter(well_nums %in% c(1, 48, 96))

########################################################################################
########################## Internal Standard Data Preparation ########################## 
########################################################################################

internal_standard_data <- internal_standard_data %>% 
  dplyr::rename(metabolite_name = "row identity (main ID)")
#Transposing the data
internal_standard_data <- internal_standard_data %>% 
                          gather(key = "plate_well", value = "value", 2:ncol(.)) %>%
                          spread(key = metabolite_name, value = value)


# Separating the first column into well, plate and patient information
internal_standard_data <- internal_standard_data %>%  
                          separate(col = plate_well, 
                                   into=c("plate_nums", "well_nums", 
                                          "Placement", "Samp", "SampID",
                                          "test_ids", "site", "timepoint"),
                                   fill = "right", convert = T)
# Dropping unrelated columns
internal_standard_data <- internal_standard_data %>% 
                          select(-c(Placement, Samp, SampID))


index_range <- 6:ncol(internal_standard_data)
# Replacing 0 values with NA
internal_standard_data <- internal_standard_data %>% 
  mutate_at(vars(index_range), funs(replace(., . == 0, NA)))

internal_standard_data <- internal_standard_data %>% dplyr::rename(test_id = test_ids) %>% 
  mutate(test_id = as.integer(test_id))

####################################################################
########################## Final touches! ########################## 
####################################################################

## Removing the metabolites that are missing all their data
tmp_data <- CPET_data 
cur_nrow <- tmp_data %>% nrow
CPET_metabolites_percent_missing <- 
  tmp_data %>% select(starts_with('mzid')) %>% 
  is.na() %>% colSums() %>% 
  ( function(x) {tibble(mzid = names(x), percent_missing = x/cur_nrow * 100)} )

all_missing_mzid <- CPET_metabolites_percent_missing %>% filter(percent_missing == 100) %>% 
  select(mzid) %>% as.matrix 

if ( (all_missing_mzid %>% nrow) > 0 ){
  CPET_data %<>% select(-one_of(all_missing_mzid %>% as.vector()))
}



# test <- CPET_data
# CPET_data <- test[,1:6006]
# index_range <- 6:ncol(CPET_data)
# 
# df = CPET_data
# metabolite_indices = index_range
# impute_missing = IMPUTE_MISSING
# do_log_transform = LOG_TRANSFORM
# do_scale = SCALE_AND_CENTER
# do_center = SCALE_AND_CENTER
# dichotomize_threshold = dichotomize_threshold
# impute_normal = IMPUTE_NORMAL

index_range <- 6:ncol(CPET_data)
CPET_data <- normalize_data(df = CPET_data, metabolite_indices = index_range,
                            impute_missing = IMPUTE_MISSING,
                            do_log_transform = LOG_TRANSFORM, 
                            do_scale = SCALE_AND_CENTER, do_center = SCALE_AND_CENTER,
                            dichotomize_threshold = dichotomize_threshold,
                            impute_normal = IMPUTE_NORMAL)

index_range <- 6:ncol(patient_data)
patient_data <- normalize_data(patient_data, index_range,
                               impute_missing = IMPUTE_MISSING,
                               do_log_transform = LOG_TRANSFORM, 
                               do_scale = SCALE_AND_CENTER, do_center = SCALE_AND_CENTER,
                               dichotomize_threshold = dichotomize_threshold,
                               impute_normal = IMPUTE_NORMAL)

index_range <- 3:ncol(plates_pooled_plasma)
plates_pooled_plasma <- normalize_data(plates_pooled_plasma, index_range,
                                       impute_missing = IMPUTE_MISSING,
                                       do_log_transform = LOG_TRANSFORM, 
                                       do_scale = SCALE_AND_CENTER, do_center = SCALE_AND_CENTER,
                                       dichotomize_threshold = dichotomize_threshold,
                                       impute_normal = IMPUTE_NORMAL)

index_range <- 6:ncol(internal_standard_data)
internal_standard_data <- normalize_data(internal_standard_data, index_range,
                                         impute_missing = IMPUTE_MISSING,
                                         do_log_transform = F, 
                                         do_scale = F, do_center = F,
                                         dichotomize_threshold = dichotomize_threshold,
                                         impute_normal = IMPUTE_NORMAL)

################################################################################
########################## Extracting Eicosanoids IDs ##########################
################################################################################

if (USE_UPDATED_DATASET){
  tmp <- input_data %>% select(Identity, MZ, RT) %>% mutate(Identity = Identity %>% str_replace("\xa0", "")) %>% 
    filter( (Identity != "") & !(str_detect(Identity, "Int_Std_")) )
  
  eic_ids <- tmp %>% filter(str_detect(Identity, "EIC_") | str_detect(Identity, "Eicosanoid_"))
  eicosanoids_ids  <-  tibble(mzid = my_join_columns(eic_ids$MZ, eic_ids$RT), compound_id = eic_ids$Identity )
  
  eicosanoids_ids %<>% group_by(compound_id) %>% dplyr::mutate(counts = n(), num = row_number()) %>% ungroup() %>% 
    mutate(new_ids = sprintf("%s_duplicate_%d/%d", compound_id, num, counts), 
           compound_id = if_else(counts > 1, new_ids, compound_id)) %>% select(-counts, -num, -new_ids)
  
  other_ids <- tmp %>% filter(!str_detect(Identity, "EIC_") & !str_detect(Identity, "Eicosanoid_") )
  other_named_ids  <-  tibble(mzid = my_join_columns(other_ids$MZ, other_ids$RT), compound_id = other_ids$Identity )
  

  other_named_ids %<>% group_by(compound_id) %>% dplyr::mutate(counts = n(), num = row_number()) %>% ungroup() %>% 
    mutate(new_ids = sprintf("%s_duplicate_%d/%d", compound_id, num, counts), 
           compound_id = if_else(counts > 1, new_ids, compound_id)) %>% select(-counts, -num, -new_ids)
}else{
  eicosanoids_ids <- my_join_columns(eicosanoids_data$`m/z`, eicosanoids_data$RT) %>% 
    tibble(mzid = ., compound_id = eicosanoids_data$`Compound IDs`)  
}



# Natural log transforming ntprobnp
if (LOG_TRANSFORM_NTPROBNP){
  pheno_data %<>% mutate(log_ntprobnp = ntprobnp %>% log(base = exp(1)) )
}

outname <- 'processed_data'

if (LOG_TRANSFORM){
  outname %<>% paste0('_logTransformed')
}

if (IMPUTE_MISSING){
  if (IMPUTE_NORMAL){
    outname %<>% paste0('_normal_imputed')
  }else{
  outname %<>% paste0('_constant_imputed')
  }
}

if(SCALE_AND_CENTER){
  outname %<>% paste0('_scaled')
}

if (USE_UPDATED_DATASET){
  outname %<>% paste0('_new_dataset')
}

if(ADD_EXTRA_ANNOTATIONS){
  outname %<>% paste0('_extraAnnotations')
}

outname <-  file.path(data_path, paste0(outname, '.bin'))

# Saving the workspace
if (NOT_TRANSFORM){
  
  CPET_data_unTransformed <- CPET_data
  patient_data_unTransformed <- patient_data
  plates_pooled_plasma_unTransformed <- plates_pooled_plasma
  meta_df_unTransformed <- meta_df
  internal_standard_data_unTransformed <- internal_standard_data
  metabolites_percent_missing_unTransformed <- metabolites_percent_missing
  patients_percent_missing_unTransformed <- patients_percent_missing
  pheno_data_unTransformed <- pheno_data
  
  rm(CPET_data, patient_data, plates_pooled_plasma, meta_df, internal_standard_data, 
     metabolites_percent_missing, pheno_data)
  
  save(CPET_data_unTransformed, patient_data_unTransformed, plates_pooled_plasma_unTransformed, 
       meta_df_unTransformed, internal_standard_data_unTransformed, 
       eicosanoids_ids, other_named_ids, metabolites_percent_missing_unTransformed, patients_percent_missing_unTransformed, pheno_data_unTransformed,
       all_params,
       file = outname
  )
} else{
  

  save(CPET_data, patient_data, plates_pooled_plasma, meta_df, internal_standard_data, 
       eicosanoids_ids, other_named_ids, metabolites_percent_missing, patients_percent_missing, pheno_data, 
       all_params,
       file = outname
  )
}



