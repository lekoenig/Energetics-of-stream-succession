
## Project: Semi-mechanistic modeling of river metabolism recovery following storm disturbances
## Script: Filter Powell Center sites
## LE Koenig
## last updated June 2021

## The objective of this script is to filter sites from the Powell Center database of stream metabolism time series to use in the stream successional energetics project.


## Load packages: 
library(dplyr)         # general data cleaning and manipulation
library(tidyr)         # general data cleaning and manipulation
library(purrr)         # general data cleaning and manipulation
library(lubridate)     # format timestamps
library(dataRetrieval) # interface with NWIS
library(sbtools)       # interface with Science Base

source("./R/Analysis_Functions.R")

# define where to save filtered Powell Center sites:
in_save_dir <- "./data/in/Appling_PC_data"
out_save_dir <- "./data/out"

##===================================================================##
##                  Read in Powell Center data                       ##
##===================================================================##

# Read in NWIS sites from Powell Center data set (Appling et al. 2018):
pc_sites <- download_site_data(in_save_dir)

# Read in Powell Center data (metabolism estimates and predictors) - exports file called 'daily_predictions.csv':
metab_all <- download_metab_est(in_save_dir) 

# Read in Powell Center data (model diagnostics) - exports file called 'diagnostics.csv':
diagnostics <- download_model_diagnostics(in_save_dir)

##===================================================================##
##                  Filter Powell Center sites                       ##
##===================================================================##

# Function to calculate pearson correlation coefficient between ER and K600 for each df within the list:
ER_K_corr <- function(df){
  out <- split(x=df,f = df$site_name) %>%
    map(~ cor(x=.$ER,y=.$K600,method="pearson")) %>%
    data.frame %>% gather(.,site_name,cor.coef)
}

# Function to calculate median growing season GPP for each df within the list:
est_med_seasonal_GPP <- function(df){
  out <- split(x=df,f = df$site_name) %>%
    map(~ grow.season.GPP(.)) %>%
    data.frame %>% gather(.,site_name,GPP_grow_median)
}

# Filter data based on defined criteria:
metab_filter <- metab_all %>%
                # 1. select only those sites where we have high confidence in the model:
                  left_join(.,diagnostics[,c("site","model_confidence")],by=c("site_name"="site")) %>%
                  filter(.,model_confidence=="H") %>%
                # 2. select sites where the correlation between ER and K600 < 0.4:
                  left_join(.,ER_K_corr(.),by="site_name") %>% 
                  filter(abs(cor.coef) < 0.4) %>%
                # 3. Filter for high median GPP during growing season (120-d continuous period of maximum productivity):
                  left_join(.,est_med_seasonal_GPP(.),by="site_name") %>%
                  filter(GPP_grow_median > 3) 
                  
# Inspect filtered correlation coefficients and growing season GPP:
quantile(unique(metab_filter$cor.coef))
quantile(unique(metab_filter$GPP_grow_median))

# Find out whether sites are co-located with turbidity sensors:
metab_filtered_sites <- bind_rows(lapply(unique(metab_filter$site_name),check_turbidity_data)) %>%
                        mutate(turb_interval = interval(start=turbidity_start,end = turbidity_end)) %>%
                        # create new columns that indicate metab dates:
                        left_join(.,bind_rows(lapply(.$site_name,find_metab_dates)),by="site_name") %>%
                        mutate(metab_interval = interval(start=metab_start_date,end=metab_end_date)) %>%
                        # check whether the metabolism and turbidity date ranges overlap:
                        mutate(check_turbidity_range = int_overlaps(turb_interval,metab_interval)) %>%
                        # filter sites that overlap with turbidity data (at least 180 days):
                        filter(.,turbidity_total_days > 180 & check_turbidity_range == "TRUE") %>% 
                        # bind with site info from Science Base and select columns:
                        left_join(.,pc_sites,by="site_name") %>%
                        select(site_name,long_name,nhdplus_id,drain_area_km2,lat,lon,coord_datum,alt,alt_datum,site_type,dvqcoefs.c,dvqcoefs.f,
                               dvqcoefs.a,dvqcoefs.b,dvqcoefs.k,dvqcoefs.m,struct.canal_flag,struct.dam_flag,turbidity_data,turbidity_start,
                               turbidity_end,turbidity_total_days,metab_start_date,metab_end_date,metab_days,check_turbidity_range)

##===================================================================##
##                      Export filtered data                         ##
##===================================================================##

# Export filtered sites:
write.csv(metab_filtered_sites,paste(out_save_dir,"/site_data_filtered.csv",sep=""),row.names = FALSE)

# Export filtered metabolism predictions
metab_filtered_preds <- filter(metab_all,metab_all$site_name %in% metab_filtered_sites$site_name)
write.csv(metab_filtered_preds,paste(out_save_dir,"/daily_predictions_filtered.csv",sep=""),row.names = FALSE)

##===================================================================##
##              Export model fits for filtered sites                 ##
##===================================================================##
    
# Read in Powell Center data (model fits):
fit_ids <- get_modfit_ids()
fit_ids_sites_filtered <- fit_ids %>% filter(site_name %in% metab_filtered_sites$site_name)
fits <- lapply(fit_ids_sites_filtered[,"file_id"],download_model_fits)

# Export fits:
names(fits) <- fit_ids_sites_filtered$site_name
saveRDS(fits,"./data/out/model_fits_filtered.rds")
    