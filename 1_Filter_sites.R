
## Project: Semi-mechanistic modeling of river metabolism recovery following storm disturbances
## Script: Filter Powell Center sites
## LE Koenig
## last updated August 2021

## The objective of this script is to filter sites from the Powell Center database of stream metabolism time series to use in the stream successional energetics project.

## Load packages: 
library(Appling356MetabolicRegimes) # interface with Appling data set on ScienceBase (https://github.com/lekoenig/Appling356MetabolicRegimes)
library(dplyr)         # general data cleaning and manipulation
library(tidyr)         # general data cleaning and manipulation
library(purrr)         # general data cleaning and manipulation
library(ggplot2)       # make plots
library(lubridate)     # format timestamps
library(dataRetrieval) # interface with NWIS
library(sbtools)       # interface with Science Base
library(maps)          # access some basemap data
library(sf)            # geospatial data manipulation 

source("./R/Analysis_Functions.R")

# Path to save filtered Powell Center data:
in_save_dir <- "./data/Appling_PC_data/in/"
out_save_dir <- "./data/Appling_PC_data/out/"


##===================================================================##
##                  Read in Powell Center data                       ##
##===================================================================##

# Read in site information from Appling et al. 2018 data set:
pc_sites <- download_site_data(in_save_dir)

# Read in daily metabolism estimates and predictor variables:
metab_all <- download_metabolism_estimates(in_save_dir) 

# Read model diagnostics:
diagnostics <- download_model_diagnostics(in_save_dir)


##===================================================================##
##                    Map Powell Center sites                        ##
##===================================================================##

site_info <- metab_all %>% group_by(site_name) %>% filter(!is.na(GPP)) %>% summarize(n = n())

state <- map_data("state") 
usa <- map_data("usa") %>% st_as_sf(.,coords=c("long","lat"),crs=4326)
pc_sites_sp <- left_join(site_info,pc_sites,by="site_name") %>%
               st_as_sf(.,coords=c("lon","lat"),crs=4269) %>% st_transform(.,4326) %>%
               st_filter(st_as_sfc(st_bbox(usa)),
                         .predicate=st_within) %>%
               mutate(lon = st_coordinates(.)[,1],
                      lat = st_coordinates(.)[,2])

ggplot() + geom_polygon(data=state,aes(x=long,y=lat,group=group),fill=NA,color="gray50",size=.25)+
           geom_point(data=pc_sites_sp,aes(x=lon,y=lat,color=n))+
           viridis::scale_color_viridis(name="Daily metabolism \nestimates",option = "D")+
           coord_map("albers",lat0=30,lat1=40) +
           theme_bw() + 
           theme(legend.position = c(0.9, 0.2),legend.text = element_text(size=10),legend.title=element_text(size=11),
                 axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
                 axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
                 panel.border = element_blank(),panel.grid=element_blank()) 
  

##===================================================================##
##                  Filter Powell Center sites                       ##
##===================================================================##

# Function to calculate pearson correlation coefficient between ER and K600 for each df within the list:
ER_K_corr <- function(df){
  out <- split(x=df,f = df$site_name) %>%
    purrr::map(~ cor(x=.$ER,y=.$K600,method="pearson")) %>%
    data.frame %>% gather(.,site_name,cor.coef)
}

# Function to calculate median growing season GPP for each df within the list:
est_med_seasonal_GPP <- function(df){
  out <- split(x=df,f = df$site_name) %>%
    purrr::map(~ grow.season.GPP(.)) %>%
    data.frame %>% gather(.,site_name,GPP_grow_median)
}

# Filter data based on defined criteria:
metab_filter <- metab_all %>%
                  mutate(year = lubridate::year(date)) %>%
                # 1. Select only those sites where we have high confidence in the model:
                  left_join(.,diagnostics[,c("site","model_confidence")],by=c("site_name"="site")) %>%
                  filter(.,model_confidence=="H") %>%
                # 2. Select sites where the correlation between ER and K600 < 0.4:
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
                        # Create new columns that indicate metab dates:
                        left_join(.,bind_rows(lapply(.$site_name,find_metab_dates)),by="site_name") %>%
                        mutate(metab_interval = interval(start=metab_start_date,end=metab_end_date)) %>%
                        # Check whether the metabolism and turbidity date ranges overlap:
                        mutate(check_turbidity_range = int_overlaps(turb_interval,metab_interval)) %>%
                        # Filter sites that overlap with turbidity data (at least 180 days):
                        filter(.,turbidity_total_days > 180 & check_turbidity_range == "TRUE") %>% 
                        # Bind with site info from ScienceBase and select columns:
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
fits <- list()
for(i in seq_along(metab_filtered_sites$site_name)){
  fit_list <- download_model_outputs(sitename = metab_filtered_sites$site_name[i],overwrite_file = TRUE)
  #fit_list <- list(fit)
  names(fit_list) <- metab_filtered_sites$site_name[i]
  fits <- c(fits,fit_list)
  print(i)
}

# Create a list containing one data frame for each site within metab_filtered_sites:
cols_select <- c("date","GPP_daily_sd","ER_daily_sd","K600_daily_sd","valid_day","warnings","errors")

obs_err_ls <- lapply(fits,function(x) {
  x %>% lapply(.,"[[","daily") %>% lapply(.,"[",cols_select) %>% bind_rows(.id="res") %>%
    filter(valid_day=="TRUE") %>% 
    mutate(date2 = as.Date(date)) %>%
    arrange(date2) %>%
    select(date2,res,GPP_daily_sd,ER_daily_sd,K600_daily_sd,valid_day,warnings,errors) %>%
    rename(date = date2)
})

# Export fits:
saveRDS(obs_err_ls,"./data/Appling_PC_data/out/model_fits_filtered.rds")
    