# Project: Semi-mechanistic modeling of river metabolism recovery following storm disturbances
# Script: Save storm sub
# LE Koenig
# last updated August 2021

# Load packages:
library(dplyr)         # general data cleaning and manipulation

source("./R/Analysis_Functions.R")

##===================================================================##
##     Focus on one river as a reference for data simulations        ##
##===================================================================##

metab_data <- read.csv("./data/Appling_PC_data/out/daily_predictions_filtered.csv",header=TRUE)

# Select South Branch Potomac, WV as a test case:
pot_data <- metab_data[which(metab_data$site_name=="nwis_01608500"),] %>%
            # convert light units from W m-2 to umol m-2 s-1, and relativize:
            mutate(light_umolm2s = shortwave*0.21739,
            light_rel = min_max_norm(light_umolm2s),
            date = as.Date(as.POSIXct(as.character(date), format="%Y-%m-%d")))

# Pull observation error from model fits: 
obs_err <- readRDS("./data/Appling_PC_data/out/model_fits_filtered.rds") %>% 
           .[[which(names(.)==pot_data$site_name[1])]] %>%
           select(date,GPP_daily_sd,ER_daily_sd)

# Select one storm for South Branch Potomac as an example/point of reference:
pot_data$storm_id <- NA
pot_data$storm_id[c(which(pot_data$date=="2012-01-13"):which(pot_data$date=="2012-04-04"))] <- 1
pot_data$storm_id[c(which(pot_data$date=="2012-03-26"):which(pot_data$date=="2012-05-03"))] <- 2

storm <- filter(pot_data,storm_id==2) %>% 
         left_join(.,obs_err[,c("date","GPP_daily_sd","ER_daily_sd")],
                   by="date")

# For now, filter out days where GPP < 0:
#storm <- storm %>% 
#         filter(storm$GPP>0) %>%
#         mutate(GPP = zoo::na.approx(GPP),
#                GPP_sd = zoo::na.approx(GPP_daily_sd),
#                light_rel = min_max_norm(light_umolm2s),
#                light_rel_smooth = zoo::rollmean(light_rel,k=3,fill=NA,align="left"))
#storm$light_rel_smooth[which(is.na(storm$light_rel_smooth))] <- rnorm(length(which(is.na(storm$light_rel_smooth))),mean(tail(storm$light_rel_smooth,n = 7),na.rm=TRUE),0.02)

# Save subset:
write.csv(storm,"./data/Appling_PC_data/out/SouthBrPotomac_storm_sub.csv",row.names=FALSE)



