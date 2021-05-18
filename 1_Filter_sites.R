
## Project: Semi-mechanistic modeling of river metabolism recovery following storm disturbances
## Script: Filter sites
## LE Koenig
## last updated 30 January 2020

## The objective of this script is to filter sites from the Powell Center database of stream metabolism time series to use in the stream successional energetics project.


## Load packages: 
library(dplyr)       # general data cleaning and manipulation
library(ggplot2)     # create plots
library(lubridate)   # format timestamps
library(sf)          # working with spatial data
library(mapview)     # plot spatial data
library(lutz)        # look up time zones

## require mda.streams package (see Appling et al. 2018: https://www.nature.com/articles/sdata2018292)
#install.packages("mda.streams", dependencies = TRUE, 
#                 repos = c("https://owi.usgs.gov/R","https://cran.rstudio.com"))
library(mda.streams)

source("./R/Analysis_Functions.R")

##===================================================================##
##                  Read in Powell Center data                       ##
##===================================================================##

# Read in NWIS sites from Powell Center data set (Appling et al. 2018; data downloaded 30 January 2020):
  #pc.dat <- read.table("./data/PC_data/daily_predictions.tsv",header=T,stringsAsFactors = F)
  #pc.dat$date <- as.POSIXct(as.character(pc.dat$date), format="%Y-%m-%d")
  pc.sites <- read.table("./data/PC_data/site_data.tsv",header=T,stringsAsFactors = F)

# Subset Powell Center data:
  #dat.sub <- pc.dat[,c("site_name","date","GPP","GPP.lower","GPP.upper", "GPP.Rhat",
  #                    "ER","ER.lower","ER.upper","K600","K600.lower","K600.upper",
  #                    "temp.water","discharge","shortwave","velocity")]

# Read in Powell Center dataset from Maite (contains all sites):
  dailytot <- readRDS("./data/PC_data/dailytot.rds")
  dailytot$date <- as.POSIXct(as.character(dailytot$date), format="%Y-%m-%d")
  dailytot$doy <- as.numeric(strftime(dailytot$date, format = "%j"))
  dailytot$year <- lubridate::year(dailytot$date)

##===================================================================##
##                  Filter Powell Center sites                       ##
##===================================================================##

# 1. Select only those sites where we have high confidence in the model:
  dailytot_high <- dailytot[which(dailytot$confidence=="H"),]  # results in 254 sites

# 2. Select sites where the correlation between ER and K600 < 0.4:
  
  # Calculate ER-K600 pearson correlation coefficient for each site within the high confidence subset:
  sitelist <- split(dailytot_high, dailytot_high$site)
  
  sitelist_ERKcorr <- data.frame(site = names(sitelist),
                                 corr_ERK = NA)
  sitelist_ERKcorr$corr_ERK <- as.numeric(do.call("rbind",lapply(sitelist,function(x){
                                cor(na.omit(subset(x,select=c(ER,K600))),method="pearson")[2,1]
  })))
  
  # Plot distribution of ER-K600 correlation coefficients:
  hist(sitelist_ERKcorr$corr_ERK)
  
  # Choose sites where ER-K600 correlation < 0.4 (results in 179 sites):
  dailytot_ERKcorr_good <- dailytot_high[which(dailytot_high$site %in% sitelist_ERKcorr$site[which(abs(sitelist_ERKcorr$corr_ERK) < 0.4)]),]

# 3. Filter sites by high productivity
  
  sitelist2 <- split(dailytot_ERKcorr_good, dailytot_ERKcorr_good$site)
  
  # Approach A: Median GPP across the whole time series:
    sitelist_medianGPP <- data.frame(site = names(sitelist2),
                                     median_GPP = NA)
    sitelist_medianGPP$median_GPP <- as.numeric(do.call("rbind",lapply(sitelist2,function(x){
                                     median(x$GPP,na.rm=T)
    })))
    
    # Choose sites where median GPP greater than or equal to 3 g/m2/d (results in 46 sites):
    dailytot_GPP <- dailytot_ERKcorr_good[which(dailytot_ERKcorr_good$site %in% sitelist_medianGPP$site[which(sitelist_medianGPP$median_GPP >= 3)]),]
  
  
  # Approach B: If some sites are productive during the growing season, but "dead" during the winter, median GPP would be biased.
  # Median GPP during growing season (120-d continuous period of maximum productivity)

    sitelist_seasonal_medianGPP <- data.frame(site = names(sitelist2),
                                   median_GPP = NA)
    sitelist_seasonal_medianGPP$median_GPP <- as.numeric(do.call("rbind",lapply(sitelist2,try(Grow.season.GPP))))
    
    # Choose sites where median seasonal GPP greater than or equal to 3 g/m2/d (results in 86 sites):
    dailytot_GPP_seasonal <- dailytot_ERKcorr_good[which(dailytot_ERKcorr_good$site %in% sitelist_seasonal_medianGPP$site[which(sitelist_seasonal_medianGPP$median_GPP >= 3)]),]
    
    
# 4. Find sites that are co-located with turbidity sensors?
   
    library(dataRetrieval)
    
    parameterCd <- "63680"  # In situ turbidity
    # see also: 
    # https://or.water.usgs.gov/grapher/fnu.html
    # using turbidity methods EPA 180.1 or ISO7027 NTU or FNU
    parameterINFO <- readNWISpCode(parameterCd)

    sites <- unique(dailytot$site)

    site_data <- data.frame(sites = rep(NA,length(sites)),
                           name = rep(NA,length(sites)),
                           lat = rep(NA,length(sites)),
                           lon = rep(NA,length(sites)),
                           datum = rep(NA,length(sites)),
                           DrArea_km2 = rep(NA,length(sites)),
                           turbidity_data = rep(NA,length(sites)))
    
    for(i in 1:length(sites)){
      site.no <- substring(sites[i],6)
      site.info <- readNWISsite(site.no)
      
      site_data$sites[i] <- sites[i]
      site_data$name[i] <- site.info$station_nm
      site_data$lat[i] <- site.info$dec_lat_va
      site_data$lon[i] <- site.info$dec_long_va
      site_data$datum[i] <- site.info$coord_datum_cd
      site_data$DrArea_km2[i] <- site.info$drain_area_va*2.58999
      site_data$turbidity_data[i] <- check.turbidity.data(nwis_site = site.no)
      
      print(i)
    } 
    
    #write.csv(site_data,"./output/PC_site_data.csv",row.names = FALSE) # 209/356 sites have in situ turbidity data (not taking into account data overlap/length of record)
    
   # How many of the filtered sites (from above) have turbidity data?
    # Approach A above, 16 sites
    Filtered_sites_medianGPP <- data.frame(site = unique(dailytot_GPP$site))
    Filtered_sites_medianGPP <- left_join(Filtered_sites_medianGPP,site_data,by=c("site"="sites"))
    Filtered_sites_medianGPP <- Filtered_sites_medianGPP[which(Filtered_sites_medianGPP$turbidity_data=="YES"),]
    
    # Approach B above, X sites, 38 sites
    Filtered_sites_seasonalGPP <- data.frame(site = unique(dailytot_GPP_seasonal$site))
    Filtered_sites_seasonalGPP <- left_join(Filtered_sites_seasonalGPP,site_data,by=c("site"="sites"))
    Filtered_sites_seasonalGPP <- Filtered_sites_seasonalGPP[which(Filtered_sites_seasonalGPP$turbidity_data=="YES"),]
    
  # How many turbidity days are represented?
    for(i in 1:length(Filtered_sites_seasonalGPP$site)){
    
      metab.sub <- dailytot_GPP_seasonal[which(dailytot_GPP_seasonal$site==Filtered_sites_seasonalGPP$site[i]),]
      Filtered_sites_seasonalGPP$metab_start[i] <- as.character(min(metab.sub$date))
      Filtered_sites_seasonalGPP$metab_end[i] <- as.character(max(metab.sub$date))
      Filtered_sites_seasonalGPP$metab_days[i] <- length(unique(metab.sub$date))
      
      turb.pcodes <- c("63675","63676","63677","63678","63679","63680","63681","63682","72213","63683","63684")
      site.no <- substring(Filtered_sites_seasonalGPP$site[i],6)
      site.data <- whatNWISdata(siteNumber=site.no,service="uv")
      turbidity.data <- site.data[which(site.data$parm_cd %in% turb.pcodes),]
      Filtered_sites_seasonalGPP$turbidity_pcode[i] <- turbidity.data$parm_cd[1]
      Filtered_sites_seasonalGPP$turbidity_start[i] <- as.character(min(turbidity.data$begin_date))
      Filtered_sites_seasonalGPP$turbidity_end[i] <- as.character(max(turbidity.data$end_date))
      Filtered_sites_seasonalGPP$turbidity_days[i] <- sum(turbidity.data$count_nu)
      print(i)
    }
    
  # Do the metabolism and turbidity date ranges overlap?
    Filtered_sites_seasonalGPP$check_range <- (Filtered_sites_seasonalGPP$metab_start <= Filtered_sites_seasonalGPP$turbidity_end) & (Filtered_sites_seasonalGPP$metab_end <= Filtered_sites_seasonalGPP$turbidity_start)
        
  # How many of our filtered sites have > 180 days of turbidity data?
    sites <- Filtered_sites_seasonalGPP[which(Filtered_sites_seasonalGPP$check_range=="FALSE" & Filtered_sites_seasonalGPP$turbidity_days>180),]
    dim(sites)

  # Plot geographic location of points:
    sites.sf.NAD83sub <- sites[which(sites$datum=="NAD83"),] %>% 
                      sf::st_as_sf(.,coords=c("lon","lat"),crs=4269) %>%
                      sf::st_transform(.,4326)
    sites.sf.NAD27sub <- sites[which(sites$datum=="NAD27"),] %>% 
                      sf::st_as_sf(.,coords=c("lon","lat"),crs=4267) %>%
                      sf::st_transform(.,4326)
    sites.sf <- rbind(sites.sf.NAD83sub,sites.sf.NAD27sub)
    
  # Add a column that indicates the site's time zone:
    sites.sf$tz <- tz_lookup(sites.sf, crs = 4326, method = "accurate", warn = TRUE)
    
  # Map filtered sites:
    mapview(sites.sf)
    
  # Convert sf object to a data frame and export:
    sites2 <- sites.sf %>%
              mutate(Lat = st_coordinates(sites.sf)[,2],
                     Lon = st_coordinates(sites.sf)[,1],
                     datum = "WGS84") %>%
              st_drop_geometry(.)
      
    write.csv(sites2,"./output/data/PC_site_data_filtered.csv",row.names = FALSE)
    
    
    
    
    
    
    