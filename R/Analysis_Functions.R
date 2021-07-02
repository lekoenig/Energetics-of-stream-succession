
## Project: Semi-mechanistic modeling of river metabolism recovery following storm disturbances
## Script: Analysis Functions
## LE Koenig
## last updated June 2021


##===================================================================##
##             READ IN FILTERED POWELL CENTER DATA                   ##
##===================================================================##

load_filtered_PC_data <- function(){
  
  # Read in list of filtered sites:
  filtered.sites <- read.csv("./data/out/site_data_filtered.csv",header=TRUE,stringsAsFactors = FALSE)
  
  # Read in Powell Center dataset from Maite (contains all sites):
  dailytot <- readRDS("./data/in/Appling_data_manual/dailytot.rds") %>% 
    rename("DATE"="date") %>%
    mutate(date = as.POSIXct(as.character(DATE),format="%Y-%m-%d"),
           doy = as.numeric(strftime(date, format = "%j")),
           year = year(date)) %>%
    select(-DATE) %>%
    # filter Powell Center data set based on list of filtered sites above:
    filter(.,site %in% unique(filtered.sites$site)) %>%
    # create a new column for storm number:
    mutate(storm_id = NA)
  
  # Split filtered dataset into a list based on site name:
  metab.dat <- split(dailytot, dailytot$site)
}


##===================================================================##
##                  DOWNLOAD POWELL CENTER DATA                      ##
##===================================================================##

## 1) Function to download metabolism site data:
download_site_data <- function(save_dir){
  
  # ScienceBase item id for Appling river metabolism data repo (https://www.sciencebase.gov/catalog/item/59bff507e4b091459a5e0982)
  # 1. Site data:
  sb_item_id <- "59bff64be4b091459a5e098b"
  
  dl_dir <- tempdir()
  
  # Download zip folders containing site data:
  site_data <- item_file_download(sb_item_id, dest_dir=dl_dir)
  
  # Read file (metadata: https://www.sciencebase.gov/catalog/file/get/59eb9c0ae4b0026a55ffe389?f=__disk__d6%2F07%2Fbb%2Fd607bb041a2005311ff1318d695ed79907f439cf&transform=1&allowOpen=true)
  site_data_df <- read.csv(paste(dl_dir,"/site_data.tsv",sep=""), sep = "\t",header=TRUE)
  
  # Clean up files and save site data in save_dir:
  files <- list.files(dl_dir, full.names = TRUE) 
  invisible(file.remove(files))
  
  write.csv(site_data,paste(save_dir,"/site_data.csv",sep=""))
  
  return(site_data_df)
}


## 2A) Function to get model output item id's:
get_modfit_ids <- function(site){
  
  # ScienceBase item id for Appling river metabolism data repo (https://www.sciencebase.gov/catalog/item/59bff507e4b091459a5e0982):
  # 6. Model outputs:
  sb_item_id <- "59eb9ba6e4b0026a55ffe37f"
  
  # Check whether item exists on ScienceBase:
  #identifier_exists(sb_item_id)
  sb_id_children <- item_list_children(sb_item_id,limit=356)
  if(length(sb_id_children)<356){print("Warning: Fewer than expected time series files found on ScienceBase. Check connection and/or try re-running download function")}
  
  # Return model fits (names, urls, and ScienceBase item id's):
  file_names <- sapply(sb_id_children, function(item) item$title)
  file_urls <- sapply(sb_id_children,function(item) item$link$url)
  file_ids <- sapply(sb_id_children,function(item) item$id)

  # return ids:
  out <- data.frame(site_name = gsub("_fits","",x = file_names), file_id = file_ids)
  return(out)
}

## 2B) Function to download model outputs from ScienceBase:
download_model_fits <- function(file_id){
  
  dl_dir <- tempdir()
  
  # Download zip folders containing model fits:
  model_zip <- item_file_download(file_id, dest_dir=dl_dir)
  site_name <- stringr::str_sub(stringr::str_extract(model_zip[1], "nwis_\\s*(.*?)\\_"),start=1,end = -2)
  mod_name <- substr(model_zip[1],
                     start = stringr::str_locate(model_zip, "nwis")[1],
                     stop=stringr::str_locate(model_zip,".zip")[1]-1)
  
  # Unpack zips:
  unzip(zipfile = model_zip[1], exdir = dl_dir)
  
  # Read file
  mod_ts <- read.csv(paste(dl_dir,"/daily.tsv",sep=""), sep = "\t",header=TRUE) %>%
            rename("DATE"="date") %>% 
            mutate(site = site_name,
                   fit = mod_name,
                   date = as.POSIXct(as.character(DATE),format="%Y-%m-%d"),
                   SiteYr = year(date))
  
  # Clean up files and save model fits in save_dir:
  files <- list.files(dl_dir, full.names = TRUE) 
  invisible(file.remove(files))
  #read.csv(mod_ts,paste(save_dir,site_name,"_fit.csv",sep=""))
  return(mod_ts)
}

## 3) Function to download metabolism estimates and predictors:
download_model_diagnostics <- function(save_dir){
  
  # ScienceBase item id for Appling river metabolism data repo (https://www.sciencebase.gov/catalog/item/59bff507e4b091459a5e0982)
  # 7. Model diagnostics:
  sb_item_id <- "59eb9bafe4b0026a55ffe382"
  
  dl_dir <- tempdir()
  
  # Download zip folders containing model fits:
  model_diagn_zip <- item_file_download(sb_item_id, dest_dir=dl_dir)
  
  # Unpack zips:
  unzip(zipfile = model_diagn_zip[grep(".zip",model_diagn_zip)], exdir = dl_dir)
  
  # Read file (metadata: https://www.sciencebase.gov/catalog/file/get/59eb9c0ae4b0026a55ffe389?f=__disk__d6%2F07%2Fbb%2Fd607bb041a2005311ff1318d695ed79907f439cf&transform=1&allowOpen=true)
  model_diagnostics <- read.csv(paste(dl_dir,"/diagnostics.tsv",sep=""), sep = "\t",header=TRUE) 

  # Clean up files and save model fits in save_dir:
  files <- list.files(dl_dir, full.names = TRUE) 
  invisible(file.remove(files))
  
  write.csv(model_diagnostics,paste(save_dir,"/diagnostics.csv",sep=""))
  
  return(model_diagnostics)
}


## 4) Function to download metabolism estimates and predictors:
download_metab_est <- function(save_dir){
  
  # ScienceBase item id for Appling river metabolism data repo (https://www.sciencebase.gov/catalog/item/59bff507e4b091459a5e0982)
  # 8. Metabolism estimates and predictors:
  sb_item_id <- "59eb9c0ae4b0026a55ffe389"
  
  dl_dir <- tempdir()
  
  # Download zip folders containing model fits:
  model_zip <- item_file_download(sb_item_id, dest_dir=dl_dir)

  # Unpack zips:
  unzip(zipfile = model_zip[grep(".zip",model_zip)], exdir = dl_dir)
  
  # Read file (metadata: https://www.sciencebase.gov/catalog/file/get/59eb9c0ae4b0026a55ffe389?f=__disk__d6%2F07%2Fbb%2Fd607bb041a2005311ff1318d695ed79907f439cf&transform=1&allowOpen=true)
  metab_est <- read.csv(paste(dl_dir,"/daily_predictions.tsv",sep=""), sep = "\t",header=TRUE) %>%
               rename("DATE"="date") %>%
               mutate(date = as.POSIXct(as.character(DATE), format="%Y-%m-%d"),
                      doy = as.numeric(strftime(date, format = "%j")),
                      year = lubridate::year(date)) %>%
               select(-"DATE")
          
  # Clean up files and save model fits in save_dir:
  files <- list.files(dl_dir, full.names = TRUE) 
  invisible(file.remove(files))
  
  write.csv(metab_est,paste(save_dir,"/daily_predictions.csv",sep=""))
  
  return(metab_est)
}


##===================================================================##
##                ESTIMATE MEDIAN (SEASONAL) GPP                     ##
##===================================================================##

  grow.season.GPP <- function(data){
    
    # Subset years with sufficient observations (only calculating median GPP for years with 120 obervations or more):
    years <- data %>% group_by(year) %>% tally(!is.na(GPP))
    datasub <- data[which(data$year %in% years$year[years$n>90]),]
    
    if(length(datasub$site)>0){
      datasub.ls <- split(datasub,datasub$year)
      MedianGPP_yr <- list()
      
      for(i in 1:length(unique(datasub$year))){
        
        # Identify 120-day window of maximum productivity
        datasub_yr <- datasub.ls[[i]] %>% mutate(sumGPP120 = zoo::rollapply(data=GPP,width=91,FUN=sum,align="center",fill=NA,na.rm=T))
        window <- datasub_yr[c((which.max(datasub_yr$sumGPP120)-45):(which.max(datasub_yr$sumGPP120)+45)),]
        
        # Calculate median productivity during 120-day windows identified above:
        MedianGPP_window <- median(window$GPP,na.rm=T)
        
        MedianGPP_yr[[i]] <- MedianGPP_window
        
      }
      
      out <- median(do.call("rbind",MedianGPP_yr))} else {
        out <- NA}
    
    return(out)
  }

    
##===================================================================##
##               FIND METABOLISM START AND END DATES                 ##
##===================================================================##
  
find_metab_dates <- function(site){
  df <- metab_filter %>% filter(site_name == site)
  out <- data.frame(
    site_name = site,
    metab_start_date = as.character(min(df$date)),
    metab_end_date = as.character(max(df$date)),
    metab_days = length(unique(df$date))
  )
}  
  

##===================================================================##
##       CHECK WHETHER NWIS SITE HAS IN SITU TURBIDITY DATA          ##
##===================================================================##
  
  check_turbidity_data <- function(nwis_site){
    
    # usgs parameterCd <- "63680"  # In situ turbidity
    # see also: https://or.water.usgs.gov/grapher/fnu.html
    # using turbidity methods EPA 180.1 or ISO7027 NTU or FNU

    turb_pcodes <- c("63675","63676","63677","63678","63679","63680","63681","63682","72213","63683","63684")
    
    site_no <- substring(nwis_site,6,100)
    site_data <- whatNWISdata(siteNumber=site_no,service="uv")
    turbidity_data <- turb_pcodes %in% site_data$parm_cd
    turbidity_flag <- ifelse("TRUE" %in% turbidity_data,"YES","NO.DATA")
    
    if(turbidity_flag == "YES"){
    turbidity_data <- site_data[which(site_data$parm_cd %in% turb_pcodes),] %>%
                      select(site_no,parm_cd,begin_date,end_date,count_nu) %>%
                      summarize(turbidity_start = as.character(min(begin_date)),
                                turbidity_end = as.character(max(end_date)),
                                turbidity_days = sum(count_nu))
    }
    site_info <- readNWISsite(site_no)
    
    out <- data.frame(
      site_name = nwis_site,
      drain_area_km2 = site_info$drain_area_va*2.58999, # convert mi^2 to km^2
      turbidity_data = turbidity_flag,
      turbidity_start = ifelse(turbidity_flag == "YES",turbidity_data$turbidity_start,NA),
      turbidity_end = ifelse(turbidity_flag == "YES",turbidity_data$turbidity_end,NA),
      turbidity_total_days = ifelse(turbidity_flag == "YES",turbidity_data$turbidity_days,NA)
    )
    
    return(out)
  }
  

##===================================================================##
##          STANDARDIZE TIME SERIES INTO RELATIVE VALUES             ##
##===================================================================##
  
min_max_norm <- function(x) {
    (x - min(x,na.rm=T)) / (max(x,na.rm=T) - min(x,na.rm=T))
}


##===================================================================##
##             PLOT PARAMETER POSTERIOR DISTRIBUTIONS                ##
##===================================================================##

post.plot <- function(param_post_df){
    post_plot <- as.data.frame(param_post_df) %>%
    ggplot() + geom_density(aes(x=param_post_df),fill="blue",alpha=.3) +
    theme_cowplot() + theme(axis.title=element_text(size=11),axis.text=element_text(size=10),
    plot.title = element_text(size=12))

return(post_plot)
    
}

