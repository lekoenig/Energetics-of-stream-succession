
## Project: Semi-mechanistic modeling of river metabolism recovery following storm disturbances
## Script: Analysis Functions
## LE Koenig
## last updated June 2021


##===================================================================##
##             READ IN FILTERED POWELL CENTER DATA                   ##
##===================================================================##

load_filtered_PC_data <- function(){
  
  # Read in list of filtered sites:
  filtered.sites <- read.csv("./data/out/PC_site_data_filtered.csv",header=TRUE,stringsAsFactors = FALSE)
  
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

## 1) Function to get model output item id's:
get_modfit_ids <- function(){
  
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
}

## 2) Function to download model outputs from ScienceBase:
download_model_fits <- function(file_id){
  
  dl_dir <- tempdir()
  
  # Download zip folders containing model fits:
  model_zip <- item_file_download(file_id, dest_dir=dl_dir)
  site_name <- substr(model_zip[1],
                      start = stringr::str_locate(model_zip, "nwis")[1],
                      stop = stringr::str_locate(model_zip, "min_fit.zip")[1]-4)
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
                   SiteYr = year(date)) %>%
            select(site,fit,date,SiteYr,GPP_daily_mean,GPP_daily_2.5pct,GPP_daily_97.5pct,GPP_daily_Rhat,
                   ER_daily_mean,ER_daily_2.5pct,ER_daily_97.5pct,ER_daily_Rhat,K600_daily_mean,
                   K600_daily_2.5pct,K600_daily_97.5pct,K600_daily_Rhat)
  
  # Clean up files and save model fits in save_dir:
  files <- list.files(dl_dir, full.names = TRUE) 
  invisible(file.remove(files))
  #write.csv(mod_ts,paste(save_dir,site_name,"_fit.csv",sep=""))
  return(mod_ts)
}


## 3) Function to download metabolism estimates and predictors:
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
               rename("light_Wm2" = "shortwave",
                      "discharge_m3s" = "discharge",
                      "velocity_ms" = "velocity")

  # Clean up files and save model fits in save_dir:
  files <- list.files(dl_dir, full.names = TRUE) 
  invisible(file.remove(files))
  
  write.csv(metab_est,paste(save_dir,"daily_predictions.csv",sep=""))
  
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
##       CHECK WHETHER NWIS SITE HAS IN SITU TURBIDITY DATA          ##
##===================================================================##
  
  
  check.turbidity.data <- function(nwis_site){
    
    #turbidity.data <- readNWISuv(site.no, "63675", 
    #                             "2007-01-01", "2018-12-31")
    #turbidity.flag <- ifelse(length(turbidity.data)>0,"YES","NO.DATA")
    
    turb.pcodes <- c("63675","63676","63677","63678","63679","63680","63681","63682","72213","63683","63684")
  
    site.data <- whatNWISdata(siteNumber=site.no,service="uv")
    turbidity.data <- turb.pcodes %in% site.data$parm_cd
    turbidity.flag <- ifelse("TRUE" %in% turbidity.data,"YES","NO.DATA")
    return(turbidity.flag)
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
