## Project: Semi-mechanistic modeling of river metabolism recovery following storm disturbances
## Script: Analysis Functions
## LE Koenig
## last updated August 2021


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



