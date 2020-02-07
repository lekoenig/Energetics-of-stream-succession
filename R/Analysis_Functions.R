
## Analysis Functions used for Energetics of stream succession R project



##===================================================================##
##                Estimate median (seasonal) GPP                     ##
##===================================================================##

  Grow.season.GPP <- function(data){
    
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
##       CHEcK WHETHER NWIS SITE HAS IN SITU TURBIDITY DATA          ##
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
  


  
  
  
  
