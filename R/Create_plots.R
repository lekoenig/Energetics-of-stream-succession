## Project: Semi-mechanistic modeling of river metabolism recovery following storm disturbances
## Script: Create plots
## LE Koenig
## last updated August 2021

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


##===================================================================##
##                       PLOT SIMULATED DATA                         ##
##===================================================================##

sim.plot <- function(df){
  
  plot <- df %>% 
          ggplot() + 
          geom_line(aes(x=time,y=GPPsim),color="darkgreen",alpha=.5) + 
          geom_point(aes(x=time,y=GPPsim),color="darkgreen") + 
          labs(x=expression(Days~since~storm),y=expression(Simulated~GPP~(g~O[2]~m^-2~d^-1))) + 
          theme_bw() + theme(panel.grid = element_blank())
  return(plot)
  
}





