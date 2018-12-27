rm(list = ls()) # clear global environment
graphics.off() # close all graphics
library(pacman) # needs to be installed first
p_load(R.matlab, plotly, extrafont, grDevices, gridExtra,
       dplyr, stringr, tidyverse, utils, reshape2,
       anomalize, MVN, forecast, fractal,
       ecp, dfphase1, 
       MALDIquant, TSclust,
       knitr, kableExtra, dplyr) 

load(file = "./Data/RGenerated/NormMedianFilteredData.RData")
changepoints_mph1_mf <- list()

for (i in 1:15) {
  df <- get(paste0("subject",i,"_medianf_norm"))
  window.size <- 2
  rows.to.keep <- nrow(df[,2:4]) %/% window.size
  remainder <- nrow(df[,2:4]) %% window.size
  df.MultivariateChangepoint<- df[1:(nrow(df)-remainder),
                                  2:4]
  dataholder.array <- t(df.MultivariateChangepoint) %>% 
    array(c(3, window.size, rows.to.keep))
  y <- mphase1(dataholder.array,
               isolated=FALSE, plot = FALSE, alpha=0.0001)
  mcp <- (df[unlist(y$alasso[2]),1]*window.size) %>% round(digits = 0)
  mcp <- mcp[order(mcp)]
  mcp <- (mcp - window.size) + 1 # to make the change point at the begining of the period and not end (to be similar to the ecp method)
  changepoints_mph1_mf[[i]] <- y
  
  cat('####',paste0("P", i), "{-}",'\n')
  df_transformed <- melt(get(paste0("subject",i,"_medianf_norm")),
                         id.vars = "percent.time.from.start",
                         measure.vars=c("scaled.stride.len",
                                        "scaled.stride.height",
                                        "stride.duration"
                         )
  ) # ggplot data needs to be tall
  
  assign(paste0("g",i),
         ggplot(data = df_transformed,
                aes(x=percent.time.from.start, y=value, group=variable,
                    color=variable,shape=variable)) + 
           geom_line() + theme_bw() + 
           ggtitle(paste("dfphase1 package Changepoints for Participant",i)) + 
           theme(legend.position="none",
                 #axis.text.x=element_text(angle=90,hjust=1),
                 plot.title = element_text(hjust = 0.5)) +
           facet_wrap(~variable,nrow=3, scales = "free_y") +
           geom_vline(xintercept= mcp)
  )
  
  print(get(paste0("g",i)))
  
  cat('\n')
  
  cat(paste0('<source> <p> Based on the <b>dfphase1 package</b>, the number of changepoints for the cusums of the median filtered data for participant ', i,
             ' is equal to: ',
             length(mcp),
             '. These changepoints are located at: ',
             paste(mcp, collapse = ", "),
             '. The results from the analysis above can be found at: <a href="https://github.com/fmegahed/fatigue-changepoint/blob/master/Data/RGenerated/ChangePointsMedianFilterMphase1.RData">ChangePointsMedianFilterMphase1.RData</a>. </p> </source>')
  )
  cat('\n \n')
}
save(changepoints_mph1_mf,
     file = "./Data/RGenerated/ChangePointsMedianFilterMphase1.RData")