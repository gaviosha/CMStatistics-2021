
#==========================================
#
# Analysis of NO2 concentraition in Spain
#
#==========================================


## Helpers
##


library(magrittr)

ts.plot.nice <- function(xx, main, ylab)
{
  #' Plot polutant time series nicely
  #'
  #'@param xx data frame, from \data
  
  plot(xx$Concentration, type = 'l', 
       xaxt = 'n', 
       xlab = 'Date', 
       ylab = ylab, 
       main = main)
  
  axis(1, at = seq(1, length(xx$Concentration), length.out = 20), 
       label = format.Date(xx$Date, "%Y-%m")[seq(1, length(xx$Date), length.out = 20)]
  )
}

## Madrid
##

Madrid <- read.csv("data/NO2_ES_Madrid_2020.csv") %>% na.exclude()
rnsp.Madrid <- nsp_robust(Madrid$Concentration, method = "runs", deg = 1)
ts.plot.nice(Madrid, "Madrid", "NO2")
draw_rects(rnsp.Madrid, yrange = c(-200,200))


## Barcelona
##

Barcelona <- read.csv("data/NO2_ES_Barcelona_2020.csv") %>% na.exclude()
rnsp.Barcelona <- nsp_robust(Barcelona$Concentration, method = "runs", deg = 1)
ts.plot.nice(Barcelona, "Barcelona", "NO2")
draw_rects(rnsp.Barcelona, yrange = c(-200,200))


## Valencia
##

Valencia <- read.csv("data/NO2_ES_Valencia_2020.csv") %>% na.exclude()
rnsp.Valencia <- nsp_robust(Valencia$Concentration, method = "runs", deg = 1)
ts.plot.nice(Valencia, "Valencia", "NO2")
draw_rects(rnsp.Valencia, yrange = c(-200,200))


## Sevilla
##

Sevilla <- read.csv("data/NO2_ES_Sevilla_2020.csv") %>% na.exclude()
rnsp.Sevilla <- nsp_robust(Sevilla$Concentration, method = "runs", deg = 1)
ts.plot.nice(Sevilla, "Sevilla", "NO2")
draw_rects(rnsp.Sevilla, yrange = c(-200,200))



