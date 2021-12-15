
#====================================
#
# Examples of method's performance
#
#====================================

blocks <- list(
  name = "blocks",
  signal = c(rep(0,205),rep(14.64,62),rep(-3.66,41),rep(7.32,164),rep(-7.32,40)), 
  simulate = function() c(rep(0,205),rep(14.64,62),rep(-3.66,41),rep(7.32,164),rep(-7.32,40)) + rnorm(512, sd = 7)
)

waves <- list(
  name = "waves",
  signal = c((1:150) * (2**-3), (150:1) * (2**-3), (1:150) * (2**-3), (150:1) * (2**-3)), 
  simulate = function() c((1:150) * (2**-3), (150:1) * (2**-3), (1:150) * (2**-3), (150:1) * (2**-3)) + rnorm(600, sd = 3)
)


#-------------------------------
#
# Piecewise-constant plots
#
#-------------------------------

set.seed(42)

blocks$simulate() -> ee 


## plot signal only


blocks$signal %>% plot(
  type = "l",
  xlab = "time",
  ylab = " ",
  lwd = 3,
  ylim = c(min(ee)-1,max(ee+1))
)


## plot signal and noise

plot(ee ~ c(1:length(ee)), 
     pch = 20,
     xlab = "time", 
     ylab = " ",
     col = rgb(red = 0, green = 0, blue = 0, alpha = 0.25)
     )

lines(
  blocks$signal,
  lwd = 3,
  )


## add detection intervals to plot

rnsp.obj <- nsp_robust(ee, method = "runs")

draw_rects(rnsp.obj, yrange = c(-100,100))


#--------------------------------
#
# Piecewise-linear plots
#
#--------------------------------

waves$simulate() -> ee 


## plot signal only


waves$signal %>% plot(
  type = "l",
  xlab = "time",
  ylab = " ",
  lwd = 3,
  ylim = c(min(ee)-1,max(ee+1))
)


## plot signal and noise

plot(ee ~ c(1:length(ee)), 
     pch = 20,
     xlab = "time", 
     ylab = " ",
     col = rgb(red = 0, green = 0, blue = 0, alpha = 0.25)
)

lines(
  waves$signal,
  lwd = 3,
)


## add detection intervals to plot

rnsp.obj <- nsp_robust(ee, method = "runs", deg = 1)

draw_rects(rnsp.obj, yrange = c(-100,100))

