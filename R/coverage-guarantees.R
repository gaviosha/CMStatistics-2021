
#========================================
#
# Coverage Guarantees Simulation
#
#========================================

library(magrittr)

#---------------------------------
# 
# Models from RNSP paper
#
#---------------------------------


model.list <- list(
  plain.gauss = function() rnorm(100), 
  plain.poisson = function() as.numeric(rpois(200, 1)),
  heterogeneous.gauss = function() c(rep(1, 100), rep(8, 50), rep(1, 100)) * rnorm(250),
  symmetric.bernoulli = function() as.numeric(rbinom(200, 1, .5)),
  plain.cauchy = function() rcauchy(100),
  mix.2 = function() rpois(200, 5) + rnorm(200)/30
)

## Illustrative plots
##

set.seed(42)

model.list$plain.gauss() %>%
  plot(
    type = "b",
    pch = 20,
    xlab = "", 
    ylab = "",
    col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5)
  )

model.list$plain.poisson() %>% 
  plot(
    type = "b",
    pch = 20,
    xlab = "", 
    ylab = "",
    col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5)
  )

model.list$heterogeneous.gauss() %>% 
  plot(
    type = "b",
    pch = 20,
    xlab = "", 
    ylab = "",
    col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5)
  )

model.list$symmetric.bernoulli() %>% 
  plot(
    type = "b",
    pch = 20,
    xlab = "", 
    ylab = "",
    col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5)
  )

model.list$plain.cauchy() %>% 
  plot(
    type = "b",
    pch = 20,
    xlab = "", 
    ylab = "",
    col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5)
  )

model.list$mix.2() %>% 
  plot(
    type = "b",
    pch = 20,
    xlab = "", 
    ylab = "",
    col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5)
  )

## Load methods
## 

library(mqs)


## Numerical study
##



n.reps <- 100

pb <- txtProgressBar(min = 1, max = n.reps, style = 3)

coverage.probabilities <- matrix(0, ncol = 3, nrow = length(model.list))
rownames(coverage.probabilities) <- names(model.list)
colnames(coverage.probabilities) <- c("multiresolution", "runs", "mqs")


for (ii in 1:length(model.list))
{
  
  multiresolution.detections <- numeric(n.reps)
  runs.detections.old <- numeric(n.reps)
  mqs.detections <- numeric(n.reps)
  
  
  cat("\r", names(model.list)[ii])
  
  for (jj in 1:n.reps)
  {
    ee <- model.list[[ii]]()
    multiresolution.detections[jj] <- nrow(nsp_robust(ee, method = "multiresolution")$intervals) == 0
    runs.detections[jj] <- nrow(nsp_robust(ee, method = "runs")$intervals) == 0
    mqs.detections[jj] <- nrow(mqse(ee, alpha = .1)$confInt) == 0
    
    setTxtProgressBar(pb, jj)
  }
  
  coverage.probabilities[ii,1] <- 100 * mean(multiresolution.detections)
  coverage.probabilities[ii,2] <- 100 * mean(runs.detections)
  coverage.probabilities[ii,3] <- 100 * mean(mqs.detections)
  
}