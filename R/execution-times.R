
#========================
#
# Execution times
#
#========================

library(mqs)
library(microbenchmark)

set.seed(42)


### Get execution times
###

nn.seq <- c(100,500,1000)

nsp.mr <- c()
nsp.runs <- c()
mqs <- c()

for (nn in nn.seq)
{
  ee <- rnorm(nn)
  
  nsp.mr <- c(nsp.mr, microbenchmark::microbenchmark(nsp_robust(ee, method = "multiresolution"), times = 1)$time / 10**9)
  
  nsp.runs <- c(nsp.runs, microbenchmark::microbenchmark(nsp_robust(ee, method = "runs"), times = 1)$time / 10**9)
  
  mqs <- c(mqs, microbenchmark::microbenchmark(mqse(ee, alpha = .1), times = 1)$time / 10**9)
  
  print(nn)
}



### Plot execution times
###

barplot(
  data.frame(
    nsp.runs,
    rep(0,3),
    rep(0,3)
  ) %>% 
    as.matrix() %>% 
    t(), 
  beside = TRUE,
  ylab = "execution time (seconds)",
  names.arg = c(100, 500, 1000), 
)

legend(1,1.5,
       c("runs","Fryzlewicz (2021)", "Jula Vanegas et al. (2021)"),
       fill = gray.colors(3),
       cex = 0.8)

barplot(
  data.frame(
    nsp.runs,
    nsp.mr,
    rep(0,3)
  ) %>% 
    as.matrix() %>% 
    t(), 
  beside = TRUE,
  ylab = "execution time (seconds)",
  names.arg = c(100, 500, 1000)
)

legend(1,50,
       c("runs","Fryzlewicz (2021)", "Jula Vanegas et al. (2021)"),
       fill = gray.colors(3),
       cex = 0.8)

barplot(
  data.frame(
    nsp.runs,
    nsp.mr,
    mqs
  ) %>% 
    as.matrix() %>% 
    t(), 
  beside = TRUE,
  ylab = "execution time (seconds)",
  names.arg = c(100, 500, 1000)
  )

legend(1,400,
       c("runs","Fryzlewicz (2021)", "Jula Vanegas et al. (2021)"),
       fill = gray.colors(3),
       cex = 0.8)

barplot(
  data.frame(
    nsp.runs,
    nsp.mr,
    mqs
  ) %>% 
    as.matrix() %>% 
    t() %>% 
    log(), 
  beside = TRUE,
  ylab = "execution time (log-seconds)",
  names.arg = c(100, 500, 1000)
)

legend(1,6,
       c("runs","Fryzlewicz (2021)", "Jula Vanegas et al. (2021)"),
       fill = gray.colors(3),
       cex = 0.8)
