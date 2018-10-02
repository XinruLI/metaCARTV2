rm(list = ls())
source("BumpXvalidV2.R")
sysTime <- Sys.time()
print(sysTime)
set.seed(sysTime)
Metrics0 <- Metrics1 <- Metrics2 <- matrix(NA, ncol = 3, nrow = 100)
for (i in 1:1){
  x1 <- sample(0:1, replace = T, 100)
  x2 <- sample(0:1, replace = T, 100)
  x3 <- runif(100)
  x4 <- runif(100)
  #----                        Model C                       ----#
  y.grand <- (x1  == 1)* 0.5
  dat <- data.frame(SimDataForBump(y.grand, n.ave = 80, tau = 0.1),
                    x1, x2, x3, x4)
  res0 <- REmrt(g ~ x1 + x2 + x3 + x4, vi = vi, data = dat,
                maxL  = 8, n.fold = 10, c = 1)
  Metrics0[i, 1] <- length(res0$n)
  Metrics0[i, 2] <- sum(c("x1") %in% res0$mod)
  if (is.null(res0$mod)) {
    Metrics0[i, 3] <- 0
  } else {
    Metrics0[i, 3] <- sum(!res0$mod %in% c("x1"))
  }
  
  res1 <- BumpNestingXvalid(g ~ x1 + x2 + x3 + x4, vi = vi, data = dat,
                            maxL  = 8, n.fold = 10, B = 25, c = 1)
  Metrics1[i, 1] <- length(res1$n)
  Metrics1[i, 2] <- sum(c("x1") %in% res1$mod)
  Metrics1[i, 3] <- sum(!res1$mod %in% c("x1"))
  
  res2 <- XcrossNestingBump(g ~ x1 + x2 + x3 + x4, vi = vi, data = dat,
                            maxL  = 8, n.fold = 10, B = 25, c = 1)
  res2 <- rbind(c(0, res0$cv.res[1, ]), res2)
  maxL <- which(res2[,2]<= (min(res2[,2]) + res2[,3][which.min(res2[,3])]))[1] - 1
  if (maxL == 0) {
    Metrics2[i, 1] <- 1
    Metrics2[i, 2] <- 0
    Metrics2[i, 3] <- 0
  } else {
    mf <- res0$data
    resB2 <- bump(mf, maxL = maxL, minsplit = 6, minbucket = 3, 
                  cp = 0.01, lookahead = F, 25) 
    mods <- unique(resB2$tree$mod[!is.na(resB2$tree$mod)])
    Metrics2[i, 1] <- length(resB2$g)
    Metrics2[i, 2] <- sum(c("x1") %in% mods)
    Metrics2[i, 3] <- sum(! mods %in% c("x1"))
  }
}
#save(Metrics0, Metrics1, Metrics2, file = "BumpXvalidModelBPilot")
print(Sys.time())
Sys.time() - sysTime
