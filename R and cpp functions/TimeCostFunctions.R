load("mf2")
library(Rcpp)
source("Xvalid_all_cpp.R")
source("REmrt0.R")
# source("REmrt_GS_update.R")
library("profvis")
res1 <- REmrt.fit0(mf, maxL = 120, minbucket=3, minsplit=4, delQ=0.00001, lookahead=F)
res2 <- REmrt_GS_cpp(mf, maxL = 120, minbucket=3, minsplit=4, delQ=0.00001, lookahead=F)
all(abs(res1$tree$Qb - res2$tree$Qb) < 1e-10)

system.time(REmrt.fit0(mf, maxL = 120, minbucket=3, minsplit=4, delQ=0.00001, lookahead=F))
system.time(REmrt_GS_cpp(mf, maxL = 120, minbucket=3, minsplit=4, delQ=0.00001, lookahead=F))


profvis(expr = {
  REmrt_GS_cpp(mf, maxL = 120, minbucket=3, minsplit=2, delQ=0.00001, lookahead=F) 
}, interval = 0.01, prof_output = "ice-prof")
profvis(expr = { for (i in 1:10){
  Xvalid_all(REmrt_GS_cpp, mf, 20, 10, 3, 4, 0.001, FALSE) }
}, interval = 0.01, prof_output = "ice-prof")
system.time(for (i in 1:10){
  Xvalid_all(REmrt_GS_cpp, mf, 10, 10, 3, 4, 0.001, FALSE) 
  })
