load("mf")
source("REmrt_SSS.R")
source("REmrt0.R")
source("REmrt_GS_update.R")
library("profvis")
profvis(expr = {
  REmrt_GS(mf, maxL = 10, minbucket=3, minsplit=2, delQ=0.00001, lookahead=F) 
}, interval = 0.01, prof_output = "ice-prof")
res <- REmrt.SSS0(mf, maxL = 5, minbucket=3, minsplit=2, delQ=0.00001, lookahead=F) 
res$tree
res0 <- REmrt.fit0(mf, maxL = 5, minbucket=3, minsplit=2, delQ=0.00001, lookahead=F) 
res0$tree
res$tree

system.time(
  for (i in 1:100) re.cutoff(mf$efk, mf$`(vi)`, mf$x1, cnode = res0$node.split[,5], pl = "9")
)
system.time(
  for (i in 1:100) GS.fit0(mf$efk, mf$`(vi)`, mf$x1, cnode = res0$node.split[,5], pl = "9")
)
