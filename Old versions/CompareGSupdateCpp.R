source("GSupdate_Cpp.R")
source("GSupdate.R")
load("mf2")
g <- mf$efk
vi <- mf$`(vi)`
cnode <- sample(1:4, length(g), replace = T)
inx.s <- cnode == 1
x <- runif(sum(inx.s))
res1 <- re.cutoff(g, vi, x, inx.s, cnode, 1)
res2 <- re.cutoff.cpp(g, vi, x, inx.s, cnode, 1)
res1
res2
system.time( for (i in 1:1000){
  re.cutoff(g, vi, x, inx.s, cnode, 1)
  }) #2.954
system.time(for (i in 1:1000){
  re.cutoff.cpp(g, vi, x, inx.s, cnode, 1)
  }) #0.173
