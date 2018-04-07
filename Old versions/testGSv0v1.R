library(Rcpp)
source("GSupdate_Cpp.R")
source("GS_V1.R")
load("mf2")
g <- mf$efk
vi <- mf$`(vi)`
cnode <- sample(1:20, length(g), replace = T)
inx.s <- cnode == 1
x <- runif(sum(inx.s))
res1 <- re.cutoff_cpp(g, vi, x, inx.s, cnode, 1)
res2 <- re.cutoff.cpp(g, vi, x, inx.s, cnode, 1)
cnode <- sample(1, length(g), replace = T)
inx.s <- cnode == 1
 re.cutoff_cpp(g, vi, x, inx.s, cnode, 1)
 re.cutoff.root(g, vi, x, 1)
all(abs(res1 - res2) < 1e-10)
system.time( for (i in 1:1000){
  re.cutoff_cpp(g, vi, x, inx.s, cnode, 1)
})
system.time( for (i in 1:1000){
  re.cutoff.cpp(g, vi, x, inx.s, cnode, 1)
})
