# load("mf")
# source("REmrtGS_cpp.R")
sourceCpp("partition.cpp")
# library(Rcpp)
# fit <- REmrt_GS_cpp(mf, maxL = 10, minbucket=3, minsplit=2, delQ=0.00001, lookahead=F)
# x <- fit
# newdata <- fit$data
# res <- REmrt.prednode_cpp(x, newdata)
REmrt.prednode_cpp <- function(x, newdata) {
  tt <- terms(x$data)
  new.ms <- model.frame(delete.response(tt), newdata)
  # the moderators in the new data
  old.ms <- model.frame(delete.response(tt), x$data)
  # the moderators in the old data
  tree <- x[["tree"]]
  nameM <- 1:ncol(new.ms)
  names(nameM) <- colnames(new.ms)
  inxM <- nameM[tree$mod]
  boolName <- sapply(old.ms, is.numeric) # tell if a moderator is numeric
  names(boolName) <- colnames(old.ms)
  boolNumeric <- boolName[tree$mod]
  partition(tree, new.ms, boolNumeric, inxM, x$cpt, old.ms)
}

# system.time( for (i in 1:100){
#   REmrt.prednode_cpp(x, newdata)
#   })
# system.time(for (i in 1:100){
#   REmrt.prednode(x, newdata)
#   })
# 
# 
# 
# 
# 
# 
#   tt <- terms(x$data)
#   new.ms <- model.frame(delete.response(tt), newdata)
#   old.ms <- model.frame(delete.response(tt), x$data)
#   tree <- x[["tree"]]
#   nameM <- 1:ncol(new.ms)
#   names(nameM) <- colnames(new.ms)
#   inxM <- nameM[tree$mod]
#   boolName <- sapply(old.ms, is.numeric)
#   names(boolName) <- colnames(old.ms)
#   boolNumeric <- boolName[tree$mod]
#   sourceCpp("partition.cpp")
#   set.seed(10)
#   K = 20
#   x1 <- rnorm(K)
#   x2 <- sample(letters[1:7], K, replace = TRUE)
#   x3 <- sample(0:1, K, replace = TRUE) 
#   x4 <- runif(K)
#   x5 <- sample(1:4, K, replace = TRUE)
#   newdata <- data.frame(x1=x1, x2=x2, x3= x3 , x4 =x4 ,x5=x5)
#   new.ms <- model.frame(delete.response(tt), newdata)
#   
#   
# 
#     res1 <- REmrt.prednode_cpp(x, newdata)
#     res2 <- REmrt.prednode(x, newdata)
#     for (i in 1:ncol(res1)){
#     print(all.equal(res1[,i], res2[,i]))
#       }
# 
#   
#   