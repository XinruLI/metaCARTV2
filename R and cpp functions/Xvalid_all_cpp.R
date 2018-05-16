library(Rcpp)
source("REmrtGSV1_cpp.R")
source("Pred_Node_cpp.R")
#source("REmrt0.R")
#source("Xvalid_for_all.R")
sourceCpp("testxvalid.cpp")
Xvalid_all <- function(func, mf, maxL, n.fold, minbucket, minsplit, delQ, lookahead){
  N <- nrow(mf)
  if (n.fold > N | n.fold <=1 ) stop("n.fold is not acceptable")
  if (maxL > N) {
    warning("The maximum number of split is too large")
    maxL <- N-1
  }
  
  pred <- matrix(NA, nrow = N, ncol = maxL+1)
  inx <- sample(1:N)
  inx.xvalid <- c(round(seq(from = 1, by= N/n.fold, length.out = n.fold)), N+1)
  for (i in 1:n.fold){
    inx.test <- inx[inx.xvalid[i]:(inx.xvalid[i+1]-1)]
    test <- mf[inx.test, ]
    train <- mf[-inx.test, ]
    fit.train <- do.call(func, list(mf = train, maxL, minbucket, minsplit, delQ, lookahead))
    yi.train <- model.response(fit.train$data)
    vi.train <- c(t(fit.train$data["(vi)"]))
    tau2 <- fit.train$tree$tau2
    nsplt <- nrow(fit.train$tree)
    train.y <- ComputeY(fit.train$node.split, yi.train, vi.train,tau2) 
    test.nodes_cpp <- REmrt.prednode_cpp(fit.train, test)
    test.y <- PredY(train.y, test.nodes_cpp)
    if (any(is.na(test.y))) {
      inx.NA <- which(is.na(test.y), arr.ind=TRUE)
      test.y <- ReplaceNA(inx.NA, test.y, yi.train, vi.train, tau2)
    }
    pred[inx.test, 1:nsplt] <- test.y
  }
  y <- model.response(mf)
  if (!is.null(dim(pred))) {
    x.error <- apply(pred,2, function(x) sum((y-x)^2)/sum((y-mean(y))^2))
    sdx.error <- apply(pred, 2, function(x)  sqrt(sum(((y-x)^2-mean((y-x)^2))^2))/sum((y-mean(y))^2))
  } else {
    x.error <- sum((y-pred)^2)/sum((y-mean(y))^2)
    sdx.error <- sqrt(sum(((y-pred)^2-mean((y-pred)^2))^2))/sum((y-mean(y))^2)
  }
  cbind(x.error, sdx.error)
}


# load("mf")
# set.seed(1)
# res1 <- Xvalid_all(REmrt_GS_cpp, mf, 10, 10, 3, 4, 0.001, FALSE)
# set.seed(1)
# res2 <- REmrt.xvalid(REmrt_GS_cpp, mf, 10, 10, 3, 4, 0.001, FALSE)
# all.equal(res1, res2)
# system.time(Xvalid_all(REmrt_GS_cpp, mf, 10, 10, 3, 4, 0.001, FALSE))
# # #  0.690 
#  system.time(REmrt.xvalid(REmrt.fit0, mf, 10, 10, 3, 4, 0.001, FALSE))
# # #  8.328
# # 
# n.fold = 10
# maxL = 10
# func <- REmrt_GS_cpp
# N <- nrow(mf)
# minbucket = 3
# minsplit = 4
# delQ = 0.0001
# lookahead = FALSE
# if (n.fold > N | n.fold <=1 ) stop("n.fold is not acceptable")
# if (maxL > N) {
#   warning("The maximum number of split is too large")
#   maxL <- N-1
# }
# pred <- matrix(NA, nrow = N, ncol = maxL+1)
# inx <- sample(1:N)
# inx.xvalid <- c(round(seq(from = 1, by= N/n.fold, length.out = n.fold)), N+1)
# i=1
# inx.test <- inx[inx.xvalid[i]:(inx.xvalid[i+1]-1)]
# test <- mf[inx.test, ]
# train <- mf[-inx.test, ]
# fit.train <- do.call(func, list(mf = train, maxL, minbucket, minsplit, delQ, lookahead))
# yi.train <- model.response(fit.train$data)
# vi.train <- c(t(fit.train$data["(vi)"]))
# nodes <- fit.train$node.split
# 
# nsplt <- nrow(fit.train$tree)
# train.wy <- sapply(1:nsplt, function(x) yi.train/(vi.train+fit.train$tree[,2][x]))
# train.wts <-  sapply(1:nsplt, function(x) 1/(vi.train+fit.train$tree[,2][x]))
# train.pred <- sapply(1:nsplt, function(x) tapply(train.wy[, x], fit.train$node.split[, x], sum)/tapply(train.wts[, x], fit.train$node.split[, x], sum))
# 
# 
# train.y <- ComputeY(nodes, yi.train, vi.train,tau2) 
# test.nodes_cpp <- REmrt.prednode_cpp(fit.train, test)
# 
# test.nodes_cpp[1:5, 1:5] <- NA
# test.nodes_cpp
# 
# test.y
# test.pred <- lapply(1:nsplt, function(x) train.pred[[x]][as.character(test.nodes_cpp[ ,x])])
# if (any(is.na(test.y))) {
#   inx.NA <- which(is.na(test.y), arr.ind=TRUE)
# }
# 
# test.y.rmna <- ReplaceNA(inx.NA, test.y, yi.train, vi.train, tau2)
# test.y.rmna
# test.pred.rmna <- sapply(1:nsplt, function(x){test.pred[[x]][is.na(test.pred[[x]])] <- sum(train.wy[, x])/sum(train.wts[, x]);test.pred[[x]]})
# test.y
# for (i in 1:nsplt){
#   if(all(abs(train.y[[i]] - train.pred[[i]]) < 1e-12))  {print(i)}
# }
# for (i in 1:nsplt){
#   if(all(abs(test.y[,i] - test.pred[[i]]) < 1e-12))  {print(i)}
# }
# 
# for (i in 1:nsplt){
#   if(all(abs(test.y.rmna[,i] - test.pred.rmna [,i]) < 1e-12))  {print(i)}
# }
# 
# for (i in 1:ncol(test.nodes_cpp)){
#   print(i)
#        print(all.equal(test.nodes[,i], test.nodes_cpp[,i]))
# }
# 
# 
