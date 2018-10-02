load("BumpingTestData")
library(metacartv2)
compute_subgrp_es <- function(fit){
  # given an initial tree model 
  # compute the subgroup effect sizes
  # add fit$g to the initial tree
  allNodes <- metacartv2::prednode_cpp(fit, fit$data)
  TNodes <- allNodes[, ncol(allNodes)] 
  depth <- nrow(fit$tree)
  tau2 <- fit$tree$tau2[depth]
  vi.star <- fit$data$`(vi)` + tau2
  y <- model.response(fit$data)
  subnodes <- fit$node.split[[ncol(fit$node.split)]]
  wy.star <- y/vi.star
  n <- tapply(y, subnodes, length)
  g <- tapply(wy.star, subnodes, sum)/tapply(1/vi.star, subnodes, sum)
  fit$g <- g
  fit
}
bump <- function(data, maxL, B = 50) {
  # maxL: the number of splits
  # data: model.data.frame
  # return: the model resulted from bumping
  mse <- numeric(B)
  fitB <- list()
  K <- nrow(data)
  # put in the original data
  fit.o <- metacartv2::REmrt_GS_(data, maxL = maxL, minsplit = 5, minbucket = 3, cp = 0.001,lookahead = F )
  fit.o <- compute_subgrp_es(fit.o)
  nodes <- metacartv2::prednode_cpp(fit.o, data)[,nrow(fit.o$tree)]
  y <- fit.o$g[as.character(nodes)]
  names(y) <- NULL
  mse[1] <- mean((y-model.response(data))^2)
  fitB[[1]] <- fit.o
  for (b in 2:B) {
    inx.b <- sample(1:K, K, replace = T)
    data.b <- data[inx.b, ]
    fit.b <- metacartv2::REmrt_GS_(data.b, maxL = maxL, minsplit = 5, minbucket = 3, cp = 0.001, lookahead = F)
    fit.b <- compute_subgrp_es(fit.b)
    nodes <- metacartv2::prednode_cpp(fit.b, data)[,nrow(fit.b$tree)]
    y <- as.vector(fit.b$g[as.character(nodes)])
    mse[b] <- mean((y-model.response(data))^2)
    fitB[[b]] <- fit.b
  }
  inxmin <- which.min(mse)
  fitB[[inxmin]]
}
# bumpling nested in Cross-validation

XerrorBump <- function(data, maxL, n.fold = 10, B = 50){
  N <- nrow(data)
  inx <- sample(1:N) # randomly shuffle the rows
  names(inx) <- as.numeric(cut(1:N, b = 10))
  predy <- numeric(N)
  for (i in 1:n.fold) {
    inx.test <- inx[which(names(inx) == i)]
    train <- data[-inx.test, ]
    test <- data[inx.test, ]
    train.fit <- bump(data = train, maxL = maxL, B)
    nodes.test <- metacartv2::prednode_cpp(train.fit, test)[,nrow(train.fit$tree)]
    y.test <- train.fit$g[as.character(nodes.test)]
    predy[inx.test] <- y.test
  }
  y <- model.response(data)
  c(sum((y-predy)^2)/sum((y-mean(y))^2),
    sqrt(sum(((y-predy)^2-mean((y-predy)^2))^2))/sum((y-mean(y))^2))
  
}
XcrossNestingBump <- function(data, maxL, n.fold, B){
  res <- sapply(1:maxL, function(x) XerrorBump(data = data, 
                                                       x, n.fold = n.fold, 
                                                       B = B))
  res <- cbind(1:maxL, t(res))
  colnames(res) <-c("no.split", "xerror", "sd.xerror")
  res

}
xres <- XcrossNestingBump(datBump, 5, 10, 25)
library(dplyr)
mutate(as.data.frame(xres),
       within1SE = xerror <= min(xerror) + sd.xerror[which.min(xerror)])
L  = 2
bump(datBump, L, 50)$tree

# cross-validation nested in bumping
formula <- as.formula(efk ~ x1 + x2 + x3 + x4)
data <- datBump
data$vi <- datBump$`(vi)`
maxL  = 5
n.fold = 10
B = 25
BumpNestingXvalid <- function(formula, vi, data, maxL, n.fold, B){
  mse <- numeric(B)
  fitB <- list()
  K <- nrow(data)
  # original data
  fit.o <- metacartv2::REmrt(formula = formula, vi = vi, n.fold = n.fold, 
                    data = data, maxL = maxL)
  y = predict(fit.o, data)[[1]]
  mse[1] <- mean((y-model.response(data))^2)
  fit.o$data <- NULL
  fitB[[1]] <- fit.o
  for (b in 2:B) {
    inx.b <- sample(1:K, K, replace = T)
    data.b <- data[inx.b, ]
    fit.b <- metacartv2::REmrt(formula = formula, vi = vi, n.fold = n.fold, 
                               data = data.b, maxL = maxL)
    y <- predict(fit.b, data)[[1]]
    mse[b] <- mean((y-model.response(data))^2)
    fit.b$data <- NULL
    fitB[[b]] <- fit.b
  }
  inxmin <- which.min(mse)
  fitB[[inxmin]]
}

datBump$vi <- datBump$`(vi)`
res <- BumpNestingXvalid(efk ~ x1 + x2 + x3 + x4, vi = vi, data = datBump,
                         maxL  = 5, n.fold = 10, B = 25)
res
