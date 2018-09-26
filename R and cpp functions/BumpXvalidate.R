load("BumpingTestData")
compute_subgrp_es <- function(fit){
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
  fit.o <- metacartv2::REmrt_GS_cpp(data, maxL = maxL, minsplit = 5, minbucket = 3, delQ = 0.001 )
  fit.o <- compute_subgrp_es(fit.o)
  nodes <- metacartv2::prednode_cpp(fit.o, data)[,nrow(fit.o$tree)]
  y <- fit.o$g[as.character(nodes)]
  names(y) <- NULL
  mse[1] <- mean((y-model.response(data))^2)
  fitB[[1]] <- fit.o
  for (b in 2:B) {
    inx.b <- sample(1:K, K, replace = T)
    data.b <- data[inx.b, ]
    fit.b <- metacartv2::REmrt_GS_cpp(data.b, maxL = maxL, minsplit = 5, minbucket = 3, delQ = 0.001 )
    fit.b <- compute_subgrp_es(fit.b)
    nodes <- metacartv2::prednode_cpp(fit.b, data)[,nrow(fit.b$tree)]
    y <- as.vector(fit.b$g[as.character(nodes)])
    mse[b] <- mean((y-model.response(data))^2)
    fitB[[b]] <- fit.b
  }
  inxmin <- which.min(mse)
  fitB[[inxmin]]
}
# Cross-validation
data = datBump
maxL  = 2
n.fold = 10
B = 20

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
CrossBump <- function(data, maxL, n.fold, B){
  res <- sapply(1:maxL, function(x) CrossValidWithBump(data = data, 
                                                       x, n.fold = n.fold, 
                                                       B = B))
  res <- cbind(1:maxL, t(res))
  colnames(res) <-c("no.split", "xerror", "sd.xerror")
  res

}
xres <- CrossBump(datBump, 8, 10, 30)
transmute(as.data.frame(xres), plussd =  xerror + sd.xerror)
L  = 5
bump(datBump, 5, 50)$tree
