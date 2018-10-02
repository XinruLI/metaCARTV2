library(metacartv2)
SimDataForBump <- function(y.grand, n.ave, tau){
  K <- length(y.grand)
  NumGen <- function(n, k){
    # a function to generate the within-study sample sizes based on the average sample size
    #
    # Arguments:
    #   n: the average sample size, needs to be a scalar
    #   k: the number of sudies, needs to be a scalar
    #
    # Returns:
    # a numeric vector of length k
    n <- n
    k <- k
    samp <- rnorm(k, n, n/3)
    samp <- round(samp, 0)  # make the generated number integer
    samp[samp < 10] <- 10  # avoid negative, 0 or too small sample size that is not realistic
    return(samp)
  }
  n <- NumGen(n.ave, K)
  y.true <- rnorm(rep(1,K), y.grand, tau)  # sample the true effect size for single study
  
  cmk <- 1 - 3/(8*n-9)  #  approximation for
  #gamma(n - 1) / (sqrt(n - 1) * gamma((2 * n - 3) / 2))
  g <- vi <- numeric(K)  # generate the effect size from a non-central t-distribution
  for (k in 1:K) {
    samp <- rt(1, df=2*n[k]-2, ncp=y.true[k]*sqrt(n[k]/2))
    g[k] <- samp*cmk[k]/sqrt(n[k]/2)
    vi[k] <- cmk[k]^2*(2*n[k]-2)*
      (1+n[k]*y.true[k]^2/2)/((2*n[k]-4)*n[k]/2)-y.true[k]^2  # note the vi is the sampling variance but not the sampling error
  }
  data.frame(g = g, vi = vi)
  
}
BumpNestingXvalid <- function(formula, vi, data, maxL, 
                              c = 1, minbucket = 3L, cp = 0.01,
                              n.fold, B){
  mse <- numeric(B)
  fitB <- list()
  K <- nrow(data)
  # numNodes <- numeric(B)
  
  # original data
  fit.o <- metacartv2::REmrt(formula = formula, vi = vi, n.fold = n.fold, 
                             data = data, maxL = maxL, c = c, 
                             minbucket = minbucket, cp = cp)
  y = predict(fit.o, data)[[1]]
  mse[1] <- mean((y-model.response(fit.o$data))^2)
  fitB[[1]] <- fit.o
  #numNodes[[1]] <- length(fit.o$n)
  for (b in 2:B) {
    inx.b <- sample(1:K, K, replace = T)
    data.b <- data[inx.b, ]
    fit.b <- metacartv2::REmrt(formula = formula, vi = vi, n.fold = n.fold, 
                               data = data.b, maxL = maxL, c = c, 
                               minbucket = minbucket, cp = cp)
    y <- predict(fit.b, data)[[1]]
    mse[b] <- mean((y-model.response(fit.o$data))^2)
    fit.b$data <- NULL
    fitB[[b]] <- fit.b
    #numNodes[[b]] <- length(fit.b$n)
  }
  inxmin <- which.min(mse)
  fitB[[inxmin]]
}

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
bump <- function(mf, maxL, minsplit, minbucket, 
                 cp,lookahead, B) {
  # maxL: the number of splits
  # data: model.data.frame
  # return: the model resulted from bumping
  mse <- numeric(B)
  fitB <- list()
  K <- nrow(mf)
  # put in the original data
  fit.o <- metacartv2::REmrt_GS_(mf, maxL = maxL, minsplit = 5, minbucket = 3, cp = 0.001,lookahead = F )
  fit.o <- compute_subgrp_es(fit.o)
  nodes <- metacartv2::prednode_cpp(fit.o, mf)[,nrow(fit.o$tree)]
  y <- fit.o$g[as.character(nodes)]
  names(y) <- NULL
  mse[1] <- mean((y-model.response(mf))^2)
  fitB[[1]] <- fit.o
  for (b in 2:B) {
    inx.b <- sample(1:K, K, replace = T)
    mf.b <- mf[inx.b, ]
    fit.b <- metacartv2::REmrt_GS_(mf.b, maxL = maxL, minsplit = 5, minbucket = 3, cp = 0.001, lookahead = F)
    fit.b <- compute_subgrp_es(fit.b)
    nodes <- metacartv2::prednode_cpp(fit.b, mf)[,nrow(fit.b$tree)]
    y <- as.vector(fit.b$g[as.character(nodes)])
    mse[b] <- mean((y-model.response(mf))^2)
    fitB[[b]] <- fit.b
  }
  inxmin <- which.min(mse)
  fitB[[inxmin]]
}
XerrorBump <- function(mf, maxL, minsplit, minbucket, 
                       cp,lookahead,
                       B, n.fold){
  N <- nrow(mf)
  inx <- sample(1:N) # randomly shuffle the rows
  names(inx) <- as.numeric(cut(1:N, b = 10))
  predy <- numeric(N)
  for (i in 1:n.fold) {
    inx.test <- inx[which(names(inx) == i)]
    train <- mf[-inx.test, ]
    test <- mf[inx.test, ]
    train.fit <- bump(mf = train, maxL = maxL, minsplit = minsplit, 
                      minbucket = minbucket, 
                      cp = cp,lookahead = lookahead,
                      B = B)
    nodes.test <- metacartv2::prednode_cpp(train.fit, test)[,nrow(train.fit$tree)]
    y.test <- train.fit$g[as.character(nodes.test)]
    predy[inx.test] <- y.test
  }
  y <- model.response(mf)
  c(sum((y-predy)^2)/sum((y-mean(y))^2),
    sqrt(sum(((y-predy)^2-mean((y-predy)^2))^2))/sum((y-mean(y))^2))
  
}
XcrossNestingBump <- function(formula, vi, data, maxL, minsplit = 6, minbucket = 3, 
                              cp = 0.01, lookahead = FALSE, n.fold = 10, B = 25){
  Call <- match.call()
  indx <- match(c("formula", "data", "vi"),
                names(Call), nomatch = 0L)
  temp <- Call[c(1L, indx)]
  temp[[1L]] <- quote(stats::model.frame)
  mf <- eval.parent(temp)
  res <- sapply(1:maxL, function(x) XerrorBump(mf = mf, maxL = maxL, minsplit = minsplit, 
                                               minbucket = minbucket, 
                                               cp = cp,lookahead = lookahead,
                                               B = B, n.fold = n.fold))
  res <- cbind(1:maxL, t(res))
  
  colnames(res) <-c("no.split", "xerror", "sd.xerror")
  res
  
}



# 
# x1 <- sample(0:1, replace = T, 100)
# x2 <- sample(0:1, replace = T, 100)
# x3 <- runif(100)
# x4 <- runif(100)
# y.grand <- ((x1 + x2) == 1)* 0.5
# dat <- data.frame(SimDataForBump(y.grand, n.ave = 80, tau = 0.1),
#            x1, x2, x3, x4)
# res <- BumpNestingXvalid(g ~ x1 + x2 + x3 + x4, vi = vi, data = dat,
#                   maxL  = 8, n.fold = 10, B = 25, c = 1)
# 
# res <- XcrossNestingBump(g ~ x1 + x2 + x3 + x4, vi = vi, data = dat,
#                          maxL  = 8, n.fold = 10, B = 25, c = 1)
# mutate(as.data.frame(res),
#        within1SE = xerror <= min(xerror) + sd.xerror[which.min(xerror)])
# L = 2
# 
# mf <- REmrt(g ~ x1 + x2 + x3 + x4, vi = vi, data = dat,
#             maxL  = 1, n.fold = 2)$data
# bump(mf, maxL = maxL, minsplit = 6, minbucket = 3, 
#                  cp = 0.01, lookahead = F, 25) 
