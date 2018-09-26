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
set.seed(55)
x1 <- sample(0:1, replace = T, 100)
x2 <- sample(0:1, replace = T, 100)
x3 <- runif(100)
x4 <- runif(100)
# two-way interaction
y.grand <- 0.5 * (x2 > 0.5 & x1 < 0.5) + 0.5 * (x1 > 0.5 & x2 < 0.5)
tau <- 0.25
K <- 100
n <- NumGen(160, K)
y.true <- rnorm(rep(1,K), y.grand, tau)  # sample the true effect size for single study

cmk <- 1 - 3/(8*n-9)  #  approximation for
#gamma(n - 1) / (sqrt(n - 1) * gamma((2 * n - 3) / 2))
efk <- vark <- numeric(K)  # generate the effect size from a non-central t-distribution
for (k in 1:K) {
  samp <- rt(1, df=2*n[k]-2, ncp=y.true[k]*sqrt(n[k]/2))
  efk[k] <- samp*cmk[k]/sqrt(n[k]/2)
  vark[k] <- cmk[k]^2*(2*n[k]-2)*
    (1+n[k]*y.true[k]^2/2)/((2*n[k]-4)*n[k]/2)-y.true[k]^2  # note the vark is the sampling variance but not the sampling error
}

data.meta <- data.frame(efk, vark, x1, x2, x3, x4)
res <-metacartv2::REmrt(efk ~ x1 + x2 + x3 + x4, vi = vark, data = data.meta)
res
res2 <- metacartv2::REmrt(efk ~ x1 + x2 + x3 + x4, vi = vark, data = data.meta, lookahead = F)
datBump <- res2$data
save(datBump, file = "BumpingTestData")
rm(list = setdiff(ls(), list("res2")))
# bumping
bump <- function(B = 50, data, maxL) {
  mse <- numeric(B)
  fitB <- list()
  K <- nrow(data)
  for (b in 1:B) {
    inx.b <- sample(1:K, K, replace = T)
    data.b <- data[inx.b, ]
    fit.b <- metacartv2::REmrt_GS_cpp(data.b, maxL = maxL, minsplit = 5, minbucket = 3, delQ = 0.001 )
    allNodes <- metacartv2::prednode_cpp(fit.b, data)
    TNodes <- allNodes[, ncol(allNodes)] 
    depth <- nrow(fit.b$tree)
    tau2 <- fit.b$tree$tau2[depth]
    vi.star <- data.b$`(vi)` + tau2
    y <- data.b$efk
    subnodes <- c(t(fit.b$node.split[ncol(fit.b$node.split)]))
    wy.star <- y/vi.star
    n <- tapply(y, subnodes, length)
    g <- tapply(wy.star, subnodes, sum)/tapply(1/vi.star, subnodes, sum)
    fit.b$g <- g
    class(fit.b) <- "REmrt"
    predy.b <- predict(fit.b, data)$pred.g
    mse[b] <- mean((predy.b - data$efk)^2)
    mse.oob[b] <- mean((predy.b[-inx.b] - data$efk[-inx.b])^2)
    fitB[[b]] <- fit.b$tree
  }
  inxmin <- which.min(mse)
  list(min.mse = mse[inxmin], 
       mse.oob = mse.oob[inxmin], 
       model = fitB[[inxmin]])
}
bumpRes <- lapply(1:10, function(x) bump(50, res2$data, x))
mse <- sapply(bumpRes, function(x) x[[1]])
plot(1:10, mse, xlab = "# of splits", type = "b", main = "bumping")
mse.oob <- sapply(bumpRes, function(x) x[[2]])
plot(1:10, mse.oob, xlab = "# of splits", type = "b", main = "bumping")

bump.oob <- function(B = 50, data, maxL) {
  mse <- numeric(B)
  fitB <- list()
  K <- nrow(data)
  len <- numeric(B)
  for (b in 1:B) {
    inx.b <- sample(1:K, K, replace = T)
    data.b <- data[inx.b, ]
    fit.b <- metacartv2::REmrt_GS_cpp(data.b, maxL = maxL, minsplit = 5, minbucket = 3, delQ = 0.001 )
    allNodes <- metacartv2::prednode_cpp(fit.b, data)
    TNodes <- allNodes[, ncol(allNodes)] 
    depth <- nrow(fit.b$tree)
    tau2 <- fit.b$tree$tau2[depth]
    vi.star <- data.b$`(vi)` + tau2
    y <- data.b$efk
    subnodes <- c(t(fit.b$node.split[ncol(fit.b$node.split)]))
    wy.star <- y/vi.star
    n <- tapply(y, subnodes, length)
    g <- tapply(wy.star, subnodes, sum)/tapply(1/vi.star, subnodes, sum)
    fit.b$g <- g
    class(fit.b) <- "REmrt"
    predy.b <- predict(fit.b, data)$pred.g
    #mse[b] <- mean((predy.b - data.meta$efk)^2)
    # print(b)
    # print(mean((predy.b - data.meta$efk)^2))
    # print(mean((predy.b[-inx.b] - data.meta$efk[-inx.b])^2))
    # print(length(predy.b[-inx.b]))
    mse[b] <- mean((predy.b[-inx.b] - data$efk[-inx.b])^2)
    fitB[[b]] <- fit.b$tree
    len[b] <- length(predy.b[-inx.b])
  }
  inxmin <- which.min(mse)
  list(min.mse = min(mse), model = fitB[[inxmin]], len = len[inxmin])
}
bumpRes.oob <- lapply(1:10, function(x) bump.oob(50, res2$data, x))
mse.oob <- sapply(bumpRes.oob, function(x) x[[1]])
plot(1:10, mse.oob, xlab = "# of splits", type = "b", main = "bumping with oob")


for (b in 1:B) {
  inx.b <- sample(1:K, K, replace = T)
  data.b <- res2$data[inx.b, ]
  fit.b <- metacartv2::REmrt_GS_cpp(data.b, maxL = 3, minsplit = 5, minbucket = 3, delQ = 0.001 )
  allNodes <- metacartv2::prednode_cpp(fit.b, data.meta)
  TNodes <- allNodes[, ncol(allNodes)] 
  depth <- nrow(fit.b$tree)
  tau2 <- fit.b$tree$tau2[depth]
  vi.star <- data.b$`(vi)` + tau2
  y <- data.b$efk
  subnodes <- c(t(fit.b$node.split[ncol(fit.b$node.split)]))
  wy.star <- y/vi.star
  n <- tapply(y, subnodes, length)
  g <- tapply(wy.star, subnodes, sum)/tapply(1/vi.star, subnodes, sum)
  fit.b$g <- g
  class(fit.b) <- "REmrt"
  predy.b <- predict(fit.b, data.meta)$pred.g
  mse[b] <- mean((predy.b - data.meta$efk)^2)
  fitB[[b]] <- fit.b$tree
}
inxmin <- which.min(mse)
fitB[[inxmin]]
first.split = ifelse(sapply(fitB, function(x) x$mod[2]) %in% c("x1", "x2"), 2, 1)
plot(1:B, mse, type = "p", col = first.split )
lines(1:B, mse, col = "lightgrey")
mean(first.split == 2)



x1 <- sample(0:1, replace = T, 100)
x2 <- sample(0:1, replace = T, 100)
x3 <- runif(100)
x4 <- sample(0:1, replace = T, 100)
x5 <- runif(100)
x6 <- runif(100)
# three-way interaction
y.grand <- 0.5 * (x2 > 0.5 & x1 < 0.5 & x3 < 0.5) + 0.5 * (x1 > 0.5 & x2 < 0.5 & x3 > 0.5) 
tau <- 0.25
K <- 100
n <- NumGen(160, K)
y.true <- rnorm(rep(1,K), y.grand, tau)  # sample the true effect size for single study

cmk <- 1 - 3/(8*n-9)  #  approximation for
#gamma(n - 1) / (sqrt(n - 1) * gamma((2 * n - 3) / 2))
efk <- vark <- numeric(K)  # generate the effect size from a non-central t-distribution
for (k in 1:K) {
  samp <- rt(1, df=2*n[k]-2, ncp=y.true[k]*sqrt(n[k]/2))
  efk[k] <- samp*cmk[k]/sqrt(n[k]/2)
  vark[k] <- cmk[k]^2*(2*n[k]-2)*
    (1+n[k]*y.true[k]^2/2)/((2*n[k]-4)*n[k]/2)-y.true[k]^2  # note the vark is the sampling variance but not the sampling error
}

data.meta <- data.frame(efk, vark, x1, x2, x3, x4, x5 , x6)
res <-metacartv2::REmrt(efk ~ x1 + x2 + x3 + x4 + x5 + x6, vi = vark, data = data.meta, maxL = 10)
res$tree
res2 <- metacartv2::REmrt(efk ~ x1 + x2 + x3 + x4 + x5 + x6, vi = vark, data = data.meta, lookahead = T, maxL = 10)
res2$tree
# bumping
B <- 500
mse <- numeric(B)
fitB <- list()
for (b in 1:B) {
  inx.b <- sample(1:K, K, replace = T)
  data.b <- res2$data[inx.b, ]
  fit.b <- metacartv2::REmrt_GS_cpp(data.b, maxL = 10, minsplit = 5, minbucket = 3, delQ = 0.001 )
  allNodes <- metacartv2::prednode_cpp(fit.b, data.meta)
  TNodes <- allNodes[, ncol(allNodes)] 
  depth <- nrow(fit.b$tree)
  tau2 <- fit.b$tree$tau2[depth]
  vi.star <- data.b$`(vi)` + tau2
  y <- data.b$efk
  subnodes <- c(t(fit.b$node.split[ncol(fit.b$node.split)]))
  wy.star <- y/vi.star
  n <- tapply(y, subnodes, length)
  g <- tapply(wy.star, subnodes, sum)/tapply(1/vi.star, subnodes, sum)
  fit.b$g <- g
  class(fit.b) <- "REmrt"
  predy.b <- predict(fit.b, data.meta)$pred.g
  mse[b] <- mean((predy.b - data.meta$efk)^2)
  fitB[[b]] <- fit.b$tree
}
inx.min <- which.min(mse)
fitB[[inx.min]]
first.split = ifelse(sapply(fitB, function(x) x$mod[2]) %in% c("x1", "x2", "x3"), 2, 1)
plot(1:B, mse, type = "p", col = first.split )
lines(1:B, mse, col = "lightgrey")
mean(first.split == 2)
