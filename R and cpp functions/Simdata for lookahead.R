set.seed(1637)
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
K <- 60
x1 <- sample(0:1, replace = T,K)
x2 <- sample(0:1, replace = T,K)
x3 <- runif(K)
x4 <- sample(0:1, replace = T,K)

y.grand <- 0.5 * (x2 > 0.5 & x1 < 0.5) + 0.5 * (x1 > 0.5 & x2 < 0.5)
tau <- 0.1
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
x1 <- factor(x1)
x2 <- factor(x2)
x4 <- factor(x4)
dat.balanced <- data.frame(g = efk, vi = vark, x1, x2, x3, x4)
set.seed(1)
res <-metacartv2::REmrt(efk ~ x1 + x2 + x3 + x4, vi = vark, data = data.meta, c = 1,
                        maxL = 5L, minsplit = 6L, cp = 1e-04,
                        minbucket = 3L, xval = 10, lookahead = FALSE)
metacartv2::REmrt_GS_cpp2(res$data, maxL = 5L, minsplit = 6L, minbucket = 3L,delQ = 0,lookahead = FALSE)$tree
res2 <- metacartv2::REmrt(efk ~ x1 + x2 + x3 + x4, vi = vark, data = data.meta, lookahead = T,
                          c = 1,
                          maxL = 5L, minsplit = 6L, cp = 1e-04,
                          minbucket = 3L, xval = 10)
plot(res2)

inx <- sample(1:K)
inx.xvalid <- c(round(seq(from = 1, by= K/60, length.out = 60)), K+1)

first.split1 <- first.split2 <- NULL
for (i in 1:60){
  inx.test <- inx[inx.xvalid[i]:(inx.xvalid[i+1]-1)]
  train <- res$data[-inx.test, ]
  temp1 <- metacartv2::REmrt_GS_cpp2(train, maxL = 5L, minsplit = 6L, minbucket = 3L,delQ = 0,lookahead = FALSE)
  first.split1 <- c(first.split1, temp1$tree$mod[2])
  temp2 <- metacartv2::REmrt_GS_cpp2(train, maxL = 5L, minsplit = 6L, minbucket = 3L,delQ = 0,lookahead = T)
  first.split2 <- c(first.split2, temp2$tree$mod[2])

}

