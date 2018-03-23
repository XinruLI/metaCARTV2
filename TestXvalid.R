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
SimData <- function(beta, mods, formula,
                    K, n, tau){
  # a function to simulate meta-analysis data set with moderators  
  # corrected for non-central t
  #
  # Arguments :
  #      beta : a numercia vector which contains the overall effect size
  #             and other coefficients. The length should equal to the 
  #             column number of model matrix
  #      mods : the dataframe of moderators.
  #   formula : a formula object to specify the relationship between the 
  #             effect size and moderators  
  #         K : a scalar which is the number of studies
  #        n  : a numerical vector of length K, the  within-study sample
  #             size
  #       tau : a scalar, the effect sizes variance
  #     
  #
  # Returns   :
  # A data frame which is simulated by the given parameters
  mods <- mods
  beta <- beta
  formula <- as.formula(formula)
  K <- K
  n <- n
  tau <- tau
  # step 1
  X <- model.matrix(formula, data = mods)
  # step 2
  delta <- X%*%beta  # compute the average effect size corresponding to model and the coefficients
  # step 3
  deltak <- rnorm(rep(1,K), delta, tau)  # sample the true effect size for single study
  cmk <- 1 - 3/(8*n-9)  #  approximation for
  #gamma(n - 1) / (sqrt(n - 1) * gamma((2 * n - 3) / 2))
  efk <- vark <- numeric(K)  # generate the effect size from a non-central t-distribution
  for (k in 1:K) {
    samp <- rt(1, df=2*n[k]-2, ncp=deltak[k]*sqrt(n[k]/2))
    efk[k] <- samp*cmk[k]/sqrt(n[k]/2)
    vark[k] <- cmk[k]^2*(2*n[k]-2)*
      (1+n[k]*deltak[k]^2/2)/((2*n[k]-4)*n[k]/2)-deltak[k]^2  # note the vark is the sampling variance but not the sampling error
  }
  dat   <- data.frame(trail=1:K, efk=efk, vark=vark)
  return(dat)
}
PreData <- function(formula, data, vi, c = 1, maxL = 5L, minsplit = 2L, delQ = 0.001, minbucket = 3, n.fold = 10, lookahead = FALSE, ...){
  Call <- match.call()
  indx <- match(c("formula", "data", "vi"),
                names(Call), nomatch = 0L)
  if (indx[1] == 0L)
    stop("a 'formula' argument is required")
  if (indx[3] == 0L)
    stop("The sampling variances need to be specified")
  if (!is.logical(lookahead))
    stop("The 'lookahead' argument needs to be a logical value")
  if (maxL < 2 & (lookahead == TRUE) )
    stop("The 'maxL' should be at least 2 when applying look-ahead strategy")
  temp <- Call[c(1L, indx)]
  temp[[1L]] <- quote(stats::model.frame)
  mf <- eval.parent(temp)
  mf
}
set.seed(207)
K <- 120
n <- 80

source("REmrt_GS_update.R")
source("Xvalid_for_all.R")
# The same for both K = 20/80
gamma <- 0
i = 1
  x1 <- rnorm(K)
  x2 <- sample(letters[1:5], K, replace = TRUE)
  x3 <- sample(0:1, K, replace = TRUE) 
  x4 <- runif(K)
  x5 <- sample(1:4, K, replace = TRUE)
  mods <- data.frame(m1 = (x1 > -0.92)*1, m2 = (x2 %in% c("a","c","e"))*1,
                     m3 = x3, x4 = x4, x5 = x5)
  dat.sim1 <- SimData(c(0,0.8), mods, "~ m1:m2:m3", K, NumGen(n,K), 0)
  dat.sim1 <- data.frame(dat.sim1, x1, x2, x3, x4 , x5)
  mf <- PreData(efk~x1+x2+x3+x4+x5, dat.sim1, vi = vark)
  res.up <- REmrt_GS(mf, maxL = 5, minbucket=5, minsplit=6, delQ=5, lookahead=F)
  res.up$tree
  REmrt.xvalid(REmrt_GS, mf, n.fold = 10, maxL = 10, minbucket=5, minsplit=6, delQ=1, lookahead=F)
