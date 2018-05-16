# THE SSS FUNCTION
expit <- function(x) (tanh(x/2)+1)/2  	# CORRECT & VERIFIED
# THE SSS FUNCTION TO COMPUTE QB WITH GLOBAL TAU2 
REbetQ.ord <- function(c, a = 50, x, g, vi, cnode, pl) {
  # Arguments:
  # c: the split point
  # a: the scale parameter
  # x: the moderator
  # g: a numeric vector of effect sizes
  # vi: a numeric vertor of sampling variance
  # Return: the opposite number of the between-subgroup Q
  #         the opposite sign is for the optimize function
  g.split <- g[cnode == pl ]
  vi.split <- vi[cnode == pl ]
  inx.left <- cnode != pl
  theta.r <- sum(g.split/vi.split*expit(a*(x-c)))/sum(expit(a*(x-c))/vi.split)
  theta.l <- sum(g.split/vi.split*(1-expit(a*(x-c))))/sum((1-expit(a*(x-c)))/vi.split)
  q.left <- tapply(g[inx.left]^2/vi[inx.left], cnode[inx.left], sum) -
    tapply(g[inx.left]/vi[inx.left], cnode[inx.left], function(x) (sum(x))^2)/
    tapply(1/vi[inx.left], cnode[inx.left], sum)
  q.r <- sum((g.split-theta.r)^2/vi.split*expit(a*(x-c))) 
  q.l <-sum((g.split-theta.l)^2/vi.split* (1- expit(a*(x-c))))
  df <- length(g) - length(unique(cnode)) - 1
  c.rl <- sum(1/vi) - sum(expit(a*(x-c))/vi.split^2)/sum(expit(a*(x-c))/vi.split) -
    sum((1-expit(a*(x-c)))/vi.split^2)/sum((1-expit(a*(x-c)))/vi.split) - 
    sum(tapply(1/vi[inx.left]^2, cnode[inx.left], sum)/tapply(1/vi[inx.left], cnode[inx.left], sum))
  tau2 <- max(c(0, (q.r + q.l + sum(q.left) - df)/ c.rl))
  wts <- 1/(vi + tau2)
  wts.split <- wts[cnode == pl]
  g.mu <- sum(g*wts)/sum(wts)
  theta.rstar <- sum(g.split*wts.split*expit(a*(x-c)))/sum(expit(a*(x-c))*wts.split)
  theta.lstar <- sum(g.split*wts.split*(1-expit(a*(x-c))))/sum((1-expit(a*(x-c)))*wts.split)
  q.within <- sum((g.split-theta.rstar)^2*wts.split*expit(a*(x-c))) +
    sum((g.split-theta.lstar)^2*wts.split* (1- expit(a*(x-c)))) +
    sum(tapply(g[inx.left]^2 * wts[inx.left], cnode[inx.left], sum) - 
          tapply(g[inx.left] * wts[inx.left], cnode[inx.left], function(x) (sum(x))^2)/
          tapply(wts[inx.left], cnode[inx.left], sum))
  q.total <- sum((g-g.mu)^2*wts)
  return(q.within-q.total)
}
rebetQ<- function(yi, vi, mods){
  wts = 1/vi
  wy = wts*yi
  wy2 = wts * yi^2
  Q <- tapply(wy2, mods, sum) - tapply(wy, mods, function(x) (sum(x))^2)/tapply(wts, mods, sum)
  df <- tapply(wy, mods, length)-1
  C <- tapply(wts, mods, sum) - tapply(wts, mods, function(x) sum(x^2))/ tapply(wts, mods, sum)
  tau2 <- (sum(Q) - sum(df))/sum(C)
  tau2 <- max(0, tau2)
  wstar = 1/(vi+tau2)
  wystar = wstar*yi
  wy2star = wstar*yi^2
  Qstar <- tapply(wy2star, mods, sum) - tapply(wystar, mods, function(x) (sum(x))^2)/tapply(wstar, mods, sum)
  Qstar.total <- sum(wy2star) - (sum(wystar))^2/sum(wstar)
  Qbet <- Qstar.total - sum(Qstar)
  if (is.na(Qbet)) {
    Qbet <- Inf
  }
  return(c(Qbet, tau2))
  
}
# FIND THE SPLIT POINT BY MAXIMIZING QB
bestcut.ord <- function(x, g, vi, a=NULL,
                        cnode, pl,
                        alpha.endcut=.02,
                        fe = FALSE,
                        multi.start=T, n.starts=3) {
  # Arguments:
  # a: the scale parameter
  # x: the moderator
  # g: a numeric vector of effect sizes
  # vi: a numeric vertor of sampling variance
  # alpha.endcut: to define the range of the split point
  # mult.start: whether to use multiple starts for the optimize function
  # n.starts: the number of multiple starts
  n <- length(x)
  # FIX a
  if (is.null(a)) a <- sqrt(n)
  # FINDING THE SEARCH RANGE TO AVOID ENDCUT PREFERENCE PROBLEM
  sigma <- sd(x); mu <- mean(x)
  x <- scale(x)  # IMPORTANT TO STANDARDIZE x IN ORDER TO APPLY A CONSTANT a
  LB <- quantile(x, probs = alpha.endcut); UB <- quantile(x, probs =1-alpha.endcut); 
  if (fe == TRUE) { 
    obj.func <- FEbetQ.ord
  } else {
    obj.func <- REbetQ.ord
  }
  if (multi.start==T) {
    B <- seq(LB, UB, length.out=n.starts)
    Q.min <- 1e15
    for (b in 2:n.starts) {
      OPT <- optimize(obj.func, lower=LB, upper=UB, maximum=F, 
                      a = a, g = g, x = x, vi = vi,  cnode = cnode, pl = pl)
      if (OPT$objective < Q.min) {
        Q.min <- OPT$objective
        cstar <- OPT$minimum
      }
    }
  } else {
    cstar <- optimize(obj.func, lower=LB, upper=UB, maximum=F, 
                      a = a, g = g, x = x, vi = vi,  cnode = cnode, pl = pl)$minimum
  }
  cstar <- cstar*sigma + mu # TRANSFORM BACK
  return(cstar)
}
# THE FUNCTION TO FIT RE META-CART USING SSS STRATEGY
REmrt.SSS0 <- function(mf, maxL, minbucket, minsplit, delQ, lookahead){
  y <- model.response(mf)
  vi <- c(t(mf["(vi)"]))
  mods.names <-  labels(terms(mf))
  mods <- mf[mods.names]
  nmod <- ncol(mods)
  cpt <- list()
  nodemark <- data.frame(node = rep(1, nrow(mf)))
  res <- data.frame(Qb = rebetQ(y, vi, mods = nodemark)[1],
                    tau2 = rebetQ(y, vi, mods = nodemark)[2],
                    split = NA, mod = NA, pleaf = NA)
  delta.Q <- Inf
  
  for (i in 1:maxL) {
    Dev<- -Inf
    tres <- NULL
    cnode <- nodemark[ ,i]
    len.node <- tapply(vi, cnode, length)
    nodes <- names(len.node) [len.node >= minsplit]
    if (length(nodes) == 0) break
    for (j in 1:length(nodes)) {
      pleaf.inx <- cnode == as.numeric(nodes[j])
      for (k in 1:nmod){
        xk <- mods[pleaf.inx, k]
        c.splits <- unique(xk)
        if (length(c.splits) < 2) { # MODERATORS CANNOT BE SPLIT
          next
        } else {
          #----------------------DICHOTOMOUS MODERATOR-------------------------#
          if (length(c.splits) == 2) { 
            c.star <- sort(xk)[1]
            cnode.tmp <- cnode
            cnode.tmp[pleaf.inx] <- ifelse( xk == c.star, 2*i, 2*i+1)
            if (min(table(cnode.tmp[pleaf.inx])) <= minbucket) {
              Dev.new <- -Inf
            } else {
              temp <- rebetQ(y, vi, mods = as.factor(cnode.tmp))
              Dev.new <- temp[1]
            }
            if (Dev.new > Dev) {
              Dev <- temp[1]
              tcpt <- c.star
              msplit <- paste(mods.names[k], "=", tcpt, collapse = " ")
              tres <- data.frame(Qb = temp[1], tau2 = temp[2],
                                 split = msplit, mod = mods.names[k], pleaf = as.numeric(nodes[j]))
              tnode <- cnode.tmp
            }
          } else {
            #--------------------------ORDINAL MODERATOR-------------------------#             
            if (is.numeric(xk) == TRUE) {
              pl <- nodes[j]
              c.star <- bestcut.ord(x = xk, a = 50, g= y, vi=vi, 
                                    cnode = cnode, pl = pl)
              cnode.tmp <- cnode
              cnode.tmp[pleaf.inx] <- ifelse( xk <= c.star, 2*i, 2*i+1)
              if (min(table(cnode.tmp[pleaf.inx])) <= minbucket) {
                Dev.new <- -Inf
              } else {
                temp <- rebetQ(y, vi, mods = as.factor(cnode.tmp))
                Dev.new <- temp[1]
              }
              if (Dev.new > Dev) {
                Dev <- temp[1]
                tcpt <- c.star
                msplit <- paste(mods.names[k], "<=", tcpt, collapse = " ")
                tres <- data.frame(Qb = temp[1], tau2 = temp[2],
                                   split = msplit, mod = mods.names[k], pleaf = as.numeric(nodes[j]))
                tnode <- cnode.tmp
              }
            } else {
              
              #-----------------------CATEGORICAL MODERATOR-------------------------#               
              
              pl <- nodes[j]
              xk.rank <- rank(tapply(y[pleaf.inx], xk, mean))
              xk.ord <- xk.rank[as.character(xk)]
              c.star <- bestcut.ord(x = xk.ord, a = 50, g= y, vi=vi, 
                                    cnode = cnode, pl = pl)
              cnode.tmp <- cnode
              cnode.tmp[pleaf.inx] <- ifelse( xk.ord <= c.star, 2*i, 2*i+1)
              if (min(table(cnode.tmp[pleaf.inx])) <= minbucket) {
                Dev.new <- -Inf
              } else {
                temp <- rebetQ(y, vi, mods = as.factor(cnode.tmp))
                Dev.new <- temp[1]
              }
              if (Dev.new > Dev) {
                Dev <- temp[1]
                tcpt <- names(xk.rank[xk.rank <= c.star])
                msplit <- paste(mods.names[k], "=", paste(tcpt, collapse = "/"), collapse = " ")
                tres <- data.frame(Qb = temp[1], tau2 = temp[2],
                                   split = msplit, mod = mods.names[k], pleaf = as.numeric(nodes[j]))
                tnode <- cnode.tmp
              }
              
            }
          }
        }
      }
    }
    if (is.null(tres)) {
      break
    } else {
      delta.Q <- abs(tres$Qb[1] - res$Qb[i])
      if (delta.Q < delQ) break
      new.node <- tnode
      nodemark <- cbind(nodemark, new.node)
      res <- rbind(res, tres)
      cpt[[i]] <- tcpt
    }
  }
  
  
  
  list(tree = res, node.split = nodemark, cpt = cpt, data = mf)
  
  
}

