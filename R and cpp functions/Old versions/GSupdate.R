# GREEDY SEARCH WITH UPDATE SCHEME FOR RE MODEL
# TO DO:
# FIX MINSPLIT, MINBUCKET PROBLEM
# FROM THE ROOT NODE
re.cutoff.root <- function(g, vi, x, minbucket) {
  n <- length(g)
  x.sort <- sort(x)
  c.split <- (x.sort[-1] - x.sort[-n]) != 0
  if (minbucket > 1) {
    c.split[c(1:(minbucket-1), (n-minbucket+1):(n-1))] <- FALSE
  }
  
  if (sum(c.split) == 0) {
    return(NULL)
  } else {
    g.sort <- g[order(x)]
    vi.sort <- vi[order(x)]
    df <- length(g) - 2
    wy <- g.sort / vi.sort
    wy2 <- g.sort^2 / vi.sort
    wts <- 1/vi.sort
    wts2 <- wts^2
    cwy <- cumsum(wy[-n])
    cwy2 <- cumsum(wy2[-n])
    cwts <- cumsum(wts[-n])
    cwts2 <- cumsum(wts2[-n])
    Qrl <- sum(wy2) - cwy^2/cwts - (sum(wy) - cwy)^2/(sum(wts) - cwts)
    C.tau2 <- sum(1/vi) - cwts2/cwts - (sum(wts2) - cwts2)/(sum(wts) - cwts) 
    tau2 <- pmax(0, (Qrl - df)/C.tau2)
    # NOT POSSIBLE TO USE CUMSUM, TRY USING MATRIX
    w.star <- 1/ (vi.sort[-n] %*% t(rep(1,n-1)) + rep(1, n-1) %*% t(tau2))
    swts <- colSums(w.star) + 1/(vi.sort[n] + tau2) # half time cost compared to sapply
    w.temp <- w.star
    w.temp[lower.tri(w.temp)] <- 0
    cwts0 <- colSums(w.temp) 
    wy.star <- (g.sort[-n] %*% t(rep(1,n-1))) * w.star
    swy <- colSums(wy.star) + g.sort[n]/(vi.sort[n] + tau2) 
    wy.temp <- wy.star
    wy.temp[lower.tri(wy.temp)] <- 0
    cwy0 <- colSums(wy.temp)
    qb.star <-  (swy - cwy0)^2/(swts - cwts0) +
      cwy0^2/cwts0 - swy^2/swts 
    inx.star <- which(qb.star ==  max(qb.star[c.split]))
    cstar <- x.sort[inx.star]
    res <- c(cstar, qb.star[inx.star], tau2[inx.star])
    res
  }
} 
# AT LEAST ONE NON-PARENT LEAVES
re.cutoff <- function(g, vi, x, inx.s, cnode, minbucket) {
  n <- sum(inx.s)
  x.sort <- sort(x)
  c.split <- (x.sort[-1] - x.sort[-n]) != 0
  if (minbucket > 1) {
    c.split[c(1:(minbucket-1), (n-minbucket+1):(n-1))] <- FALSE
  }
  if (sum(c.split) == 0) {
    return(NULL)
  } else {
    g.sort <- g[inx.s][order(x)]
    vi.sort <- vi[inx.s][order(x)]
    inx.left <- inx.s == FALSE
    g.left <- g[inx.left]
    vi.left <- vi[inx.left]
    cnode.left <- cnode[inx.left]
    q.left <- sum(tapply(g.left^2/vi.left, cnode.left, sum) -
                    tapply(g.left/vi.left, cnode.left, function(x) (sum(x))^2)/
                    tapply(1/vi.left, cnode.left, sum))
    C.left <- tapply(1/vi.left^2, cnode.left, sum)/tapply(1/vi.left, cnode.left, sum)
    df <- length(g) - length(unique(cnode))-1
    wy <- g.sort / vi.sort
    wy2 <- g.sort^2 / vi.sort
    wts <- 1/vi.sort
    wts2 <- wts^2
    cwy <- cumsum(wy[-n])
    cwy2 <- cumsum(wy2[-n])
    cwts <- cumsum(wts[-n])
    cwts2 <- cumsum(wts2[-n])
    Qrl <- sum(wy2) - cwy^2/cwts - (sum(wy) - cwy)^2/(sum(wts) - cwts)
    # q.l <- cwy2 - cwy^2/cwts # VERIFIED
    # q.r <- sum(wy2)-cwy2 - (sum(wy) - cwy)^2/(sum(wts) - cwts) # VERIFIED
    C.tau2 <- sum(1/vi) - cwts2/cwts - (sum(wts2) - cwts2)/(sum(wts) - cwts) - sum(C.left)
    tau2 <- pmax(0, (Qrl + q.left - df)/C.tau2) 
    # NOT POSSIBLE TO USE CUMSUM, TRY USING MATRIX
    w.star <- 1/ (vi.sort[-n] %*% t(rep(1,n-1)) + rep(1, n-1) %*% t(tau2))
    swts <- colSums(w.star) + 1/(vi.sort[n] + tau2) # half time cost compared to sapply
    w.temp <- w.star
    w.temp[lower.tri(w.temp)] <- 0
    cwts0 <- colSums(w.temp) 
    wy.star <- (g.sort[-n] %*% t(rep(1,n-1))) * w.star
    swy <- colSums(wy.star) + g.sort[n]/(vi.sort[n] + tau2) 
    wy.temp <- wy.star
    wy.temp[lower.tri(wy.temp)] <- 0
    cwy0 <- colSums(wy.temp)
    w.left <- 1/ (vi.left %*% t(rep(1, length(tau2))) + rep(1, sum(inx.left)) %*% t(tau2))
    cnode.left <- cnode.left
    # SLOW w.test <- aggregate(x=data.frame(w.left), by = list(cnode.left), sum)
    w.left.bynode <- split(data.frame(w.left), f = cnode.left)
    sw.left <- sapply(w.left.bynode, colSums, USE.NAMES = FALSE)
    wy.left <- g.left %*% t(rep(1,n-1)) * w.left
    wy.left.bynode <- split(data.frame(wy.left), f = cnode.left)
    swy.left <- sapply(wy.left.bynode, colSums, USE.NAMES = FALSE)
    sswy <- swy + rowSums(swy.left)
    ssw <- swts + rowSums(sw.left)
    qb.star <-  (swy - cwy0)^2/(swts - cwts0) + rowSums(swy.left^2/sw.left) +
      cwy0^2/cwts0 - sswy^2/ssw 
    inx.star <- which(qb.star ==  max(qb.star[c.split]))
    cstar <- x.sort[inx.star]
    res <- c(cstar, qb.star[inx.star], tau2[inx.star])
    res
  }
  
} 
#re.cutoff(mf$efk, mf$`(vi)`, mf$x1, cnode = res0$node.split[,5], pl = "9")

GS.fit0 <- function(g, vi, x, cnode, pl){
  mod.order <- unique(x)
  cpoints <- sort(mod.order)
  pleaf.inx <- cnode == pl
  cnode.test <- cnode
  Dev <- -Inf
  for (g in 1:length(cpoints)) {
    cnode.test[pleaf.inx] <- ifelse( x[pleaf.inx] <= cpoints[g], 2*i, 2*i+1)
    temp <- rebetQ(g, vi, mods = as.factor(cnode.test))
    Dev.new <- temp[1]
    if (Dev.new > Dev) {
      Dev <- temp[1]
      tcpt <- cpoints[g]
      tau2 <- temp[2]
    }
  }
  c(tcpt, tau2, Dev)
  
  
}
#GS.fit0(mf$efk, mf$`(vi)`, mf$x1, cnode = res0$node.split[,5], pl = "9")
  