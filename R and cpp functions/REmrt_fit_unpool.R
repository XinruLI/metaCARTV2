# A FUNCTION TO SEARCH FOR THE BEST SPLIT POINT
# BASED ON MAXIMIZING THE UNPOOLED BETWEEN-SUBGROUPS Q.
# ESTIMATES ARE COMPUTED LOCALLY RATHER THAN GLOBALLY.

# Note: R is nearly as fast as using Rcpp.

re.find.cutoff.unpool <- function(g, vi, x, minbucket) {
  n <- length(x)
  xinx.order <- order(x)
  x.sort = x[xinx.order]
  # CHECK IF THERE IS POSSIBLE SPLIT POINTS
  c.split <- (x.sort[-1] - x.sort[-n]) != 0
  # CONSTRICT THE MINIMAL NUMBER OF STUDIES IN CHILD NODES
  if (minbucket > 1) {
    c.split[c(1:(minbucket-1), (n-minbucket+1):(n-1))] <- FALSE
  }
  if (all(c.split == FALSE)) {
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
    qb.star <- (swy - cwy0)^2/(swts - cwts0) +
      cwy0^2/cwts0 - swy^2/swts 
    inx.star <- which(qb.star ==  max(qb.star[c.split]))
    cstar <- x.sort[inx.star]
    res <- c(cstar, qb.star[inx.star], tau2[inx.star])
    res
  }
} 

REmrt_GS_unpool <- function(mf, maxL, minbucket, minsplit, lookahead){
  y <- model.response(mf)
  vi <- c(t(mf["(vi)"]))
  mods.names <-  labels(terms(mf))
  mods <- mf[mods.names]
  nmod <- ncol(mods)
  cpt <- list()
  nodemark <- data.frame(node = rep(1, nrow(mf)))
  res.Qb = 0
  res.tau2 = (sum(y^2/vi) - (sum(y/vi))^2/sum(1/vi)   - length(y)+1)/(sum(1/vi)-sum(1/vi^2)/sum(1/vi)) #VERIFIED
  res.split = NA
  res.mod = NA
  res.pleaf = NA
  
  for (i in 1) {
    Dev<- -Inf
    TQb <- Ttau2 <- Tsplit <- Tmod <- Tpleaf <- NULL
    cnode <- nodemark[ ,i]
    len.node <- tapply(vi, cnode, length)
    nodes <- names(len.node) [len.node >= minsplit]
    for (pl in nodes) {
      pleaf.inx <- cnode == pl
      for (k in 1:nmod) {
        xk <- mods[pleaf.inx, k]
        c.splits <- unique(xk)
        if (length(c.splits) < 2) next
        if (is.numeric(xk)) {
          # NUMERIC VARIABLE
          temp <- re.find.cutoff.unpool(y[pleaf.inx], vi[pleaf.inx], xk, minbucket)
          if (is.null(temp)) {
            Dev.new <- -Inf
          } else {
            Dev.new <- temp[2]
          }
          if (Dev.new > Dev) {
            Dev <- temp[2]
            c.star <- temp[1]
            msplit <- paste(mods.names[k], "<=", c.star, collapse = " ")
            TQb = temp[2]
            Ttau2 = temp[3]
            Tsplit = msplit
            Tmod = mods.names[k]
            Tpleaf = as.numeric(pl)
            new.node <- cnode
            new.node[pleaf.inx] <- ifelse( xk <= c.star, 2*i, 2*i+1)
          }
        } else {
          xk.rank <- rank(tapply(y[pleaf.inx], xk, mean))
          xk.ordinal <- xk.rank[as.character(xk)]
          temp <- re.find.cutoff.unpool(y[pleaf.inx], vi[pleaf.inx], xk.ordinal, minbucket)
          if (is.null(temp)) {
            Dev.new <- -Inf
          } else {
            Dev.new <- temp[2]
          }
          if (Dev.new > Dev) {
            Dev <- temp[2]
            c.star <- names(xk.rank[xk.rank <= temp[1]])
            msplit <- paste(mods.names[k], "=", paste(c.star, collapse = "/"), collapse = " ")
            TQb = temp[2]
            Ttau2 = temp[3]
            Tsplit = msplit
            Tmod = mods.names[k]
            Tpleaf = as.numeric(pl)
            new.node <- cnode
            new.node[pleaf.inx] <- ifelse( xk.ordinal <= temp[1], 2*i, 2*i+1)
          }
        }
      }
    }
    
  }
  while(i <= maxL) {
    Dev <- -Inf
    nodemark <- cbind(nodemark, new.node)
    res.Qb <- c(res.Qb, TQb)
    res.tau2 <- c(res.tau2, Ttau2)
    res.split <- c(res.split, Tsplit)
    res.mod <- c(res.mod,Tmod)
    res.pleaf <- c(res.pleaf, Tpleaf)
    cpt[[i]] <- c.star
    i <- i+1
    TQb <- Ttau2 <- Tsplit <- Tmod <- Tpleaf <- NULL
    cnode <- nodemark[ ,i]
    len.node <- tapply(vi, cnode, length)
    nodes <- names(len.node) [len.node >= minsplit]
    for (pl in nodes) {
      pleaf.inx <- cnode == pl
      for (k in 1:nmod) {
        xk <- mods[pleaf.inx, k]
        c.splits <- unique(xk)
        if (length(c.splits) < 2) next
        if (is.numeric(xk)) {
          # NUMERIC VARIABLE
          temp <- re.find.cutoff.unpool(y[pleaf.inx], vi[pleaf.inx], xk, minbucket)
          if (length(temp) == 0) {
            Dev.new <- -Inf
          } else {
            Dev.new <- temp[2]
          }
          if (Dev.new > Dev) {
            Dev <- temp[2]
            c.star <- temp[1]
            msplit <- paste(mods.names[k], "<=", c.star, collapse = " ")
            TQb = temp[2]
            Ttau2 = temp[3]
            Tsplit = msplit
            Tmod = mods.names[k]
            Tpleaf = as.numeric(pl)
            new.node <- cnode
            new.node[pleaf.inx] <- ifelse( xk <= c.star, 2*i, 2*i+1)
          }
        } else {
          xk.rank <- rank(tapply(y[pleaf.inx], xk, mean))
          xk.ordinal <- xk.rank[as.character(xk)]
          temp <- re.find.cutoff.unpool(y[pleaf.inx], vi[pleaf.inx], xk.ordinal, minbucket)
          if (length(temp) == 0) {
            Dev.new <- -Inf
          } else {
            Dev.new <- temp[2]
          }
          if (Dev.new > Dev) {
            Dev <- temp[2]
            c.star <- names(xk.rank[xk.rank <= temp[1]])
            msplit <- paste(mods.names[k], "=", paste(c.star, collapse = "/"), collapse = " ")
            TQb = temp[2]
            Ttau2 = temp[3]
            Tsplit = msplit
            Tmod = mods.names[k]
            Tpleaf = as.numeric(pl)
            new.node <- cnode
            new.node[pleaf.inx] <- ifelse( xk.ordinal <= temp[1], 2*i, 2*i+1)
          }
        }
      }
    }
    
  }
  list(tree = data.frame(Qb = res.Qb, tau2 = res.tau2, split = res.split,
                         mod = res.mod, pleaf = res.pleaf),
       node.split = nodemark, cpt = cpt, data = mf)
  
}
# load("mf")
# res_unpool <- REmrt_GS_unpool(mf, maxL = 10, 1, 2, F)
# res_unpool$tree
# table(res_unpool$node.split[,4])
