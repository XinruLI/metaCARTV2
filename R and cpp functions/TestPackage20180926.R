library(metacartv2)
load("/Users/xinruli/Documents/Meta-CART\ Projects/Meta-CART\ R-Package/The\ Package\ v2/R\ and\ cpp\ functions/mf")
mf$vi <- mf[,5]
res1 <- REmrt(efk~x1+x2+x3, data = mf, vi = vi, minsplit = 2, minbucket = 1, cp = 0.001, lookahead = T)
plot(res1)
res <- REmrt_GS_(mf,minsplit = 2, maxL = 5, minbucket = 1, cp = 0.001, lookahead = F)$tree
class(res$split)
while(delta.Q >= cp & i < maxL) {
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
        temp <- re.cutoff_cpp(y, vi, xk, pleaf.inx, cnode, minbucket)
        if (is.null(temp)) {
          Dev.new <- -Inf
        } else {
          Dev.new <- temp[2]
        }
        if (Dev.new > Dev) {
          Dev <- temp[2]
          c.star <- temp[1]
          msplit <- paste(mods.names[k], "<", c.star, collapse = " ")
          TQb = temp[2]
          Ttau2 = temp[3]
          Tsplit = msplit
          Tmod = mods.names[k]
          Tpleaf = as.numeric(pl)
          new.node <- cnode
          new.node[pleaf.inx] <- ifelse( xk < c.star, 2*i, 2*i+1)
        }
      } else {
        xk.rank <- rank(tapply(y[pleaf.inx], xk, mean))
        xk.ordinal <- xk.rank[as.character(xk)]
        temp <-  re.cutoff_cpp(y, vi, xk.ordinal, pleaf.inx, cnode, minbucket)
        if (is.null(temp)) {
          Dev.new <- -Inf
        } else {
          Dev.new <- temp[2]
        }
        if (Dev.new > Dev) {
          Dev <- temp[2]
          c.star <- names(xk.rank[xk.rank < temp[1]])
          msplit <- paste(mods.names[k], "=", paste(c.star, collapse = "/"), collapse = " ")
          TQb = temp[2]
          Ttau2 = temp[3]
          Tsplit = msplit
          Tmod = mods.names[k]
          Tpleaf = as.numeric(pl)
          new.node <- cnode
          new.node[pleaf.inx] <- ifelse( xk.ordinal < temp[1], 2*i, 2*i+1)
        }
      }
    }
  }
  nodemark <- cbind(nodemark, new.node)
  res.Qb <- c(res.Qb, TQb)
  res.tau2 <- c(res.tau2, Ttau2)
  res.split <- c(res.split, Tsplit)
  res.mod <- c(res.mod,Tmod)
  res.pleaf <- c(res.pleaf, Tpleaf)
  cpt[[i]] <- c.star
  if (is.null(TQb)) { 
    delta.Q <- -Inf
  } else {
    delta.Q <- abs(TQb - res.Qb[i])
  }
}