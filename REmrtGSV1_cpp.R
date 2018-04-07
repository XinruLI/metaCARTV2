source("GSupdate_Cpp.R")
# TO DO:
# ADD CHECK minsplit SHOULD BE LARGER THAN minbucket
# THE minsplit should be larger than 2?
REmrt_GS_cpp <- function(mf, maxL, minbucket, minsplit, delQ, lookahead){
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
  delta.Q <- Inf
  
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
          temp <- re.cutoff.cpp(y, vi, xk, pleaf.inx, cnode, minbucket)
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
          temp <- re.cutoff.cpp(y, vi, xk.ordinal, pleaf.inx, cnode, minbucket)
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
    if (is.null(TQb)) { 
      delta.Q <- -Inf
    } else {
      delta.Q <- abs(TQb - res.Qb[i])
    }
  }
  while(delta.Q >= delQ & i <= maxL) {
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
          temp <- re.cutoff.cpp(y, vi, xk, pleaf.inx, cnode, minbucket)
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
          temp <-  re.cutoff.cpp(y, vi, xk.ordinal, pleaf.inx, cnode, minbucket)
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
    if (is.null(TQb)) { 
      delta.Q <- -Inf
    } else {
      delta.Q <- abs(TQb - res.Qb[i])
    }
    
  }
  list(tree = data.frame(Qb = res.Qb, tau2 = res.tau2, split = res.split,
                         mod = res.mod, pleaf = res.pleaf),
       node.split = nodemark, cpt = cpt, data = mf)
  
}

# res.up <- REmrt_GS(mf, maxL = 10, minbucket=3, minsplit=2, delQ=0.00001, lookahead=F)
# system.time(res.up <- REmrt_GS(mf, maxL = 10, minbucket=3, minsplit=2, delQ=0.00001, lookahead=F) )
# # 0.672
# source("REmrt0.R")
# system.time(
# res0 <- REmrt.fit0(mf, maxL = 10, minbucket=3, minsplit=2, delQ=0.00001, lookahead=F) 
# )
# 0.917 
