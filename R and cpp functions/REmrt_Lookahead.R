set.seed(1456)

source("/Users/xinruli/Documents/Visit\ to\ El\ Paso/Permutation/Permutation_funcs.R")
K = 160
n = 160
x1 <- sample(0:1, K, replace = TRUE)
x2 <- sample(1:5, K, replace = TRUE) # nominal
x3 <- sample(1:6, K, replace = TRUE)  # binary
x4 <- rep("A", K)
mods <- data.frame(x1 = x1, x2 = x2, x3 = x3, x4=x4)
mods.form <- data.frame(x = (x2 %in% c(2,4)) * (x3>3))
dat.sim1 <- SimData(c(0, 0.8), mods.form, "~ x", K, NumGen(n,K), 0.1)
dat <- data.frame(dat.sim1, mods)
dat$x2 <- letters[x2]
mf <- metacartv2::REmrt(efk ~ x1+ x2 + x3 + x4, vi = vark, data = dat)$data
mf
rm(list = setdiff(ls(), list("mf")))
minbucket = 3
lookahead = T
minsplit = 6


REmrt_GS_cpp2 <- function(mf, maxL, minbucket, minsplit, delQ, lookahead){
  #===================  Error message  ======================#
  if (minbucket >= minsplit) {
    stop("minbucket should be smaller than minsplit")
  }
  if (!is.logical(lookahead)) {
    stop("lookahead should be TRUE or FALSE")
  }
  find_triplet <- function(xk, leaf.labels, y, vi, minbucket, minsplit) {
    # given a moderator xk, and a vector of leaf membership
    # return the combination of parent leaf and the splitting point 
    leaves <- as.numeric(names(table(leaf.labels) >= minsplit))
    if (is.numeric(xk)) {
      xk.rank <- NULL
      tempQ <-  sapply(leaves, function(x)
        metacartv2:::re.cutoff_cpp(y, vi, xk[leaf.labels == x], leaf.labels == x, leaf.labels, minbucket)
      )
    } else {
      xk.rank <- lapply(leaves,function(x) rank(tapply(y[leaf.labels == x], xk[leaf.labels == x], mean)))
      names(xk.rank) <- leaves
      tempQ <- sapply(leaves, function(x) 
        metacartv2:::re.cutoff_cpp(y, vi, xk.rank[[as.character(x)]][as.character(xk[leaf.labels == (x)])], 
                                   leaf.labels == x, leaf.labels, minbucket) )
    }
    if(class(tempQ) == "list") {
      if (all(sapply(tempQ, is.null))) {
        list(pleaf = NA, cstar = c(NA, -Inf, NA), rank = NULL)
      } else { #### PROBLOMATIC
        pleaf <- which.max(sapply(tempQ, function(x) if (is.null(x)) -Inf else x[2]))
        list(pleaf = leaves[pleaf], cstar = tempQ[[pleaf]], rank = xk.rank[[pleaf]])
      }
      
    } else {
      pleaf <- which.max(tempQ[2,])
      list(pleaf = leaves[pleaf], 
           cstar = tempQ[ ,pleaf],
           rank = xk.rank[[pleaf]])
    }
  }
  
  y <- model.response(mf)
  vi <- c(t(mf["(vi)"]))
  mods.names <-  labels(terms(mf))
  mods <- mf[mods.names]
  nmod <- ncol(mods)
  cpt <- list()
  nodemark <- data.frame(rep(1, nrow(mf)))
  res.Qb = 0
  res.tau2 = (sum(y^2/vi) - (sum(y/vi))^2/sum(1/vi)   - length(y)+1)/(sum(1/vi)-sum(1/vi^2)/sum(1/vi)) #VERIFIED
  res.split = NA
  res.mod = NA
  res.pleaf = NA
  delta.Q <- Inf
  if (lookahead & maxL > 1) {
    make_first_split <- function(xk, y) {
      # A function to list all possible split points for the first split
      # Return : a list with returned child leaves (child.leaf), 
      #         levels of the moderator (rank), 
      #         and if the moderator is numeric (is.num)
      if (is.numeric(xk)) {
        split.points <- sort(unique(xk))
        res <- list(child.leaf = sapply(split.points[-1], function(x) ifelse(xk < x, 2, 3)),
                    rank1 = split.points, is.num = TRUE)
      } else {
        split.points <- sort(rank(tapply(y, xk, mean)))
        res <- list(child.leaf = sapply(split.points[-1], function(x) ifelse(split.points[as.character(xk)] < x, 2, 3)),
                    rank1 = split.points, is.num = FALSE)
      }
    }
    find_second_split <- function(xk, first.splits, y, vi, minbucket) {
      # Given a moderator and all possible first splits
      # Return the optimal combination of first two splits, and the corresponding Q-between
      # split1: the first element is the index of the first moderator
      #         the second element is the split point of the first moderator
      # split2: the first element is the selected parent leaf for the 2nd split
      #         the second element contains the split point of the 2nd split, Q, the tau-square
      Q.temp <- -Inf
      Csplit.temp <- list(Q = Q.temp)
      
      for (i in 1:length(first.splits)) {
        temp.first.split <- first.splits[[i]]$child.leaf
        if (length(temp.first.split) == 0) next
        res <- lapply(1:ncol(temp.first.split), function(x) find_triplet(xk, temp.first.split[ ,x], y, vi, minbucket))
        tempQs <- sapply(1:ncol(temp.first.split), function(x) {if (is.null(res[[x]])) {-Inf} else {res[[x]]$cstar[2]}})
        Q.max <- max(tempQs)
        Csplit1.max <- c(i, (first.splits[[i]]$rank1[which.max(tempQs)] + first.splits[[i]]$rank1[which.max(tempQs)+1])/2, which.max(tempQs))
        Csplit2.max <- res[[which.max(tempQs)]]
        if (Q.max > Q.temp) {
          Q.temp <- Q.max
          Csplit.temp <- list(Q = Q.temp,
                              split1 = Csplit1.max, 
                              split2 = Csplit2.max)
        }
        
      }
      Csplit.temp
      
    }
    compute_rebetQ <- function(yi, vi, mods){
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
    first.splits <- lapply(mods, function(x) make_first_split(x, y))
    second.splits <- lapply(1:nmod, function(x) find_second_split(mods[ ,x], first.splits, y, vi, minbucket))
    tempQs <- sapply(1:length(second.splits), function(x) second.splits[[x]]$Q)
    x2.inx <- which.max(tempQs)
    x1.inx <- second.splits[[x2.inx]]$split1[1]
    res.mod <- c(NA, mods.names[c(x1.inx, x2.inx)])
    res.pleaf <- c(NA, 1, second.splits[[x2.inx]]$split2$pleaf)
    if (first.splits[[x1.inx]]$is.num) {
      cstar1 <- second.splits[[x2.inx]]$split1[2]
      msplit1 <- paste(mods.names[x1.inx], "<", cstar1, collapse = " ")
    } else {
      cstar1 <- names(first.splits[[x1.inx]]$rank[first.splits[[x1.inx]]$rank < second.splits[[x2.inx]]$split1[2]])
      msplit1 <- paste(mods.names[x1.inx], "=", paste(cstar1, collapse = "/"), collapse = " ")
    }
    node1 <- first.splits[[x1.inx]]$child.leaf[ ,second.splits[[x2.inx]]$split1[3]]
    names(node1) <- NULL
    pleaf.inx <- node1 == second.splits[[x2.inx]]$split2$pleaf
    node2 <- node1
    if (is.null(second.splits[[x2.inx]]$split2$rank)) {
      cstar2 <- second.splits[[x2.inx]]$split2$cstar[1]
      node2[pleaf.inx] <- ifelse(mods[pleaf.inx, x2.inx] < cstar2, 4, 5)
      msplit2 <- paste(mods.names[x2.inx], "<", cstar2 , collapse = " ")
    } else {
      cstar2 <- names(second.splits[[x2.inx]]$split2$rank[second.splits[[x2.inx]]$split2$rank < second.splits[[x2.inx]]$split2$cstar[1]])
      node2[pleaf.inx] <- ifelse(mods[pleaf.inx, x2.inx] %in% cstar2, 4, 5)
      msplit2 <- paste(mods.names[x2.inx], "=", paste(cstar2, collapse = "/"), collapse = " ")
    }
    nodemark <- cbind(nodemark, node1, node2)
    cpt[[1]] <- cstar1
    cpt[[2]] <- cstar2
    Q.split1 <- compute_rebetQ(y, vi, node1)
    res.Qb <- c(0, Q.split1[1], second.splits[[x2.inx]]$Q)
    res.tau2 <- c(res.tau2, Q.split1[2], second.splits[[x2.inx]]$split2$cstar[3])
    res.split <- c(NA, msplit1, msplit2)
    Dev <- res.Qb[3]
    delta.Q <- res.Qb[3] - res.Qb[2]
    i = 2
    while(delta.Q >= delQ & i < maxL) {
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
            temp <-  re.cutoff_cpp(y, vi, xk.ordinal, pleaf.inx, cnode, minbucket)
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
  } else{
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
            temp <- re.cutoff_cpp(y, vi, xk, pleaf.inx, cnode, minbucket)
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
            temp <- re.cutoff_cpp(y, vi, xk.ordinal, pleaf.inx, cnode, minbucket)
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
            temp <- re.cutoff_cpp(y, vi, xk, pleaf.inx, cnode, minbucket)
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
            temp <-  re.cutoff_cpp(y, vi, xk.ordinal, pleaf.inx, cnode, minbucket)
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
  }
  list(tree = data.frame(Qb = res.Qb, tau2 = res.tau2, split = as.character(res.split),
                         mod = res.mod, pleaf = res.pleaf),
       node.split = nodemark, cpt = cpt, data = mf)
  
}





y <- model.response(mf)
vi <- c(t(mf["(vi)"]))
mods.names <-  labels(terms(mf))
mods <- mf[mods.names]
nmods <- ncol(mods)
make_first_split <- function(xk, y) {
  if (is.numeric(xk)) {
    split.points <- sort(unique(xk))
    res <- list(child.leaf = sapply(split.points[-1], function(x) ifelse(xk < x, 2, 3)),
                rank1 = split.points, is.num = TRUE)
  } else {
    split.points <- sort(rank(tapply(y, xk, mean)))
    res <- list(child.leaf = sapply(split.points[-1], function(x) ifelse(split.points[as.character(xk)] < x, 2, 3)),
                rank1 = split.points, is.num = FALSE)
  }
}




first.splits <- lapply(mods, function(x) make_first_split(x, y))
# second splits#
temp_func1 <- function(xk, first.splits, y, vi, minbucket) {
  # Given a moderator and all possible first splits
  # Return the optimal combination of first two splits, and the corresponding Q-between
  # split1: the first element is the index of the first moderator
  #         the second element is the split point of the first moderator
  # split2: the first element is the selected parent leaf for the 2nd split
  #         the second element contains the split point of the 2nd split, Q, the tau-square
  find_triplet <- function(xk, leaf.labels, y, vi, minbucket) {
    # given a moderator xk, and a vector of leaf membership
    # return the combination of parent leaf and the splitting point 
    if (is.numeric(xk)) {
      xk.rank <- NULL
      tempQ <-  sapply(2:3, function(x)
        metacartv2:::re.cutoff_cpp(y, vi, xk[leaf.labels == x], leaf.labels == x, leaf.labels, minbucket)
      )
    } else {
      xk.rank <- lapply(2:3,function(x) rank(tapply(y[leaf.labels == x], xk[leaf.labels == x], mean)))
      
      tempQ <- sapply(2:3, function(x) 
        metacartv2:::re.cutoff_cpp(y, vi, xk.rank[[x-1]][as.character(xk[leaf.labels == (x)])], 
                                   leaf.labels == x, leaf.labels, minbucket) )
    }
    if(class(tempQ) == "list") {
      if (all(sapply(tempQ, is.null))) {
        list(pleaf = NA, cstar = c(NA, -Inf, NA), rank = NULL)
      } else {
        pleaf <- which(sapply(tempQ, is.null) == FALSE)
        list(pleaf = pleaf + 1, cstar = tempQ[[pleaf]], rank = xk.rank[[pleaf]])
      }
      
    } else {
      pleaf <- which.max(tempQ[2,])
      list(pleaf = pleaf + 1, 
           cstar = tempQ[ ,pleaf],
           rank = xk.rank[[pleaf]])
    }
  }
  Q.temp <- -Inf
  Csplit.temp <- list(Q = Q.temp)
 for (i in 1:length(first.splits)) {
   temp.first.split <- first.splits[[i]]$child.leaf
   res <- lapply(1:ncol(temp.first.split), function(x) find_triplet(xk, temp.first.split[ ,x], y, vi, minbucket))
   tempQs <- sapply(1:ncol(temp.first.split), function(x) {if (is.null(res[[x]])) {-Inf} else {res[[x]]$cstar[2]}})
   Q.max <- max(tempQs)
   Csplit1.max <- c(i, (first.splits[[i]]$rank1[which.max(tempQs)] + first.splits[[i]]$rank1[which.max(tempQs)+1])/2)
   Csplit2.max <- res[[which.max(tempQs)]]
   if (Q.max > Q.temp) {
     Q.temp <- Q.max
     Csplit.temp <- list(Q = Q.temp,
                         split1 = Csplit1.max, 
                         split2 = Csplit2.max)
   }
   
 }
  Csplit.temp
  
}

xk <- mods[ ,3]
leaf.labels <- node2


find_triplet_old <- function(xk, leaf.labels, y, vi, minbucket) {
  # given a moderator xk, and a vector of leaf membership
  # return the combination of parent leaf and the splitting point 
  if (is.numeric(xk)) {
    xk.rank <- NULL
    tempQ <-  sapply(2:3, function(x)
      metacartv2:::re.cutoff_cpp(y, vi, xk[leaf.labels == x], leaf.labels == x, leaf.labels, minbucket)
    )
  } else {
    xk.rank <- lapply(2:3,function(x) rank(tapply(y[leaf.labels == x], xk[leaf.labels == x], mean)))
    
    tempQ <- sapply(2:3, function(x) 
      metacartv2:::re.cutoff_cpp(y, vi, xk.rank[[x-1]][as.character(xk[leaf.labels == (x)])], 
                                 leaf.labels == x, leaf.labels, minbucket) )
  }
  if(class(tempQ) == "list") {
    if (all(sapply(tempQ, is.null))) {
      list(pleaf = NA, cstar = c(NA, -Inf, NA), rank = NULL)
    } else {
      pleaf <- which(sapply(tempQ, is.null) == FALSE)
      list(pleaf = pleaf + 1, cstar = tempQ[[pleaf]], rank = xk.rank[[pleaf]])
    }
    
  } else {
    pleaf <- which.max(tempQ[2,])
    list(pleaf = pleaf + 1, 
         cstar = tempQ[ ,pleaf],
         rank = xk.rank[[pleaf]])
  }
}

res <- REmrt_GS_cpp2(mf,  maxL = 5, minbucket = 5, minsplit = 6, delQ = 0.001, lookahead = T)
res2 <- REmrt_GS_cpp2(mf,  maxL = 5, minbucket = 5, minsplit = 6, delQ = 0.001, lookahead = F)
