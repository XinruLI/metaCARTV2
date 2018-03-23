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

REmrt.fit0 <- function(mf, maxL, minbucket, minsplit, delQ, lookahead){
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
  
  if (lookahead == TRUE) {
    # first two split to examine the combination
    i <- 1
    Dev<- -Inf
    cnode <- nodemark[ ,i]
    len.node <- tapply(vi, cnode, length)
    nodes <- names(len.node) [len.node >= minsplit]
    tres <- NULL
    if (length(nodes) == 0) break  # if none of the nodes has the length bigger than minsplit, then stop
    for (j in 1:length(nodes)){  # search through nodes
      pleaf.inx <- cnode == as.numeric(nodes[j])
      tres1 <- tres2 <- NULL
      for (k1 in 1:(nmod-1)){
        chosenmod1 <- mods[pleaf.inx, k1]
        if(sapply(mods[mods.names[k1]], is.numeric) == FALSE) {  # if the moderator is non-numeric
          mod.order1 <- rank(tapply(y[pleaf.inx],chosenmod1,mean))  # sort the nominal variable by group means
          cmod.ordinal1 <- mod.order1[as.character(chosenmod1)]
        } else {  # for nummeric variables
          mod.order1 <- unique(chosenmod1)
          cmod.ordinal1 <- chosenmod1
        }
        cpoints1 <- sort(mod.order1)
        if (length(cpoints1) >= 2) {
          for (g1 in 1:(length(cpoints1)-1)) {  # le goes to left, > goes to right
            cnode.test1 <- cnode
            cnode.test1[pleaf.inx] <- ifelse( cmod.ordinal1 <= cpoints1[g1], 2*i, 2*i+1)
            nodes.test2 <- names(table(cnode.test1)>= minsplit)
            if (length(nodes.test2) == 0) break
            for (j2 in 1:length(nodes.test2)){  # search through all possible children leaves
              pleaf.inx2 <- cnode.test1 == as.numeric(nodes.test2[j2])
              for (k2 in 1:nmod){  # search through the combinations of moderators
                if (k2 ==  k1) next  # jump over the same moderator
                chosenmod2 <- mods[pleaf.inx2, k2]
                if(sapply(mods[mods.names[k2]], is.numeric) == FALSE) {
                  mod.order2 <- rank(tapply(y[pleaf.inx2],chosenmod2,mean))
                  cmod.ordinal2 <- mod.order2[as.character(chosenmod2)]
                } else {
                  mod.order2 <- unique(chosenmod2)
                  cmod.ordinal2 <- chosenmod2}
                cpoints2 <- sort(mod.order2)
                
                if (length(cpoints2) >= 2) {
                  for (g2 in 1:(length(cpoints2)-1)) {
                    cnode.test2 <- cnode.test1
                    cnode.test2[pleaf.inx2] <- ifelse( cmod.ordinal2 <= cpoints2[g2], 2*2, 2*2+1)
                    if (min(table(cnode.test2)) < minbucket) {
                      Dev.new <- -Inf
                    } else {
                      temp <- rebetQ(y, vi, mods = as.factor(cnode.test2))
                      Dev.new <- temp[1]
                    }
                    if (Dev.new > Dev) {
                      Dev <- Dev.new
                      temp0 <- rebetQ(y, vi, mods = as.factor(cnode.test1))
                      Dev1 <- temp0[1]
                      tau2.1 <- temp0[2]
                      
                      if(sapply(mods[mods.names[k1]], is.numeric) == FALSE) {
                        tcpt1 <- names(mod.order1[mod.order1 <= cpoints1[g1]])
                        msplit1 <- paste(mods.names[k1], "=", paste(names(mod.order1[mod.order1 <= cpoints1[g1]]), collapse = "/"), collapse = " ")
                      } else {
                        tcpt1 <- cpoints1[g1]
                        msplit1 <- paste(mods.names[k1], "<=", tcpt1, collapse = " ")
                      }
                      
                      tres1 <- data.frame(Qb = Dev1, tau2 =tau2.1,
                                          split = msplit1, mod = mods.names[k1], pleaf = as.numeric(nodes[j]))
                      newnode1 <- cnode.test1
                      Dev2 <- Dev.new
                      tau2.2 <- temp[2]
                      if(sapply(mods[mods.names[k2]], is.numeric) == FALSE) {
                        tcpt2 <- names(mod.order2[mod.order2 <= cpoints2[g2]])
                        msplit2 <- paste(mods.names[k2], "=", paste(names(mod.order2[mod.order2 <= cpoints2[g2]]), collapse = "/"), collapse = " ")
                      } else {
                        tcpt2 <- cpoints2[g2]
                        msplit2 <- paste(mods.names[k2], "<=", tcpt2, collapse = " ")
                      }
                      
                      tres2 <- data.frame(Qb = Dev2, tau2 =tau2.2,
                                          split = msplit2, mod = mods.names[k2], pleaf = as.numeric(nodes.test2[j2]))
                      newnode2 <- cnode.test2
                    }
                  }
                }
              }
            }
            
          }
        } # the first split: cnode.test1
        
        
        
        
        
      }
    }
    
    if (!is.null(tres2)) {
      delta.Q <- abs(tres2$Qb)
      if (delta.Q < delQ) break
      nodemark <- cbind(nodemark, newnode1, newnode2)
      res <- rbind(res, tres1, tres2)
      cpt[[1]] <- tcpt1
      cpt[[2]] <- tcpt2
    }
    
    # continue here for the 3rd split
    
    if (maxL > 2 & (nrow(res)==3) ) {
      for (i in 3:maxL) {
        Dev<- -Inf
        tres <- NULL
        cnode <- nodemark[ ,i]
        len.node <- tapply(vi, cnode, length)
        nodes <- names(len.node) [len.node >= minsplit]
        if (length(nodes) == 0) break
        for (j in 1:length(nodes)) {
          pleaf.inx <- cnode == as.numeric(nodes[j])
          for (k in 1:nmod){
            chosenmod <- mods[pleaf.inx, k]
            if(sapply(mods[mods.names[k]], is.numeric) == FALSE) {
              mod.order <- rank(tapply(y[pleaf.inx],chosenmod,mean))
              cmod.ordinal <- mod.order[as.character(chosenmod)]
            } else {
              mod.order <- unique(chosenmod)
              cmod.ordinal <- chosenmod}
            cpoints <- sort(mod.order)
            if (length(cpoints) >= 2) {
              for (g in 1:(length(cpoints)-1)) {
                cnode.test <- cnode
                cnode.test[pleaf.inx] <- ifelse( cmod.ordinal <= cpoints[g], 2*i, 2*i+1)
                if (min(table(cnode.test[pleaf.inx])) < minbucket) {
                  Dev.new <- -Inf
                } else {
                  temp <- rebetQ(y, vi, mods = as.factor(cnode.test))
                  Dev.new <- temp[1]
                }
                if (Dev.new > Dev) {
                  Dev <- temp[1]
                  if(sapply(mods[mods.names[k]], is.numeric) == FALSE) {
                    tcpt <- names(mod.order[mod.order <= cpoints[g]])
                    msplit <- paste(mods.names[k], "=", paste(names(mod.order[mod.order <= cpoints[g]]), collapse = "/"), collapse = " ")
                  } else {
                    tcpt <- cpoints[g]
                    msplit <- paste(mods.names[k], "<=", tcpt, collapse = " ")
                  }
                  tres <- data.frame(Qb = temp[1], tau2 = temp[2],
                                     split = msplit, mod = mods.names[k], pleaf = as.numeric(nodes[j]))
                  tnode <- cnode.test
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
    }
    
    
  } else {
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
          chosenmod <- mods[pleaf.inx, k]
          if(sapply(mods[mods.names[k]], is.numeric) == FALSE) {
            mod.order <- rank(tapply(y[pleaf.inx],chosenmod,mean))
            cmod.ordinal <- mod.order[as.character(chosenmod)]
          } else {
            mod.order <- unique(chosenmod)
            cmod.ordinal <- chosenmod}
          cpoints <- sort(mod.order)
          if (length(cpoints) >= 2) {
            for (g in 1:(length(cpoints)-1)) {
              cnode.test <- cnode
              cnode.test[pleaf.inx] <- ifelse( cmod.ordinal <= cpoints[g], 2*i, 2*i+1)
              if (min(table(cnode.test[pleaf.inx])) < minbucket) {
                Dev.new <- -Inf
              } else {
                temp <- rebetQ(y, vi, mods = as.factor(cnode.test))
                Dev.new <- temp[1]
              }
              if (Dev.new > Dev) {
                Dev <- temp[1]
                if(sapply(mods[mods.names[k]], is.numeric) == FALSE) {
                  tcpt <- names(mod.order[mod.order <= cpoints[g]])
                  msplit <- paste(mods.names[k], "=", paste(names(mod.order[mod.order <= cpoints[g]]), collapse = "/"), collapse = " ")
                } else {
                  tcpt <- cpoints[g]
                  msplit <- paste(mods.names[k], "<=", tcpt, collapse = " ")
                }
                tres <- data.frame(Qb = temp[1], tau2 = temp[2],
                                   split = msplit, mod = mods.names[k], pleaf = as.numeric(nodes[j]))
                tnode <- cnode.test
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
  }
  
  
  
  list(tree = res, node.split = nodemark, cpt = cpt, data = mf)
  
  
}
