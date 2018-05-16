library(Rcpp)
source("REmrtGSV1_cpp.R")
source("Xvalid_all_cpp.R")

prune.REmrt <- function(x, size){
  x$tree$CP <- -c(NA, diff(x$tree$tau2))
  tree <- x$tree
  nodes <- x$node.split
  temp <- nodes[,ncol(nodes)]
  id <- unique(temp)
  # subT.rm <- NULL
  new <- tree
  old.nodes <- 1:(2*nrow(tree) - 1)
  new.nodes <- old.nodes
  while(length(id) > size){

    id2 <- id[id%%2 == 0L]
    id2.inx <- (id2 + 1) %in% id
    id.rm <- id2[id2.inx]
    pl <- tree$pleaf[id.rm/2 +1L]
    if (length(id.rm) > 1) {
      remove <- id.rm[-which.max(tree$CP[match(pl, tree$pleaf, 0L)])]
    } else {remove <- id.rm}

    for (nn in remove){
      temp[temp == nn] <- tree$pleaf[nn/2 +1L]
      temp[temp == (nn+1)] <- tree$pleaf[nn/2 +1L]
  #     subT.rm <- rbind(subT.rm,
  #                   tree[match(tree$pleaf[nn/2 +1L], tree$pleaf, 0L), ])
      new <- new[-match(tree$pleaf[nn/2 +1L], new$pleaf, 0L), ]
  #     
      new.nodes[match(c(nn, nn+1), old.nodes, 0L)] <- NA
      new.nodes[old.nodes > nn] <- new.nodes[old.nodes > nn] - 2
    }
    id <- unique(temp)
  }
  names(new.nodes) <- old.nodes
  new$pleaf <- new.nodes[as.character(new$pleaf)]
  new$tau2 <- new$Qb <- new$CP <- NULL
  temp <- new.nodes[as.character(temp)]
  return(list(new = new, lable.changed = new.nodes, Tnodes = temp))
  
}
grow.REmrt <- function(mf, x, size, minbucket=5, minsplit=6) {
  new <- x[[1]]
  Tnodes <- x[[3]]
  y <- model.response(mf)
  vi <- c(t(mf["(vi)"]))
  mods.names <-  labels(terms(mf))
  mods <- mf[mods.names]
  nmod <- ncol(mods)
  add <- NULL
  for (i in nrow(new):(size-1)){
    Dev <- -Inf
    cnode <- Tnodes
    len.node <- tapply(vi, cnode, length)
    nodes <- names(len.node) [len.node >= minsplit]
    for (pl in unique(nodes)) {
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
    add <- rbind(add, data.frame(split = Tsplit,
                                 mod = Tmod, pleaf = Tpleaf))
    Tnodes <- new.node
  }
  return(list(rbind(new, add), Tnodes, Dev))

  
  
  
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

 load("mf")
 res <- REmrt_GS_cpp(mf, maxL = 15, minbucket=5, minsplit=6, delQ=1, lookahead=F)
x <- prune.REmrt(res, 8)
x2 <- grow.REmrt(mf, x, 8,  minbucket=5, minsplit=6)
rebetQ(mf$efk, mf$`(vi)`, x2[[2]])

 
 
# res$tree$CP <- -c(NA, diff(res$tree$tau2))
# tree <- res$tree
# nodes <- res$node.split
# size.pr = 6
# temp <- nodes[,ncol(nodes)]
# id <- unique(temp)
# subT.rm <- NULL
# new <- tree
# old.nodes <- 1:27
# new.nodes <- old.nodes
# while(length(id) > size.pr){
#   
#   id2 <- id[id%%2 == 0L]
#   id2.inx <- (id2 + 1) %in% id
#   id.rm <- id2[id2.inx]
#   pl <- tree$pleaf[id.rm/2 +1L]
#   if (length(id.rm) > 1) {
#     remove <- id.rm[-which.max(tree$CP[match(pl, tree$pleaf, 0L)])]
#   } else {remove <- id.rm}
#   
#   for (nn in remove){
#     temp[temp == nn] <- tree$pleaf[nn/2 +1L]
#     temp[temp == (nn+1)] <- tree$pleaf[nn/2 +1L]
#     subT.rm <- rbind(subT.rm,
#                   tree[match(tree$pleaf[nn/2 +1L], tree$pleaf, 0L), ])
#     new <- new[-match(tree$pleaf[nn/2 +1L], new$pleaf, 0L), ]
#     
#     new.nodes[match(c(nn, nn+1), old.nodes, 0L)] <- NA
#     new.nodes[old.nodes > nn] <- new.nodes[old.nodes > nn] - 2
#   }
#   id <- unique(temp)
#   
# }
# table(temp)
# subT.rm
# new
# names(new.nodes) = old.nodes
# new.nodes
# rebetQ(mf$efk, mf$`(vi)`, temp)[1]
# rebetQ(mf$efk, mf$`(vi)`, nodes[,6])[1]
# tree[6, ]
# table(nodes[,6])
# 
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