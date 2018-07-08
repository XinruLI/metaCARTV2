# 2018-05-28
# done: add minbucket & minsplit
# To do: add yr, yl, ser, sel in the cutoff function
#        add node number, remove msplit (since the data.frame is slow)
#        add track of parents 
load("/Users/xinruli/Documents/Meta-CART\ Projects/Meta-CART\ R-Package/The\ Package\ v2/R\ and\ cpp\ functions/mf")

FE.GS.cutoff <- function(yi, vi, xk, minbucket){
  n <- length(yi)
  xk.order <- order(xk)
  xk <- xk[xk.order]
  c.split <- (xk[-1] - xk[-n]) != 0
  if (minbucket > 1) {
    c.split[c(1:(minbucket-1), (n-minbucket+1):(n-1))] <- FALSE
  }
  
  if (sum(c.split) == 0) {
    return(NULL)
  } else {
    yi <- yi[xk.order]
    vi <- vi[xk.order]
    
    wy <- yi / vi
    wy2 <- yi^2 / vi
    wts <- 1/vi
    wts2 <- wts^2
    cwy <- cumsum(wy[-n])
    cwy2 <- cumsum(wy2[-n])
    cwts <- cumsum(wts[-n])
    cwts2 <- cumsum(wts2[-n])
    Ql <- cwy2 - cwy^2/cwts
    Qr <- sum(wy2) - cwy2 - (sum(wy) - cwy)^2/(sum(wts) - cwts)
    Qrl <- Ql + Qr
    inx.star <- which(Qrl ==  min(Qrl[c.split]))
    res <- c((xk[inx.star]+xk[inx.star+1])/2, Ql[inx.star], Qr[inx.star],
             cwy[inx.star]/cwts[inx.star], (sum(wy)-cwy[inx.star])/(sum(wts)-cwts[inx.star]),
             sqrt(1/cwts[inx.star]), sqrt(1/(sum(wts)-cwts[inx.star])))
    names(res) <- c("c.star", "Ql", "Qr", "d.l", "d.r", "se.l", "se.r")
    return(res)
  }
}
# Check the function
# FE.GS.cutoff(mf$efk, mf$`(vi)`, mf$x3, 3)
# library(metafor)
# rma(yi = efk, vi = `(vi)`, mods =  ~ x3, data = mf, method = "FE")
# rma(yi = efk, vi = `(vi)`, subset =  x3 == 1, data = mf, method = "FE")
# rma(yi = efk, vi = `(vi)`, subset =  x3 == 0, data = mf, method = "FE")
FE.split <- function(y, vi, mods, mod.types, node, minsplit = 2L, minbucket = 1L){
  splt.var <- node$splt.candidates
  inx.pleaf <- node$inx.leaf
  if (length(inx.pleaf) < minsplit | length(splt.var) == 0) {
    res <- list(split = NULL, cpt = NULL, nodes = node)
    } else {
      Q <- node$Q
      Dev <- -Inf
      y <- y[inx.pleaf]
      vi <- vi[inx.pleaf]
      for (k in splt.var){
        xk <- mods[inx.pleaf, k]
        c.splits <- unique(xk)
        if(length(c.splits) == 1) {
          splt.var <- setdiff(splt.var, k)
          next
        } else {
          if (mod.types[k] == TRUE) {
            xk.ord <- xk
          } else {
            xk.rank <- rank(tapply(y, xk, mean))
            xk.ord <- xk.rank[as.character(xk)]
          }
          temp <- FE.GS.cutoff(y, vi, xk.ord, minbucket)
          if (is.null(temp)) {
            Dev.new <- -Inf
          } else {
            Dev.new <- Q - sum(temp[2:3])
          }
          if (Dev.new > Dev) {
            Dev <- Dev.new
            if (mod.types[k] == TRUE) {
              c.star <- temp[1]
              msplit <- paste(k, "<", c.star, collapse = " ")
            } else {
              c.star <- names(xk.rank[xk.rank <= temp[1]])
              msplit <- paste(k, "=", paste(c.star, collapse = "/"), collapse = " ")
              
            }
            Q.left <- temp[2]
            Q.right <- temp[3]
            mod = k
            inx.left <- inx.pleaf[xk.ord < temp[1]]
            inx.right <- inx.pleaf[xk.ord >= temp[1]]
          }
        }
      }
      if(is.finite(Dev)) {
        split <- data.frame(Dev, mod, msplit)
        cpt <- list(c.star)
        child.left <- list(splt.candidates = mods.names, 
                           Q = Q.left,
                           inx.leaf = inx.left)
        child.right <- list(splt.candidates = mods.names, 
                            Q = Q.right,
                            inx.leaf = inx.right)
        offspring <- list(child.left, child.right)
        res <- list(split = split, cpt = cpt, nodes = offspring)
        
      } else (res <- list(split = NULL, cpt = NULL, nodes = node))
    
    }

  res
  
}

y <- model.response(mf)
vi <- c(t(mf["(vi)"]))
mods.names <-  labels(terms(mf))
mods <- mf[mods.names]
cpt <- NULL
tree.frame <- NULL
offsprings <- NULL
mod.types <- sapply(mods, is.numeric) # If a moderator is numeric
new.parents <- list(list(splt.candidates = mods.names, 
                         Q = sum(y^2/vi) - (sum(y/vi))^2/sum(1/vi),
                         inx.leaf = 1:length(y)))
maxD = 3

for (j in 1:maxD) {
  for (i in 1:length(new.parents)) {
    children <- FE.split(y, vi, mods, mod.types, new.parents[[i]])
    tree.frame <- rbind(tree.frame, children[["split"]])
    cpt <- do.call(c, list(cpt, children$cpt))
    offsprings <- do.call(c, list(offsprings, children$nodes))
    print(2^(j-1)+i-1)
  }
  new.parents <- offsprings
  offsprings <- NULL
  names(new.parents) <- (2^j) : (2^(j+1)-1)
}
length(new.parents)
