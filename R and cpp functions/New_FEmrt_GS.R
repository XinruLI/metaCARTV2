# a FE to split the nodes based on depth
# To do list: to add minsplit. and do something else other than storing the indication for every child leaves
load("/Users/xinruli/Documents/Meta-CART\ Projects/Meta-CART\ R-Package/The\ Package\ v2/R\ and\ cpp\ functions/mf")
Find_Split_Point_FE <- function(yi, vi, xk, minbucket){
  # given the moderator, the effecti size and the sampling variance
  # find the best split point
  # To do: add minbucekt and minsplit option
  n <- length(yi)
  xk.order <- order(xk)
  xk <- xk[xk.order]
  yi <- yi[xk.order]
  vi <- vi[xk.order]
  c.split <- (xk[-1] - xk[-n]) != 0
  if (minbucket > 1) {
    c.split[c(1:(minbucket-1), (n-minbucket+1):(n-1))] <- FALSE
  }
  if (all(c.split == FALSE)) {
    return(NULL)
  } else {
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
    res <- c((xk[inx.star]+xk[inx.star+1])/2, Ql[inx.star], Qr[inx.star])
    names(res) <- c("c.star", "Ql", "Qr")
    return(res)
  }
  
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
  inx.star <- which(Qrl ==  min(Qrl))
  res <- c((xk[inx.star]+xk[inx.star+1])/2, Ql[inx.star], Qr[inx.star])
  names(res) <- c("c.star", "Ql", "Qr")
  return(res)
  
}
FE.split <- function(y, vi, mods, mod.types, node, minbucekt){
  # Given the effect size, the sampling variance, the data frame of 
  # all moderators, whether the moderators are numeric, and a list 
  # containing information about the parent leaves
  # node is a list with each elements containing information of one
  # parent leaf. So each element should be a list contains three values
  # splt.candidates: the name of the candidate moderators
  # Q: the within subgroup Q
  # inx.pleaf: a logical vector indicate if a study belongs to this leaf
  Q <- node$Q
  splt.var <- node$splt.candidates
  inx.pleaf <- node$inx.leaf
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
      temp <- Find_Split_Point_FE(y, vi, xk.ord, minbucekt)
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
        inx.left <- inx.right <- inx.pleaf
        inx.left[inx.pleaf] <-  xk.ord < temp[1]
        inx.right[inx.pleaf] <- xk.ord >= temp[1]
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
    
  } else (res <- node)
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
             inx.leaf = rep(TRUE, length(y))))

maxD = 3

for (j in 1:maxD) {
  for (i in 1:length(new.parents)) {
  children <- FE.split(y, vi, mods, mod.types, new.parents[[i]], minbucekt = 3)
  tree.frame <- rbind(tree.frame, children[["split"]])
  cpt <- do.call(c, list(cpt, children$cpt))
  offsprings <- do.call(c, list(offsprings, children$nodes))
  }
  new.parents <- offsprings
  offsprings <- NULL
}


for (i in 1:length(new.parents)) {
  children <- FE.split(y, vi, mods, mod.types, new.parents[[i]])
  tree.frame <- rbind(tree.frame, children[["split"]])
  cpt <- do.call(c, list(cpt, children$cpt))
  offsprings <- do.call(c, list(offsprings, children$nodes))
  
}
new.parents <- offsprings
offsprings <- NULL


# Check and compare with the rpart version
library(metacart)
tree <- FEmrt(efk~x1+x2+x3, vi=vi, data = mf)




depth = 0L
Tree.temp <- list( split = NULL, nodes = list (root))

depth = 1L
cpt <- list()
Q <- Tree.temp$nodes[[1]]$Q
splt.var <- Tree.temp$nodes[[1]]$splt.candidates
inx.pleaf <- Tree.temp$nodes[[1]]$inx.leaf
Dev <- -Inf
for (k in splt.var){
  xk <- mods[pleaf.inx, k]
  c.splits <- unique(xk)
  if(length(c.splits) == 1) {
    mods.names <- setdiff(mods.names, k)
    next
  } else {
    if (mod.types[k] == TRUE) {
      xk.ord <- xk
    } else {
      xk.rank <- rank(tapply(y[pleaf.inx], xk, mean))
      xk.ord <- xk.rank[as.character(xk)]
    }
      temp <- Find_Split_Point_FE(y, vi, xk.ord, minbucket)
      if (is.null(temp)) {
        Dev.new <- -Inf
      } else {
        Dev.new <- Q - sum(temp[2:3])
      }
      if (Dev.new > Dev) {
        Dev <- Dev.new
        if (mod.types[k] == TRUE) {
          c.star <- temp[1]
          msplit <- paste(k, "<=", c.star, collapse = " ")
        } else {
          c.star <- names(xk.rank[xk.rank <= temp[1]])
          msplit <- paste(k, "=", paste(c.star, collapse = "/"), collapse = " ")
         
        }
        Q.left <- temp[2]
        Q.right <- temp[3]
        mod = k
        inx.left <- inx.right <- inx.pleaf
        inx.left[inx.pleaf] <-  xk.ord <= temp[1]
        inx.right[inx.pleaf] <- xk.ord > temp[1]
      }
  }
}
temp.split <- data.frame(Dev, mod, msplit)
cpt <- do.call(c, list(cpt, list(c.star)))

child.left <- list(splt.candidates = mods.names, 
                           Q = Q.left,
                           inx.leaf = inx.left)
child.right <- list(splt.candidates = mods.names, 
                   Q = Q.right,
                   inx.leaf = inx.right)
offspring <- list(child.left, child.right)
names(offspring) <- 2:3

66.06950+316.46880

all.equal(mf$x1 < 0.1139092, mf$x1 <= cpt[[1]])
which()