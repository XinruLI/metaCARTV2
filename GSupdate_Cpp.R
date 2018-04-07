sourceCpp("Qleft.cpp")
sourceCpp("Cleft.cpp")
sourceCpp("REQleft.cpp")
sourceCpp("reQtau2RL.cpp")
re.cutoff.cpp <- function(g, vi, x, inx.s, cnode, minbucket) {
  n <- sum(inx.s)
  ordx <- order(x)
  x.sort <- x[ordx]
  c.split <- (x.sort[-1] - x.sort[-n]) != 0
  if (minbucket > 1) {
    c.split[c(1:(minbucket-1), (n-minbucket+1):(n-1))] <- FALSE
  }
  if (all(c.split == FALSE)) {
    return(NULL)
  } else {
    g.sort <- g[inx.s][ordx]
    vi.sort <- vi[inx.s][ordx]
    inx.left <- inx.s == FALSE
    g.left <- g[inx.left]
    vi.left <- vi[inx.left]
    cnode.left <- cnode[inx.left]
    q.left <- compute_Q_left(g.left, vi.left, cnode.left, unique(cnode.left)) 
    C.left <- compute_C_left(vi.left, cnode.left, unique(cnode.left)) 
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
    C.tau2 <- sum(1/vi) - cwts2/cwts - (sum(wts2) - cwts2)/(sum(wts) - cwts) - sum(C.left)
    tau2 <- pmax(0, (Qrl + sum(q.left) - df)/C.tau2) 
    # 
    temp.unsplit <- compute_swy_tau2(g.left,vi.left, 
                                   cnode.left,tau2,
                                   unique(cnode.left))
    swy.unsplit <- temp.unsplit[,1]
    sw.unsplit <- temp.unsplit[,2]
    wy2byw.unsplit <- temp.unsplit[,3]
    
    temp.rl <- compute_rl_tau2(g.sort, vi.sort, tau2)
    swy.rl <- temp.rl[, 1]
    sw.rl <- temp.rl[, 2]
    cwy.rl <- temp.rl[, 3]
    cw.rl <- temp.rl[, 4]
    
    sswy <- swy.rl + swy.unsplit
    ssw <- sw.rl + sw.unsplit
    
    qb.star <-  (swy.rl - cwy.rl)^2/(sw.rl - cw.rl) + wy2byw.unsplit +
      cwy.rl^2/cw.rl - sswy^2/ssw 
    inx.star <- which(qb.star ==  max(qb.star[c.split]))
    cstar <- x.sort[inx.star]
    res <- c(cstar, qb.star[inx.star], tau2[inx.star])
    res
  }
  
} 

    
    
   