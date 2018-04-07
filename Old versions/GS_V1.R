# PROBLEMATIC WITH BINARY VARIABLES CANNOT BE USED!!!
# PROBLEMATIC WITH BINARY VARIABLES CANNOT BE USED!!!
# PROBLEMATIC WITH BINARY VARIABLES CANNOT BE USED!!!
sourceCpp("Compute_tau_.cpp")
sourceCpp("compute_re_Q_.cpp")
re.cutoff_cpp <- function(g, vi, x, inx.s, cnode, minbucket) {
  n <- sum(inx.s)
  xinx.order <- order(x)
  x.sort = x[xinx.order]
  SPs <- x.sort[n] - x.sort[-1]
  if (minbucket > 1) {
    SPs[c(1:(minbucket-1), (n-minbucket+1):(n-1))] <- 0
  }
  if (all(SPs == 0) ) {
    return(NULL)
  } else {

    g.sort <- g[inx.s][xinx.order]
    vi.sort <- vi[inx.s][xinx.order]
    inx.left <- inx.s == FALSE
    g.left <- g[inx.left]
    vi.left <- vi[inx.left]
    cnode.left <- cnode[inx.left]
    tau2 <- compute_tau_(g.left, vi.left, cnode.left, unique(cnode.left),
                         g.sort, vi.sort) 
    qb.star <- compute_re_Q_(g.left, vi.left, cnode.left,tau2,
                            unique(cnode.left), g.sort, vi.sort)
    inx.star <- which(qb.star ==  max(qb.star[SPs != 0]))
    cstar <- x.sort[inx.star]
    res <- c(cstar, qb.star[inx.star], tau2[inx.star])
    res
  }
  
} 




