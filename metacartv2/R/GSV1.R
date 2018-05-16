#' A function to find the split point
#'
#' @param g the effect size
#' @param vi the sampling variance
#' @param x the splitting moderator 
#' @param minbucket the minimum number of the studies in a terminal node
#' @param inx.s indicates whether a study belongs to the candidate parent leaf
#' @param cnode the terminal nodes that the studies belong to in the current tree
#' @return a vector including the split point, Q, and tau2
#' @keywords internal
re.cutoff_cpp <- function(g, vi, x, inx.s, cnode, minbucket) {
  n <- sum(inx.s)
  xinx.order <- order(x)
  x.sort = x[xinx.order]
  c.split <- (x.sort[-1] - x.sort[-n]) != 0
  if (minbucket > 1) {
    c.split[c(1:(minbucket-1), (n-minbucket+1):(n-1))] <- FALSE
  }
  if (all(c.split == FALSE)) {
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
    inx.star <- which(qb.star ==  max(qb.star[c.split]))
    cstar <- x.sort[inx.star]
    res <- c(cstar, qb.star[inx.star], tau2[inx.star])
    res
  }
  
} 




