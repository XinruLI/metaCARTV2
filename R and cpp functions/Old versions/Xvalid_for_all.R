REmrt.prednode <- function(x, newdata){
  tt <- terms(x$data)
  ms <- model.frame(delete.response(tt), newdata)
  oms <- model.frame(delete.response(tt), x$data)
  tree <- x[["tree"]]
  # if (any(sapply(ms, class) != sapply(oms, class)))
  #   stop("The type of the variables do not match")
  if(nrow(tree) < 2) {
    pred.node <- rep(1, nrow(ms))
  } else {
    tnode <- rep(1, nrow(ms))
    nodes <- tnode
    for (i in 1:(nrow(tree) - 1)){
      tinx <- which(tnode == tree[i+1, "pleaf"])
      tempm <- ms[tree[i+1, "mod"]]
      if(sapply(tempm, is.numeric) == TRUE) {
        tnode[tinx] <- ifelse(tempm[tinx,1] <= x[["cpt"]][[i]], 2*i, 2*i+1)
      } else {
        tnode[tinx] <- ifelse(tempm[tinx,1] %in% oms[,tree[i+1, "mod"]],
                              ifelse(tempm[tinx,1] %in% x[["cpt"]][[i]], 2*i, 2*i+1),
                              NA)
      }
      nodes <- cbind(nodes, tnode)
    }
    pred.node <- nodes
  }
  pred.node
}
REmrt.xvalid <- function(func, mf, maxL, n.fold, minbucket, minsplit, delQ, lookahead){
  N <- nrow(mf)
  if (n.fold > N | n.fold <=1 ) stop("n.fold is not acceptable")
  if (maxL > N) {
    warning("The maximum number of split is too large")
    maxL <- N-1
  }
  pred <- matrix(NA, nrow = N, ncol = maxL+1)
  inx <- sample(1:N)
  inx.xvalid <- c(round(seq(from = 1, by= N/n.fold, length.out = n.fold)), N+1)
  for (i in 1:n.fold){
    inx.test <- inx[inx.xvalid[i]:(inx.xvalid[i+1]-1)]
    test <- mf[inx.test, ]
    train <- mf[-inx.test, ]
    fit.train <- do.call(func, list(mf = train, maxL, minbucket, minsplit, delQ, lookahead))
    yi.train <- model.response(fit.train$data)
    vi.train <- c(t(fit.train$data["(vi)"]))
    nsplt <- nrow(fit.train$tree)
    train.wy <- sapply(1:nsplt, function(x) yi.train/(vi.train+fit.train$tree[,2][x]))
    train.wts <-  sapply(1:nsplt, function(x) 1/(vi.train+fit.train$tree[,2][x]))
    train.pred <- sapply(1:nsplt, function(x) tapply(train.wy[, x], fit.train$node.split[, x], sum)/tapply(train.wts[, x], fit.train$node.split[, x], sum))
    test.nodes <- REmrt.prednode(fit.train, test)
    if (is.null(dim(test.nodes))) {
      pred[inx.test, 1:nsplt] <- train.pred[1]
    } else {
      test.pred <- lapply(1:nsplt, function(x) train.pred[[x]][as.character(test.nodes[ ,x])])
      test.pred.rmna <- sapply(1:nsplt, function(x){test.pred[[x]][is.na(test.pred[[x]])] <- sum(train.wy[, x])/sum(train.wts[, x]);test.pred[[x]]})
      pred[inx.test, 1:nsplt] <- test.pred.rmna
    }
    
  }
  pred <- pred[ ,!is.na(colSums(pred))]
  y <- model.response(mf)
  if (!is.null(dim(pred))) {
    x.error <- apply(pred,2, function(x) sum((y-x)^2)/sum((y-mean(y))^2))
    sdx.error <- apply(pred, 2, function(x)  sqrt(sum(((y-x)^2-mean((y-x)^2))^2))/sum((y-mean(y))^2))
  } else {
    x.error <- sum((y-pred)^2)/sum((y-mean(y))^2)
    sdx.error <- sqrt(sum(((y-pred)^2-mean((y-pred)^2))^2))/sum((y-mean(y))^2)
  }
  cbind(x.error, sdx.error)
  
}

