load("mf2")
g <- mf$efk
vi <- mf$`(vi)`
cnode <- sample(1:4, length(g), replace = T)
inx.s <- cnode == 1
x <- mf$x1[inx.s]
n <- sum(inx.s)
x.sort <- sort(x)
c.split <- (x.sort[-1] - x.sort[-n]) != 0
g.sort <- g[inx.s][order(x)]
vi.sort <- vi[inx.s][order(x)]
inx.left <- inx.s == FALSE
g.left <- g[inx.left]
vi.left <- vi[inx.left]
cnode.left <- cnode[inx.left]
q.left <- sum(tapply(g.left^2/vi.left, cnode.left, sum) -
                tapply(g.left/vi.left, cnode.left, function(x) (sum(x))^2)/
                tapply(1/vi.left, cnode.left, sum))
C.left <- tapply(1/vi.left^2, cnode.left, sum)/tapply(1/vi.left, cnode.left, sum)


library(Rcpp)
sourceCpp("Qleft.cpp")
compute_Q_left(g.left, vi.left, cnode.left, unique(cnode.left)) 
sourceCpp("Cleft.cpp")
compute_C_left(vi.left, cnode.left, unique(cnode.left)) 
# 0.015 comparing to tapply 0.241 for 1000 iters


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
tau2 <- pmax(0, (Qrl + q.left - df)/C.tau2) 
w.star <- 1/ (vi.sort[-n] %*% t(rep(1,n-1)) + rep(1, n-1) %*% t(tau2))
swts <- colSums(w.star) + 1/(vi.sort[n] + tau2) # half time cost compared to sapply
w.temp <- w.star
w.temp[lower.tri(w.temp)] <- 0
cwts0 <- colSums(w.temp) 
wy.star <- (g.sort[-n] %*% t(rep(1,n-1))) * w.star
swy <- colSums(wy.star) + g.sort[n]/(vi.sort[n] + tau2) 
wy.temp <- wy.star
wy.temp[lower.tri(wy.temp)] <- 0
cwy0 <- colSums(wy.temp)
w.left <- 1/ (vi.left %*% t(rep(1, length(tau2))) + rep(1, sum(inx.left)) %*% t(tau2))
cnode.left <- cnode.left
# SLOW w.test <- aggregate(x=data.frame(w.left), by = list(cnode.left), sum)
w.left.bynode <- split(data.frame(w.left), f = cnode.left)
sw.left <- sapply(w.left.bynode, colSums, USE.NAMES = FALSE)
wy.left <- g.left %*% t(rep(1,n-1)) * w.left
wy.left.bynode <- split(data.frame(wy.left), f = cnode.left)
swy.left <- sapply(wy.left.bynode, colSums, USE.NAMES = FALSE)
wy2.w <- rowSums(swy.left^2/sw.left)

sourceCpp("REQleft.cpp")
test <- compute_swy_tau2(g.left,vi.left, 
                cnode.left,tau2,
                 unique(cnode.left))


for (i in 1:(n-1)){
  if (abs(rowSums(swy.left)[i]-test[i,1]) > 1e-12) print(i)
}
for (i in 1:(n-1)){
  if (abs(rowSums(sw.left)[i]-test[i,2]) > 1e-12) print(i)
}
for (i in 1:(n-1)){
  if (abs(wy2.w[i]-test[i,3]) > 1e-12) print(i)
}

sourceCpp("reQtau2RL.cpp")
test2 <- compute_rl_tau2(g.sort, vi.sort, tau2)
for (i in 1:(n-1)){
  if (abs(swts[i]-test2[i,2]) > 1e-12) print(i)
  if (abs(swy[i]-test2[i,1]) > 1e-12) print(i)
  }
for (i in 1:(n-1)){
  if(abs(cwy0[i] - test2[i,3]) > 1e-12) print(i)
  if(abs(cwts0[i] - test2[i,4]) > 1e-12) print(i)
}
