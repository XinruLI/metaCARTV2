source("REmrt_GS_update.R")
source("Xvalid_for_all.R")
load("mf")
REmrt.xvalid(REmrt_GS, mf, n.fold = 10, maxL = 10, minbucket=5, minsplit=6, delQ=1, lookahead=F)
maxL = 3
tree <- REmrt_GS(mf, maxL = 3, minbucket=5, minsplit=6, delQ=1, lookahead=F)
Qb <- tree$tree$Qb[4]

set.seed(13)
Qb.p <- numeric(1000)
for (p in 1:1000){
  print(p)
  inx.p <- sample(1:nrow(mf))
  mf.b <- mf
  mf.b$efk <- mf$efk[inx.p]
  mf.b$`(vi)`<- mf$`(vi)`[inx.p]
  tree.p <- REmrt_GS(mf.b, maxL = 3, minbucket=5, minsplit=6, delQ=1, lookahead=F)
  Qb.p[p] <- tree.p$tree$Qb[4]

}
mean(Qb.p > Qb)

Qb <- tree$tree$Qb[2]

set.seed(8)
Qb.p <- numeric(1000)
for (p in 1:1000){
  print(p)
  inx.p <- sample(1:nrow(mf))
  mf.b <- mf
  mf.b$efk <- mf$efk[inx.p]
  mf.b$`(vi)`<- mf$`(vi)`[inx.p]
  tree.p <- REmrt_GS(mf.b, maxL = 1, minbucket=5, minsplit=6, delQ=1, lookahead=F)
  Qb.p[p] <- tree.p$tree$Qb[2]
  
}
mean(Qb.p > Qb)
