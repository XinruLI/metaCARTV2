# Solved by restricting maxL being larger than 2
load("mf2")
#load("mf")
source("REmrt_GS_update.R")
res.up <- REmrt_GS(mf, maxL = 5, minbucket=1, minsplit=2, delQ=5, lookahead=F)
