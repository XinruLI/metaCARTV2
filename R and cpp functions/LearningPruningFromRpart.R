tree <- rpart(efk~ x1+x2+x3, data = mf, control =  rpart.control(maxdepth = 3L, cp = 0)
              )
tree$frame
prp(tree)


rt <- tree$frame$dev/120
names(rt) <- row.names(tree$frame)

a1 <- (rt[1] - sum(rt[tree$frame$var == "<leaf>"]))/7

a2 <- (rt[2] - sum(rt[as.character(c(8: 11))]))/3
a3 <- (rt[as.character(3)] - sum(rt[as.character(c(12: 15))]))/3
a4 <- (rt[as.character(4)] - sum(rt[as.character(c(8: 9))]))/1
a5 <- (rt[as.character(5)] - sum(rt[as.character(c(10: 11))]))/1
a6 <- (rt[as.character(6)] - sum(rt[as.character(c(12: 13))]))/1
a7 <- (rt[as.character(7)] - sum(rt[as.character(c(14: 15))]))/1
rank(c(a1,a2,a3,a4,a5,a6,a7)) 
# remove 89

a1 <- ((rt[1] - sum(rt[tree$frame$var == "<leaf>"])) + a4)/6

a2 <- (rt[2] - sum(rt[as.character(c(8: 11))]))/2 + a4/2
a3 <- (rt[as.character(3)] - sum(rt[as.character(c(12: 15))]))/3 
a5 <- (rt[as.character(5)] - sum(rt[as.character(c(10: 11))]))/1 
a6 <- (rt[as.character(6)] - sum(rt[as.character(c(12: 13))]))/1 
a7 <- (rt[as.character(7)] - sum(rt[as.character(c(14: 15))]))/1 
rank(c(a1,a2,a3,a5,a6,a7))

# remove 10,11

a1 <- ((rt[1] - sum(rt[tree$frame$var == "<leaf>"])) + a4 + a5)/5
a2 <- (rt[2] - sum(rt[as.character(c(8: 11))])) + a4 + a5
a3 <- (rt[as.character(3)] - sum(rt[as.character(c(12: 15))]))/3 
a6 <- (rt[as.character(6)] - sum(rt[as.character(c(12: 13))]))/1 
a7 <- (rt[as.character(7)] - sum(rt[as.character(c(14: 15))]))/1 
rank(c(a1,a2,a3,a6,a7))

# remove 12, 13

a1 <- ((rt[1] - sum(rt[tree$frame$var == "<leaf>"])) + a4 + a5 + a6)/4
a2 <- (rt[2] - sum(rt[as.character(c(8: 11))])) + a4 + a5
a3 <- (rt[as.character(3)] - sum(rt[as.character(c(12: 15))]))/2 + a6/2 
a7 <- (rt[as.character(7)] - sum(rt[as.character(c(14: 15))]))/1 
rank(c(a1,a2,a3,a7))

# remove 4,5

a1 <- ((rt[1] - sum(rt[tree$frame$var == "<leaf>"])) + a4 + a5 + a6 + a2)/3
a3 <- (rt[as.character(3)] - sum(rt[as.character(c(12: 15))]))/2 + a6/2 
a7 <- (rt[as.character(7)] - sum(rt[as.character(c(14: 15))]))/1 
rank(c(a1,a3,a7))





diff(tree$cptable[,3])
diff(tree$cptable[,3])[1]/3









mf$splt <- tree$where
summary(lm(efk~factor(splt), data = mf))

library(rpart.plot)
set.seed(1621)
x1 <- sample(0:1, 200, replace = T)
x2 <- sample(0:1, 200, replace = T)
x3 <- sample(0:1, 200, replace = T)
y <- ifelse(x1 == 0, 1, ifelse(x2 == 0, 0, ifelse(x3 == 0, 1, 0)))
tree <- rpart(y~ x1+x2+x3, method = "class",
              control =  rpart.control(maxdepth = 30L, cp = 0)
)
tree$frame
prp(tree)
table(y)


a1 <- 69/200/3
a2 <- 26/200/2
a3 <- 19/200/1

