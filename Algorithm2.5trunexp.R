x = c()
n = 2
t = 1

for (i in 1:50000){
  u = runif(n)
  if (t*max(u) <= sum(u) & sum(u) <= t){
    x = append(x,(t/sum(u))*u)
  }
}

hist(x, breaks = 50, xlim = c(0,1))