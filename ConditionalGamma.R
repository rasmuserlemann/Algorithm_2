
#Variables 
xmetrop = list()
xmetrop2 = c()
iter = 3

#The halving algorithm
half = function(u){
  e = 10^(-6)
  kk = 1
  A = 0.1
  B = 2
  while (alphahat(A,u) > 0){
    A = A/2
  }
  while(alphahat(B,u) < 0){
    B = B*2
  }
  C = (A+B)/2
  while (abs(alphahat(C,u)) > e){
    C = (A+B)/2
    if (alphahat(A,u)*alphahat(C,u) > 0){
      A = C
    }
    else {
      B = C
    }
  }
  return (C)
}

t1 = 10
t2 = -5
n = 10
MLEa = 1
MLEb = 1


#Equation to find alphahat
alphahat = function(alpha,u){
  beta = exp((sum(log(u))-t2/alpha)/n)
  value = (sum((u/beta)^alpha)-t1)
  return(value)
}

#Betahat depends on alphahat 
betahat = function(alphahat, u){
  value = exp((sum(log(u))-t2/alphahat)/n)
  return(value)
}

prior = function(alpha, beta){
  return(1)
  if (alpha > 0.5 && alpha < 1.5 && beta > 0.5 && beta < 1.5){
    value = 1
  }
  else {
    value = 0
  }
  return(value)
}


h = function(u, alpha, beta){
  pri = prior(alpha,beta)
  deriv = (1/beta)*(t1*t2-n*sum(((u/beta)^alpha)*log((u/beta)^alpha)))
  fcond = ((alpha/beta)^n)*exp((1-1/alpha)*t2)
  value = abs(fcond*pri/deriv)
  return(value)
}

chi = function(u, alpha, beta){
  return((u/beta)^alpha)
}

for (vvv in 1:5000){
xmetrop[1] = list(rgamma(n, shape = MLEa, scale = MLEb))
alphah = c()
betah = c()
alphah[1] = half(unlist(xmetrop[1]))
betah[1] = betahat(alphah[1], unlist(xmetrop[1]))
while (h(unlist(xmetrop[1]), alphah[1], betah[1])==0){
  xmetrop[1] = list(rgamma(n, shape = MLEa, scale = MLEb))
  alphah[1] = half(unlist(xmetrop[1]))
  betah[1] = betahat(alphah[1], unlist(xmetrop[1]))
}
for (i in 1:iter){
  unif2 = runif(1)
  unif1 = rgamma(n, shape=MLEa, scale=MLEb)
  alphah[i+1] = half(unif1)
  betah[i+1] = betahat(alphah[i+1], unif1)
  up = h(unif1, alphah[i+1], betah[i+1])*prod(dgamma(unlist(xmetrop[i]), shape=MLEa, scale=MLEb))
  down = h(unlist(xmetrop[i]), alphah[i], betah[i])*prod(dgamma(unif1, shape=MLEa, scale=MLEb))
  alpha = up/down
  if (unif2 <= alpha){
    xmetrop[i+1] = list(unif1)
  }
  else{
    alphah[i+1] = alphah[i]
    betah[i+1] = betah[i]
    xmetrop[i+1] = xmetrop[i]
  }
  if (i == iter){
    xmetrop2 = append(xmetrop2,chi(unlist(xmetrop[i]), alphah[i], betah[i])[1])
  }
}
}

naive = c()
iter = 100000
tt1a = c()
tt2a = c()
eps1 = 10**(-1)
eps2 = 10**(-1)
for (k in 1:iter){
  d = rgamma(n, shape=MLEa, scale=MLEb)
  tt1 = sum(d)
  tt1a = c(tt1a, tt1)
  tt2a = c(tt2a, tt2)
  tt2 = sum(log(d))
  if (abs(tt1-t1)<eps1 && abs(tt2-t2)<eps2){
    naive = c(naive, d[1])
  }
}

print(mean(tt1a))
print(mean(tt2a))

plot(ecdf(xmetrop2), title='', xlab = '', ylab='', main='', xlim=c(0,5.5), ylim = c(0,1))
par(new=TRUE)
plot(ecdf(naive), title='', xlab = '', ylab='', main='', col='red', xlim=c(0,5.5), ylim = c(0,1))
