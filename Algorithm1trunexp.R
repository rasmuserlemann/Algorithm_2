library(truncdist)

n = 2
lambda = 2
tau = 1

#Fix the sufficient statistic T(x)=2
t = 1

#Generate the initial sample for the Gibbs algorithm
x = rep.int(t/n, n)

#Gibbs algorithm for sampling from the conditional distribution, where T(x)=t
xdistr = c()
for (k in 1:25000) {
  ind = sample(1:n,2)
  a = x[ind[1]] + x[ind[2]]
  if (a <= 1){
    x[ind[1]] = runif(1, min = 0, max = a)
  }
  else {
    x[ind[1]] = runif(1, min = a-1, max = 1)
  }
  x[ind[2]] = a - x[ind[1]]
  if (k%%25 == 0) {xdistr = append(xdistr,x)}
}

#Generate 3 data points from the conditional distribution using algorithm 1. Run this 1000 times and
#estimate the distribution with those 3000 data points.
xt = c()
for (k in 1:1000){
#Draw n samples from Uni(0,1)
u = runif(n)

#Estimate the parameter lambda
f = function(lambdah){
  if (lambdah == 0){
    value = sum(tau*u)-t
  }
  else{
    value = sum(-log(1-(1-exp(-lambdah*tau))*u))-lambdah*t
  }
  return(value)
}
lambdahat = uniroot(f,c(0,10), extendInt = 'yes')$root

#Calculate a sample from the conditional distribution
xt = append(xt,-log(1-(1-exp(-lambdahat*tau))*u)/lambdahat)
}

par(mfrow=c(1,2))
plot(density(xdistr))
plot(density(xt))

ks.test(xt,xdistr)$p.value

#High p-value suggests that the initial sample x and 
#the conditional sample xt are drawn from the same distribution



