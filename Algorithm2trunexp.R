#Fix the variables
xmetrop = list()
xdistr = c()
xmetrop2 = c()
n = 2
tau = 1
t = 1

#The halving algorithm
half = function(u){
  e = 10^(-9)
  kk = 1
  A = -10
  B = 5
  C = (A+B)/2
  while (f(A,u)*f(B,u) > 0){
    A = A-k
    B = B+k
    kk = kk+1
  }
  while (abs(f(C,u)) > e){
    C = (A+B)/2
    if (f(A,u)*f(C,u) > 0){
      A = C
    }
    else {
      B = C
    }
  }
  return (B)
}

#Define a function for estimating the parameter theta
f = function(theta, u){
  value = sum(log(1+(exp(theta)-1)*u))/theta-t
  return(value)
}

#Prior density for theta
epsilon = 10^(-6)
prior = function(theta){
  if (abs(theta) > epsilon){
    return(1/(abs(theta)))
  }
  else{
    return(0)
  }
}

#Weights W_t(u)
w = function(u, thetahat){
    pri = sapply(thetahat,prior)
    deriv = (exp(thetahat)/thetahat)*sum(u/(1+(exp(thetahat)-1)*u))-t/thetahat
    value = abs(pri/deriv)
  return(value)
}

#Run algorithm 2 for 5000 times
for (j in 1:5000){

#Metropolis algorithm for drawing samples from W_t(u)f(u)
xmetrop[1] = list(runif(n))
m = 10
lambdahat = c()
lambdahat[1] = half(unlist(xmetrop[1]))
for (i in 1:m){
  unif2 = runif(1)
  unif1 = runif(n)
  lambdahat[i+1] = half(unif1)
  alpha = w(unif1, lambdahat[i+1])/w(unlist(xmetrop[i]),lambdahat[i])
  if (unif2 <= alpha){
    xmetrop[i+1] = list(unif1)
  }
  else{
    lambdahat[i+1] = lambdahat[i]
    xmetrop[i+1] = xmetrop[i]
  }
}

#Define the function chi
chi = function(u, theta){
  value = log(1+(exp(theta)-1)*u)/theta
  return(value)
}
#First m elements are for the burnin and m+1th element is chosen as the sample V. Sample is calculated by chi(V,thetahat(V, t))
xmetrop2 = append(xmetrop2,chi(unlist(xmetrop[m+1]), lambdahat[m+1]))
#Every iteration, adds n sample points to vector xmetrop2
}





#Generate initial sample for the Gibbs algorithm
x = rep.int(t/n, n)

#Gibbs algorithm for sampling from the conditional distribution, where T(x)=t.
for (k in 1:10000) {
  ind = sample(1:n,2)
  a = x[ind[1]] + x[ind[2]]
  if (a <= 1){
    x[ind[1]] = runif(1, min = 0, max = a)
  }
  else {
    x[ind[1]] = runif(1, min = a-1, max = 1)
  }
  x[ind[2]] = a - x[ind[1]]
  xdistr = append(xdistr,x)
}

#Plot a histogram of the true conditional sample and the sample given by algorithm 2
par(mfrow=c(1,2))
hist(xmetrop2, breaks = 50, xlim=c(0,1), main = "Algorithm 2")
hist(xdistr, breaks = 50, xlim=c(0,1), main = "True conditional sample")





