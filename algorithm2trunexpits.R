#Fix the variables
xmetrop2 = c()
n = 2
tau = 1
t = 1

#Define the function chi
chi = function(u, theta){
  value = log(1+(exp(theta)-1)*u)/theta
  return(value)
}

#The halving algorithm
half = function(u){
  e = 10^(-6)
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
epsilon = 10^(-3)
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

#Approximate the density by calculating w_t(u) values on a grid
step = 500
wvalues = matrix(nrow = step-1, ncol = step-1)
lambdahat = matrix(nrow = step-1, ncol = step-1)
for (v1 in 1:(step-1)){
  t1 = v1/(step)
  for (v2 in 1:(step-1)){
    t2 = v2/(step)
    t3 = c(t1,t2)
    lambdahatt = half(t3)
    wvalues[v1,v2] = w(t3,lambdahatt)
  }
}

#Integrate one of the variables out and normalize the marginal "density"
mar1pdf = c()
for (int1 in 1:(step-1)){
  mar1pdf[int1] = sum(wvalues[int1,])
}
pdf1 = mar1pdf/sum(mar1pdf)

#Calculate the marginal CDF
cdf1 = c()
for (int2 in 1:(step-1)){
  cdf1[int2] = sum(pdf1[1:int2])
}

for (b in 1:30000){
#Draw a realization from the marginal distribution
u1 = runif(1)
index1 = which.min(abs(cdf1-u1))
v1 = index1/step

#Draw a realization from the approximation of the conditional distribution
u2 = runif(1)
mar2pdf = wvalues[index1,]/sum(wvalues[index1,])
mar2cdf = c()
for (k2 in 1:step-1){
  mar2cdf[k2] = sum(mar2pdf[1:k2])
}
index2 = which.min(abs(mar2cdf-u2))
v2 = index2/step

vv = c(v1,v2)
lambdahat = half(vv)
xmetrop2 = append(xmetrop2,chi(vv,lambdahat))
}

#Plot the histogram
hist(xmetrop2, breaks = 50, xlim=c(0,1), main = "")

