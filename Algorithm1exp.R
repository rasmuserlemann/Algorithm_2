#Draw n samples from exp(lambda)
n = 100
lambda = 2
x = rexp(n,lambda)

#Sufficient statistic
t = sum(x)

#Draw n samples from exp(1)
u = rexp(n,1)

#Estimate the parameter lambda
lambdahat = sum(u)/t

#Conditional sample
xt = u/lambdahat


ks.test(xt,x)
#High p-value suggests that the initial sample x and 
#the conditional sample xt are drawn from the same distribution