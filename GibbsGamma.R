library(truncdist)

xgibbs = c()

#Sample from the gamma distribution and statistics
x = c(1,5,2,2,4)
s = sum(x)
p = prod(x)

#Variables
n = length(x)
xx = x[3:n-1]
ss = sum(xx)
pp = prod(xx)
C = s-ss
D = p/pp

sqrtpos = function(v, c){
  return(v*(1-v)^2-c)
}

#Kernel function for the conditonal density of x_n
kernel = function(v, C, D){
  return(1/(sqrt(v)*sqrt(v*(1-v)^2-4*D/(C^3))))
}

for (u in 1:1000){

epsilon = 0.0001
  
for (o in 1:n-2){
  #Step 1: Generate x* = xn
  #Metropolis algorithm to draw a sample from the kernel
  xm = c()
  m = 10
  a = uniroot(sqrtpos, c(0,1/3), extendInt = 'no', c=4*D/(C^3))$root + epsilon
  b = uniroot(sqrtpos, c(1/3,1), extendInt = 'no', c=4*D/(C^3))$root - epsilon
  xm[1] = rtrunc(1, spec = 'norm', a = a, b = b)
  for (i in 1:m){
    xm[i+1] = rtrunc(1, spec = 'norm', a = a, b = b)
    alpha = kernel(xm[i+1], C, D)/kernel(xm[i], C, D)
    u = runif(1)
    if (alpha <= u){
      xm[i+1] = xm[i]
    }
  }
  
  xn = C*xm[m+1]
  
  #Step 2: Relabel
  xtemp = c()
  xtemp[n] = xn
  
  #Step 3: Rotate the sample
  for (j in 3:n-2){
    xtemp[j] = x[j+1] 
  }
  xtemp[n-1] = xtemp[n]
  xtemp[n] = x[3]
  x[3:n] = xtemp[3:n]
  #Step 4: Recalculate variables
  xx = x[seq(3,n-1,1)]
  ss = sum(xx)
  pp = prod(xx)
  C = s - ss
  D = p/pp
  
  #Step 5
}

#Step 6: 
bb = s - sum(x[3:n])
cc = p/prod(x[3:n])
x[1] = (bb-sqrt(bb^2-4*cc))/2
x[2] = (bb+sqrt(bb^2-4*cc))/2

xgibbs = append(xgibbs, x)

}

hist(xgibbs, breaks = 100, xlab = NULL, main = NULL)
