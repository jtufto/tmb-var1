# VAR(1) in TMB
# assuming that Phi has real eigenvalues 
Phi <- matrix(c(.9,.1,-.2,.3),2,2)
eig <- eigen(Phi)
theta <- atan2(eig$vectors[2,],eig$vectors[1,])
lambda <- eig$values
vec <- matrix(c(cos(theta),sin(theta)),2,2,byrow=TRUE)
Phi <- vec %*% diag(lambda) %*% solve(vec)
Phi

rho <- -.5
sigma <- c(2,3)
Rho <- diag(2)
Rho[2,1] <- Rho[1,2] <- rho
Sigma <- diag(sigma) %*% Rho %*% diag(sigma)
Sigma

Gamma0 <- solve(diag(4)-kronecker(Phi,Phi), as.vector(Sigma))
Gamma0 <- matrix(Gamma0,2,2)
Gamma0

# Simulate some data
n <- 100
x <- matrix(NA,n,2)
x[1,] <- mvtnorm::rmvnorm(1,sigma=Gamma0)
for (t in 2:n) 
  x[t,] <- Phi %*% x[t-1,] + t(mvtnorm::rmvnorm(1,sigma=Sigma))

# Check stationary covariance against theoretical value
var(x)
Gamma0

sdy <- .1
y <- x + rnorm(2*n, sd=sdy)
#matplot(1:n, y, type="l",col=1)

library(TMB)
compile("var1.cpp")
dyn.load(dynlib("var1"))
data <- list(y=y)
parameters <- list(
  x = matrix(0,nrow(y), 2),
  logsdy = log(sdy),
  log_sigma = log(sigma),
  rho = rho,
  theta = c(0, pi/2), 
  logit_eigval = qlogis((lambda+1)/2)
)
map <- list(rho = factor(c(NA)))
obj <- MakeADFun(data,parameters,random="x")
obj$method="BFGS"
opt <- do.call(optim,obj)
rep <- sdreport(obj)
summary(rep,"report")
Phi
Sigma
