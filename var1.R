# VAR(1) in TMB
# assuming that Phi has real eigenvalues 
theta <- c(0,pi/2)
lambda <- c(.9,.3)
vec <- matrix(c(cos(theta),sin(theta)),2,2,byrow=FALSE)
Phi <- vec %*% diag(lambda) %*% solve(vec)
Phi
eigen(Phi)

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

var(x)
Gamma0

y <- x + rnorm(2*n, sd=2)
matplot(1:n, y, type="l",col=1)

library(TMB)
compile("var1.cpp")
dyn.load(dynlib("var1"))
data <- list(y=y)
burnin <- 100 
parameters <- list(
  x = matrix(0,nrow(y)+burnin,2),
  logsdy = 0,
  log_sigma = c(0,0),
  rho = 0,
  theta = c(0, pi/2),
  logit_eigval = c(0, 0)
)
obj <- MakeADFun(data,parameters)
