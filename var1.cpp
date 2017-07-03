#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;

  DATA_MATRIX(y);
  
  PARAMETER_MATRIX(x);
  PARAMETER(logsdy); Type sdy = exp(logsdy);
  
  // Sigma parameterized in terms of log of the two sd's and logit2 of rho
  PARAMETER_VECTOR(log_sigma); vector<Type> sigma = exp(log_sigma);
  PARAMETER_VECTOR(rho);
  
  // Phi parameterized in terms
  PARAMETER_VECTOR(theta);
  PARAMETER_VECTOR(logit_eigval); vector<Type> eigval = 2/(1+exp(-logit_eigval)-1); // eigenvalues
  
  matrix<Type> eigvec(2,2);
  eigvec(0,0) = cos(theta(0)); eigvec(1,0) = cos(theta(1));
  eigvec(0,1) = sin(theta(0)); eigvec(1,1) = sin(theta(1));
  
  matrix<Type> D(2,2);
  D(0,0) = eigval(0); D(0,1) = 0;
  D(1,0) = 0;         D(1,1) = eigval(1);
  
  matrix<Type> Phi = eigvec * D * eigvec.inverse();
  
  Type nll = 0;

  // Density of VAR(1) process random effects x
  vector<Type> xt = 0;
  for (int t=1; t<x.rows(); t++) {
    vector<Type> mean = Phi * xt;
    xt(0) = x(t,0);
    xt(1) = x(t,1);
    nll -= VECSCALE(UNSTRUCTURED_CORR(rho),sigma)(xt - mean);
  }
  // Conditional density of observations y given x
  for (int t=0; t<y.rows(); t++) {
    for (int i=0; i<2; i++) {
      nll -= dnorm(y(t,i), x(t+x.rows()-y.rows(),i), sdy, true);
    }
  }
  
  return nll;
}
