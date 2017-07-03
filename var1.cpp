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
  PARAMETER_VECTOR(logit_eigval); vector<Type> eigval = 2/(1+exp(-logit_eigval))-1; // eigenvalues
  
  matrix<Type> eigvec(2,2);
  eigvec(0,0) = cos(theta(0)); eigvec(0,1) = cos(theta(1));
  eigvec(1,0) = sin(theta(0)); eigvec(1,1) = sin(theta(1));
  
  matrix<Type> D(2,2);
  D(0,0) = eigval(0); D(0,1) = 0;
  D(1,0) = 0;         D(1,1) = eigval(1);
  
  matrix<Type> Phi = eigvec * D * eigvec.inverse();
  
  Type nll = 0;

  // Density of VAR(1) process random effects x
  vector<Type> x0 = x.row(0);
  nll -= sum(dnorm(x0, Type(0), Type(1), true));
  for (int t=1; t<x.rows(); t++) {
    vector<Type> xt = x.row(t);
    vector<Type> xtminus1 = x.row(t-1);
    vector<Type> Phix = Phi*xtminus1;
    vector<Type> error = xt - Phix;
    nll -= VECSCALE(UNSTRUCTURED_CORR(rho),sigma)(error);
  }
  // Conditional density of observations y given x
  for (int t=0; t<y.rows(); t++) {
    for (int i=0; i<2; i++) {
      nll -= dnorm(y(t,i), x(t+x.rows()-y.rows(),i), sdy, true);
    }
  }
  
  return nll;
}
