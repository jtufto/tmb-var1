#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;

  int d = 2; // dimension of vector process

  DATA_MATRIX(y);
  
  PARAMETER_MATRIX(x);
  PARAMETER(logsdy); Type sdy = exp(logsdy); ADREPORT(sdy);
  
  // Sigma parameterized in terms of log of the standard deviations
  PARAMETER_VECTOR(log_sigma); vector<Type> sigma = exp(log_sigma); ADREPORT(sigma);
  PARAMETER_VECTOR(rho);
  matrix<Type> Rho(d,d);
  Rho = UNSTRUCTURED_CORR(rho).cov(); ADREPORT(Rho);
  matrix<Type> Sigma(d,d);
  vector<Type> vecSigma(d*d);
  for (int i=0; i<d; i++)
    for (int j=0; j<d; j++) {
      Sigma(i,j) = sigma(i)*sigma(j)*Rho(i,j);
      vecSigma(i+j*d) = Sigma(i,j);
    }
  ADREPORT(Sigma);
  
  // Phi parameterized in terms
  PARAMETER_VECTOR(theta);
  PARAMETER_VECTOR(logit_lambda); vector<Type> lambda = 2/(1+exp(-logit_lambda))-1; ADREPORT(lambda)
  
  matrix<Type> eigvec(2,2);
  eigvec(0,0) = cos(theta(0)); eigvec(0,1) = cos(theta(1));
  eigvec(1,0) = sin(theta(0)); eigvec(1,1) = sin(theta(1));
  
  matrix<Type> D(2,2);
  D(0,0) = lambda(0); D(0,1) = 0;
  D(1,0) = 0;         D(1,1) = lambda(1);
  
  matrix<Type> Phi = eigvec * D * eigvec.inverse(); ADREPORT(Phi);
  
  Type nll = 0;

  matrix<Type> A(d*d,d*d);      // I_4 - kronecker(Phi,Phi)
  for (int i=0; i<d; i++)
    for (int j=0; j<d; j++)
      for (int k=0; k<d; k++)
        for (int l=0; l<d; l++)
          A(i*d+k, j*d+l) = -Phi(i,j)*Phi(k,l);
  for (int i=0; i<d*d; i++)
    A(i,i) += 1;
  matrix<Type> Ainv = A.inverse();

  vector<Type> vecGamma0 = Ainv*vecSigma; // vec(Stationary covariance matrix)
  matrix<Type> Gamma0(d,d);
  for (int i=0; i<d; i++)
    for (int j=0; j<d; j++)
      Gamma0(i,j) = vecGamma0(i+j*d);
  
  // Density of VAR(1) process random effects x
  nll += MVNORM(Gamma0)(x.row(0));
  for (int t=1; t<x.rows(); t++) {
    vector<Type> xt = x.row(t);
    vector<Type> xtminus1 = x.row(t-1);
    vector<Type> Phix = Phi*xtminus1;
    vector<Type> error = xt - Phix;
    nll += MVNORM(Sigma)(error); //VECSCALE(UNSTRUCTURED_CORR(rho),sigma)(error); // UNSTRUCTURED_CORR returns the negative log lik!
  }
  // Conditional density of observations y given x
  for (int t=0; t<y.rows(); t++) {
    for (int i=0; i<d; i++) {
      nll -= dnorm(y(t,i), x(t+x.rows()-y.rows(),i), sdy, true);
    }
  }
  
  return nll;
}
