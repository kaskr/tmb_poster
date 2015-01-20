#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_ARRAY(obs); /* timeSteps x stateDim */
  PARAMETER_ARRAY(X); /* State */
  PARAMETER(transf_rho);
  PARAMETER_VECTOR(logsds);
  PARAMETER_VECTOR(logsdObs);
  int timeSteps=obs.dim[1];
  int stateDim=obs.dim[0];
  // Parameter transform
  Type rho=Type(2)/(Type(1) + exp(-Type(2) * transf_rho)) - Type(1);
  vector<Type> sds=exp(logsds);
  vector<Type> sdObs=exp(logsdObs);
  // Setup the covariance matrix Sigma
  matrix<Type> cov(stateDim,stateDim);
  for(int i=0;i<stateDim;i++)
    for(int j=0;j<stateDim;j++)
      cov(i,j)=pow(rho,Type(abs(i-j)))*sds[i]*sds[j];
  // Setup the multivariate normal density
  using namespace density;
  MVNORM_t<Type> neg_log_density(cov);
  // Define likelihood
  Type ans=0;
  // Process likelihood: Draw initial state
  ans-=dnorm(vector<Type>(X.col(0)),Type(0),Type(1),1).sum();
  // Process likelihood: Draw states
  for(int i=1;i<timeSteps;i++)    
    ans+=neg_log_density(X.col(i)-X.col(i-1)); 
  // Data likelihood: Draw the data given the states
  for(int i=1;i<timeSteps;i++)
    ans-=dnorm(vector<Type>(obs.col(i)),vector<Type>(X.col(i)),sdObs,1).sum(); 
  return ans;
}
