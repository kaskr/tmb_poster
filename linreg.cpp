#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data section
  DATA_VECTOR(Y);
  DATA_VECTOR(x);
  // Parameter section
  PARAMETER(a);
  PARAMETER(b);
  PARAMETER(logSigma);
  // Procedure section
  Type nll = -sum(dnorm(Y, a + b * x, exp(logSigma), true));
  return nll;
}

