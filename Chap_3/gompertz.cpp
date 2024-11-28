#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR( log_b_t );

  // Parameters
  PARAMETER( log_d0 );
  PARAMETER( log_sigmaP );
  PARAMETER( log_sigmaM );
  PARAMETER( alpha );
  PARAMETER( rho );
  PARAMETER_VECTOR( log_d_t );

  // Objective funcction
  Type jnll = 0;
  vector<Type> log_dhat_t( log_d_t.size() );

  // Probability of random coefficients
  jnll -= dnorm( log_d_t(0), log_d0, exp(log_sigmaP), true );
  for( int t=1; t<log_b_t.size(); t++){
    log_dhat_t(t) = alpha + rho*log_d_t(t-1);
    jnll -= dnorm( log_d_t(t), log_dhat_t(t), exp(log_sigmaP), true );
  }

  // Probability of data conditional on fixed and random effect values
  for( int t=0; t<log_b_t.size(); t++){
    if( !R_IsNA(asDouble(log_b_t(t))) ){
      jnll -= dnorm( log_b_t(t), log_d_t(t), exp(log_sigmaM), true );
    }
  }

  // Reporting
  REPORT( log_dhat_t );
  ADREPORT( log_dhat_t );
  return jnll;
}
