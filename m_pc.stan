// // merge_missing is a function for imputing missing values
// functions{
//   vector merge_missing( int[] miss_indexes , vector x_obs , vector x_miss ) {
//      int M = dims(x_obs)[1];
//      int M_miss = dims(x_miss)[1];
//      vector[M] merged;
//      merged = x_obs;
//      for ( i in 1:M_miss )
//        merged[ miss_indexes[i] ] = x_miss[i];
//      return merged;
//   }
// }
data {
  int<lower=0> N;                    // number of observations, i.e. harvest episodes
  int<lower=2> K;                    // number of patch types [Kategorien]
  int<lower=1,upper=K> Y_pc[N];      // Y, choice of patch type ID
  // covariates
  int<lower=1,upper=2> season[N];    // 1-ice, 2-no ice
  int<lower=1,upper=2> gender[N];    // 1-male, 2-female
  int<lower=1> age_cat[N];           // age category
  int<lower=1> N_age_cats;           // number of age categories
  real age_cat_index[N_age_cats];    // 1 to N_age_cats; cov_exp_quad needs real
  real income[N];                    // income, standardised
  real indegree[N];                  // indegree, standardised
  real outdegree[N];                 // outdegree, standardised
  real Nhunt[N];                     // group size, number of hunters on trip, standardised
  // missing values
  int<lower=1> missn_income;          // number of missing values
  int missdex_income[missn_income];   // index of missing values in vector
}
transformed data {
  real delta = 1e-9; // add to the diagonal of the GP age covariance matrix to ensure positive definite values
}
parameters {
   matrix [2,2] WA[K];  // Kx2x2 array holding intercepts for each patch
   real bIn[K];			    // fixed effect for income
   real bDeI[K];		    // fixed effect for indegree  
   real bDeO[K];		    // fixed effect for outdegree  
   real bNh[K];			    // fixed effect for group size
   
  // Gaussian process model - age
  vector<lower=0>[K] alpha;
  vector<lower=0>[K] rho;
  vector[N_age_cats] eta; // scaling factor
   
  // impute missing values - income
  real mu_income;
  real<lower=0> sigma_income;
  vector[missn_income] income_impute;
}
transformed parameters {
  
  matrix[K, N_age_cats] fAge;
  
  // impute missing values - income
  real income_merge[N];
  income_merge = income;
  for (i in 1:missn_income) income_merge[missdex_income[i]] = income_impute[i];
  // // alternatively, use the merge_missing function (top of script):
  // real income_merge = merge_missing(missdex_income, to_vector(income), income_impute);
  
  // latent variable Gaussian process - age
  for (k in 1:K){
    matrix[N_age_cats, N_age_cats] L_cov;
    matrix[N_age_cats, N_age_cats] cov = cov_exp_quad(age_cat_index, alpha[k], rho[k]); //covariance matrix
      
    // diagonal elements
    for (n in 1:N_age_cats) {
      cov[n, n] = cov[n, n] + delta;
      }
      
    L_cov = cholesky_decompose(cov);
    fAge[k, ] = to_row_vector( L_cov * eta );
  }
}
model{
  // priors 
  for (s in 1:2){
    for (g in 1:2){
      to_vector(WA[ ,s,g]) ~ normal(0, 0.5); // patch choice intercept, by season and gender
    }
  }
  bIn ~ normal(0,0.5);    // income
  bDeI ~ normal(0,0.5);   // indegree
  bDeO ~ normal(0,0.5);   // outdegree
  bNh ~ normal(0,0.5);    // group size
  
  // priors for Gaussian process model - age
  alpha ~ normal(0, 0.5); // std_normal() is same as ~ normal(0, 1);
  rho ~ inv_gamma(5, 5);
  eta ~ normal(0, 0.5); // std_normal();
  
  // priors for missing values imputation - income
  sigma_income ~ cauchy(0, 1);
  mu_income ~ normal(0.5, 1);
  income_merge ~ normal(mu_income, sigma_income);
  
  // Likelihood functions for each patch choice model
  for ( n in 1:N ) {
    vector[K] A;
    for ( k in 1:(K-1) ) { 
      A[k] = WA[k,season[n],gender[n]]
            + fAge[k, age_cat[n]]
            + bIn[k] * income_merge[n]
            + bDeI[k] * indegree[n]
            + bDeO[k] * outdegree[n]
            + bNh[k] * Nhunt[n];
    }
    A[K] = 0;  // K-1 handling
     
    // multinomial logistic regression
    target += categorical_logit_lpmf( Y_pc[n] | A );  // same as Y_pc[n] ~ categorical_logit( A ); 
  }
}
generated quantities{
  vector[N] log_lik_pc;
   
  for ( n in 1:N ) {
    vector[K] A;
      for ( k in 1:(K-1) ) { 
        A[k] = WA[k,season[n],gender[n]]
                + fAge[k, age_cat[n]]
                + bIn[k] * income_merge[n]
                + bDeI[k] * indegree[n]
                + bDeO[k] * outdegree[n]
                + bNh[k] * Nhunt[n];
      }
      A[K] = 0;  // K-1 handling
     
    // generate the likelihood of each observation, conditional on the model
    log_lik_pc[n] = categorical_logit_lpmf( Y_pc[n] | A ); 
        
  }
}

// end

// check if code is syntactically correct
// rstan:::rstudio_stanc(paste0(here(), "/m_pc.stan"))
