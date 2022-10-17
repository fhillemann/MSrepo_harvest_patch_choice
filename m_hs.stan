data {
  int<lower=0> N;                    // number of observations, i.e. harvest episodes
  int<lower=2> K;                    // number of patch types [Kategorien]
  int<lower=0,upper=1> Y_hs[N];      // harvest success Y/N
  // covariates
  int<lower=1,upper=K> Y_pc[N];      // chosen patch type ID
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
   matrix [2,2] WS[K]; // Kx2x2 array holding intercepts
   vector[K] bIn;      // effect of income
   vector[K] bDeI;	   // effect for indegree  
   vector[K] bDeO;	   // effect for outdegree  
   vector[K] bNh;			 // effect for group size
   
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
   
  // impute missing values
  real income_merge[N];
  income_merge = income;
  for (i in 1:missn_income) income_merge[missdex_income[i]] = income_impute[i];
   
  // latent variable Gaussian process
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
  
  vector[N] S;
  
  // priors
  for (s in 1:2){
    for (g in 1:2){
      to_vector(WS[ ,s,g]) ~ normal(0, 0.5); // patch success intercept, by season and gender
    }
  }
  bIn ~ normal(0,0.5);    // income
  bDeI ~ normal(0,0.5);   // indegree
  bDeO ~ normal(0,0.5);   // outdegree
  bNh ~ normal(0,0.5);    // group size
  
  // priors for Gaussian process model - age
  alpha ~ std_normal();  // same as ~ normal(0, 1);
  rho ~ inv_gamma(5, 5);
  eta ~ std_normal();
  
  // priors for missing values imputation - income
  sigma_income ~ cauchy(0, 1);
  mu_income ~ normal(0.5, 1);
  income_merge ~ normal(mu_income, sigma_income);

  // Likelihood function for within-patch harvest success model
  for (n in 1:N) {
    S[n] = WS[Y_pc[n], season[n], gender[n]]  // intercept, patch-specific success weight
              + fAge[Y_pc[n],age_cat[n]]
              + bIn[Y_pc[n]]*income_merge[n]
              + bDeI[Y_pc[n]]*indegree[n]
              + bDeO[Y_pc[n]]*outdegree[n]
              + bNh[Y_pc[n]]*Nhunt[n];
    target += bernoulli_logit_lpmf( Y_hs[n] | S[n]); // same as bernoulli_lpmf( Y_hs[n] | inv_logit(S[n])) or Y_hs[n] ~ bernoulli_logit(S[n]);
  }
}
generated quantities {
  // calculate log likelihood for model comparisons
  // and store expected values for each observation
  vector[N] S; 
  vector[N] log_lik_hs;
   
  for (n in 1:N) {
    S[n] = WS[Y_pc[n],season[n],gender[n]]  // intercept, patch-specific success weight
            + fAge[Y_pc[n],age_cat[n]]
            + bIn[Y_pc[n]]*income_merge[n]
            + bDeI[Y_pc[n]]*indegree[n]
            + bDeO[Y_pc[n]]*outdegree[n]
            + bNh[Y_pc[n]]*Nhunt[n];
            
    log_lik_hs[n] = bernoulli_lpmf( Y_hs[n] | inv_logit(S[n]) );
  }
}

// end

// check if code is syntactically correct
// rstan:::rstudio_stanc(paste0(here(), "/m_hs.stan"))
