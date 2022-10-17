#______________________________________________________________________________________________________
#______________________________________________________________________________________________________
#
# Socio-economic variation in Inuit harvest choices and its implications for climate change adaptation
# F. Hillemann, B. A. Beheim, E. Ready
# 
# R and Stan code to simulate and analyse foraging trip data (patch choice and harvest success)
# contact: friederike_hillemann@eva.mpg.de // f.hillemann@web.de
#
#______________________________________________________________________________________________________
#______________________________________________________________________________________________________


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#                                              xxxx
# set up environment                           ####
#                                              xxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

rm(list=ls())

## packages and functions ####
library(here)
library(rstan)
library(rethinking)

# simplex transforms counts into proportions of the total
# a simplex is a vector that sums to 1 and can be interpreted as probabilities
# e.g., simplex(c(3, 1)) returns "[1] 0.75 0.25"
simplex <- function(x) x/sum(x)

# softmax normalises an input vector of real numbers into a probability distribution, 
# with probabilities proportional to the exponentials of the input
softmax <- function(x) simplex(exp(x)) # same as exp(x) / sum(exp(x))

## custom function to replace NA
na2dummy <- function(data){ 
  data[is.na(data)] <- (999)
  return(data)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx####
#                                              xxxx
# simulate data                                ####
#                                              xxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

N = 300              # N hunting trips
K = 7                # K patch types
J = 25               # J harvester
Y_pc <- rep(NA, N)   # Y response patch choice
Y_hs <- rep(NA, N)   # Y response harvest success


#__________________________________________________
## create harvest trips
dd <- data.frame( trip_id = c(1:N) ,
                  season = rbinom(N, 1, 0.6)) # ice and snow Y/N

# rename season, reflecting how data is most likely entered
dd$season <- ifelse( dd$season == 0, "summer", "winter" )


#__________________________________________________
## make some people
ppl <- NULL
for (j in 1:J) {
  j_id = j
  age = as.integer(rnorm(1, 40, 7)) 
  gender = rbinom(1, 1, 0.8) 
  income = as.integer(rnorm(1, mean=60000, sd=25000))
  indegree = rpois(1, lambda=2.5)
  outdegree = rpois(1, lambda=2.5)
  ppl = rbind(ppl, data.frame(j_id, gender, age, income, indegree, outdegree))
}

# for age and gender, make sure we have someone from each gender and age_cat
ppl$age[1:4] <- c(25, 35, 45, 55)
ppl$gender[1:2] <- c(1, 2)

# use age categories
ppl$age_cat[ppl$age %in% c( 0:30)]<- 1
ppl$age_cat[ppl$age %in% c(30:40)]<- 2
ppl$age_cat[ppl$age %in% c(40:50)]<- 3
ppl$age_cat[ppl$age %in% c(50:70)]<- 4

# rename gender (reflecting how data is most likely entered)
ppl$gender <- ifelse( ppl$gender == 0, "f", "m" )

# use standardised variables
ppl$income_s <- (ppl$income - mean(ppl$income, na.rm=TRUE) ) / sd(ppl$income, na.rm=TRUE)
ppl$indegree_s <- (ppl$indegree - mean(ppl$indegree) ) / sd(ppl$indegree)
ppl$outdegree_s <- (ppl$outdegree - mean(ppl$outdegree) ) / sd(ppl$outdegree)


#__________________________________________________
## send people harvesting
# to draw IDs, dgeom works well to reflect that few people go out a lot, most only occacionally 
# plot( dgeom( c(1:J), prob=0.2) )
j_once <- c(1:J) # making sure each hunter shows up at least once
dd$j_id <- sample(c( j_once , sample( c(1:J), size=N-length(j_once), replace=TRUE, prob=dgeom(c(1:J), 0.2))))
dd$age_cat <- ppl$age_cat[match(dd$j_id, ppl$j_id)]
dd$gender <- ppl$gender[match(dd$j_id, ppl$j_id)]
dd$income_s <- ppl$income_s[match(dd$j_id, ppl$j_id)]
dd$indegree_s <- ppl$indegree_s[match(dd$j_id, ppl$j_id)]
dd$outdegree_s <- ppl$outdegree_s[match(dd$j_id, ppl$j_id)]

# group size, standardised
Nhunt <- sample(c(1:10), size=N, replace=TRUE, prob = (dgeom(c(1:10), prob=0.3)) )
dd$Nhunt_s <- (Nhunt - mean(Nhunt)) / sd(Nhunt)


#__________________________________________________
# prepare data for Stan
# Stan needs integer IDs as input (1:x instead of charcter strings)
# e.g., if patch data were entered as "winter marine" etc, or j_id as "ID01" etc
# dd$patch_id <- as.integer(as.factor( dd$patch_cat )) 
# dd$hunter_id <- as.integer(as.factor( dd$j_id ))
dd$season_id <- ifelse( dd$season=="winter", 1 , 2 ) # 1-snow/ice, 2-ice-free
dd$gender_id = ifelse( dd$gender=="m", 1 , 2 ) # 1-male, 2-female


#__________________________________________________
## set patch choice probability ("attraction weight") intercepts
# store in an array [patch x season x gender]
# i.e., use different schedules by season and gender
# (instead of using a vector of length K, and adding gender and season as effects)
# note, some patches are winter-only, others are summer-only, some can be chosen year-round
WA <- array(NA, dim = c(K, 2, 2))
WA[ ,1,1] <- c(-2, -2, -2,  1,  1, 1.2, 0)    # season1, gender1
WA[ ,1,2] <- c(-2, -2, -2,  1,  1, 0.2, 0)    # season1, gender2
WA[ ,2,1] <- c( 1,  1,  1, -2, -2, 1.2, 0)    # season2, gender1
WA[ ,2,2] <- c( 1,  1,  1, -2, -2, 0.2, 0)    # season2, gender2


## set patch success probability intercepts
WS <- array(NA, dim = c(K, 2, 2))
WS[ ,1,1] <- logit(c(0.01, 0.01, 0.01, 0.5 , 0.6 , 0.8, 0))    # season1, gender1
WS[ ,1,2] <- logit(c(0.01, 0.01, 0.01, 0.6 , 0.7 , 0.8, 0))    # season1, gender2
WS[ ,2,1] <- logit(c(0.4 , 0.6 , 0.5 , 0.01, 0.01, 0.7, 0))    # season2, gender1
WS[ ,2,2] <- logit(c(0.6 , 0.7 , 0.5 , 0.01, 0.01, 0.7, 0))    # season2, gender2


## set other effects
# for each patch category, we set a different effects
# bIn: income, bDeI/bDeO: in-/outdegree, N hunters: group size
# bIn  <- round(rnorm(K, 0.2, 0.5), 1)
# bDeI <- round(rnorm(K, 0, 0.3), 1)
# bDeO <- round(rnorm(K, 0, 0.3), 1)
# bNh  <- round(rnorm(K, 0.2, 0.3), 1)
# bIn[K] <- bDeI[K]  <- bDeO[K] <- bNh[K] <- 0 # reference cat
# here, we use observed effect sizes
bIn  <- c(0.07, -0.04, 0.45, 0.23, -0.02, -0.49, 0.01)
bDeI <- c(-0.3, -0.1, 0.2, 0.1, -0.7, -0.6, 0)
bDeO <- c(-0.1, -0.1, -0.2, 0.2, 0.2, 0.4, 0)
bNh  <- c(0.6, 0.4, -0.1, 0.1, -0.1, -0.6, 0)



# age-category: log-odds offsets for the probability of choice
fAge <- matrix(c(0.1, 0.0, -0.2, -0.2,
                 0.4, 0.3,  0.2,  0.0,
                 0.2, 0.1, -0.1, -0.1,
                 0.2, 0.1,  0.0, -0.1,
                 0.4, 0.3,  0.1, -0.1,
                 0.3, 0.1, -0.2, -0.3,
                 0.2, 0.1,  0.0, -0.1),
               nrow=K, 
               ncol=length(unique((dd$age_cat))), 
               byrow=TRUE)



#__________________________________________________
## Y_pc: simulate patch choice data
Y_pc <- rep(NA, N)

for (n in 1:N) {
  A <- WA[ , dd$season_id[n], dd$gender_id[n]] +
        fAge[ ,dd$age_cat[n]] +
        bIn*dd$income_s[n] +
        bDeI*dd$indegree_s[n] +
        bDeO*dd$outdegree_s[n] +
        bNh*dd$Nhunt_s[n]
  Y_pc[n] <- sample(1:K, 1, prob = softmax(A))
}

dd$patch_id <- Y_pc

# >>> check: all patches can be chosen in any season ####
table(dd$patch_id, dd$season)
# not ideal that all patches can be chosen in any season
# can we weigh the effect of season more?

for (n in 1:N) {
  A <- WA[ , dd$season_id[n], dd$gender_id[n]] * 10 +
    fAge[ ,dd$age_cat[n]] +
    bIn*dd$income_s[n] +
    bDeI*dd$indegree_s[n] +
    bDeO*dd$outdegree_s[n] +
    bNh*dd$Nhunt_s[n]
  Y_pc[n] <- sample(1:K, 1, prob = softmax(A))
}

dd$patch_id <- Y_pc
table(dd$patch_id, dd$season)

#__________________________________________________
## Y_hs: simulate harvest success data

Y_hs <- rep(NA, N)

for (n in 1:N) {
  S <- WS[ , dd$season_id[n], dd$gender_id[n]] +
        fAge[ , dd$age_cat[n]] +
        bIn*dd$income_s[n] +
        bDeI*dd$indegree_s[n] +
        bDeO*dd$outdegree_s[n] + 
        bNh*dd$Nhunt_s[n]
  Y_hs[n] <- rbern(1, logistic(S[dd$patch_id[n]]))
}

dd$harvest <- Y_hs



#__________________________________________________
## missing value imputation in Stan
# realistically, we may not have income data for everyone
j_id_incomplete <- sample(J, size=sample(2,1))
ppl$income[ppl$j_id %in% j_id_incomplete] <- NA
dd$income[dd$j_id %in% ppl$j_id[is.na(ppl$income)] ] <- NA

# store count of missing values (missN) and index of missing values (missDEX)
names(which(colSums(is.na(dd)) > 0))
missn_income <- sum(is.na(dd$income))
missdex_income <- which(is.na(dd$income))

# replace NAs with 999 (or any other arbitraty number)
dd <- na2dummy(dd)




#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx####
#                                              xxxx
# fit model                                    ####
#                                              xxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

## set number of chains and iterations
ch=1
it=1000


# data list for Stan
datlist = list(N = N,                    # number of harvest episodes,
               K = K,                    # number of patches 
               Y_pc = dd$patch_id,       # patch choice ID
               Y_hs = dd$harvest,        # harvest success Y/N
               season = dd$season_id,    # 1-snow/ice, 2-ice-free 
               gender = dd$gender_id,    # 1-male, 2-female
               age_cat = dd$age_cat, 
               age_cat_index = 1:length(unique(dd$age_cat)),
               N_age_cats = length(unique(dd$age_cat)),
               income = dd$income_s,                  
               missn_income   = missn_income,     # missN (number of missing values for income)
               missdex_income = missdex_income,   # missDEX (index of missing values for income)
               indegree = dd$indegree_s, 
               outdegree = dd$outdegree_s, 
               Nhunt = dd$Nhunt_s) 


#__________________________________________________
#                                              xxxx
#* patch choice                                ####
#                                              xxxx
#__________________________________________________

dear_stan = '
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
  // // alternatively use the merge_missing function (top of script):
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
'

# run Stan model
mm_pc <- stan(model_code=dear_stan, data=datlist, chains=ch, iter=it, warmup=200)
#mm_pc <- stan(file=paste0(here(),"/m_pc.stan"), data=datlist, chains=ch, iter=it, warmup=200)

## explore sampling behavior and assess mixing across Markov chains
traceplot(mm_pc, pars="WA")
traceplot(mm_pc, pars="bIn")
traceplot(mm_pc, pars="fAge")
traceplot(mm_pc, pars=c("bDeI", "bDeO"))
traceplot(mm_pc, pars="bNh")


## compare estimated and actual probabilities of each choice
post_pc <- extract.samples(mm_pc)
pc_true <- pc_est <- pc_lb <- pc_ub <- array(NA, dim=c(N,K,2,2))

for (g in 1:length(unique(dd$gender)) ) {
  for (s in 1:length(unique(dd$season)) ) { 
    for (n in 1:N) {
      
      # calculate true prob of patch choice
      pc_true[n, ,s,g] <- softmax( WA[ ,s,g] + # intercept
                                     fAge[ ,dd$age_cat[n]] +
                                     bIn*dd$income_s[n] +
                                     bDeI*dd$indegree_s[n] +
                                     bDeO*dd$outdegree_s[n] +
                                     bNh*dd$Nhunt_s[n] )
      
      # calculate estimate of posterior prob of patch choice
      A <- post_pc$WA[ ,  , s, g] + # intercept
            post_pc$fAge[ , , dd$age_cat[n]] +
            post_pc$bIn*dd$income_s[n] +
            post_pc$bDeI*dd$indegree_s[n] +
            post_pc$bDeO*dd$outdegree_s[n] +
            post_pc$bNh*dd$Nhunt_s[n]
      pc_n <- t(apply(A, 1, softmax))
      pc_est[n, ,s,g] <- apply(pc_n, 2, mean)
      pc_lb[n, ,s,g] <- apply(pc_n, 2, rethinking::HPDI)[1,]
      pc_ub[n, ,s,g] <- apply(pc_n, 2, rethinking::HPDI)[2,]
      if (n %% 100 == 0) print(n)
    }
  }
}


plot(pc_true, pc_est, xlim = c(0, 0.5), ylim = c(0, 0.5))
abline(0, 1)
for (n in 1:N) {
  for (k in 1:K) {
    for (g in 1:length(unique(dd$gender)) ) {
      for (s in 1:length(unique(dd$season)) ) { 
        lines(c(pc_true[n,k,s,g], pc_true[n,k,s,g]), c(pc_lb[n,k,s,g], pc_ub[n,k,s,g]))
      }
    }
  }
}





#__________________________________________________
#                                              xxxx
#* harvest success                             ####
#                                              xxxx
#__________________________________________________

dear_stan='
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
'

# run Stan model
mm_hs <- stan(model_code=dear_stan, data=datlist, chains=ch, iter=it, warmup=200)
#mm_hs <- stan(file=paste0(here(),"/m_hs.stan"), data=datlist, chains=ch, iter=it, warmup=200)

## explore sampling behavior and assess mixing across Markov chains
traceplot(mm_hs, pars="bIn")
traceplot(mm_hs, pars="fAge")
traceplot(mm_hs, pars="bDeI")


# compare estimated and actual probabilities of each choice
post_hs <- extract.samples(mm_hs)
hs_true <- hs_est <- hs_lb <- hs_ub <- array(NA, dim=c(N,K,2,2))

for (g in 1:length(unique(dd$gender)) ) {
  for (s in 1:length(unique(dd$season)) ) { 
    for (n in 1:N) {
      
      # use true patch attraction intercepts to calculate true prob of patch choice
      hs_true[n, ,s,g] <- logistic ( WS[ ,s,g] +
                                      fAge[ ,dd$age_cat[n]] +
                                      bIn*dd$income_s[n] + 
                                      bDeI*dd$indegree_s[n] + 
                                      bDeO*dd$outdegree_s[n] +
                                      bNh*dd$Nhunt_s[n] )
      
      # use posterior attraction weights to calculate estimate of posterior prob of patch choice
      S <- post_hs$WS[ ,  , s, g] +
            post_hs$fAge[ , , dd$age_cat[n]] + 
            post_hs$bIn*dd$income_s[n] + 
            post_hs$bDeI*dd$indegree_s[n] +
            post_hs$bDeO*dd$outdegree_s[n] + 
            post_hs$bNh*dd$Nhunt_s[n]
      hs_n <- t(apply(S, 1, logistic))
      hs_est[n, ,s,g] <- apply(hs_n, 2, mean)
      hs_lb[n, ,s,g] <- apply(hs_n, 2, rethinking::HPDI)[1,]
      hs_ub[n, ,s,g] <- apply(hs_n, 2, rethinking::HPDI)[2,]
      if (n %% 100 == 0) print(n)
    }
  }
}


# >>> check #####

# so this doesn't make sense and it doesn't look great
# but why? too many things mixed together?!
plot(hs_true, hs_est, xlim = c(0, 1), ylim = c(0, 1))
abline(0, 1)
for (n in 1:N) {
  for (k in 1:K) {
    for (g in 1:length(unique(dd$gender)) ) {
      for (s in 1:length(unique(dd$season)) ) { 
        lines(c(hs_true[n,k,s,g], hs_true[n,k,s,g]), c(hs_lb[n,k,s,g], hs_ub[n,k,s,g]))
      }
    }
  }
}

# ... this is better, right? (note all patches are plotted although some are/should be impossible in a given season (see other "check" comment))
plot( apply(post_hs$WS[,,1,1], 2, mean), WS[,1,1] )
text( apply(post_hs$WS[,,1,1], 2, mean), WS[,1,1], cex=0.5, pos=2 )
abline(0, 1, lty = 2)

plot( apply(post_hs$WS[,,2,1], 2, mean),  WS[,2,1] )
text( apply(post_hs$WS[,,2,1], 2, mean), WS[,2,1], cex=0.5, pos=2 )
abline(0, 1, lty = 2)

# ... and then I plot true vs estimated effect of specific covariates?
# plot true vs estimated bIn
plot(bIn, apply(post_hs$bIn, 2, mean))
text(bIn, apply(post_hs$bIn, 2, mean), labels=c(1:7), cex=0.8, pos=1)
abline(0, 1, lty = 2)
bIn_lb <- apply(post_hs$bIn, 2, HPDI)[1, ]
bIn_ub <- apply(post_hs$bIn, 2, HPDI)[2, ]
for (k in seq_len(K)) {
  lines(c( bIn[k], bIn[k] ), c( bIn_lb[k], bIn_ub[k] ))
}





#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx####
#                                              xxxx
# result plots and tables                      ####
#                                              xxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

cols <- rainbow(K)
# colours used in MS are linked to season and marine/terrestrial (darker/paler tone)
# cols <- c("lavenderblush3", "darkseagreen3","gold2","powderblue","goldenrod3","steelblue","lightslategrey")

#__________________________________________________
#                                              xxxx
#* Fig.2: predictive plots                     ####
#                                              xxxx
#__________________________________________________

# example: income

xlab = "income"
# set steps of x-axis
sort(ppl$income_s)
# values might range between -2.5 and +2.5
# note, one can use observed range or one could extrapolate here
x <- seq(-2.5, 2.5, by = 0.05) 


# create a fictive hunter profile
fix_agecat <- 2
fix_indegree <- ppl$indegree_s[ppl$indegree==2][1] # use indegree_s that refers to indegree of 2
fix_outdegree <- ppl$outdegree_s[ppl$outdegree==2][1] # use outdegree_s that refers to outdegree of 2
fix_Nhunt <- sort(unique(dd$Nhunt))[1] # Nhunt ranged from 1:10, so this gives us Nhunt_s that refers to a group size of 1


# posterior probabilities of each choice and success (estimated mean, lower and upper boud)
pr_pc_est <- pr_pc_lb <- pr_pc_ub <- array(NA, dim=c(length(x), 7, 2, 2))
pr_hs_est <- pr_hs_lb <- pr_hs_ub <- array(NA, dim=c(length(x), 7, 2, 2))

for (g in 1:2) {
  for (s in 1:2) { 
    for (i in 1:length(x)) {
        A <- post_pc$W[ ,  , s, g] + 
              post_pc$fAge[,,fix_agecat] + # or post_pc$bAg*fix_age_s_40 + post_pc$bAg2*fix_age2_s_40 + 
              post_pc$bIn * x[i] + 
              post_pc$bDeI * fix_indegree +
              post_pc$bDeO * fix_outdegree +
              post_pc$bNh * fix_Nhunt
      pr_pc_n <- t(apply(A, 1, softmax)) # transform to make prob of all patches sum to 1
      pr_pc_est[i, ,s,g] <- apply(pr_pc_n, 2, mean)
      pr_pc_lb[i, ,s,g] <- apply(pr_pc_n, 2, HPDI)[1,]
      pr_pc_ub[i, ,s,g] <- apply(pr_pc_n, 2, HPDI)[2,]
      
      s_post <- 
        S <- post_hs$WS[ ,  , s, g] + 
              post_hs$fAge[,,fix_agecat] + # or post_hs$bAg*fix_age_s_40 + post_hs$bAg2*fix_age2_s_40 + 
              post_hs$bIn * x[i] + 
              post_hs$bDeI * fix_indegree +
              post_hs$bDeO * fix_outdegree +
              post_hs$bNh * fix_Nhunt
      pr_hs_n <- logistic(S) # transform into probabilities
      pr_hs_est[i, ,s,g] <- apply(pr_hs_n, 2, mean)
      pr_hs_lb[i, ,s,g] <- apply(pr_hs_n, 2, HPDI)[1,]
      pr_hs_ub[i, ,s,g] <- apply(pr_hs_n, 2, HPDI)[2,]
    }
  }
} 


# plot (and save the 4-pages output as PDF on desktop)
# pdf(file=paste0(file.path(path.expand('~'),'Desktop'), "/test.pdf"), width=5, height=8)
for (g in 1:2) {
  for (s in 1:2) {
    
    title_s <- dd$season[dd$season_id == s][1]
    title_g <- dd$gender[dd$gender_id == g][1]
    
    par(mfrow=c(K, 2), mar=c(2,4,1,4))
    
    for (k in 1:K) { 
      plot(1, 1, type = "n", xlim=c(min(x), max(x)), ylim=c(0,1), las=1,
           xlab = xlab, ylab = "pr(choose)", cex.axis=0.9,
           main = paste( title_s, title_g, sep=", "), cex.main=0.6, xaxt="n")
      axis(1, at = c(-1.5, 0, 1.5))
      points(x, pr_pc_est[ ,k,s,g], type = "l", col=cols[k])
      polygon(c(x, rev(x)), c(pr_pc_lb[ ,k,s,g], rev(pr_pc_ub[ ,k,s,g])), col=col.alpha(cols[k], 0.2), border=NA)
      
      plot(1, 1, type = "n", xlim=c(min(x), max(x)), ylim=c(0,1), las=1,
           xlab = xlab, ylab = "pr(success)", cex.axis=0.9, 
           main = paste( title_s, title_g, sep=", "), cex.main=0.6, xaxt="n")
      axis(1, at = c(-1.5, 0, 1.5))
      points(x, pr_hs_est[ ,k,s,g], type = "l", col=cols[k])
      polygon(c(x, rev(x)), c(pr_hs_lb[ ,k,s,g], rev(pr_hs_ub[ ,k,s,g])), col=col.alpha(cols[k], 0.2), border=NA)
    }
  }
}
# dev.off() # produce pdf




#__________________________________________________
#                                              xxxx
#* Fig.3: Prob PC ~ Prob HS                    ####
#                                              xxxx
#__________________________________________________

# create a fictive hunter profile
fix_agecat <- 2
fix_income <- median(ppl$income_s, na.rm = TRUE) # use median of income (standardised)
fix_indegree <- ppl$indegree_s[ppl$indegree==2][1] # use indegree_s that refers to indegree of 2
fix_outdegree <- ppl$outdegree_s[ppl$outdegree==2][1] # use outdegree_s that refers to outdegree of 2
fix_Nhunt <- sort(unique(dd$Nhunt))[1] # Nhunt ranged from 1:10, so this gives us Nhunt_s that refers to a group size of 1
  
# posterior probabilities of each choice and success (estimated mean, lower and upper boud)
pc_est <- pc_lb <- pc_ub <- array(NA, dim=c(K, 2, 2))
hs_est <- hs_lb <- hs_ub <- array(NA, dim=c(K, 2, 2))

for (g in 1:2 ) {
  for (s in 1:2 ) { 
    
    A <- post_pc$WA[ , ,s,g] + 
          post_pc$fAge[ , ,fix_agecat] + 
          post_pc$bIn*fix_income +
          post_pc$bDeI*fix_indegree +
          post_pc$bDeO*fix_outdegree +
          post_pc$bNh*fix_Nhunt
    
    pc_n <- t(apply(A, 1, softmax)) # transform to make prob of all patches sum to 1
    pc_est[ ,s,g] <- apply(pc_n, 2, mean)
    pc_lb[ ,s,g] <- apply(pc_n, 2, HPDI)[1,]
    pc_ub[ ,s,g] <- apply(pc_n, 2, HPDI)[2,]
    
    
    S <- post_hs$WS[ , ,s,g] + 
        post_hs$fAge[ , ,fix_agecat] +
        post_hs$bIn*fix_income +
        post_hs$bDeI*fix_indegree +
        post_hs$bDeO*fix_outdegree +
        post_hs$bNh*fix_Nhunt
    
    hs_n <- logistic(S) # transform into probabilities
    hs_est[ ,s,g] <- apply(hs_n, 2, mean)
    hs_lb[ ,s,g] <- apply(hs_n, 2, HPDI)[1,]
    hs_ub[ ,s,g] <- apply(hs_n, 2, HPDI)[2,]
    
  }
}

par(mfrow=c(2,2), mar=c(4.1, 4.1, 2.1, 2.1))
for (g in 1:2 ) {
  for (s in 1:2) {
    
    title_s <- dd$season[dd$season_id == s][1]
    title_g <- dd$gender[dd$gender_id == g][1]
    
    plot(pc_est[ ,s,g], hs_est[ ,s,g],
         xlim=c(0,0.6), ylim=c(0,1), las=1,
         xlab="prob of choosing", ylab="prob of success",
         pch=3, lwd=3, col=cols,
         main = paste( title_s, title_g, sep=", "))
     for (k in 1:K) {
       lines( x=c(pc_lb[k,s,g], pc_ub[k,s,g]), y=c(hs_est[k,s,g], hs_est[k,s,g]), col=cols[k], lty=3)
       lines( x=c(pc_est[k,s,g], pc_est[k,s,g]), y=c(hs_lb[k,s,g], hs_ub[k,s,g]), col=cols[k], lty=3)
     }
  }#s
}#g
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 4.1))



#__________________________________________________
#                                              xxxx
#* supplementatry Figure: effect sizes         ####
#                                              xxxx
#__________________________________________________

plot(mm_hs, pars="bIn", 
     show_density = FALSE,
     point_est = "mean",
     ci_level = 0.89,
     outer_level = 0.89,
     fill_color = cols, 
     outline_color = cols,
     est_color = cols)




#__________________________________________________
#                                              xxxx
#* supplementatry Table: effect sizes          ####
#                                              xxxx
#__________________________________________________

rethinking::precis( as.data.frame(mm_pc), depth=2 )
# or
print(mm_pc, probs = c(0.025, 0.5, 0.975))
print(mm_pc, probs = c(0.025, 0.5, 0.975), pars="WA")
print(mm_pc, probs = c(0.025, 0.5, 0.975), pars="bIn")
print(mm_pc, probs = c(0.025, 0.5, 0.975), pars="fAge")
print(mm_pc, probs = c(0.025, 0.5, 0.975), pars="bDeI")
print(mm_pc, probs = c(0.025, 0.5, 0.975), pars="bDeO")
print(mm_pc, probs = c(0.025, 0.5, 0.975), pars="bNh")

