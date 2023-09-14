#______________________________________________________________________________________________________
#______________________________________________________________________________________________________
#
# R and Stan code to simulate and analyse foraging trip data (patch choice and harvest success)
# 
# Repository available and maintained on Github and Data Dryad:
# https://github.com/fhillemann/MSrepo_harvest_patch_choice.git
# https://doi.org/10.5061/dryad.k3j9kd5dv 
#
# Accompanying manuscript:
# Socio-economic predictors of Inuit hunting choices and their implications for climate change adaptation
# F. Hillemann, B. A. Beheim, E. Ready
# Phil. Trans. R. Soc. B 2023, 378: 20220395. 
# https://doi.org/10.1098/rstb.2022.0395
#
#______________________________________________________________________________________________________
#______________________________________________________________________________________________________


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#                                              xxxx
# reproducibility notes                        ####
#                                              xxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

rm(list=ls())

# all analyses reported in MS were conducted using R version 4.2.1 and Stan 2.21.0
R.Version()$version.string
stan_version()

seed <- 186
set.seed(seed) # produces the output file data_seed186

# seed <- 152
# set.seed(152) # produces the output file data_seed152


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#                                              xxxx
# set up environment                           ####
#                                              xxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

## packages and functions ####
library(rstan) # instructions for downloading and installing RStan: https://mc-stan.org/users/interfaces/rstan.html
library(rethinking) # alternatively, install via devtools::install_github("rmcelreath/rethinking@slim")
if (!require(here)) install.packages(here); library(here) # to locate loading (stan) files; or use .Rproj
if (!require(xtable)) install.packages(xtable); library(xtable) # to export LaTeX-ready tables


options(warnPartialMatchDollar=TRUE) 
# squawks if you refer to a variable in a data frame that only partially matches 


# simplex transforms counts into proportions of the total
# a simplex is a vector that sums to 1 and can be interpreted as probabilities
# e.g., simplex(c(3, 1)) returns "[1] 0.75 0.25"
simplex <- function(x) x/sum(x)

# softmax normalises an input vector of real numbers into a probability distribution, 
# with probabilities proportional to the exponentials of the input
softmax <- function(x) simplex(exp(x)) # same as exp(x) / sum(exp(x))

# custom function to replace NA
na2dummy <- function(data){ 
  data[is.na(data)] <- (999)
  return(data)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx####
#                                              xxxx
# load or simulate data                        ####
#                                              xxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


#__________________________________________________
## load example data
dd <- read.csv("dataset_seed186.csv")
# or
# dd <- read.csv("dataset_seed152.csv")

# NOW, jump to line 285, section "fit models"
# ALTERNATIVELY, use code below to simulate data


#__________________________________________________
#__________________________________________________
## simulate data (and modify the parameters)

N = 300              # N hunting trips
K = 7                # K patch types
J = 25               # J harvester
Y_pc <- rep(NA, N)   # Y response patch choice
Y_hs <- rep(NA, N)   # Y response harvest success


#__________________________________________________
## create harvest trips
dd <- data.frame( trip_id = c(1:N) ,
                  season_id = rbinom(N, 1, 0.4) + 1) # 1-snow/ice, 2-ice-free

# rename season, reflecting how data is most likely entered
dd$season <- ifelse( dd$season_id == 1, "ice", "ice-free")

# create hunt group size + standardise
dd$Nhunt <- sample(c(1:10), size=N, replace=TRUE, prob = (dgeom(c(1:10), prob=0.3)) )
dd$Nhunt_s <- (dd$Nhunt - mean(dd$Nhunt)) / sd(dd$Nhunt)



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
## send people harvesting, i.e. draw focal IDs
dd$j_id <- c(sample( c(1:J), size=N, replace=TRUE, prob=rep(1/J, J) ) ) 
dd$age_cat <- ppl$age_cat[match(dd$j_id, ppl$j_id)]
dd$gender <- ppl$gender[match(dd$j_id, ppl$j_id)]
dd$gender_id = ifelse( dd$gender=="m", 1 , 2 ) # 1-male, 2-female
dd$income_s <- ppl$income_s[match(dd$j_id, ppl$j_id)]
dd$indegree_s <- ppl$indegree_s[match(dd$j_id, ppl$j_id)]
dd$outdegree_s <- ppl$outdegree_s[match(dd$j_id, ppl$j_id)]



#__________________________________________________
## set patch choice probability ("attraction weight") intercepts
# store in an array [patch x season x gender]
# i.e., use different schedules by season and gender
# (instead of using a vector of length K, and adding gender and season as effects)
# note, some patches are winter-only, others are summer-only, some can be chosen year-round
# the centered probability distribution for patch choice for season i and gender j is
# softmax(WA[,i,j])
WA <- array(NA, dim = c(K, 2, 2))
WA[ ,1,1] <- c(-200, -200, -200, 1,      1,    1.2, 0)    # season1, gender1
WA[ ,1,2] <- c(-200, -200, -200, 1,      1,    0.2, 0)    # season1, gender2
WA[ ,2,1] <- c(   1,    1,    1, 1.2, -200, -200,   0)    # season2, gender1
WA[ ,2,2] <- c(   1,    1,    1, 0.2, -200, -200,   0)    # season2, gender2


## set patch success probability intercepts (as log-odds)
WS <- array(NA, dim = c(K, 2, 2))
WS[ ,1,1] <- logit(c(0,    0,    0,    0.1 , 0.6 , 0.8 , 0.2))    # season1, gender1
WS[ ,1,2] <- logit(c(0,    0,    0,    0.6 , 0.7 , 0.8 , 0.25))   # season1, gender2
WS[ ,2,1] <- logit(c(0.3 , 0.6 , 0.5 , 0.2 , 0 ,   0   , 0.8))    # season2, gender1
WS[ ,2,2] <- logit(c(0.6 , 0.7 , 0.6 , 0.1 , 0 ,   0   , 0.3))    # season2, gender2


## set other effects
# for each patch category, we set a different effects
# bIn: income, bDeI/bDeO: in-/outdegree, N hunters: group size
# here, we use observed effect sizes
bIn  <- c( 0.07, -0.04,  0.45, 0.23, -0.02, -0.49, 0.01)
bDeI <- c(-0.3 , -0.1 ,  0.2 , 0.1 , -0.7 , -0.6 , 0   )
bDeO <- c(-0.1 , -0.1 , -0.2 , 0.2 ,  0.2 ,  0.4 , 0   )
bNh  <- c( 0.6 ,  0.4 , -0.1 , 0.1 , -0.1 , -0.6 , 0   )
# alternatively, use:
# bIn  <- round(rnorm(K, 0.2, 0.5), 1)
# bDeI <- round(rnorm(K, 0  , 0.3), 1)
# bDeO <- round(rnorm(K, 0  , 0.3), 1)
# bNh  <- round(rnorm(K, 0.2, 0.3), 1)



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
        fAge[ ,dd$age_cat[n]] +  # alternatively can use a linear age model here
        bIn*dd$income_s[n] +
        bDeI*dd$indegree_s[n] +
        bDeO*dd$outdegree_s[n] +
        bNh*dd$Nhunt_s[n]
  Y_pc[n] <- sample(1:K, 1, prob = softmax(A))
}

dd$patch_id <- Y_pc

#__________________________________________________
## Y_hs: simulate harvest success data

Y_hs <- rep(NA, N)

for (n in 1:N) {
  S <- WS[ , dd$season_id[n], dd$gender_id[n]] +
        fAge[ , dd$age_cat[n]] +  # alternatively can use a linear age model here
        bIn*dd$income_s[n] +
        bDeI*dd$indegree_s[n] +
        bDeO*dd$outdegree_s[n] + 
        bNh*dd$Nhunt_s[n]
  Y_hs[n] <- rbern(1, logistic(S[dd$patch_id[n]]))
}

dd$harvest <- Y_hs



#__________________________________________________
## missing values

# realistically, we may not have all information available for everyone
# randomly choose some individuals with incomplete data (here: income)

dd$income_s_complete <- dd$income_s

j_id_incomplete <- sample(J, size=sample(2,1))
ppl$income_s[ppl$j_id %in% j_id_incomplete] <- NA
dd$income_s[dd$j_id %in% ppl$j_id[is.na(ppl$income_s)] ] <- NA

# for missing value imputation in Stan, we want to store
# missN - count of missing values, and
# missDEX - index of missing values
names(which(colSums(is.na(dd)) > 0))
missn_income <- sum(is.na(dd$income_s))
missdex_income <- which(is.na(dd$income_s))

# replace NAs with 999 (or any other arbitraty number)
dd <- na2dummy(dd)


#__________________________________________________
# optional: save data
# write.csv(dd, file = paste0("dataset_seed", seed, ".csv"))



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx####
#                                              xxxx
# fit models                                   ####
#                                              xxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

## set number of chains and iterations (values used in the analysis presented in the MS: ch=4, it=6000)
ch=4
it=2000


#__________________________________________________
# prepare data for Stan

# (1) Stan needs integer IDs as input instead of character strings
# e.g., if patch data were entered as "winter marine" etc, or hunters j_id as "ID01" etc
# dd$patch_id <- as.integer(as.factor( dd$patch_cat )) 
# dd$hunter_id <- as.integer(as.factor( dd$j_id ))
# similarly, we use dd$season_id (1 for snow/ice, 2 for ice-free) and dd$gender_id (1 for male, 2 for female)

# (2) Stan needs data in list format
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

# show Stan model code in console
writeLines(readLines("m_pc.stan"))

# run Stan model
mm_pc <- stan(file="m_pc.stan", data=datlist, chains=ch, iter=it, warmup=200)

## explore sampling behavior and assess mixing across Markov chains
traceplot(mm_pc, pars="WA")
traceplot(mm_pc, pars="bIn")
traceplot(mm_pc, pars="fAge")
traceplot(mm_pc, pars=c("bDeI", "bDeO"))
traceplot(mm_pc, pars="bNh")


## visually compare estimated and true (simulated) probabilities of each choice
post_pc <- extract.samples(mm_pc)
pc_true <- pc_est <- pc_lb <- pc_ub <- matrix(NA, nrow = N, ncol = K)

for (n in 1:N) {  
  # calculate true prob of patch choice
  pc_true[n,] <- softmax( WA[ , dd$season_id[n], dd$gender_id[n]] + # intercept
                           fAge[ , dd$age_cat[n]] + 
                           bIn *dd$income_s_complete[n] +
                           bDeI*dd$indegree_s[n] +
                           bDeO*dd$outdegree_s[n] +
                           bNh *dd$Nhunt_s[n] )
  
  # calculate estimate of posterior prob of patch choice
  A <- post_pc$WA[ , , dd$season_id[n], dd$gender_id[n]] + # intercept
        post_pc$fAge[ , , dd$age_cat[n]] +
        post_pc$bIn *post_pc$income_merge[,n] +
        post_pc$bDeI*dd$indegree_s[n] +
        post_pc$bDeO*dd$outdegree_s[n] +
        post_pc$bNh *dd$Nhunt_s[n]
  pc_n <- t(apply(A, 1, softmax))
  pc_est[n, ] <- apply(pc_n, 2, mean)
  pc_lb[n, ] <- apply(pc_n, 2, rethinking::HPDI)[1,]
  pc_ub[n, ] <- apply(pc_n, 2, rethinking::HPDI)[2,]
  if (n %% 100 == 0) print(n)
}


## inspect performance
plot(pc_true, pc_est, xlim = c(0, 0.8), ylim = c(0, 0.8),
  xlab = "true probability of choice", ylab = "estimated probability of choice",
  main = "patch choice model, simulated parameter recovery")
for (n in 1:N) {
  for (k in 1:K) {
    lines(c(pc_true[n,k], pc_true[n,k]), c(pc_lb[n,k], pc_ub[n,k]), col = col.alpha("dodgerblue", 0.2))
  }
}
abline(0, 1, lty = 2)


#__________________________________________________
#                                              xxxx
#* harvest success                             ####
#                                              xxxx
#__________________________________________________

# show Stan model code in console
writeLines(readLines("m_hs.stan"))

# run Stan model
mm_hs <- stan(file="m_hs.stan", data=datlist, chains=ch, iter=it, warmup=200)

## explore sampling behavior and assess mixing across Markov chains
traceplot(mm_hs, pars="bIn")
traceplot(mm_hs, pars="fAge")
traceplot(mm_hs, pars="bDeI")


# visually compare estimated and true (simulated) probabilities of harvest success
post_hs <- extract.samples(mm_hs)
hs_true <- hs_est <- hs_lb <- hs_ub <- rep(NA, length.out = N)

for (n in 1:N) {

  # use true patch attraction intercepts to calculate true prob of harvest success
  hs_true[n] <- logistic ( WS[dd$patch_id[n], dd$season_id[n], dd$gender_id[n]] +
                            fAge[dd$patch_id[n], dd$age_cat[n]] +
                            bIn[dd$patch_id[n]] * dd$income_s_complete[n] + 
                            bDeI[dd$patch_id[n]] * dd$indegree_s[n] + 
                            bDeO[dd$patch_id[n]] * dd$outdegree_s[n] +
                            bNh[dd$patch_id[n]] * dd$Nhunt_s[n] )
  
  # use posterior attraction weights to calculate estimate of posterior prob of harvest success
  S <- post_hs$WS[ , dd$patch_id[n], dd$season_id[n], dd$gender_id[n]] +
        post_hs$fAge[ , dd$patch_id[n], dd$age_cat[n]] + 
        post_hs$bIn[ , dd$patch_id[n]] * post_hs$income_merge[ , n] + 
        post_hs$bDeI[ , dd$patch_id[n]] * dd$indegree_s[n] +
        post_hs$bDeO[ , dd$patch_id[n]] * dd$outdegree_s[n] + 
        post_hs$bNh[ , dd$patch_id[n]] * dd$Nhunt_s[n]
  pr_hs_n <- logistic(S)
  hs_est[n] <- mean(pr_hs_n)
  hs_lb[n] <- rethinking::HPDI(pr_hs_n)[1]
  hs_ub[n] <- rethinking::HPDI(pr_hs_n)[2]
  if (n %% 100 == 0) print(n)
}


## inspect performance
plot(hs_true, hs_est, xlim = c(0, 1), ylim = c(0, 1),
  xlab = "true probability of success", ylab = "estimated probability of success",
  main = "harvest success model, simulated parameter recovery")
for (n in 1:N) lines(c(hs_true[n], hs_true[n]), c(hs_lb[n], hs_ub[n]), col = col.alpha("blue", 0.2))
abline(0, 1, lty = 2)





#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx####
#                                              xxxx
# result plots and tables                      ####
#                                              xxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

cols <- rainbow(K)

# colours used in MS are linked to season, marine/terrestrial (darker/paler tone)
cols <- c("lavender blush 3" , "yellow green", "gold2", "powderblue", "goldenrod3", "steelblue", "lightslategrey")

#__________________________________________________
#                                              xxxx
#* Fig.2: predictive plots                     ####
#                                              xxxx
#__________________________________________________

# example: income
xlab = "income"

## set steps of x-axis
# either used use observed range of income or extrapolate
min(dd$income_s_complete)
max(dd$income_s_complete)
x <- seq(-2.5, 2.5, by = 0.05) 

# create a fictive hunter profile
fix_agecat <- 2
fix_indegree <- ppl$indegree_s[ppl$indegree==2][1] # use indegree_s that refers to indegree of 2
fix_outdegree <- ppl$outdegree_s[ppl$outdegree==2][1] # use outdegree_s that refers to outdegree of 2
fix_Nhunt <- sort(unique(dd$Nhunt_s))[1] # Nhunt ranged from 1:10, so this gives us Nhunt_s that refers to a group size of 1


# posterior probabilities of each choice and success (estimated mean, lower and upper bound)
pr_pc_est <- pr_pc_lb <- pr_pc_ub <- array(NA, dim=c(length(x), 7, 2, 2))
pr_hs_est <- pr_hs_lb <- pr_hs_ub <- array(NA, dim=c(length(x), 7, 2, 2))

for (g in 1:2) {
  for (s in 1:2) { 
    for (i in 1:length(x)) {
      A <- post_pc$WA[ ,  , s, g] + 
            post_pc$fAge[ , , fix_agecat] +
            post_pc$bIn * x[i] + # loop through the range of income values set above
            post_pc$bDeI * fix_indegree +
            post_pc$bDeO * fix_outdegree +
            post_pc$bNh * fix_Nhunt
      pr_pc_n <- t(apply(A, 1, softmax)) # transform to make prob of all patches sum to 1
      pr_pc_est[i, ,s,g] <- apply(pr_pc_n, 2, mean)
      pr_pc_lb[i, ,s,g] <- apply(pr_pc_n, 2, HPDI)[1,]
      pr_pc_ub[i, ,s,g] <- apply(pr_pc_n, 2, HPDI)[2,]
      
      S <- post_hs$WS[ ,  , s, g] + 
            post_hs$fAge[ , , fix_agecat] +
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


# plot (and save the 4-pages output as a single PDF)
pdf(file="figure2.pdf", width=5, height=8)
for (g in 1:2) {
  for (s in 1:2) {
    
    title_s <- dd$season[dd$season_id == s][1]
    title_g <- dd$gender[dd$gender_id == g][1]
    
    par(mfrow=c(K, 2), mar=c(2,4,1,4))
    
    for (k in 1:K) { 
      plot(1, 1, type = "n", xlim=c(min(x), max(x)), ylim=c(0,1), las=1,
           xlab = xlab, ylab = "pr(choose)", cex.axis=0.9, xaxt="n",
           main = paste( title_s, title_g, sep=", "), cex.main=0.6)
      axis(1, at = c(-1.5, 0, 1.5))
      points(x, pr_pc_est[ ,k,s,g], type = "l", col=cols[k])
      polygon(c(x, rev(x)), c(pr_pc_lb[ ,k,s,g], rev(pr_pc_ub[ ,k,s,g])), col=col.alpha(cols[k], 0.2), border=NA)
      
      plot(1, 1, type = "n", xlim=c(min(x), max(x)), ylim=c(0,1), las=1,
           xlab = xlab, ylab = "pr(success)", cex.axis=0.9, xaxt="n",
           main = paste( title_s, title_g, sep=", "), cex.main=0.6)
      axis(1, at = c(-1.5, 0, 1.5))
      points(x, pr_hs_est[ ,k,s,g], type = "l", col=cols[k])
      polygon(c(x, rev(x)), c(pr_hs_lb[ ,k,s,g], rev(pr_hs_ub[ ,k,s,g])), col=col.alpha(cols[k], 0.2), border=NA)
    }
  }
}
dev.off() # produce pdf




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
fix_Nhunt <- sort(unique(dd$Nhunt_s))[1] # Nhunt ranged from 1:10, so this gives us Nhunt_s that refers to a group size of 1
  
# posterior probabilities of each choice and success (estimated mean, lower and upper bound)
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

# plot (and save the 4-panel output as PDF)
pdf(file="figure3.pdf", width=5, height=5)
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
dev.off() # produce pdf
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 4.1))



#__________________________________________________
#                                              xxxx
#* supplementary Figure: effect sizes         ####
#                                              xxxx
#__________________________________________________

# alternatively, use rstan::stan_plot (rstan plots behave like ggplot objects)
plot(mm_pc, pars="bIn", 
     show_density = FALSE,
     point_est = "mean",
     ci_level = 0.89,
     outer_level = 0.89,
     fill_color = cols, 
     outline_color = cols,
     est_color = cols)



#__________________________________________________
#                                              xxxx
#* supplementary Table: effect sizes          ####
#                                              xxxx
#__________________________________________________

tt <- as.data.frame(rethinking::precis(mm_pc, prob=0.89, digits=2, depth=3,
                                       pars=c("WA", "bIn", "fAge", "bDeI", "bDeO", "bNh")))

# alternatively, use
print(mm_pc)

tt <- as.data.frame(summary(mm_pc, 
                            pars = c("WA", "bIn", "fAge", "bDeI", "bDeO", "bNh"),
                            probs = c(0.055, 0.945))$summary)


# reformat table
pars <- row.names(tt)
tt <- as.data.frame(lapply(tt, round, digits = 3))
tt$pars <- pars

tt$confound <- sapply(strsplit(tt$pars, "[", fixed=TRUE), getElement, 1)
tt$confound[tt$confound %in% "WA"]  <- "intercept"   # WA[k,s,g]
tt$confound[tt$confound %in% "bIn"] <- "income"      # bIn[k]
tt$confound[tt$confound %in% "fAge"]<- "age"         # fAge[k,agecat]
tt$confound[tt$confound %in% "bDeI"]<- "in-degree"
tt$confound[tt$confound %in% "bDeO"]<- "out-degree"
tt$confound[tt$confound %in% "bNh"] <- "N hunters"

tt$patch <- sapply(strsplit(tt$pars, "[", fixed=TRUE), getElement, 2)
tt$patch <- lapply(tt$patch, substr, 1,1)
tt$patch_cat[tt$patch %in% "1"]<- "incidental"
tt$patch_cat[tt$patch %in% "2"]<- "inland spring"
tt$patch_cat[tt$patch %in% "3"]<- "inland summer"
tt$patch_cat[tt$patch %in% "4"]<- "inland winter"
tt$patch_cat[tt$patch %in% "5"]<- "marine summer"
tt$patch_cat[tt$patch %in% "6"]<- "marine winter"
tt$patch_cat[tt$patch %in% "7"]<- "tidal"

tt <- tt[ , c("pars", "confound", "patch_cat", "mean", "sd", "X5.5.", "X94.5.")]
colnames(tt) <- c("pars", "Counfound", "Patch Category", "Mean", "SD", "5.5%", "94.5%")
tt$Counfound[tt$Counfound=="intercept"] <- paste("intercept", rep(c("wi, m", "wi, f", "su, m", "su, f"), 7))
tt$Counfound[tt$Counfound=="age"] <- rep(c("age <30", "age 30-40", "age 40-50", "age 50+"), 7)
tt$pars <- NULL


# export table for LaTeX
print(xtable(tt, type = "latex"), include.rownames=FALSE, file = "TableS1_effectsizes.tex")

