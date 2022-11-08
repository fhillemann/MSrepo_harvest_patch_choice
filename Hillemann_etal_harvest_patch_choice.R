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


options(warnPartialMatchDollar=TRUE) 
# squawks if you refer to a variable in a data frame that only partially matches 

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

set.seed(1)

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
# to draw IDs, dgeom works well to reflect that few people go out a lot, most only occasionally 
# plot( dgeom( c(1:J), prob=0.2) )
j_once <- c(1:J) # making sure each hunter shows up at least once
dd$j_id <- sample(c( j_once , sample( c(1:J), size=N-length(j_once), replace=TRUE, prob=dgeom(c(1:J), 0.2))))
dd$age_cat <- ppl$age_cat[match(dd$j_id, ppl$j_id)]
dd$gender <- ppl$gender[match(dd$j_id, ppl$j_id)]
dd$income_s <- ppl$income_s[match(dd$j_id, ppl$j_id)]
dd$indegree_s <- ppl$indegree_s[match(dd$j_id, ppl$j_id)]
dd$outdegree_s <- ppl$outdegree_s[match(dd$j_id, ppl$j_id)]

# create hunt group size + standardise
dd$Nhunt <- sample(c(1:10), size=N, replace=TRUE, prob = (dgeom(c(1:10), prob=0.3)) )
dd$Nhunt_s <- (dd$Nhunt - mean(dd$Nhunt)) / sd(dd$Nhunt)


#__________________________________________________
# prepare data for Stan
# Stan needs integer IDs as input (1:x instead of character strings)
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
# the centered probability distribution for patch choice for season i and gender j is
# softmax(WA[,i,j])
WA <- array(NA, dim = c(K, 2, 2))
WA[ ,1,1] <- c(-200, -200, -200,  1,  1, 1.2, 0)    # season1, gender1
WA[ ,1,2] <- c(-200, -200, -200,  1,  1, 0.2, 0)    # season1, gender2
WA[ ,2,1] <- c( 1,  1,  1, 1.2, -200, -200, 0)    # season2, gender1
WA[ ,2,2] <- c( 1,  1,  1, 0.2, -200, -200, 0)    # season2, gender2


## set patch success probability intercepts (as log-odds)
WS <- array(NA, dim = c(K, 2, 2))
WS[ ,1,1] <- logit(c(0, 0, 0, 0.1 , 0.6 , 0.8, 0.2))    # season1, gender1
WS[ ,1,2] <- logit(c(0, 0, 0, 0.6 , 0.7 , 0.8, 0.25))    # season1, gender2
WS[ ,2,1] <- logit(c(0.3 , 0.6 , 0.5 , 0.2, 0, 0, 0.8))    # season2, gender1
WS[ ,2,2] <- logit(c(0.6 , 0.7 , 0.6 , 0.1, 0, 0, 0.3))    # season2, gender2


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
dd$income_s[dd$j_id %in% ppl$j_id[is.na(ppl$income)] ] <- NA

# store count of missing values (missN) and index of missing values (missDEX)
names(which(colSums(is.na(dd)) > 0))
missn_income <- sum(is.na(dd$income_s))
missdex_income <- which(is.na(dd$income_s))

# replace NAs with 999 (or any other arbitraty number)
dd <- na2dummy(dd)




#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx####
#                                              xxxx
# fit models                                   ####
#                                              xxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

## set number of chains and iterations
ch=3
it=2000


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

# run Stan model
mm_pc <- stan(file=paste0(here(),"/m_pc.stan"), data=datlist, chains=ch, iter=it, warmup=200)

## explore sampling behavior and assess mixing across Markov chains
traceplot(mm_pc, pars="WA")
traceplot(mm_pc, pars="bIn")
traceplot(mm_pc, pars="fAge")
traceplot(mm_pc, pars=c("bDeI", "bDeO"))
traceplot(mm_pc, pars="bNh")


## compare estimated and actual probabilities of each choice
post_pc <- extract.samples(mm_pc)
pc_true <- pc_est <- pc_lb <- pc_ub <- matrix(NA, nrow = N, ncol = K)

for (n in 1:N) {  
  # calculate true prob of patch choice
  pc_true[n,] <- softmax( WA[,dd$season_id[n],dd$gender_id[n]] + # intercept
                                  fAge[ ,dd$age_cat[n]] +
                                  bIn*dd$income_s[n] +
                                  bDeI*dd$indegree_s[n] +
                                  bDeO*dd$outdegree_s[n] +
                                  bNh*dd$Nhunt_s[n] )
  
  # calculate estimate of posterior prob of patch choice
  A <- post_pc$WA[ , , dd$season_id[n], dd$gender_id[n]] + # intercept
        post_pc$fAge[ , , dd$age_cat[n]] +
        post_pc$bIn*dd$income_s[n] +
        post_pc$bDeI*dd$indegree_s[n] +
        post_pc$bDeO*dd$outdegree_s[n] +
        post_pc$bNh*dd$Nhunt_s[n]
  pc_n <- t(apply(A, 1, softmax))
  pc_est[n, ] <- apply(pc_n, 2, mean)
  pc_lb[n, ] <- apply(pc_n, 2, rethinking::HPDI)[1,]
  pc_ub[n, ] <- apply(pc_n, 2, rethinking::HPDI)[2,]
  if (n %% 100 == 0) print(n)
}

plot(pc_true, pc_est, xlim = c(0, 0.8), ylim = c(0, 0.8),
  xlab = "true probability of choice", ylab = "estimated probability of choice",
  main = "patch choice model, simulated parameter recovery")
abline(0, 1)
for (n in 1:N) {
  for (k in 1:K) {
    lines(c(pc_true[n,k], pc_true[n,k]), c(pc_lb[n,k], pc_ub[n,k]), col = gray(0.5, 0.5))
  }
}


#__________________________________________________
#                                              xxxx
#* harvest success                             ####
#                                              xxxx
#__________________________________________________

# run Stan model
mm_hs <- stan(file=paste0(here(),"/m_hs.stan"), data=datlist, chains=ch, iter=it, warmup=200)

## explore sampling behavior and assess mixing across Markov chains
traceplot(mm_hs, pars="bIn")
traceplot(mm_hs, pars="fAge")
traceplot(mm_hs, pars="bDeI")


# compare estimated and actual probabilities of each choice
post_hs <- extract.samples(mm_hs)
hs_true <- hs_est <- hs_lb <- hs_ub <- rep(NA, length.out = N)


for (n in 1:N) {

  # use true patch attraction intercepts to calculate true prob of harvest success
  hs_true[n] <- logistic ( WS[dd$patch_id[n], dd$season_id[n], dd$gender_id[n]] +
                                  fAge[dd$patch_id[n], dd$age_cat[n]] +
                                  bIn[dd$patch_id[n]]*dd$income_s[n] + 
                                  bDeI[dd$patch_id[n]]*dd$indegree_s[n] + 
                                  bDeO[dd$patch_id[n]]*dd$outdegree_s[n] +
                                  bNh[dd$patch_id[n]]*dd$Nhunt_s[n] )
  
  # use posterior attraction weights to calculate estimate of posterior prob of harvest success
  S <- post_hs$WS[ , dd$patch_id[n], dd$season_id[n], dd$gender_id[n]] +
        post_hs$fAge[ , dd$patch_id[n], dd$age_cat[n]] + 
        post_hs$bIn[, dd$patch_id[n]]*dd$income_s[n] + 
        post_hs$bDeI[, dd$patch_id[n]]*dd$indegree_s[n] +
        post_hs$bDeO[, dd$patch_id[n]]*dd$outdegree_s[n] + 
        post_hs$bNh[, dd$patch_id[n]]*dd$Nhunt_s[n]
  pr_hs_n <- logistic(S)
  hs_est[n] <- mean(pr_hs_n)
  hs_lb[n] <- rethinking::HPDI(pr_hs_n)[1]
  hs_ub[n] <- rethinking::HPDI(pr_hs_n)[2]
  if (n %% 100 == 0) print(n)
}

plot(hs_true, hs_est, xlim = c(0, 1), ylim = c(0, 1),
  xlab = "true probability of success", ylab = "estimated probability of success",
  main = "harvest success model, simulated parameter recovery")
abline(0, 1, lty = 2)
for (i in 1:N) lines(c(hs_true[i], hs_true[i]), c(hs_lb[i], hs_ub[i]), col = col.alpha("dodgerblue", 0.2))


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
fix_Nhunt <- sort(unique(dd$Nhunt_s))[1] # Nhunt ranged from 1:10, so this gives us Nhunt_s that refers to a group size of 1


# posterior probabilities of each choice and success (estimated mean, lower and upper bound)
pr_pc_est <- pr_pc_lb <- pr_pc_ub <- array(NA, dim=c(length(x), 7, 2, 2))
pr_hs_est <- pr_hs_lb <- pr_hs_ub <- array(NA, dim=c(length(x), 7, 2, 2))

for (g in 1:2) {
  for (s in 1:2) { 
    for (i in 1:length(x)) {
      A <- post_pc$WA[ ,  , s, g] + 
            post_pc$fAge[,,fix_agecat] + # alternatively can use a linear age model here
            post_pc$bIn * x[i] + 
            post_pc$bDeI * fix_indegree +
            post_pc$bDeO * fix_outdegree +
            post_pc$bNh * fix_Nhunt
      pr_pc_n <- t(apply(A, 1, softmax)) # transform to make prob of all patches sum to 1
      pr_pc_est[i, ,s,g] <- apply(pr_pc_n, 2, mean)
      pr_pc_lb[i, ,s,g] <- apply(pr_pc_n, 2, HPDI)[1,]
      pr_pc_ub[i, ,s,g] <- apply(pr_pc_n, 2, HPDI)[2,]
      
      S <- post_hs$WS[ ,  , s, g] + 
            post_hs$fAge[,,fix_agecat] + # alternatively can use a linear age model here
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


# plot (and save the 4-pages output as PDF)
pdf(file=paste0(file.path(here(), "/figure2.pdf")), width=5, height=8)
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
#* supplementary Figure: effect sizes         ####
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
#* supplementary Table: effect sizes          ####
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

