# Transdimensional MCMC Bayesian estimation and comparison of nested logistic choice models

rm(list = ls())

# Load necessary libraries
library(dplyr)
library(R2jags)

# For reproducibility purposes
seed <- 123
set.seed(seed)

# Load dataset
setwd("DATA_LOCATION")
data <- read.csv("DATA_FILENAME.csv")

# Get rid of missing data
data <- data %>% filter(SUBJECT != "NaN", !is.na(RT))

# Transform CONTEXT (condition) into factor
data$CONTEXT <- factor(data$CONTEXT, levels=c(0,1), labels=c("other","self"))

# Extract condition-specific data: Choose Self and Choose Other
data.self  <- data %>% filter(CONTEXT == "self")
data.other <- data %>% filter(CONTEXT == "other")

# Number of observations per condition
Nobs.self  <- length(data.self$ACCEPT)
Nobs.other <- length(data.other$ACCEPT)

# Number of subjects per condition
Nsubj.self  <- max(data.self$SUBJECT)
Nsubj.other <- max(data.other$SUBJECT)

# Extract and rename variables of interest for convenience
# Subject ID for each observation
subj.self  <- data.self$SUBJECT
subj.other <- data.other$SUBJECT

# Task-relevant Rating for each observation
rating.self.rel  <- scale(data.self$RATING_SELF_AVG,   center = TRUE, scale = FALSE)[,1]
rating.other.rel <- scale(data.other$RATING_OTHER_AVG, center = TRUE, scale = FALSE)[,1]

# Task-irrelevant Rating for each observation
rating.self.irr  <- scale(data.self$RATING_OTHER_AVG, center = TRUE, scale = FALSE)[,1]
rating.other.irr <- scale(data.other$RATING_SELF_AVG, center = TRUE, scale = FALSE)[,1]

# Task-relevant Probability for each observation
prob.self.rel  <- data.self$P_SELF
prob.other.rel <- data.other$P_OTHER

# Task-irrelevant Probability for each observation
prob.self.irr  <- data.self$P_OTHER
prob.other.irr <- data.other$P_SELF

# Choice for each observation
choice.self  <- data.self$ACCEPT
choice.other <- data.other$ACCEPT

# Choice of prior for hierarchical model probabilites (Dirichlet)
dp <- c(1,1,1,1)

# Specify input data to JAGS model
model.data.self  <- c("Nobs.self", "Nsubj.self", "subj.self", "choice.self", "rating.self.rel", "rating.self.irr", "prob.self.rel", "prob.self.irr", "dp")
model.data.other <- c("Nobs.other", "Nsubj.other", "subj.other", "choice.other", "rating.other.rel", "rating.other.irr", "prob.other.rel", "prob.other.irr", "dp")

# Specify JAGS parameters
n.chains <- 4        # Number of chains
n.iter   <- 250000     # Number of samples
n.burnin <- 50000     # Number of burn-in samples
n.thin   <- 5        # Thinning parameter

# Initialize MCMC chains (same for both conditions)
inits1 <- list("alphagmean"=0, "beta1gmean"=0, "beta2gmean"=0, "beta3gmean"=0,
               "alphagprec"=0.001, "beta1gprec"=0.001, "beta2gprec"=0.001, "beta3gprec"=0.001,
               "pmod"=c(0.9, 1/30, 1/30, 1/30))
inits2 <- list("alphagmean"=0, "beta1gmean"=0, "beta2gmean"=0, "beta3gmean"=0,
               "alphagprec"=0.001, "beta1gprec"=0.001, "beta2gprec"=0.001, "beta3gprec"=0.001,
               "pmod"=c(1/30, 0.9, 1/30, 1/30))
inits3 <- list("alphagmean"=0, "beta1gmean"=0, "beta2gmean"=0, "beta3gmean"=0,
               "alphagprec"=0.001, "beta1gprec"=0.001, "beta2gprec"=0.001, "beta3gprec"=0.001,
               "pmod"=c(1/30, 1/30, 0.9, 1/30))
inits4 <- list("alphagmean"=0, "beta1gmean"=0, "beta2gmean"=0, "beta3gmean"=0,
               "alphagprec"=0.001, "beta1gprec"=0.001, "beta2gprec"=0.001, "beta3gprec"=0.001,
               "pmod"=c(1/30, 1/30, 1/30, 0.9))
model.inits <- list(inits1, inits2, inits3, inits4)

# Specify which parameters to monitor and store
save.all = FALSE
model.params.prob <- c("pmod")
model.params.all  <- c("alpha", "beta1", "beta2", "beta3",
                       "alphagmean", "beta1gmean", "beta2gmean", "beta3gmean",
                       "alphagprec", "beta1gprec", "beta2gprec", "beta3gprec",
                       "mod", "pmod")
if (save.all) {
  model.params <- model.params.all
} else {
  model.params <- model.params.prob
  }

# Model specification
# CHOOSE SELF condition
nested.logistic.model.self <- function() {

    # Model observations
    for (obs in 1:Nobs.self) {
      
      # Binary choice data (accept = 1 / reject = 0) ~ Bernoulli
      choice.self[obs] ~ dbern(p.accept[obs])
      
      # Define p.accept for nested models using model indices
      logit(p.accept[obs]) <- alpha[subj.self[obs]] +
                              in.mod.dv[obs] * beta1[subj.self[obs]] * rating.self.rel[obs] * prob.self.rel[obs] +
                              in.mod.dv.irr[obs] * beta2[subj.self[obs]] * rating.self.irr[obs] * prob.self.irr[obs] +
                              in.mod.abs.diff[obs] * beta3[subj.self[obs]] * abs((rating.self.rel[obs] * prob.self.rel[obs]) - (rating.self.irr[obs] * prob.self.irr[obs]))
      
      mod1[obs] <- mod[subj.self[obs]] == 1
      mod2[obs] <- mod[subj.self[obs]] == 2
      mod3[obs] <- mod[subj.self[obs]] == 3
      mod4[obs] <- mod[subj.self[obs]] == 4
      
      in.mod.dv[obs]       <- mod2[obs] + mod3[obs] + mod4[obs]
      in.mod.dv.irr[obs]   <- mod3[obs] + mod4[obs]
      in.mod.abs.diff[obs] <- mod4[obs]
      
    }
    
    # Specify prior for lower-level parameters (individual subjects)
    for (subj in 1:Nsubj.self) {
      
      alpha[subj] ~ dnorm(alphagmean, alphagprec)
      beta1[subj] ~ dnorm(beta1gmean, beta1gprec)
      beta2[subj] ~ dnorm(beta2gmean, beta2gprec)
      beta3[subj] ~ dnorm(beta3gmean, beta3gprec)
      
      mod[subj]  ~ dcat(pmod[1:4])
      
    }
    
    # Specify hierarchical priors (population-level parameters)
    # Means
    alphagmean ~ dnorm(0, 0.001)
    beta1gmean ~ dnorm(0, 0.001)
    beta2gmean ~ dnorm(0, 0.001)
    beta3gmean ~ dnorm(0, 0.001)
    # Precisions
    alphagprec ~ dgamma(.1, .1)
    beta1gprec ~ dgamma(.1, .1)
    beta2gprec ~ dgamma(.1, .1)
    beta3gprec ~ dgamma(.1, .1)
    # Model probabilities
    pmod[1:4] ~ ddirch(dp[1:4])
  
}

# CHOOSE OTHER condition
nested.logistic.model.other <- function() {
  
  # Model observations
  for (obs in 1:Nobs.other) {
    
    # Binary choice data (accept = 1 / reject = 0) ~ Bernoulli
    choice.other[obs] ~ dbern(p.accept[obs])
    
    # Define p.accept for nested models using model indices
    logit(p.accept[obs]) <- alpha[subj.other[obs]] +
                            in.mod.dv[obs] * beta1[subj.other[obs]] * rating.other.rel[obs] * prob.other.rel[obs] +
                            in.mod.dv.irr[obs] * beta2[subj.other[obs]] * rating.other.irr[obs] * prob.other.irr[obs] +
                            in.mod.abs.diff[obs] * beta3[subj.other[obs]] * abs((rating.other.rel[obs] * prob.other.rel[obs]) - (rating.other.irr[obs] * prob.other.irr[obs]))
    
    mod1[obs] <- mod[subj.other[obs]] == 1
    mod2[obs] <- mod[subj.other[obs]] == 2
    mod3[obs] <- mod[subj.other[obs]] == 3
    mod4[obs] <- mod[subj.other[obs]] == 4
    
    in.mod.dv[obs]       <- mod2[obs] + mod3[obs] + mod4[obs]
    in.mod.dv.irr[obs]   <- mod3[obs] + mod4[obs]
    in.mod.abs.diff[obs] <- mod4[obs]
    
  }
  
  # Specify prior for lower-level parameters (individual subjects)
  for (subj in 1:Nsubj.other) {
    
    alpha[subj] ~ dnorm(alphagmean, alphagprec)
    beta1[subj] ~ dnorm(beta1gmean, beta1gprec)
    beta2[subj] ~ dnorm(beta2gmean, beta2gprec)
    beta3[subj] ~ dnorm(beta3gmean, beta3gprec)
    
    mod[subj]  ~ dcat(pmod[1:4])
    
  }
  
  # Specify hierarchical priors (population-level parameters)
  # Means
  alphagmean ~ dnorm(0, 0.001)
  beta1gmean ~ dnorm(0, 0.001)
  beta2gmean ~ dnorm(0, 0.001)
  beta3gmean ~ dnorm(0, 0.001)
  # Precisions
  alphagprec ~ dgamma(.1, .1)
  beta1gprec ~ dgamma(.1, .1)
  beta2gprec ~ dgamma(.1, .1)
  beta3gprec ~ dgamma(.1, .1)
  # Model probabilities
  pmod[1:4] ~ ddirch(dp[1:4])
  
}

# Run JAGS and update until convergence
model.self  <- jags(data = model.data.self, inits = model.inits, parameters.to.save = model.params,
                    n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin,
                    model.file = nested.logistic.model.self, jags.seed = seed)
     
model.other <- jags(data = model.data.other, inits = model.inits, parameters.to.save = model.params,
                    n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin,
                    model.file = nested.logistic.model.other, jags.seed = seed)     

model.self.upd  <- autojags(model.self,  n.iter = 250000, n.thin = 10, Rhat = 1.01, n.update = 2, progress.bar = "text")
model.other.upd <- autojags(model.other, n.iter = 250000, n.thin = 10, Rhat = 1.01, n.update = 2, progress.bar = "text")
