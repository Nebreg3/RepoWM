library(mixtools)
library(ggplot2)
library(gridExtra)

# Load custom negative log-likelihood function
source("llh_covars.R")

# Load main data
# load("data_main.RData")
# load("Data/max_llh.RData")

# Covariates matrix (Time, Age, Gender, Interaction, Sine, Cosine)
t <- rep(seq(1, 96), 4) / 96
pr3$interaction <- as.numeric(pr3$age) * as.numeric(pr3$sexe)
covars <- cbind(t, pr3$age, pr3$sexe, pr3$interaction, sin(2 * pi * t / 3), cos(2 * pi * t / 3))

# Initial values via mixtools (Direct estimation for MLE start)
w0 <- 0.7
q0 <- 0.5
prova <- regmixEM(pr3$incid, covars, lambda=c(w0, 1 - w0), 
                  beta=matrix(c(rep(mean(pr3$incid), 2), rep(0, 12)), ncol=2), 
                  sigma=c(sd(pr3$incid), sd(pr3$incid) / q0), 
                  k=2, addintercept=TRUE, epsilon=1e-16, maxit=10000)

# Linear regression for initial covariate values
linmod <- lm(pr3$incid ~ covars[, 1:6])

# Maximum likelihood estimation for `w` and `q`
max.llh <- nlm(f=llh, p=c(log(prova$lambda[1] / (1 - prova$lambda[1])), 
                          -0.5, linmod$coefficients, prova$sigma[1], 
                          log((prova$beta[1,1] / prova$beta[1,2]) / (1 - prova$beta[1,1] / prova$beta[1,2]))), 
              data=pr3$incid, covars=covars, hessian=TRUE)
q <- exp(max.llh$estimate[11]) / (1 + exp(max.llh$estimate[11]))

# Estimate `y` with `w`, covariate weights, and transformations
y_est <- numeric(384)
w <- numeric(384)

for (i in 1:384) {
  j <- ifelse((i %% 96) / 96 == 0, 1, (i %% 96) / 96)
  w[i] <- exp(max.llh$estimate[1] + max.llh$estimate[2] * j) / 
          (1 + exp(max.llh$estimate[1] + max.llh$estimate[2] * j))
  cov_weights <- max.llh$estimate[3:10]
  m <- sum(cov_weights * c(1, j, pr3$age[i], pr3$sexe[i], pr3$interaction[i], covars[i, 5:6]))
  y_est[i] <- w[i] * m + (1 - w[i]) * m / q
}

# Function definitions (unchanged unless stated)
calc_posterior <- function(incid, cov_effects, weight, covariates) {
  mean_group1 <- cov_effects[3] + sum(cov_effects[4:9] * covariates[1:6])
  mean_group2 <- mean_group1 / q
  sd_group1 <- max.llh$estimate[10]
  sd_group2 <- sd_group1 / q
  prob_group1 <- weight * dnorm(incid, mean=mean_group1, sd=sd_group1)
  prob_group2 <- (1 - weight) * dnorm(incid, mean=mean_group2, sd=sd_group2)
  posterior_prob <- prob_group1 / (prob_group1 + prob_group2)
  return(posterior_prob)
}

posterior_probabilities <- function(data, covariates) {
  w <- numeric(nrow(data))
  post <- matrix(0, nrow=nrow(data), ncol=2)
  
  for (i in 1:nrow(data)) {
    cov_effects <- max.llh$estimate[1:10] * covariates[i, ]
    w[i] <- exp(cov_effects[1] + cov_effects[2] * i / 96) / 
            (1 + exp(cov_effects[1] + cov_effects[2] * i / 96))
    post[i, 1] <- calc_posterior(data$incid[i], cov_effects, w[i], covariates[i, ])
    post[i, 2] <- 1 - post[i, 1]
  }
  return(post)
}

# Main demographic analysis
process_demographic_group <- function(data, covariates, q_value, start_year, end_year, ylim_range, main_title) {
  post <- posterior_probabilities(data, covariates)
  xrec <- reconstruct_incid(data$incid, post, q_value)
  calc_reconstruction_error(data$incid, xrec)
  plot_time_series(data$incid, xrec, start_year, end_year, ylim_range, main_title)
}

process_demographic_group(pr3[pr3$sexe == 0 & pr3$age == 0, ], covars, q, 2009, 2016, c(9, 32), "Women 15-29 years old")
process_demographic_group(pr3[pr3$sexe == 0 & pr3$age == 1, ], covars, q, 2009, 2016, c(1, 8), "Women 30-94 years old")
process_demographic_group(pr3[pr3$sexe == 1 & pr3$age == 0, ], covars, q, 2009, 2016, c(4, 32), "Men 15-29 years old")
process_demographic_group(pr3[pr3$sexe == 1 & pr3$age == 1, ], covars, q, 2009, 2016, c(1, 12), "Men 30-94 years old")

print("Done")
