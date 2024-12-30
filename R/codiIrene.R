library(mixtools)
library(ggplot2)
library(gridExtra)

# Load custom log-likelihood function
source("R/llh_covars.R")

# Load main dataset
load("Data/data_main.RData")

# Covariates matrix (Time, Age, Gender, Interaction, Sine, Cosine)
t <- rep(seq(1, 96), 4) / 96
pr3$interaction <- as.numeric(pr3$age) * as.numeric(pr3$sexe)
covars <- cbind(t, pr3$age, pr3$sexe, pr3$interaction, sin(2 * pi * t / 3), cos(2 * pi * t / 3))

# Initial values for mixtools
w0 <- 0.7
q0 <- 0.5
prova <- regmixEM(pr3$incid, covars, lambda=c(w0, 1 - w0), 
                  beta=matrix(c(rep(mean(pr3$incid), 2), rep(0, 12)), ncol=2), 
                  sigma=c(sd(pr3$incid), sd(pr3$incid) / q0), 
                  k=2, addintercept=TRUE, epsilon=1e-16, maxit=10000)

# Linear regression for initial values
linmod <- lm(pr3$incid ~ covars[, 1:6])

# Maximum likelihood estimation for `w` and `q`
max.llh <- nlm(f=llh, p=c(log(prova$lambda[1] / (1 - prova$lambda[1])), 
                          -0.5, linmod$coefficients, prova$sigma[1], 
                          log((prova$beta[1, 1] / prova$beta[1, 2]) / (1 - prova$beta[1, 1] / prova$beta[1, 2]))), 
              data=pr3$incid, covars=covars, hessian=TRUE)

# Calculate q as a function of time
c <- 0  # Adjust this value if needed
q_t <- numeric(384)
for (i in 1:384) {
  q_t[i] <- exp(c + max.llh$estimate[11] * t[i]) / (1 + exp(c + max.llh$estimate[11] * t[i]))
}

# Function to calculate posterior probabilities
posterior_probabilities <- function(data, covariates, max_llh, q_t) {
  n <- nrow(data)  # Number of rows in the subset
  post <- matrix(nrow=n, ncol=2)
  w <- numeric(n)
  
  # Ensure covariates have the expected size
  if (n > nrow(covariates)) {
    stop("The number of rows in 'data' exceeds that in 'covariates'. Please check your data.")
  }
  
  for (i in 1:n) {
    # Calculate w
    w[i] <- exp(max.llh$estimate[1] + max.llh$estimate[2] * t[i]) /
            (1 + exp(max.llh$estimate[1] + max.llh$estimate[2] * t[i]))
    
    # Safeguard against out-of-bounds column access
    cov_subset <- covariates[i, 1:6]  # Adjust if covariates have fewer columns
    mean_val <- max.llh$estimate[3] + sum(max.llh$estimate[4:9] * cov_subset)
    
    # Verify if q_t[i] is valid
    if (i > length(q_t)) {
      stop("Index 'i' exceeds the length of 'q_t'.")
    }
    
    sd_val <- max.llh$estimate[10]
    post[i, 1] <- w[i] * dnorm(data$incid[i], mean=mean_val, sd=sd_val)
    post[i, 2] <- (1 - w[i]) * dnorm(data$incid[i], mean=mean_val / q_t[i], sd=sd_val / q_t[i])
  }
  
  # Normalize the posterior probabilities
  post <- post / rowSums(post)
  return(post)
}

# Function to reconstruct incidence
reconstruct_incid <- function(incid, post, q_t) {
  xrec <- numeric(length(incid))
  for (i in 1:length(incid)) {
    xrec[i] <- ifelse(post[i, 2] > 0.5, incid[i], incid[i] / q_t[i])
  }
  return(xrec)
}

# Function to plot time series
plot_time_series <- function(incid, xrec, start_year, end_year, ylim_range, main_title) {
  ts.plot(ts(incid, start=c(start_year, 1), end=c(end_year, 12), freq=12), ylim=ylim_range, ylab="Incidence x 100,000", main=main_title)
  lines(seq(start_year, end_year + 0.99, 1/12), xrec, col="red", lty=2)
  legend("topright", legend=c("Observed", "Reconstructed"), col=c("black", "red"), lty=c(1, 2))
}

# Process demographic group
process_demographic_group <- function(data, covariates, q_t, start_year, end_year, ylim_range, main_title) {
  post <- posterior_probabilities(data, covariates, max_llh, q_t)
  xrec <- reconstruct_incid(data$incid, post, q_t)
  plot_time_series(data$incid, xrec, start_year, end_year, ylim_range, main_title)
}

# Analyze by demographic group
process_demographic_group(pr3[pr3$sexe == 0 & pr3$age == 0, ], covars, q_t, 2009, 2016, c(9, 32), "Women 15-29 years old")
process_demographic_group(pr3[pr3$sexe == 0 & pr3$age == 1, ], covars, q_t, 2009, 2016, c(1, 8), "Women 30-94 years old")
process_demographic_group(pr3[pr3$sexe == 1 & pr3$age == 0, ], covars, q_t, 2009, 2016, c(4, 32), "Men 15-29 years old")
process_demographic_group(pr3[pr3$sexe == 1 & pr3$age == 1, ], covars, q_t, 2009, 2016, c(1, 12), "Men 30-94 years old")

print("Done")
