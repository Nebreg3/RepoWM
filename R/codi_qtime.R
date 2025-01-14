library(mixtools)
library(ggplot2)
library(gridExtra)

# Load custom negative log-likelihood function
source("R/llh_covars_qtime.R")

# Load main data
load("R/data_main.RData")
# load("Data/max_llh.RData")

# Covariates matrix (Time, Age, Sex, Interaction, Sine, Cosine)
t <- rep(seq(1, 96), 4) / 96
pr3$interaction <- as.numeric(pr3$age) * as.numeric(pr3$sexe)
covars <- cbind(t, pr3$age, pr3$sexe, pr3$interaction, sin(2 * pi * t / 3), cos(2 * pi * t / 3))

# Initial values via mixtools (Direct estimation for MLE start)
w0 <- 0.7
q0 <- 0.5
prova <- regmixEM(pr3$incid, covars, lambda = c(w0, 1 - w0), 
          beta = matrix(c(rep(mean(pr3$incid), 2), rep(0, 12)), ncol = 2), 
          sigma = c(sd(pr3$incid), sd(pr3$incid) / q0), 
          k = 2, addintercept = TRUE, epsilon = 1e-16, maxit = 10000)

# Linear regression for initial covariate values
linmod <- lm(pr3$incid ~ covars[, 1:6])

# Maximum likelihood estimation for `w` and time-dependent `q`
# Perform nonlinear minimization of the llh function
max.llh <- nlm(
  f = llh, 
  p = c(
  log(prova$lambda[1] / (1 - prova$lambda[1])),  # Logit transformation of prova$lambda[1]
  -0.5,  # Fixed initial value for w slope
  linmod$coefficients,  # Coefficients from the linear model
  prova$sigma[1],  # First element of prova$sigma
  0,  # Initial value for q intercept (gamma_0)
  0   # Initial value for q slope (gamma_1)
  ), 
  data = pr3$incid,  # Data for the llh function
  covars = covars,  # Covariates for the llh function
  hessian = TRUE  # Request Hessian matrix
)

# Extract dynamic q as a function of time
q_time <- function(j) {
  exp(max.llh$estimate[11] + max.llh$estimate[12] * j) /
  (1 + exp(max.llh$estimate[11] + max.llh$estimate[12] * j))
}

# Estimate `y` with `w`, dynamic `q`, covariate weights, and transformations
y_est <- numeric(384)
w <- numeric(384)

for (i in 1:384) {
  j <- ifelse((i %% 96) / 96 == 0, 1, (i %% 96) / 96)
  w[i] <- exp(max.llh$estimate[1] + max.llh$estimate[2] * j) / 
      (1 + exp(max.llh$estimate[1] + max.llh$estimate[2] * j))
  q_t <- q_time(j)
  cov_weights <- max.llh$estimate[3:10]
  m <- sum(cov_weights * c(1, j, pr3$age[i], pr3$sexe[i], pr3$interaction[i], covars[i, 5:6]))
  y_est[i] <- w[i] * m + (1 - w[i]) * m / q_t
}

# Calculate posterior probabilities for a given incidence and covariate effects
#' @param incid Numeric value representing the incidence.
#' @param cov_effects Numeric vector of coefficients for covariates.
#' @param weight Numeric value representing the weight of the first group.
#' @param covariates Numeric vector of covariates for the individual.
#' @return Numeric value of the posterior probability for the first group.
calc_posterior <- function(incid, cov_effects, weight, covariates, time) {
  q_t <- q_time(time)
  print(q_t)
  mean_group1 <- cov_effects[3] + sum(cov_effects[4:9] * covariates[1:6])
  mean_group2 <- mean_group1 / q_t
  sd_group1 <- max.llh$estimate[10]
  sd_group2 <- sd_group1 / q_t
  
  prob_group1 <- weight * dnorm(incid, mean = mean_group1, sd = sd_group1)
  prob_group2 <- (1 - weight) * dnorm(incid, mean = mean_group2, sd = sd_group2)
  
  posterior_prob <- prob_group1 / (prob_group1 + prob_group2)
  return(posterior_prob)
}

# Calculate posterior probabilities for all data points in a dataset
#' @param data Data frame containing the incidence values.
#' @param covariates Matrix of covariates for the dataset.
#' @param max_llh A list containing the maximum likelihood estimates for the parameters.
#' @param q Numeric value representing the underreporting factor.
#' @return Matrix of posterior probabilities for both groups (group 1 and group 2).
posterior_probabilities <- function(data, covariates, max_llh) {
  n <- 96
  post <- matrix(nrow = n, ncol = 2)
  w <- vector()
  
  for (i in 1:n) {
  time <- i / n
  w[i] <- exp(max_llh$estimate[1] + max_llh$estimate[2] * time) / 
      (1 + exp(max_llh$estimate[1] + max_llh$estimate[2] * time))
  q_t <- q_time(time)
  
  mean_val <- max_llh$estimate[3] + 
    max_llh$estimate[4] * time + 
    max_llh$estimate[5] * covariates[i, 2] + 
    max_llh$estimate[6] * covariates[i, 3] + 
    max_llh$estimate[7] * covariates[i, 4] +
    max_llh$estimate[8] * covariates[i, 5] + 
    max_llh$estimate[9] * covariates[i, 6]
  
  sd_val <- max.llh$estimate[10]
  
  numerator1 <- w[i] * dnorm(data$incid[i], mean = mean_val, sd = sd_val)
  numerator2 <- (1 - w[i]) * dnorm(data$incid[i], mean = mean_val / q_t, sd = sd_val / q_t)
  denominator <- numerator1 + numerator2
  print(data$incid[i])
  print(mean_val / q_t)
  print(sd_val / q_t)
  cat("\n\n")
  post[i, 1] <- numerator1 / denominator
  post[i, 2] <- numerator2 / denominator
  }
  
  return(post)
}

# Reconstruct incidences based on posterior probabilities
#' @param incid Numeric vector of observed incidence values.
#' @param post Matrix of posterior probabilities.
#' @param q Numeric value representing the underreporting factor.
#' @return Numeric vector of reconstructed incidences.
reconstruct_incid <- function(incid, post, q_time_func) {
  xrec <- numeric(length(incid))
  for (i in 1:length(incid)) {
  time <- i / length(incid)  
  q_t <- q_time_func(time)  
  xrec[i] <- ifelse(post[i, 2] > 0.5, incid[i], incid[i] / q_t)
  }
  return(xrec)
}

# Plot time series of observed and reconstructed incidences
#' @param incid Numeric vector of observed incidence values.
#' @param xrec Numeric vector of reconstructed incidences.
#' @param start_year Integer representing the start year of the time series.
#' @param end_year Integer representing the end year of the time series.
#' @param ylim_range Numeric vector defining the y-axis range for the plot.
#' @param main_title Character string representing the title of the plot.
plot_time_series <- function(incid, xrec, start_year, end_year, ylim_range, main_title) {
  ts.plot(ts(incid, start = c(start_year, 1), end = c(end_year, 12), freq = 12), ylim = ylim_range, ylab = "Incidence x 100,000", main = main_title)
  lines(seq(start_year, end_year + 0.99, 1 / 12), xrec, col = "red", lty = 2)
  legend("topright", legend = c("Observed", "Reconstructed"), col = c("black", "red"), lty = c(1, 2))
}

# Process and analyze incidence data for a specific demographic group
#' @param data Data frame containing the demographic group's incidence data.
#' @param covariates Matrix of covariates for the demographic group.
#' @param start_year Integer representing the start year of the time series.
#' @param end_year Integer representing the end year of the time series.
#' @param ylim_range Numeric vector defining the y-axis range for the plot.
#' @param main_title Character string representing the title of the plot and output file.
process_demographic_group <- function(data, covariates, start_year, end_year, ylim_range, main_title) {
  cat("Processing demographic group\n")
  post <- posterior_probabilities(data, covariates, max.llh)
  cat("Posterior probabilities calculated\n")
  xrec <- reconstruct_incid(data$incid, post, q_time)
  cat("Incidence reconstructed\n")
  plot_time_series(data$incid, xrec, start_year, end_year, ylim_range, main_title)
  cat("Time series plotted\n")
  
  # Load the paper result
  file_name <- gsub(" ", "-", main_title)
  xrec_paper <- readRDS(paste0("R/data/", file_name, ".RDS"))

  # Plot the percentage difference
  percentage_diff <- (xrec_paper - xrec) / xrec_paper * 100
  p <- ts.plot(ts(percentage_diff, start = c(start_year, 1), end = c(end_year, 12), freq = 12), ylim = c(-5, 5), ylab = "Percentage Difference", main = paste0("Percentage Difference: ", main_title))
}

process_demographic_group(pr3[pr3$sexe == 0 & pr3$age == 0, ], covars[covars[, 2] == 0 & covars[, 3] == 0, ], 2009, 2016, c(9, 32), "Women 15-29 years old")
process_demographic_group(pr3[pr3$sexe == 0 & pr3$age == 1, ], covars[covars[, 2] == 1 & covars[, 3] == 0, ], 2009, 2016, c(1, 8), "Women 30-94 years old")
process_demographic_group(pr3[pr3$sexe == 1 & pr3$age == 0, ], covars[covars[, 2] == 0 & covars[, 3] == 1, ], 2009, 2016, c(4, 32), "Men 15-29 years old")
process_demographic_group(pr3[pr3$sexe == 1 & pr3$age == 1, ], covars[covars[, 2] == 1 & covars[, 3] == 1, ], 2009, 2016, c(1, 12), "Men 30-94 years old")
