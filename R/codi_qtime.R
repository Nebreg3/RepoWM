library(mixtools)
library(ggplot2)
library(gridExtra)
library(hrbrthemes)
library(ggpubr)

# Load custom negative log-likelihood function
source("R/llh_covars_qtime.R")

# Load main data
load("R/data_main.RData")

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
  # ts.plot(ts(incid, start = c(start_year, 1), end = c(end_year, 12), freq = 12), ylim = ylim_range, ylab = "Incidence x 100,000", main = main_title)
  # lines(seq(start_year, end_year + 0.99, 1 / 12), xrec, col = "red", lty = 2)
  # legend("topright", legend = c("Observed", "Reconstructed"), col = c("black", "red"), lty = c(1, 2))
  
  p <- ggplot(data = data.frame(incid = incid, xrec = xrec, time = seq(start_year, end_year + 0.99, 1 / 12))) +
    geom_line(aes(x = time, y = incid), color = "black") +
    geom_line(aes(x = time, y = xrec), color = "red", linetype = "dashed") +
    labs(x = "Year", y = "Incidence x 100,000", title = main_title) +
    scale_x_continuous(breaks = seq(start_year, end_year, 1), 
                       labels = scales::label_number(accuracy = 1)) +
    theme_ipsum_ps()

  show(p)
  
  return(p)
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
  p1 <- plot_time_series(data$incid, xrec, start_year, end_year, ylim_range, main_title)
  cat("Time series plotted\n")
  
  # Load the paper result
  file_name <- gsub(" ", "-", main_title)
  xrec_paper <- readRDS(paste0("R/data/", file_name, ".RDS"))

  # Plot the percentage difference
  percentage_diff <- (xrec_paper - xrec) / xrec_paper * 100

  # ts.plot(ts(percentage_diff, start = c(start_year, 1), end = c(end_year, 12), freq = 12), ylim = c(-5, 5), ylab = "Percentage Difference", main = paste0("Percentage Difference: ", main_title))
  
  p2 <- ggplot(data = data.frame(percentage_diff = percentage_diff, time = seq(start_year, end_year + 0.99, 1 / 12))) +
    geom_line(aes(x = time, y = percentage_diff)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = "Year", y = "Percentage Difference", title = main_title) +
    scale_x_continuous(breaks = seq(start_year, end_year, 1), 
                       labels = scales::label_number(accuracy = 1)) +
    theme_ipsum_ps()
  show(p2)
  
  return(p1)
}

#' Compute global residuals for the model
#' @param max_llh A list containing the maximum likelihood estimates for the parameters.
#' @param pr3 Data frame containing the incidence data.
#' @param covars Matrix of covariates for the dataset.
#' @param q_time Function representing the dynamic underreporting factor.
perform_residual_analysis <- function(max_llh, pr3, covars, q_time) {
  # Initialize vectors
  y_est <- vector()
  w     <- vector()
  
  # Loop to calculate w and y_est
  for (i in 1:384) {
    j <- (i %% 96) / 96
    if (j == 0) j <- 1
    
    w[i] <- exp(max_llh$estimate[1] + max_llh$estimate[2] * j) /
      (1 + exp(max_llh$estimate[1] + max_llh$estimate[2] * j))
    
    m <- max_llh$estimate[3] + max_llh$estimate[4] * j + 
      max_llh$estimate[5] * pr3$age[i] +
      max_llh$estimate[6] * pr3$sexe[i] +
      max_llh$estimate[7] * pr3$inter[i] +
      max_llh$estimate[8] * covars[i, 5] +
      max_llh$estimate[9] * covars[i, 6]
    
    y_est[i] <- w[i] * m + (1 - w[i]) * m / q_time(j)
  }
  
  # Aggregate y_est
  y_est_agg_temp <- data.frame(
    t = rep(seq(1:96), 4),
    sexe = c(rep(0, 192), rep(1, 192)),
    edat = c(rep(0, 96), rep(1, 96), rep(0, 96), rep(1, 96)),
    y_est
  )
  
  tw <- sum(unique(pr3$Pob)) 
  y_est_agg <- aggregate(y_est_agg_temp$y_est, by = list(y_est_agg_temp$t), FUN = sum)
  y_agg <- aggregate(pr3$incid, by = list(pr3$mes_any_problema), FUN = sum)
  
  # Normalize y_est_agg and y_agg by year
  for (year in 2009:2016) {
    indices <- ((year - 2009) * 12 + 1):((year - 2008) * 12)
    norm_factor <- sum(unique(pr3$Pob[pr3$Year == year])) / tw
    y_est_agg$x[indices] <- y_est_agg$x[indices] * norm_factor
    y_agg$x[indices] <- y_agg$x[indices] * norm_factor
  }
  
  # Calculate residuals
  resid <- y_est_agg$x - y_agg$x
  
  # Compute ACF and PACF
  bacf <- acf(resid, lag.max = 10, plot = FALSE)
  bacfdf <- with(bacf[1:10], data.frame(lag, acf))
  ciline <- qnorm((1 - 0.95) / 2) / sqrt(length(y_agg$x))
  
  ACF <- ggplot(data = bacfdf, mapping = aes(x = as.integer(lag), y = acf)) +
    geom_hline(aes(yintercept = 0)) +
    geom_hline(aes(yintercept = ciline), linetype = 2) +
    geom_hline(aes(yintercept = -ciline), linetype = 2) +
    geom_segment(mapping = aes(xend = lag, yend = 0)) +
    labs(x = "Lag", y = "", title = "ACF") +
    theme_ipsum_ps()
  
  bacf <- pacf(resid, lag.max = 10, plot = FALSE)
  bacfdf <- with(bacf, data.frame(lag, acf))
  
  PACF <- ggplot(data = bacfdf, mapping = aes(x = as.integer(lag), y = acf)) +
    geom_hline(aes(yintercept = 0)) +
    geom_hline(aes(yintercept = ciline), linetype = 2) +
    geom_hline(aes(yintercept = -ciline), linetype = 2) +
    geom_segment(mapping = aes(xend = lag, yend = 0)) +
    labs(x = "Lag", y = "", title = "PACF") +
    theme_function()
  
  # Combine ACF and PACF plots
  ggarrange(ACF, PACF, ncol = 2, nrow = 1)
}

p1 <- process_demographic_group(pr3[pr3$sexe == 0 & pr3$age == 0, ], covars[covars[, 2] == 0 & covars[, 3] == 0, ], 2009, 2016, c(9, 32), "Women 15-29 years old")
p2 <- process_demographic_group(pr3[pr3$sexe == 0 & pr3$age == 1, ], covars[covars[, 2] == 1 & covars[, 3] == 0, ], 2009, 2016, c(1, 8), "Women 30-94 years old")
p3 <- process_demographic_group(pr3[pr3$sexe == 1 & pr3$age == 0, ], covars[covars[, 2] == 0 & covars[, 3] == 1, ], 2009, 2016, c(4, 32), "Men 15-29 years old")
p4 <- process_demographic_group(pr3[pr3$sexe == 1 & pr3$age == 1, ], covars[covars[, 2] == 1 & covars[, 3] == 1, ], 2009, 2016, c(1, 12), "Men 30-94 years old")

# compute_global_residues(max.llh, pr3, covars, q)

ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)

perform_residual_analysis(max.llh, pr3, covars, q_time)

