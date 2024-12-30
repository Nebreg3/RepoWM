llh <- function(pars, data, covars) 
{
  loglik <- 0
  for (i in 1:length(data))
  {
    # Calculate time-dependent w using pars[1] (intercept) and pars[2] (slope)
    logit.w <- pars[1] + pars[2] * covars[i, 1]
    w       <- exp(logit.w) / (1 + exp(logit.w))
    
    # Calculate time-dependent q using pars[11] (intercept) and pars[12] (slope)
    logit.q <- pars[11] + pars[12] * covars[i, 1]
    q       <- exp(logit.q) / (1 + exp(logit.q))
    
    # Calculate the log-likelihood contribution for observation i
    loglik <- loglik + log(
      w * dnorm(
        data[i], 
        mean = (pars[3] + pars[4] * covars[i, 1] + pars[5] * covars[i, 2] + 
                pars[6] * covars[i, 3] + pars[7] * covars[i, 4] + 
                pars[8] * covars[i, 5] + pars[9] * covars[i, 6]), 
        sd = pars[10]
      ) + 
      (1 - w) * dnorm(
        data[i], 
        mean = (pars[3] + pars[4] * covars[i, 1] + pars[5] * covars[i, 2] + 
                pars[6] * covars[i, 3] + pars[7] * covars[i, 4] + 
                pars[8] * covars[i, 5] + pars[9] * covars[i, 6]) / q, 
        sd = pars[10] / q
      )
    )
  }
  
  # Return the negative log-likelihood
  return((-1) * loglik)
}
