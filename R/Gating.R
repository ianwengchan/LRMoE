## Gating Function: Multiple Logit
#' Computes the logit regression weights in log.
#'
#' @param x An N * P covariate matrix, where N is sample size. The first column MUST be 1.
#' @param alpha A g * P matrix. Logit regression coefficients.
#' @return The log of \code{tcrossprod(x,alpha)} with dim N * g, with each row normalized by \code{rowLogSumExps}
#'
#' @keywords internal
#'
#' @export GateLogit
GateLogit = function(x, alpha)
{
  gate.body=tcrossprod(x,alpha)
  return(sweep(gate.body, 1, rowLogSumExps(gate.body), FUN = "-", check.margin = FALSE))
}

## Gating Function with Random Effects: Multiple Logit
#' Computes the logit regression weights in log.
#'
#' @param x An N * P covariate matrix, where N is sample size. The first column MUST be 1.
#' @param alpha A g * P matrix. Logit regression coefficients.
#' @param w An N * L random effect matrix, where N is sample size, L is number of random effects.
#' @param beta A g * L matrix. Random effect coefficients.
#' @return The log of \code{tcrossprod(x,alpha)+tcrossprod(w,beta)} with dim N * g, with each row normalized by \code{rowLogSumExps}
#'
#' @keywords internal
#'
#' @export GateLogitRandom
GateLogitRandom = function(x, alpha, w, beta)
{
  gate.body=tcrossprod(x,alpha)+tcrossprod(w,beta)
  return(sweep(gate.body, 1, rowLogSumExps(gate.body), FUN = "-", check.margin = FALSE))
}


LogitGatingEval = function(X, alpha, beta, t, ww.list)
{
  W = ProduceW(t, ww.list)
  gate.body=tcrossprod(X,alpha)+tcrossprod(W,beta)
  return(sweep(gate.body, 1, rowLogSumExps(gate.body), FUN = "-", check.margin = FALSE))
}


LogitGatingSim = function(X, alpha, beta, t, ww.mu.list, ww.Sigma.list)
{
  ww.list.sample <- as.list(rep(NA, length(ww.mu.list)))
  for (l in 1:length(ww.mu.list)){
    ww.list.sample[[l]] <- mvtnorm::rmvnorm(1, mean = ww.mu.list[[l]], sigma = ww.Sigma.list[[l]])
  }
  W <- ProduceW(t, ww.list.sample)
  gate.body=tcrossprod(X,alpha)+tcrossprod(W,beta)
  return(sweep(gate.body, 1, rowLogSumExps(gate.body), FUN = "-", check.margin = FALSE))
}


DiagMvNormal_KL = function(mu, Sigma)
{
  return(0.5*(-sum(log(Sigma)) + sum(Sigma) + sum(mu^2) - length(mu)))
}


ww_KL = function(ww.mu.list, ww.Sigma.list)
{
  temp = 0
  for (l in 1:length(ww.mu.list)){
    temp = temp + DiagMvNormal_KL(ww.mu.list[[l]], ww.Sigma.list[[l]])
  }
  return(temp)
}
