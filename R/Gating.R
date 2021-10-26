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
