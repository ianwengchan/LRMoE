## With Random Effects, M-Step for alpha
#' ECM: M-Step for logit regression coefficients \code{alpha}.
#'
#' @param X A N*P matrix of numerical covariates.
#' @param alpha A g*P matrix of old logit regression coefficients.
#' @param W A N*L matrix of random effects (from last iteration).
#' @param beta A g*L matrix of random effect coefficients (from last iteration).
#' @param comp.zkz.e.list An object returned by \code{EMEzkz}.
#' @param alpha.iter.max Numeric: maximum number of iterations.
#' @param penalty TRUE/FALSE, which indicates whether penalty is applied.
#' @param hyper.alpha A numeric of penalty applied to \code{alpha}.
#'
#' @return \code{alpha.new} Updated logit regression coefficients.
#'
#' @importFrom matrixStats rowLogSumExps
#'
#' @keywords internal
#'
#' # @export EMMalpha.random
EMMalpha.random = function(X, alpha, W, beta, comp.zkz.e.list,
                           alpha.iter.max, penalty, hyper.alpha)
{
  comp.zpzk = XPlusYColTimesZ(comp.zkz.e.list$z.e.obs, comp.zkz.e.list$z.e.lat, comp.zkz.e.list$k.e)
  # comp.zkz.e.list$z.e.obs + sweep(comp.zkz.e.list$z.e.lat, 1, comp.zkz.e.list$k.e, FUN = "*", check.margin = FALSE)

  sample.size.n = nrow(X)
  n.covar.p = ncol(X)
  n.comp = nrow(alpha)
  iter=array(0, dim = c(n.comp, 1))
  alpha.new = alpha
  alpha.old = alpha - Inf
  comp.zpzk.marg = apply(comp.zpzk, 1, sum)

  for (j in 1:(n.comp-1)) # The last component's alpha's are always kept at zero (reference category).
  {
    while ((iter[j]<=alpha.iter.max)&(sum((alpha.old[j,]-alpha.new[j,])^2)>10^(-8))) # Stopping criteria: (alpha.iter.max) iterations, or small difference
    {
      alpha.old[j,]=alpha.new[j,]
      gate.body=tcrossprod(X,alpha.new)+tcrossprod(W,beta)
      pp = exp(gate.body-rowLogSumExps(gate.body))
      qqj = exp(rowLogSumExps(array(gate.body[,-j],dim=c(sample.size.n,n.comp-1)))-rowLogSumExps(gate.body))
      dQ = EMalphadQ(X, comp.zpzk[,j], comp.zpzk.marg, pp[,j]) - if(penalty){alpha.new[j,]/hyper.alpha^2} else{0}
      # apply(sweep(X,1,comp.zpzk[,j]-comp.zpzk.marg*exp(gate.body[,j]-rowLogSumExps(gate.body)),FUN="*",check.margin=FALSE),2,sum)-if(penalty){alpha.new[j,]/hyper.alpha^2} else{0}
      dQ2 = EMalphadQ2(X, comp.zpzk.marg, pp[,j], qqj) - if(penalty){diag(1/hyper.alpha^2,nrow = n.covar.p, ncol = n.covar.p)} else{diag(10^(-7),nrow = n.covar.p, ncol = n.covar.p)}
      # -crossprod(sweep(X,1,comp.zpzk.marg*exp(rowLogSumExps(array(gate.body[,-j],dim=c(sample.size.n,n.comp-1)))+gate.body[,j]-2*rowLogSumExps(gate.body)),FUN="*",check.margin=FALSE),X)-if(penalty){diag(1/hyper.alpha^2,nrow = n.covar.p, ncol = n.covar.p)} else{diag(10^(-7),nrow = n.covar.p, ncol = n.covar.p)}

      alpha.new[j,]=alpha.new[j,] + crossprod(dQ, chol2inv(chol(-dQ2))) # -crossprod(dQ,solve(dQ2))
      iter[j] = iter[j]+1
    }
  }

  return(alpha.new)
}

## With Random Effects, M-Step for beta
#' ECM: M-Step for random effect coefficients \code{beta}.
#'
#' @param X A N*P matrix of numerical covariates.
#' @param alpha A g*P matrix of logit regression coefficients (from last iteration).
#' @param W A N*L matrix of random effects (from last iteration).
#' @param beta A g*L matrix of old random effect coefficients.
#' @param comp.zkz.e.list An object returned by \code{EMEzkz}.
#' @param beta.iter.max Numeric: maximum number of iterations.
#' @param penalty TRUE/FALSE, which indicates whether penalty is applied.
#' @param hyper.beta A numeric of penalty applied to \code{beta}.
#'
#' @return \code{beta.new} Updated random effect coefficients.
#'
#' @importFrom matrixStats rowLogSumExps
#'
#' @keywords internal
#'
#' # @export EMMbeta.random
EMMbeta.random = function(X, alpha, W, beta, comp.zkz.e.list,
                          beta.iter.max, penalty, hyper.beta)
{
  comp.zpzk = XPlusYColTimesZ(comp.zkz.e.list$z.e.obs, comp.zkz.e.list$z.e.lat, comp.zkz.e.list$k.e)
  # comp.zkz.e.list$z.e.obs + sweep(comp.zkz.e.list$z.e.lat, 1, comp.zkz.e.list$k.e, FUN = "*", check.margin = FALSE)

  sample.size.n = nrow(X)
  n.rand.l = ncol(W)
  n.comp = nrow(beta)
  iter=array(0, dim = c(n.comp, 1))
  beta.new = beta
  beta.old = beta - Inf
  comp.zpzk.marg = apply(comp.zpzk, 1, sum)

  for (j in 1:(n.comp-1)) # The last component's beta's are always kept at zero (reference category).
  {
    while ((iter[j]<=beta.iter.max)&(sum((beta.old[j,]-beta.new[j,])^2)>10^(-8))) # Stopping criteria: (beta.iter.max) iterations, or small difference
    {
      beta.old[j,]=beta.new[j,]
      gate.body=tcrossprod(X,alpha)+tcrossprod(W,beta.new)
      pp = exp(gate.body-rowLogSumExps(gate.body))
      qqj = exp(rowLogSumExps(array(gate.body[,-j],dim=c(sample.size.n,n.comp-1)))-rowLogSumExps(gate.body))
      dQ = EMalphadQ(W, comp.zpzk[,j], comp.zpzk.marg, pp[,j]) - if(penalty){beta.new[j,]/hyper.beta^2} else{0}
      # apply(sweep(X,1,comp.zpzk[,j]-comp.zpzk.marg*exp(gate.body[,j]-rowLogSumExps(gate.body)),FUN="*",check.margin=FALSE),2,sum)-if(penalty){alpha.new[j,]/hyper.alpha^2} else{0}
      dQ2 = EMalphadQ2(W, comp.zpzk.marg, pp[,j], qqj) - if(penalty){diag(1/hyper.beta^2,nrow = n.rand.l, ncol = n.rand.l)} else{diag(10^(-7),nrow = n.rand.l, ncol = n.rand.l)}
      # -crossprod(sweep(X,1,comp.zpzk.marg*exp(rowLogSumExps(array(gate.body[,-j],dim=c(sample.size.n,n.comp-1)))+gate.body[,j]-2*rowLogSumExps(gate.body)),FUN="*",check.margin=FALSE),X)-if(penalty){diag(1/hyper.alpha^2,nrow = n.covar.p, ncol = n.covar.p)} else{diag(10^(-7),nrow = n.covar.p, ncol = n.covar.p)}

      beta.new[j,]=beta.new[j,] + crossprod(dQ, chol2inv(chol(-dQ2))) # -crossprod(dQ,solve(dQ2))
      iter[j] = iter[j]+1
    }
  }

  return(beta.new)
}


#' @keywords internal
EM_E_z_obs <- function(gate_expert_ll_comp, gate_expert_ll) {
  return(exp(XColMinusY(gate_expert_ll_comp, gate_expert_ll)))
}

#' @keywords internal
EM_E_z_lat <- function(gate_expert_tn_bar_comp, gate_expert_tn_bar) {
  tmp = exp(XColMinusY(gate_expert_tn_bar_comp, gate_expert_tn_bar))
  tmp[is.na(tmp)] = 1/ncol(gate_expert_tn_bar_comp)
  return(tmp)
}

#' @keywords internal
EM_E_k <- function(gate_expert_tn) {
  return(expm1(-gate_expert_tn))
}

#' @keywords internal
EM_E_z_zero_obs <- function(yl, p_old, gate_expert_ll_pos_comp){
  tmp = ifelse(yl==0, p_old/(p_old + (1-p_old)*exp(gate_expert_ll_pos_comp)), 0.0)
  return(tmp)
}

#' @keywords internal
EM_E_z_zero_lat <- function(tl, p_old, gate_expert_tn_bar_pos_comp){
  tmp = ifelse(tl>0, p_old/(p_old + (1-p_old)*exp(gate_expert_tn_bar_pos_comp)), 0.0)
  return(tmp)
}

#' @keywords internal
EM_M_zero <- function(z_zero_e_obs, z_pos_e_obs, z_zero_e_lat, z_pos_e_lat, k_e) {
  num = sum(z_zero_e_obs + (z_zero_e_lat * k_e))
  denom = num + sum(z_pos_e_obs + (z_pos_e_lat * k_e))
  return(num/denom)
}


