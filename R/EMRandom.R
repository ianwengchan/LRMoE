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

#' ## Single Update for beta and ww jointly
#' #' Single update in M-step for random effect coefficients \code{beta} and realization of random effects \code{ww}.
#' #' Specific for j-th latent class and l-th random effect (suppressing subscripts below).
#' #'
#' #' @param betajl A scalar, random effect coefficient beta of the j-th latent class and l-th random effect.
#' #' @param wl A vector of length Sl, realization of the l-th random effect (with Sl clusters).
#' #' @param dbeta A scalar, dQ/dbeta.
#' #' @param dw A vector of length Sl, dQ/dw.
#' #' @param dbeta2 A scalar, d2Q/dbeta2.
#' #' @param dwdbeta A vector of length Sl, d2Q/dwdbeta.
#' #' @param dw2 A Sl*Sl matrix, d2Q/dw2.
#' #'
#' #' @return \code{betajl.update} Updated random effect coefficient betajl.
#' #' @return \code{wl.update} Updated realization of random effects wl.
#' #'
#' #' @importFrom matrixStats rowLogSumExps
#' #'
#' #' @keywords internal
#' #'
#' #' # @export EMMupdatebetaww.random
#' EMMupdatebetaww.random = function(betajl, wl, dbeta, dw, dbeta2, dwdbeta, dw2)
#' {
#'   # Entries of the inverse of Hessian
#'   inverse = solve(dw2 - 1/dbeta2*tcrossprod(dwdbeta)) # = lower_right
#'   Binverse = t(dwdbeta)%*%inverse # = B%*%inverse
#'
#'   upper_left = 1/dbeta2 + (1/dbeta2)^2 * Binverse%*%dwdbeta
#'   upper_right = -1/dbeta2 * Binverse
#'
#'   betajl.update = betajl - (upper_left * dbeta + upper_right%*%dw)
#'   wl.update = wl - (dbeta*t(upper_right) + inverse%*%dw)
#'
#'   return(list(betajl.update = betajl.update, wl.update = wl.update))
#' }
#'
#' ## With Random Effects, M-Step for beta and ww jointly
#' #' ECM: M-Step for random effect coefficients \code{beta} and realization of random effects \code{ww}.
#' #'
#' #' @param X A N*P matrix of numerical covariates.
#' #' @param alpha A g*P matrix of logit regression coefficients (from last iteration).
#' #' @param t A list of L matrices, where matrix l is a N*Sl matrix of 0 and 1's,
#' #'          indicating the known clustering of the observations with respect to the l-th random effect.
#' #' @param ww A list of L vectors of old realization of random effects.
#' #'           The l-th vector is of length Sl, representing the Sl clusters of the l-th random effect.
#' #' @param beta A g*L matrix of old random effect coefficients.
#' #' @param comp.zkz.e.list An object returned by \code{EMEzkz}.
#' #' @param beta.iter.max Numeric: maximum number of iterations.
#' #' @param penalty TRUE/FALSE, which indicates whether penalty is applied.
#' #' @param hyper.beta A numeric of penalty applied to \code{beta}.
#' #'
#' #' @return \code{beta.new} Updated random effect coefficients.
#' #' @return \code{ww.new} Updated realization of random effects.
#' #'
#' #' @importFrom matrixStats rowLogSumExps
#' #'
#' #' @keywords internal
#' #'
#' #' # @export EMMbetaww.random
#' EMMbetaww.random = function(X, alpha, t, ww, beta, comp.zkz.e.list,
#'                           beta.iter.max, penalty, hyper.beta)
#' {
#'   comp.zpzk = XPlusYColTimesZ(comp.zkz.e.list$z.e.obs, comp.zkz.e.list$z.e.lat, comp.zkz.e.list$k.e)
#'   # comp.zkz.e.list$z.e.obs + sweep(comp.zkz.e.list$z.e.lat, 1, comp.zkz.e.list$k.e, FUN = "*", check.margin = FALSE)
#'
#'   sample.size.n = nrow(X)
#'   n.rand.l = length(ww)
#'   n.comp = nrow(beta)
#'   # iter=array(0, dim = c(n.comp, 1))
#'   beta.new = beta
#'   beta.old = beta - Inf
#'   ww.new = ww
#'   # ww.old = ww
#'   comp.zpzk.marg = apply(comp.zpzk, 1, sum)
#'
#'   for (l in 1:n.rand.l)
#'   {
#'     iter=array(0, dim = c(n.comp, 1))
#'     tl = t[[l]] # retrieve the corresponding t
#'     # wl = ww.new[[l]] # retrieve the corresponding ww
#'     # ww.old[[l]] <- ww.new[[l]] - Inf
#'
#'     for (j in 1:(n.comp-1)) # The last component's beta's are always kept at zero (reference category).
#'     {
#'       while ((iter[j]<=beta.iter.max)&(sum((beta.old[j,l]-beta.new[j,l])^2)>10^(-8))) # Stopping criteria: (beta.iter.max) iterations, or small difference
#'       {
#'         beta.old[j,l] = beta.new[j,l]
#'         W.new = ProduceW(t, ww.new)
#'         gate.body = tcrossprod(X,alpha)+tcrossprod(W.new,beta.new)
#'         pp = exp(gate.body-rowLogSumExps(gate.body))
#'         qqj = exp(rowLogSumExps(array(gate.body[,-j],dim=c(sample.size.n,n.comp-1)))-rowLogSumExps(gate.body))
#'
#'         dbeta = EMalphadQ(W.new, comp.zpzk[,j], comp.zpzk.marg, pp[,j])[l] # ignore penalty for now
#'         dbeta2 = EMalphadQ2(W.new, comp.zpzk.marg, pp[,j], qqj)[l,l] # ignore penalty for now
#'
#'         dw = EMwwdQ(comp.zpzk, pp, beta.new[,l], tl, ww.new[[l]]) # ww has no penalty; its prior has been included in the M-step
#'         dw2 = EMwwdQ2(pp, beta.new[,l], tl)
#'
#'         wbetal = pp%*%beta.new[,l]
#'         dwdbeta = EMwwbetadQ2(comp.zpzk[,j], pp[,j], pp, beta.new[j,l], wbetal, ww.new[[l]], tl)
#'
#'         # Use square root beta to ensure positive updates
#'         sqrtbeta = sqrt(beta.new[j,l])
#'         dsqrtbeta = dbeta*(2*sqrtbeta)
#'         dsqrtbeta2 = dbeta*2 + (2*sqrtbeta)^2 * dbeta2
#'         dwdsqrtbeta = dwdbeta*(2*sqrtbeta)
#'
#'         # single_update = EMMupdatebetaww.random(sqrtbeta, ww.new[[l]][-1], dsqrtbeta, dw[-1], dsqrtbeta2, dwdsqrtbeta[-1], dw2[-1,-1])
#'         #
#'         # beta.new[j,l] = (single_update$betajl.update)^2
#'         # ww.new[[l]][-1] = single_update$wl.update
#'
#'         single_update = EMMupdatebetaww.random(sqrtbeta, ww.new[[l]], dsqrtbeta, dw, dsqrtbeta2, dwdsqrtbeta, dw2)
#'
#'         beta.new[j,l] = (single_update$betajl.update)^2
#'         ww.new[[l]] = single_update$wl.update
#'
#'         # single_update = EMMupdatebetaww.random(beta.new[j,l], ww.new[[l]], dbeta, dw, dbeta2, dwdbeta, dw2)
#'         #
#'         # beta.new[j,l] = (single_update$betajl.update)
#'         # ww.new[[l]] = single_update$wl.update
#'
#'         iter[j] = iter[j] + 1
#'       }
#'     }
#'   }
#'   return(list(beta.new = beta.new, ww.new = ww.new))
#' }

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

  if (n.comp > 2)
  {
    for (j in 2:(n.comp-1)) # The last component's beta's are always kept at zero (reference category).
      # The first component's beta is estimated standard deviation
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
  }

  return(beta.new)
}

## With Random Effects, M-Step for ww
#' ECM: M-Step for realization of random effects \code{ww}.
#'
#' @param X A N*P matrix of numerical covariates.
#' @param alpha A g*P matrix of logit regression coefficients (from last iteration).
#' @param t A list of L matrices, where matrix l is a N*Sl matrix of 0 and 1's,
#'          indicating the known clustering of the observations with respect to the l-th random effect.
#' @param ww A list of L vectors of old realization of random effects.
#'           The l-th vector is of length Sl, representing the Sl clusters of the l-th random effect.
#' @param beta A g*L matrix of random effect coefficients (from last iteration).
#' @param sigma A vector of length L, where the l-th entry is the estimated standard deviation of the l-th random effect.
#' @param comp.zkz.e.list An object returned by \code{EMEzkz}.
#' @param ww.iter.max Numeric: maximum number of iterations.
#'
#' @return \code{ww.new} Updated realization of random effects.
#' @return \code{sigma.new} Updated standard deviations of random effects.
#'
#' @importFrom matrixStats rowLogSumExps
#'
#' @keywords internal
#'
#' @export EMMww.random
EMMww.random = function(X, alpha, t, ww, beta, sigma, comp.zkz.e.list, ww.iter.max)
{
  comp.zpzk = XPlusYColTimesZ(comp.zkz.e.list$z.e.obs, comp.zkz.e.list$z.e.lat, comp.zkz.e.list$k.e)
  # comp.zkz.e.list$z.e.obs + sweep(comp.zkz.e.list$z.e.lat, 1, comp.zkz.e.list$k.e, FUN = "*", check.margin = FALSE)

  sample.size.n = nrow(X)
  n.rand.l = length(ww)
  n.comp = nrow(alpha)
  iter = array(0, dim = c(n.rand.l, 1))
  ww.new = ww
  ww.old = ww
  sigma.new = sigma
  W.new = ProduceW(t, ww.new) # ProduceW from last iteration
  comp.zpzk.marg = apply(comp.zkz.e.list$z.e.obs, 1, sum)

  for (l in 1:n.rand.l)
  {
    tl = t[[l]]
    ww.old[[l]] <- ww.new[[l]] - Inf
    sigmal <- sigma.new[l]

    while ((iter[l]<=ww.iter.max)&(sum((ww.old[[l]]-ww.new[[l]])^2)>10^(-8))) # Stopping criteria: (ww.iter.max) iterations, or small difference
    {
      ww.old[[l]] = ww.new[[l]] # keeping all other wwl's unchanged; keep track of last iteration
      W.new[,l] = tl%*%ww.new[[l]] # minimal update to W: only changes the l-th column
      gate.body = tcrossprod(X,alpha) + tcrossprod(W.new,beta)
      pp = exp(gate.body-rowLogSumExps(gate.body))
      dQ = EMwwdQ(comp.zpzk, comp.zpzk.marg, pp, beta[,l], tl, ww.new[[l]], sigmal) # ww has no penalty; its prior has been included in the M-step
      dQ2 = EMwwdQ2(comp.zpzk.marg, pp, beta[,l], tl, sigmal)

      ww.new[[l]] = c(ww.new[[l]] + crossprod(dQ, chol2inv(chol(-dQ2)))) # latest
      # ww.new[[l]] = c(ww.new[[l]] - crossprod(dQ, solve(dQ2)))
      iter[l] = iter[l]+1
    }
    sigma.new[l] <- sqrt(sum(ww.new[[l]]^2) / length(ww.new[[l]]))
    # ww.new[[l]] = (ww.new[[l]] - mean(ww.new[[l]]))/sd(ww.new[[l]])
  }

  return(list(ww.new = ww.new, sigma.new = sigma.new))
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


