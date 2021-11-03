## Random Effect Mapping
#' Produce the random effect matrix W from clustering matrix t and realization of random effects ww.
#'
#' @param t A list of L matrices, where matrix l is a N * Sl matrix of 0 and 1's,
#'           indicating the known clustering of the observations with respect to the l-th random effect.
#' @param ww A list of L vectors, where the l-th vector is of length Sl, representing the Sl clusters of the l-th random effect.
#' @return An N * L matrix of random effect values, where N is sample size, L is number of random effects.
#'
#' @keywords internal
#'
#' @export ProduceW
ProduceW <- function(t, ww){
  result <- matrix(NA, nrow=nrow(t[[1]]), ncol=length(ww))
  for(j in 1:length(ww)){
    result[,j] = t[[j]]%*%ww[[j]]
  }
  return(result)
}

## Initial guess of ww
#' Simulate an initial guess of random effect realization ww from clustering matrix t.
#'
#' @param t A list of L matrices, where matrix l is a N * Sl matrix of 0 and 1's,
#'          indicating the known clustering of the observations with respect to the l-th random effect.
#' @return A list of L vectors, where the l-th vector is of length Sl, representing the Sl clusters of the l-th random effect.
#'
#' @keywords internal
#'
#' @export Initww
Initww <- function(t){
  ww.init <- as.list(rep(NA, length(t)))
  for(l in 1:length(t)){
    ww.init[[l]] <- rnorm(ncol(t[[l]]), 0, 1)
  }
  return(ww.init)
}
