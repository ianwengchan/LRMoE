## Random Effect Mapping
#' Produce the random effect matrix W from clustering matrix tl and realization of random effects wl.
#'
#' @param tl A list of L matrices, where matrix l is a N * Sl matrix of 0 and 1's,
#'           indicating the known clustering of the observations with respect to the l-th random effect.
#' @param wl A list of L vectors, where the l-th vector is of length Sl, representing the Sl clusters of the l-th random effect.
#' @return An N * L matrix of random effect values, where N is sample size, L is number of random effects.
#'
#' @keywords internal
#'
#' @export ProduceW
ProduceW <- function(tl, wl){
  result <- matrix(NA, nrow=nrow(tl[[1]]), ncol=length(wl))
  for(j in 1:length(wl)){
    result[,j] = tl[[j]]%*%wl[[j]]
  }
  return(result)
}
