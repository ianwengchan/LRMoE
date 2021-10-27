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
