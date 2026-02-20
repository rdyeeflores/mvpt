
#' Calculate Matrix A (Internal)
#'
#' For a fitted SEM, computes the A information matrix from Vuong (1989).
#'
#' @param SEMfitted A lavaan-fitted SEM
#' @return Matrix A
#' @keywords internal
calc_A <- function(SEMfitted){
  n <- nobs(SEMfitted)
  tmpvc <- n * vcov(SEMfitted, remove.duplicated = TRUE)
  chol2inv(chol(tmpvc))
}

#' Calculate Score Contribution Matrix (Internal)
#'
#' For a fitted SEM, computes score contribution matrix (sc). 
#'
#' @param SEMfitted A lavaan-fitted SEM
#' @return SC Matrix 
#' @keywords internal
calc_sc <- function(SEMfitted){
  sandwich::estfun(SEMfitted, remove.duplicated=TRUE)
}

#' Calculate Matrix B (Internal)
#'
#' For pairs of SEMS, compute the sample adjusted cross-product B information matrix from Vuong (1989).
#'
#' @param sc1 Mdoel 1 SC Matrix
#' @param sc2 Model 2 SC Matrix
#' @param n Number of observations
#' @return Matrix B
#' @keywords internal
calc_B <- function(sc1, sc2, n){
  crossprod(sc1, sc2)/n
}
