#' Calculate Vuong-Wald Core MVP Test Statistic (Internal)
#'
#' For a list of fitted SEMs and a shared path, a multivariate sampling distribution is used to compare path values between these SEMs.
#'
#' @param SEMfitted_list A list of lavaan-fitted SEMs
#' @param path A lavaan syntax path
#' @return List of lists and more test components
#' @keywords internal
VW_core <- function(SEMfitted_list, path){
  
  ## Number of observations and models
  n <- nobs(SEMfitted_list[[1]])
  M <- length(SEMfitted_list)
  
  ## WARNING: Models with negative variances
  all_pos_param <- vector(length=M)
  for (i in 1:M) {
    all_pos_param[i] <- lavInspect(SEMfitted_list[[i]], "post.check") 
  }
  if(all(all_pos_param)==FALSE){
    warning("WARNING: Negative variance computed. Interpret MVP test results with caution.")
  }
  
  ## Gathering needed components per model by using helper functions 
  A <- vector(mode = "list", length=M)
  sc  <- vector(mode = "list", length=M)
  for (i in 1:M) {
    A[[i]]  <- calc_A(SEMfitted_list[[i]])
    sc[[i]] <- calc_sc(SEMfitted_list[[i]])
  }
  
  ## Pasting smaller matrices together into SIGMA
  prev_cols <- NULL
  prev_rows <- NULL
  for (i in 1:M) { ## first row to last
    for (j in 1:M) { ## first col to last
      curr_cols <- solve(A[[i]]) %*% calc_B(sc[[i]], sc[[j]], n) %*% solve(A[[j]])
      prev_cols <- cbind(prev_cols, curr_cols)
    }
    curr_rows <- prev_cols
    prev_rows <- rbind(prev_rows, curr_rows)
    prev_cols <- NULL
  }
  SIGMA <- prev_rows
  
  ## Computing a Wald CHI-square test 
  THETA  <- NULL
  for (i in 1:M) {THETA <- c(THETA, coef(SEMfitted_list[[i]]))}
  path <- gsub(" ", "", path)
  H      <- clubSandwich::constrain_equal(as.character(path), THETA, reg_ex = TRUE) 
  df     <- nrow(H)
  HSH <- H %*% (SIGMA/n) %*% t(H)
  lambda <- 1e-8   # adjust as needed to help regularize matrix inversion
  HSH_regu <- HSH + lambda * diag(nrow(HSH))
  CHI.sq <- t(H %*% THETA) %*% solve(HSH_regu) %*% (H %*% THETA) 
  p.val  <- pchisq(CHI.sq, df, lower.tail=FALSE)
  
  ## Collecting shared path values from each model in a new vector
  sharedparamvals <- THETA[names(THETA) == path]
  for (i in seq_along(sharedparamvals)) {
    names(sharedparamvals)[i] <- paste0("M", seq(length(sharedparamvals)))[i]
  }
  
  ## OUTPUT: List of test components
  list(M=M, n=n, path=path, sharedparamvals=sharedparamvals, CHI.sq=CHI.sq, p.val=p.val, SIGMA=SIGMA, H=H, THETA=THETA)
  
}

