#' SEM Auto-fitting (Internal)
#'
#' Takes a lavaan-ready list of models and fits all with same basic settings, then produces a fitted list. Not passing sem() estimator and missing arguments yet (maximum likelihood default).
#'
#' @param subMEC_lavaan_ready A list of pre-fitted SEMs in lavaan syntax
#' @param data A data frame to fit the given SEM 
#' @return A lavaan-fitted model list
#' @keywords internal
auto_sem <- function(subMEC_lavaan_ready, data, missing = "ml", estimator = "ML"){
  
  # STOP: Enforce missing and estimator defaults
  if (tolower(missing) != "ml" || tolower(estimator) != "ml") {
    stop("The default missing = 'ml' and estimator = 'ML' arguments cannot currently be changed.")
  }
  
  fit_list <- list()
  for (i in 1:length(subMEC_lavaan_ready)) {
    fit_list[[i]] <- sem(model = subMEC_lavaan_ready[[i]], 
                         data = data, 
                         fixed.x = FALSE, 
                         missing = missing,
                         estimator = estimator) ## avoiding drops bc sandwich::estfun()
  }
  fit_list  
}

