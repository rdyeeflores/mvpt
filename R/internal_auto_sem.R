#' SEM Auto-fitting (Internal)
#'
#' Takes a lavaan-ready list of models and fits all with same basic settings, then produces a fitted list. Not passing sem() estimator and missing arguments yet (maximum likelihood default).
#'
#' @param fam_lavaan_ready A list of pre-fitted SEMs in lavaan syntax
#' @param data A data frame to fit the given SEM 
#' @return A lavaan-fitted model list
#' @keywords internal
auto_sem <- function(fam_lavaan_ready, data, missing = "ml", estimator = "ML"){
  
  # STOP: Enforce missing and estimator defaults
  # NOTE: Avoiding listwise deletion because of mismatch with sandwich::estfun()
  if (tolower(missing) != "ml" || tolower(estimator) != "ml") {
    stop("The default missing = 'ml' and estimator = 'ML' arguments cannot currently be changed.")
  }
  
  ## NOTE: Indexing message and warning, BUT stopping upon fit fail 
  fit_list <- list()
  failed_models <- integer()
  for (i in seq_along(fam_lavaan_ready)) {
    
    fit_list[[i]] <- withCallingHandlers(
      sem(model = fam_lavaan_ready[[i]], 
          data = data, 
          fixed.x = FALSE,
          missing = missing,
          estimator = estimator,
          check.gradient = FALSE),
      
      ## Warnings and messages indexed by model
      message = function(m) {
        message("[Model ", i, "] ", conditionMessage(m))
        invokeRestart("muffleMessage")},
      warning = function(w) {
        warning("[Model ", i, "] ", conditionMessage(w), call. = FALSE)
        invokeRestart("muffleWarning")}
    )
    
    if (!isTRUE(lavaan::lavInspect(fit_list[[i]], "converged"))) {
      failed_models <- c(failed_models, i)
    }
    
  }
  
  ## STOP: At least one model failed to converge
  if (length(failed_models) > 0) {
    stop(
      "An MVPT could not be computed because Model(s) ",
      paste(failed_models, collapse = ", "),
      " failed to converge.",
      call. = FALSE
    )
  }
  
  fit_list  
}

