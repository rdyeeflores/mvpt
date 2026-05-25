#' Print User-facing Results from mvpt() Object
#'
#' Displays a compact summary of an object returned by \code{\link{mvpt}}.
#'
#' @param x An object of class \code{mvpt}.
#' @param ... Unused.
#' @return The object \code{x}, invisibly.
#' @export
print.mvpt <- function(x, ...) { 
  
  FIGURE_list <- x[[1]]
  CORE_comp <- x[[3]]
  showplots <- x[[4]]
  reversal <- x[[5]]
  
  ## CORE components 
  M <- CORE_comp$M
  n <- CORE_comp$n
  path <- CORE_comp$path
  sharedparamvals <- round(CORE_comp$sharedparamvals, 3)
  CHI.sq  <- CORE_comp$CHI.sq
  df <- CORE_comp$df
  p.val <- CORE_comp$p.val           
  
  ## Message  
  cat(sprintf("Including the given model, %d SEMs were compared \nin an MVPT using a dataset of %d observations. \nTest results across these models were based on the \nspecified path: %s\n", M, n, path))
  cat("\n")
  cat("MVPT Results\n")
  cat("------------------\n")
  cat(sprintf("Overall chi-square = %.2f, df = %.0f        (p = %.3f)\n\n", CHI.sq, df, p.val))
  if(reversal){
    cat("Path values per model (standardized data):\n")
  }else{
    cat("Path values per model:\n")  
  }
  
  print(sharedparamvals)
  
  ## ON/OFF switch for grid figure
  if(showplots){
    
    ## FIGURE list into grid
    wrap <- patchwork::wrap_plots(FIGURE_list, ncol = 3)
    print(wrap)
    
  }
  
}

