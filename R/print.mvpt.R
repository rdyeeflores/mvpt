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
  
  ## CORE components 
  M <- CORE_comp$M
  n <- CORE_comp$n
  path <- CORE_comp$path
  sharedparamvals <- CORE_comp$sharedparamvals
  CHI.sq  <- CORE_comp$CHI.sq            
  p.val <- CORE_comp$p.val           
  
  ## Message  
  cat(sprintf("Including the given model, %d fitted SEMs were compared \nin a MVP test using a dataset of %d observations. \nTest results across these models are based on the \nshared path: %s\n", M, n, path))
  cat("\n")
  cat("MVP Test Results\n")
  cat("------------------\n")
  cat(sprintf("Overall chi-square = %.2f     (p = %.3f)\n\n", CHI.sq, p.val))
  cat("Model path values:\n")
  print(sharedparamvals)
  
  ## ON/OFF switch for grid figure
  if(showplots){
    
    ## FIGURE list into grid
    wrap <- patchwork::wrap_plots(FIGURE_list, ncol = 3)
    print(wrap)
    
  }
  
}

