#' mvpt() Follow-up: View Single Model
#' 
#' Follow-up function to see lavaan syntax, summary, and figure for a single SEM.
#' 
#' @param mvpt_output Output from mvpt().
#' @param index Model index for list of compared models. User-provided models are always indexed first.
#' @examples
#' \dontrun{
#' 
#' data("UnfairApprais")
#' 
#' lavaan_input <-
#'   "
#'   ## latent variables
#'   Rumi =~ rumi1 + rumi2
#'   Angr =~ anger1 + anger1
#'   UnApp =~ unfair1 + unfair2
#'   ## regressions
#'   Aggr ~ Angr + Rumi
#'   Rumi ~ Angr + UnApp
#'   Angr ~ UnApp
#'   "
#' 
#' mvpt_path1 <- mvpt(lavaan_input, 
#'                    path = "Rumi~UnApp", 
#'                    data = UnfairApprais, 
#'                    showplots = TRUE)
#' mvpt_path1 
#' mvptZoom(mvpt_path1, index = 3)
#' }
#' @export
mvptZoom <- function(mvpt_output, index){
  
  ## STOP: Input is not mvpt object
  if (!inherits(mvpt_output, "mvpt")) {
    stop("'mvpt_output' must be an object of class 'mvpt'.")
  }
  
  ## STOP: Index is out of bounds or non-integer
  if (!is.numeric(index) ||
      length(index) != 1 ||
      index %% 1 != 0 ||
      index < 1 ||
      index > mvpt_output$CORE_comp$M) {
    stop("'index' must be a single integer value from 1 up to the total number of compared SEMs.")
  }
  
  fam_lavaan_ready <- mvpt_output$fam_lavaan_ready
  SEMfitted_list <- mvpt_output$SEMfitted_list
  FIGURE_list <- mvpt_output$FIGURE_list
  
  cat(sprintf("\nModel %d Specification\n---------------------\n", index))
  cat(sprintf('To replicate using the lavaan package, use this model\nin sem(model, data, fixed.x = FALSE, missing = "ML")\n'))
  
  model <- unlist(fam_lavaan_ready[[index]])
  model <- gsub("\n", "\n  ", model)
  cat(sprintf('\n  "\n  %s\n  "\n', model))
    
  cat(sprintf("\nModel %d Summary\n---------------\n", index))
  print(summary(SEMfitted_list[[index]], fit.measures = TRUE, standardized=TRUE)) 
  print(FIGURE_list[[index]])

}
