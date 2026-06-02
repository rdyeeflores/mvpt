#' mvpt() Follow-up: View Single Model
#' 
#' Follow-up function to see figure and lavaan output for a single SEM.
#' 
#' @param mvpt_output Output from mvpt().
#' @param index Model index for list of compared models. Model 1 is always the given model.
#' @examples
#' \dontrun{
#' data(PoliticalDemocracy, package = "lavaan")
#' 
#' lavaan_model <- 
#'   "
#'   ## latent variables
#'   ind60 =~ x1 + x2 + x3
#'   dem60 =~ y1 + y2 + y3 + y4
#'   dem65 =~ y5 + y6 + y7 + y8
#'   ## regressions
#'   dem60 ~ ind60
#'   dem65 ~ ind60 + dem60
#'   "
#'   
#' path <- "dem60 ~ ind60"
#' 
#' mvpt_output <- mvpt(lavaan_model, path, data = PoliticalDemocracy)
#' mvptZoom(mvpt_output, index = 5)}
#' @export
mvptZoom <- function(mvpt_output, index){
  
  if (!inherits(mvpt_output, "mvpt")) {
    stop("'mvpt_output' must be an object of class 'mvpt'.")
  }
  
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
  
  ##list(fam_lavaan_ready[[index]],
  ##     summary(SEMfitted_list[[index]], standardized = TRUE, fit.measures = TRUE),
    ##   FIGURE_list[[index]]
  ##)
  
  cat(sprintf("\nModel %d Specification\n---------------------\n", index))
  #cat(unlist(fam_lavaan_ready[[index]]))
  cat(sprintf('"\n%s\n"\n', unlist(fam_lavaan_ready[[index]])))
  cat(sprintf("\nModel %d Summary\n---------------\n", index))
  print(summary(SEMfitted_list[[index]], standardized = TRUE, fit.measures = TRUE))
  print(FIGURE_list[[index]])
  
  
}
