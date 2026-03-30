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
#' path <- "dem60 ~ ind60"
#' mvpt_output <- mvpt(lavaan_model, path, data = PoliticalDemocracy)
#' mvptZoom(mvpt_output, index = 5)}
#' @export
mvptZoom <- function(mvpt_output, index){
  
  SEMfitted_list <- mvpt_output[[2]]
  FIGURE_list <- mvpt_output[[1]]
  
  list(summary(SEMfitted_list[[index]], standardized = TRUE, fit.measures = TRUE),
       FIGURE_list[[index]]
  )
  
}
