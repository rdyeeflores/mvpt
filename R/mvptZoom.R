#' mvpt() Follow-up: View Single Model
#' 
#' Follow-up function to see figure and lavaan output for a single SEM.
#' 
#' @param mvpt_output Output from mvpt().
#' @param index Model index for list of compared model, both given and auto-generated. Model 1 is always the given model.
#' @examples
#' \dontrun{
#' mvpt_output <- mvpt(lavaan_model, path, data)
#' mvptZoom(mvpt_output, index = 5)}
#' @export
mvptZoom <- function(mvpt_output, index){
  
  SEMfitted_list <- mvpt_output[[2]]
  FIGURE_list <- mvpt_output[[1]]
  
  list(summary(SEMfitted_list[[index]], standardized = TRUE),
       FIGURE_list[[index]]
  )
  
}
