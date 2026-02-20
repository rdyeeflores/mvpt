#' mvpt() Follow-up: Rank Models by Path Value 
#' 
#' Follow-up function to rank and view all models by path value.
#' 
#' @param mvpt_output Output from using mvpt() function.
#' @return An object of class...
#' @examples
#' \dontrun{
#' mvpt_output <- mvpt(lavaan_model, path, data)
#' mvptRank(mvpt_output)}
#' @export
mvptRank <- function(mvpt_output){
  
  CORE_comp <- mvpt_output[[3]]
  vals <- CORE_comp$sharedparamvals
  groups <- split(names(vals), vals)
  ranked_groups <- groups[order(as.numeric(names(groups)))]
  
  out <- unname(cbind(
    as.numeric(names(ranked_groups)),
    vapply(groups, paste, collapse = ", ", character(1))
  ))
  out
  
}
