#' mvpt() Follow-up: Rank Models by Path Value 
#' 
#' Follow-up function to rank and view all models by path value.
#' 
#' @param mvpt_output Output from mvpt().
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
