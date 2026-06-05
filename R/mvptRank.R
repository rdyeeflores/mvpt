#' mvpt() Follow-up: Rank Models by Path Value 
#' 
#' Follow-up function to rank and view all models by path value.
#' 
#' @param mvpt_output Output from mvpt().
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
#' mvptRank(mvpt_path1)
#' }
#' @export
mvptRank <- function(mvpt_output){
  
  ## STOP: Input is not mvpt object
  if (!inherits(mvpt_output, "mvpt")) {
    stop("'mvpt_output' must be an object of class 'mvpt'.")
  }
  
  CORE_comp <- mvpt_output$CORE_comp
  
  vals <- CORE_comp$sharedparamvals
  groups <- split(names(vals), vals)
  ranked_groups <- groups[order(as.numeric(names(groups)))]
  
  out <- data.frame(
    "Path values" = as.numeric(names(ranked_groups)),
    "Models" = vapply(ranked_groups, paste, collapse = ", ", character(1)),
    check.names = FALSE
  )
  
  rownames(out) <- NULL
  out
}
