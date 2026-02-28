#' Run MVP Test
#' 
#' Using the model specification of a SEM and single path within that SEM, both in \pkg{lavaan} syntax, this function auto-generates multiple other models with the same single path by only using the graphical features of the given SEM. Auto-generated models will share the same fit statistics as the given SEM (Verma & Pearl, 1991), though suggest differing relationships between variables. After auto-generation, this function then uses the supplied data to fit all models using the same settings (limited to maximum likelihood estimation for now), followed by a chi-square test across models for significant value changes in the specified path. The given model will always be indexed first and appear as "M1" in the output.
#' 
#' @param lavaan_model A SEM in lavaan syntax. 
#' @param path The path to be tested within the given SEM. This must also be in lavaan syntax (eg: Y~X).
#' @param data A data frame to fit the given SEM, and all others that may be auto-generated. 
#' @param showplots Whether to show a parameter-free figure containing all the compared SEMs, with the user-supplied model always first. This figure is for quick inspection of model specification differences. 
#' @return An object of class \code{mvpt} containing lavaan-fitted results, figures, and MVP test components. 
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
#' mvpt <- mvpt(lavaan_model, path, data = PoliticalDemocracy, showplots = TRUE)
#' mvpt}
#' @export
mvpt <- function(lavaan_model, 
                 path, 
                 data, 
                 showplots = FALSE){
  
  ## STOP: User inputs are in wrong format
  fit_try <- try(sem(lavaan_model, data=data, do.fit=FALSE), silent = TRUE)
  if (inherits(fit_try, "try-error")) {
    stop("Invalid lavaan model syntax:\n", fit_try)
  }
  if (!is.character(path) || length(strsplit(path, "~")[[1]]) != 2) {
    stop("The path must be a single character string like 'Y ~ X'.")
  }
  
  dagu <- dagu(lavaan_model, path)
  if(is.null(dagu)){ return(invisible(NULL)) }
  
  ## FIGURE_list to house all model ggdag figures
  ## NOTE: Tailoring DAG figures via ggplot instead of ggdags
  subMEC <- dagu$subMEC
  FIGURE_list <- list()
  for (i in 1:length(subMEC)) {
    td <- ggdag::tidy_dagitty(subMEC[[i]])  # node coords + edge list
    FIGURE_list[[i]] <- ggplot(td, aes(x = x, y = y, xend = xend, yend = yend)) +
      ggdag::geom_dag_edges(
        start_cap = ggraph::circle(4, "mm"),
        end_cap   = ggraph::circle(4, "mm"),
        arrow_directed = grid::arrow(length = unit(2, "mm"), type = "closed")
      ) + 
      ggdag::geom_dag_text(aes(label = name), size = 3, colour = "deepskyblue3") + 
      ggdag::theme_dag()
  }
  
  ## SEMfitted_list to house all sem() results
  subMEC_lavaan_ready <- dagu$subMEC_lavaan_ready
  SEMfitted_list <- auto_sem(subMEC_lavaan_ready, data)
  
  ## CORE test components
  CORE_comp <- VW_core(SEMfitted_list, path)
  
  ## OUTPUT: Via print() method 
  output <- list(FIGURE_list=FIGURE_list, SEMfitted_list=SEMfitted_list, CORE_comp=CORE_comp, showplots=showplots)
  class(output) <- "mvpt"
  output
  
}


