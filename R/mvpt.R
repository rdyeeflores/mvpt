#' Run MVP Test
#' 
#' Using the model specification of a SEM and single path within that SEM, both in \pkg{lavaan} syntax, this function auto-generates multiple other models with the same single path by only using the graphical features of the given SEM. Auto-generated models will share the same fit statistics as the given SEM (Verma & Pearl, 1991), though suggest differing relationships between variables. After auto-generation, this function then uses the supplied data to fit all models using the same settings (limited to maximum likelihood estimation for now), followed by a chi-square test across models for significant value changes in the specified path. The given model will always be indexed first and appear as "M1" in the output.
#' 
#' @param lavaan_model An SEM in lavaan syntax. 
#' @param path The path to be tested within the given SEM. This must also be in lavaan syntax (eg: "Y ~ X").
#' @param data A data frame to fit the given SEM, and all other SEMs that may be auto-generated. 
#' @param showplots Whether to show a parameter-free figure containing all the compared SEMs, with the user-supplied model always first (top-left). This figure is for quick inspection of model specification differences. 
#' @param reversal Whether to include auto-generated SEMs in which the given path is reversed.
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
#' mvpt_output <- mvpt(lavaan_model, path, data = PoliticalDemocracy, showplots = TRUE)
#' mvpt_output}
#' @export
mvpt <- function(lavaan_model, 
                 path, 
                 data, 
                 showplots = FALSE,
                 reversal = FALSE){
  
  ## Generalized lavaan_model for use with single model or list()
  LAV_list <- if (is.list(lavaan_model)) lavaan_model else list(lavaan_model) ## list() object no matter class
  
  ## Gathering observed variables (whether with or without LVs) 
  ## NOTE: This also catches basic syntax errors
  OVs <- unique(unlist(lapply(LAV_list, \(m) {
    v <- unique(unlist(lavaanify(m)[c("lhs","rhs")]))
    intersect(v, names(data))
  })))
  
  ## Subsetting given data to only consist of standardized observed variables
  data_std <- as.data.frame(lapply(data[OVs], \(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)))
  
  ## STOP: Model(s) generating lavaan error/warning
  for (i in seq_along(LAV_list)) {
    tryCatch(
      sem(LAV_list[[i]], data = data_std, do.fit = TRUE, fixed.x = FALSE),
      error = function(e) {
        stop(conditionMessage(e), "\nError returned by lavaan. Run lavaan_model and data in lavaan first, fix the issue reported, and only retry mvpt() once there is no error. ", call. = FALSE)},
      warning = function(w) {
        stop(conditionMessage(w), "\nWarning retuned by lavaan. Try running lavaan_model and data directly in lavaan, reviewing, and fixing before retrying mvpt().", call. = FALSE)}
      )
  }
  
  ## STOP: Model(s) include covariance parameters (update pending)
  if (any(vapply(LAV_list, function(x) grepl("~~", x), logical(1)))) {
    stop("Covariance parameters (`~~`) are not allowed in this version of mvpt().")
  }
  
  ## STOP: Model(s) include parameter labels
  ## NOTE: Using '\\*' because '*' does not work at identifying models with labels
  if (any(vapply(LAV_list, function(x) grepl("\\*", x), logical(1)))) {
    stop("Model syntax includes parameter labels, which are automatically removed by dagitty::lavaanToGraph. Before running mvpt(), please remove all labels from your lavaan-formatted model and use tilde-notation for the path argument (e.g., 'Y ~ X').")
  }
  
  ## STOP: Specified path is in not in tilde notation format
  if (!is.character(path) || length(strsplit(path, "~")[[1]]) != 2) {
    stop("Specified path must be a character string: 'Y ~ X'.")
  }
  
  ## STOP: Specified path is not in model(s)
  has_path <- function(LAV, path) {
    m <- lavaan::lavaanify(paste(LAV, collapse = "\n"), auto = FALSE, warn = FALSE)
    p <- lavaan::lavaanify(path, auto = FALSE, warn = FALSE)
    any(m$lhs == p$lhs[1] & m$op == "~" & m$rhs == p$rhs[1]) ## T/F
  }
  ok <- vapply(LAV_list, has_path, logical(1), path = path)
  if (!all(ok)) {
    stop("Specified path must be part of specified model." )
  }
  
  ## Use dagu() to get fam and fam_lavaan_ready (varies by lavaan_model class)
  dagu <- dagu(LAV_list, path, reversal)
  if(is.null(dagu)){ return(invisible(NULL)) } ## because dagu message does not exit
  
  ## FIGURE_list to house all ggdag figures (tailoring via ggplot)
  fam <- dagu$fam
  FIGURE_list <- list()
  for (i in 1:length(fam)) {
    td <- ggdag::tidy_dagitty(fam[[i]])  # node coords + edge list
    FIGURE_list[[i]] <- ggplot(td, aes(x = x, y = y, xend = xend, yend = yend)) +
      ggdag::geom_dag_edges(
        start_cap = ggraph::circle(4, "mm"),
        end_cap   = ggraph::circle(4, "mm"),
        arrow_directed = grid::arrow(length = unit(2, "mm"), type = "closed")
      ) + 
      ggdag::geom_dag_text(aes(label = name), size = 3, colour = "deepskyblue3") + 
      ggdag::theme_dag()
  }
  ##FIGURE_list[[1]] <- FIGURE_list[[1]] + 
  ##  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5))
  
  ## SEMfitted_list to house all sem() results
  fam_lavaan_ready <- dagu$fam_lavaan_ready 
  SEMfitted_list <- auto_sem(fam_lavaan_ready, data_std)
  
  ## CORE test components
  CORE_comp <- VW_core(SEMfitted_list, path, reversal)
  
  ## OUTPUT: Via print() method 
  output <- list(FIGURE_list=FIGURE_list, SEMfitted_list=SEMfitted_list, CORE_comp=CORE_comp, showplots=showplots)
  class(output) <- "mvpt"
  output
  
}


