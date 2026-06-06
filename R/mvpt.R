#' Run MVPT
#' 
#' Perform an SEM multiverse path test.
#' 
#' @details
#' \itemize{
#'   \item This function can auto-generate other SEMs by only using the graphical features of a given SEM (Verma & Pearl, 1991). Auto-generated models will share the same fit statistics, though suggest differing relationships between variables. 
#'   \item All models are fitted using maximum likelihood estimation, followed by an MVPT (chi-square) across models to determine significant value changes in the tested path. User-provided models will always be indexed first in the output ("M1"). 
#'   \item If a user already has a set of SEMs to compare, model auto-generation can be bypassed by supplying a list() of SEMs for the lavaan_input argument.
#'   }
#' 
#' @usage
#' mvpt(lavaan_input, path, data,
#'   showplots = FALSE, MEC_only = TRUE)
#' 
#' @param lavaan_input A single SEM or list() of SEMs in lavaan syntax. For a single SEM, this function will attempt to auto-generate other models. For a list(), only the listed models will be used for testing. 
#' @param path The path to be tested. This path must be in lavaan syntax (eg: "Y ~ X") and in any associated model specification.
#' @param data A data frame to fit the associated SEMs, including auto-generated SEMs. 
#' @param showplots Whether to show a parameter-free figure containing all the compared SEMs. User-provided models start at the top-left. This figure is for quick inspection of model specification differences. 
#' @param MEC_only When given a single SEM, this determines whether model auto-generation goes beyond that model's Markov equivalence class (the space of models with equal fit metrics). In this version of mvpt(), specifying FALSE allows additionally generated models with arrows directly into the outcome variable of the tested path. 
#' @return An object of class \code{mvpt} containing lavaan-fitted results, figures, and MVPT components. 
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
#' ## path test 1
#' mvpt_path1 <- mvpt(lavaan_input, 
#'                    path = "Rumi~UnApp", 
#'                    data = UnfairApprais, 
#'                    showplots = TRUE)
#' mvpt_path1 
#' 
#' ## path test 2
#' mvpt_path2 <- mvpt(lavaan_input, 
#'                    path = "Aggr~Rumi", 
#'                    data = UnfairApprais, 
#'                    showplots = TRUE)
#' mvpt_path2
#' }
#' @export
mvpt <- function(lavaan_input, 
                 path, 
                 data, 
                 showplots = FALSE,
                 MEC_only = TRUE){
  
  ## Default options
  reversal = FALSE
  
  ## Generalized lavaan_input for use with single model or list()
  LAV_list <- if (is.list(lavaan_input)) lavaan_input else list(lavaan_input) ## list() no matter class
  
  ## STOP: Model(s) generating lavaan-based error or warning
  for (i in seq_along(LAV_list)) {
    tryCatch(
      sem(LAV_list[[i]], data = data, do.fit = TRUE), ##, auto.cov.lv.x = FALSE),
      error = function(e) {
        stop(conditionMessage(e), "\nError returned by lavaan. Run models and data in lavaan first, fix the issue reported, and retry mvpt() once there is no error.", call. = FALSE)},
      warning = function(w) {
        stop(conditionMessage(w), "\nWarning retuned by lavaan. Try running models and data directly in lavaan, reviewing, and fixing before retrying mvpt().", call. = FALSE)}
      )
  }
  
  ## STOP: Model(s) include covariance parameters (update pending)
  if (any(vapply(LAV_list, function(x) grepl("~~", x), logical(1)))) {
    stop('Covariance parameters ("~~") are not allowed in the specified model input for this version of mvpt().')
  }
  
  ## STOP: Model(s) include parameter labels
  ## NOTE: Using '\\*' because '*' does not work at identifying models with labels
  if (any(vapply(LAV_list, function(x) grepl("\\*", x), logical(1)))) {
    stop('Model syntax includes parameter labels, which are not used in mvpt(). Labels are removed during conversion from lavaan to dagitty format. Please remove all labels from your any lavaan-formatted models and use tilde-notation for the any path argument (e.g., "Y ~ X").')
  }
  
  ## STOP: Specified path is in not in correct tilde-notation
  path_no_space <- gsub("\\s+", "", path)
  if (!is.character(path) ||
      length(path) != 1 ||
      !grepl("^[A-Za-z][A-Za-z0-9._]*~[A-Za-z][A-Za-z0-9._]*$", path_no_space)) {
    stop('Specified path must be a character string in tilde-notation, with only two variables and a tilde between them (e.g., "Y ~ X").')
  }
  
  ## STOP: Specified path is not in model(s)
  has_path <- function(LAV, path) {
    m <- lavaan::lavaanify(paste(LAV, collapse = "\n"), auto = FALSE, warn = FALSE)
    p <- lavaan::lavaanify(path, auto = FALSE, warn = FALSE)
    any(m$lhs == p$lhs[1] & m$op == "~" & m$rhs == p$rhs[1]) ## T/F
  }
  ok <- vapply(LAV_list, has_path, logical(1), path = path)
  if (!all(ok)) {
    stop("Specified path missing from model. Make sure this one path is part of any associated model specification." )
  }
  
  ## Subsetting to observed variables only (whether with or without LVs) 
  ## NOTE: Additionally standardizing subset if reversal = TRUE
  OVs <- unique(unlist(lapply(LAV_list, \(m) {
    v <- unique(unlist(lavaanify(m)[c("lhs","rhs")]))
    intersect(v, names(data))
  })))
  data_sub <- data[OVs]
  if(reversal) {
    data_sub <- as.data.frame(lapply(
      data_sub, \(x) (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)
  ))} 
  
  ## Using dagu() to get fam and fam_lavaan_ready (to produce figures and lavaan-fitted models)
  dagu <- dagu(LAV_list, path, MEC_only)
  if(is.null(dagu)){ return(invisible(NULL)) } ## needed because messages dont exit
  
  ## FIGURE_list to house all ggdag figures, tailoring via ggplot
  fam <- dagu$fam
  FIGURE_list <- list()
  for (i in 1:length(fam)) {
    td <- ggdag::tidy_dagitty(fam[[i]])  # node coords + edge list
    
    FIGURE_list[[i]] <- ggplot(td, aes(x = x, y = y, xend = xend, yend = yend)) +
      ggdag::geom_dag_edges(
        start_cap = ggraph::circle(4, "mm"), ## still 4
        end_cap   = ggraph::circle(6, "mm"), ## was 4
        arrow_directed = grid::arrow(length = unit(1.5, "mm"), type = "closed") ## was 2
      ) + 
      ggdag::geom_dag_text(aes(label = substr(name, 1, 5)), 
                           size = 3, 
                           colour = "deepskyblue3") +  ## still 3
      ggdag::theme_dag()
  }
  
  ## SEMfitted_list to house all sem() results
  fam_lavaan_ready <- dagu$fam_lavaan_ready 
  SEMfitted_list <- auto_sem(fam_lavaan_ready, data_sub)
  
  ## CORE test components
  CORE_comp <- vw_core(SEMfitted_list, path)
  
  ## OUTPUT: Via print() method 
  output <- list(FIGURE_list=FIGURE_list, 
                 fam_lavaan_ready=fam_lavaan_ready,
                 SEMfitted_list=SEMfitted_list, 
                 CORE_comp=CORE_comp, 
                 showplots=showplots)
  class(output) <- "mvpt"
  output
  
}


