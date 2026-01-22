#' mvpt: SEM Multiverse Path Test
#'
#' Sensitivity analysis for individual SEM path estimates 
#' across multiple automatically generated SEMs.
#'
#' @details
#' Main functions are \code{\link{mvpt}},
#' \code{\link{mvptZoom}}, and \code{\link{mvptRank}}.
#'
"_PACKAGE"

#' @importFrom lavaan sem lavInspect nobs vcov coef
#' @importFrom dagitty dagitty lavaanToGraph edges graphLayout equivalentDAGs convert
#' @importFrom ggplot2 ggplot aes
#' @importFrom stats na.omit pchisq
#' @importFrom grid unit
NULL

###############################################################################################
#### COMPONENT FUNCTIONS
##############################################################################################

#' Model Auto-generation DAG Utility (Internal)
#'
#' Takes one lavaan syntax model and path, auto-generates MEC, then gets subMEC_lavaan_ready.
#'
#' @param LAV A lavaan syntax model
#' @param path A lavaan syntax path
#' @return A model list in dagitty and another in lavaan format
#' @keywords internal
dagu <- function(LAV, path){
  
  ## Separation of LAV into regressions and LVs 
  ## NOTE: Avoids flipping arrows emanating from LVs when producing MEC
  split_syntax <- function(LAV) {
    lines <- trimws(strsplit(LAV, "\n", fixed = TRUE)[[1]])
    lines <- lines[nzchar(lines) & !grepl("^#", lines)]  # drop blanks + full-line comments
    LVs <- lines[grepl("=~", lines, fixed = TRUE)]
    reg <- lines[grepl("~",  lines) & !grepl("=~", lines, fixed = TRUE)]  # "~" but not "=~"
    list(LAV_regressions = reg, LAV_LVs = LVs)
  }
  LAV_regressions <- split_syntax(LAV)[[1]]
  LAV_LVs <- split_syntax(LAV)[[2]]
  
  ## Converting LAV_regressions to DAG format
  ## NOTE: Must manually remove dagitty-automated placement of exogenous covariances
  DAG <- lavaanToGraph(LAV_regressions)
  unidirect <- edges(DAG)[ !edges(DAG)$e == "<->", ]
  lines <- apply(unidirect, 1, function(x) paste(x[1], "->", x[2]))
  dagtext <- paste0("dag {\n", paste(lines, collapse="\n"), "\n}")
  DAG <- dagitty(dagtext)
  
  ## Giving coordinates, creating MEC
  DAG <- graphLayout(DAG) 
  MEC <- equivalentDAGs(DAG)
  
  ## STOP: Exiting if user included covariances (update PENDING)
  if (grepl("~~", LAV)) {
    stop("Covariance specifications (`~~`) are not allowed in this function.")
  }
  
  ## STOP: Exiting if user specifies an SEM or path with labels 
  ## NOTE: Using "\\*"  because "*" does not work at identifying models with labels 
  if ( grepl("\\*", LAV) ) {
    stop("Model syntax includes parameter labels, which are automatically removed by dagitty::lavaanToGraph. Please remove all labels from your model and specified path before running mvpt().")
  }
  if (!grepl("~", path)) {
    stop("For the path speficied, labeling is not allowed. Please specify using tilde-notation instead (e.g., Y ~ X).")
  }
  
  ## Turning path's "~-notion" into parts for better indexing in MEC
  parts <- strsplit(path, "~")[[1]]
  outcome_var <- trimws(parts[1])
  predictor_var <- trimws(parts[2])
  
  ## Looping with index to get subMEC (member models with same user-specified path) 
  IND <- vector()
  for (i in 1:length(MEC)) {
    IND[i] <- sum(edges(MEC[[i]])$v == predictor_var & edges(MEC[[i]])$w == outcome_var) == 1
  }
  subMEC <- MEC[IND]  
  
  ## MESSAGE: Early exit bc subMEC contains only one model; nothing to compare, nothing to compute
  if (length(subMEC) < 2) {
    message("An MVP Test could not be computed using the provided lavaan model stynax \nand path. No other models could be generated for comparison, thus halting \ncomputation of a test statistic.")
    return(invisible(NULL)) 
  }
  
  ## Turning subMEC list into lavaan syntax list, adding any previous split LVs too
  LAV_regs_multi <- list()
  subMEC_lavaan_ready <- list()
  for (i in 1:length(subMEC)) {
    LAV_regs_multi[[i]] <- convert(subMEC[[i]], to="lavaan")
    subMEC_lavaan_ready[[i]] <- paste(c(LAV_regs_multi[[i]], LAV_LVs), collapse = "\n")
  }
  
  list(subMEC=subMEC, subMEC_lavaan_ready=subMEC_lavaan_ready)
}


#' Model Auto-fitting (Internal)
#'
#' Takes lavaan ready list and fits all with same basic settings, then produces fitted list. Not passing sem() arguments yet, and assumes complete data.
#'
#' @param subMEC_lavaan_ready A list of pre-fitted SEMs in lavaan syntax
#' @param data A data frame to fit the given SEM 
#' @param na.action SEM option
#' @param subset SEM option
#' @param varcov SEM option
#' @param n SEM option
#' @return A lavaan-fitted model list
#' @keywords internal
auto_sem <- function(subMEC_lavaan_ready, data, na.action=na.omit, subset=NULL, varcov = NULL, n = NULL){
  
  ## STOP: Exiting if user specified LVs or covariances (update PENDING)
  if (anyNA(data)) {
    stop("Only complete datasets (no NAs) are allowed. Please account for this missing data before running mvpt().")
  }
  
  fit_list <- list()
  for (i in 1:length(subMEC_lavaan_ready)) {
    fit_list[[i]] <- sem(model=subMEC_lavaan_ready[[i]], data=data, fixed.x=FALSE)
  }
  fit_list  
  
}


#' Calculate Matrix A (Internal)
#'
#' For a fitted SEM, computes the A information matrix from Vuong (1989).
#'
#' @param SEMfitted A lavaan-fitted SEM
#' @return Matrix A
#' @keywords internal
calc_A <- function(SEMfitted){
  n <- nobs(SEMfitted)
  tmpvc <- n * vcov(SEMfitted, remove.duplicated = TRUE)
  chol2inv(chol(tmpvc))
}

#' Calculate Score Contribution Matrix (Internal)
#'
#' For a fitted SEM, computes score contribution matrix (sc). 
#'
#' @param SEMfitted A lavaan-fitted SEM
#' @return SC Matrix 
#' @keywords internal
calc_sc <- function(SEMfitted){
  sandwich::estfun(SEMfitted, remove.duplicated=TRUE)
}

#' Calculate Matrix B (Internal)
#'
#' For pairs of SEMS, compute the sample adjusted cross-product B information matrix from Vuong (1989).
#'
#' @param sc1 Mdoel 1 SC Matrix
#' @param sc2 Model 2 SC Matrix
#' @param n Number of observations
#' @return Matrix B
#' @keywords internal
calc_B <- function(sc1, sc2, n){
  crossprod(sc1, sc2)/n
}



#' Calculate Vuong-Wald Core MVP Test Statistic (Internal)
#'
#' For a list of fitted SEMs and shared path, a multivariate sampling distribution, is used to compare values for a path shared by these SEMs.
#'
#' @param SEMfitted_list A list of lavaan-fitted SEMs
#' @param path A lavaan syntax path
#' @return List of lists and more test components
#' @keywords internal
VW_core <- function(SEMfitted_list, path){
  
  ## Number of observations and models
  n <- nobs(SEMfitted_list[[1]])
  M <- length(SEMfitted_list)
  
  ## WARNING: Models with negative variances
  all_pos_param <- vector(length=M)
  for (i in 1:M) {
    all_pos_param[i] <- lavInspect(SEMfitted_list[[i]], "post.check") 
  }
  if(all(all_pos_param)==FALSE){
    warning("WARNING: Negative variance computed. Interpret MVP test results with caution.")
  }
  
  ## Gathering needed components per model by using helper functions 
  A <- vector(mode = "list", length=M)
  sc  <- vector(mode = "list", length=M)
  for (i in 1:M) {
    A[[i]]  <- calc_A(SEMfitted_list[[i]])
    sc[[i]] <- calc_sc(SEMfitted_list[[i]])
  }
  
  ## Pasting smaller matrices together into SIGMA
  prev_cols <- NULL
  prev_rows <- NULL
  for (i in 1:M) { ## first row to last
    for (j in 1:M) { ## first col to last
      curr_cols <- solve(A[[i]]) %*% calc_B(sc[[i]], sc[[j]], n) %*% solve(A[[j]])
      prev_cols <- cbind(prev_cols, curr_cols)
    }
    curr_rows <- prev_cols
    prev_rows <- rbind(prev_rows, curr_rows)
    prev_cols <- NULL
  }
  SIGMA <- prev_rows
  
  ## Computing a Wald CHI-square test 
  THETA  <- NULL
  for (i in 1:M) {THETA <- c(THETA, coef(SEMfitted_list[[i]]))}
  H      <- clubSandwich::constrain_equal(as.character(path), THETA, reg_ex = TRUE) 
  df     <- nrow(H)
  HSH <- H %*% (SIGMA/n) %*% t(H)
  lambda <- 1e-8   # adjust as needed to help regularize matrix inversion
  HSH_regu <- HSH + lambda * diag(nrow(HSH))
  CHI.sq <- t(H %*% THETA) %*% solve(HSH_regu) %*% (H %*% THETA) 
  p.val  <- pchisq(CHI.sq, df, lower.tail=FALSE)
  
  ## Collecting shared path values from each model in a new vector
  sharedparamvals <- THETA[names(THETA) == path]
  for (i in seq_along(sharedparamvals)) {
    names(sharedparamvals)[i] <- paste0("M", seq(length(sharedparamvals)))[i]
  }

  ## OUTPUT: List of test components
  list(M=M, n=n, path=path, sharedparamvals=sharedparamvals, CHI.sq=CHI.sq, p.val=p.val, SIGMA=SIGMA, H=H, THETA=THETA)
  
}




#########################################################################
#### MAIN FUNCTIONS
#########################################################################

#' Run MVP Test
#' 
#' Given a model and path in lavaan syntax, plus data frame, this function auto-generates other models with the same path, fits each one, and tests for path values differences.
#' 
#' @param lavaan_model A pre-fitted SEM in lavaan syntax. 
#' @param path The path to be tested with the given SEM. This must also be in lavaan syntax (eg: Y~X).
#' @param data A data frame to fit the given SEM, and all others that may be auto-generated. 
#' @param showplots Whether to show a figure containing all the compared SEMs (The given model is always first). This figure is free of path values and meant for quick inspection of model specification differences. 
#' @param na.action lavaan option
#' @param subset lavaan option
#' @param varcov lavaan option
#' @param n lavaan option
#' @return An object of class \code{mvpt} containing lavaan-fitted results, figures, and MVP test components. 
#' @examples
#' \dontrun{
#' data(PoliticalDemocracy, package = "lavaan")
#' lavaan_model <- 
#'   "
#'   ## regressions
#'   dem60 ~ ind60
#'   dem65 ~ ind60 + dem60
#'   ## latent variables
#'   ind60 =~ x1 + x2 + x3
#'   dem60 =~ y1 + a*y2 + b*y3 + c*y4
#'   dem65 =~ y5 + a*y6 + b*y7 + c*y8
#'   "
#' path <- "dem60 ~ ind60"
#' data <- PoliticalDemocracy
#' output <- mvpt(lavaan_model, path, data)
#' output}
#' @export
mvpt <- function(lavaan_model, path, data, 
                 showplots=FALSE,
                 na.action=na.omit, subset=NULL, varcov = NULL, n = NULL){
  
  ## STOP: User inputs are in wrong format
  fit_try <- try(sem(lavaan_model, data=data, do.fit=FALSE), silent = TRUE)
  if (inherits(fit_try, "try-error")) {
    stop("Invalid lavaan model syntax:\n", fit_try)
  }
  if (!is.character(path) || length(strsplit(path, "~")[[1]]) != 2) {
    stop("The path must be a single character string like 'y ~ x'.")
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


#' View Single Model
#' 
#' Follow-up function to see figure and lavaan results for a single model.
#' 
#' @param mvpt_output Output from using mvpt() function.
#' @param index Model index for list of compared model, both given and auto-generated. Model 1 is always the given model.
#' @return An object of class...
#' @examples
#' \dontrun{
#' mvpt_output <- mvpt(lavaan_model, path, data)
#' mvptZoom(mvpt_output, index = 5)}
#' @export
mvptZoom <- function(mvpt_output, index){
  
  FIGURE_list <- mvpt_output[[1]]
  SEMfitted_list <- mvpt_output[[2]]
  
  list(FIGURE_list[[index]], 
       summary(SEMfitted_list[[index]], standardized = TRUE))
  
}



#' Rank Models by Path Value 
#' 
#' Follow-up function to rank and view models with largest path values (updated rank by absolute value pending).
#' 
#' @param mvpt_output Output from using mvpt() function.
#' @param top Number of top models to include in report.  
#' @return An object of class...
#' @examples
#' \dontrun{
#' mvpt_output <- mvpt(lavaan_model, path, data)
#' mvptRank(mvpt_output, top = 3)}
#' @export
mvptRank <- function(mvpt_output, top){
  
  CORE_comp <- mvpt_output[[3]]
  x <- CORE_comp$sharedparamvals
  x_top <- sort(x, decreasing = TRUE)[seq_len(top)]
  x_top
  
}


#' Print an mvpt() object
#'
#' Displays a compact summary of an object returned by \code{\link{mvpt}}.
#'
#' @param x An object of class \code{mvpt}.
#' @param ... Unused.
#' @return The object \code{x}, invisibly.
#' @export
print.mvpt <- function(x, ...) { 
  
  FIGURE_list <- x[[1]]
  CORE_comp <- x[[3]]
  showplots <- x[[4]]
  
  ## CORE components 
  M <- CORE_comp$M
  n <- CORE_comp$n
  path <- CORE_comp$path
  sharedparamvals <- CORE_comp$sharedparamvals
  CHI.sq  <- CORE_comp$CHI.sq            
  p.val <- CORE_comp$p.val           
  
  ## Message  
  cat(sprintf("Including the given model, %d fitted SEMs were compared in a MVP test \nusing a dataset of %d observations. Test results across these models \nare based on the shared path: %s\n", M, n, path))
  cat("\n")
  cat("MVP Test Results\n")
  cat("------------------\n")
  cat(sprintf("Overall chi-square = %.2f     (p = %.3f)\n\n", CHI.sq, p.val))
  cat("Model path values:\n")
  print(sharedparamvals)
  
  ## ON/OFF switch for grid figure
  if(showplots){
    
    ## FIGURE list into grid
    wrap <- patchwork::wrap_plots(FIGURE_list, ncol = 3)
    print(wrap)
    
  }
  
}



