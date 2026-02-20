#' DAG Utility for Model Auto-generation (Internal)
#'
#' Takes one model and path in lavaan syntax, auto-generates the MEC, then gets subMEC_lavaan_ready.
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
  dagtext <- paste0("dag {\n", paste(lines, collapse = "\n"), "\n}")
  DAG <- dagitty(dagtext)
  
  ## Giving coordinates, creating MEC
  DAG <- graphLayout(DAG) 
  MEC <- equivalentDAGs(DAG)
  
  ## STOP: Exiting if user included covariances (update PENDING)
  if (grepl("~~", LAV)) {
    stop("Covariance parameters (`~~`) are not allowed in this version of mvpt().")
  }
  
  ## STOP: Exiting if user specifies an SEM with labels (no need to repeat for path) 
  ## NOTE: Using "\\*"  because "*" does not work at identifying models with labels 
  if ( grepl("\\*", LAV) ) {
    stop("Model syntax includes parameter labels, which are automatically removed by dagitty::lavaanToGraph. Before running mvpt(), please remove all labels from your lavaan-formatted model and use tilde-notation for the path argument (e.g., Y ~ X).")
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
