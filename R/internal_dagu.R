#' DAG Utility for Model Auto-Generation (Internal)
#'
#' Takes one model and path in lavaan syntax, auto-generates the MEC, then gets subMEC_lavaan_ready.
#'
#' @param LAV_list One or more lavaan syntax models in a list 
#' @param path A lavaan syntax path
#' @param reversal Option to allow reversal of given path
#' @return A model list in dagitty and another in lavaan format
#' @keywords internal
dagu <- function(LAV_list, path, reversal = FALSE){
  
  ## Function separating LAV into regressions and LVs, with former turned to DAGs for FIGURE_list
  ## NOTE: Avoids flipping arrows emanating from LVs when producing MEC
  split_syntax <- function(LAV) {
    lines <- trimws(strsplit(LAV, "\n", fixed = TRUE)[[1]])
    lines <- lines[nzchar(lines) & !grepl("^#", lines)]  ## drop blanks + full-line comments
    LVs <- lines[grepl("=~", lines, fixed = TRUE)]
    regs <- lines[grepl("~",  lines) & !grepl("=~", lines, fixed = TRUE)]  ## "~" but not "=~"
    list(regs=regs, LVs=LVs)
  }
  
  ## Function converting LAV_regressions to DAG format
  ## NOTE: Must manually remove dagitty-automated placement of exogenous covariances
  regs2DAG <- function(regs){
    DAG <- lavaanToGraph(regs)
    unidirect <- edges(DAG)[ !edges(DAG)$e == "<->", ]
    lines <- apply(unidirect, 1, function(x) paste(x[1], "->", x[2]))
    dagtext <- paste0("dag {\n", paste(lines, collapse = "\n"), "\n}")
    DAG <- dagitty(dagtext)
    ## Giving coordinates
    DAG <- graphLayout(DAG) 
    DAG
  
  }
  
  ## Taking LAV_list to make a list of DAGs from each LAV's set of regression equations
  split_syntax_list <- lapply(LAV_list, split_syntax)
  regs_list <- lapply(split_syntax_list, function(x) x[[1]])
  regs2DAG_list <- lapply(regs_list, regs2DAG)
  
  ## Using regs2DAG_list to get fam and fam_lavaan_ready
  ## NOTE: Fork assumes regs2DAG_list with one user-given model ready for auto-generation
  if(length(regs2DAG_list) == 1){

    ## YES auto-generation 
    DAG <- regs2DAG_list[[1]]
    
    ## Getting MEC and indexing to get non-reversed MEC subgroup
    ## NOTE: Turning path's "~-notion" into parts for better indexing
    MEC <- equivalentDAGs(DAG)
    parts <- strsplit(path, "~")[[1]] 
    outcome_var <- trimws(parts[1])
    predictor_var <- trimws(parts[2])
    IND <- vector()
    for (i in 1:length(MEC)) {
      IND[i] <- sum(edges(MEC[[i]])$v == predictor_var & edges(MEC[[i]])$w == outcome_var) == 1
    }
    subMEC <- MEC[IND]  

    ## Family dependent on reversal choice
    ## NOTE: Not using fam <- MEC so the non-reversed models appear first
    if(reversal){
      fam <- c(subMEC, MEC[!IND])
    }else{
      fam <- subMEC
    }
    
    ## Turning fam list into lavaan syntax list, and adding the same previously split-off LVs (if any)
    fam_lavaan_ready <- list()
    for (i in 1:length(fam)) {
      fam_lavaan_ready[[i]] <- paste(c(convert(fam[[i]], to="lavaan"), split_syntax_list[[1]]$LVs), collapse = "\n")
    }
    
  }else{
    
    ## NO auto-generation (but still getting regs2DAG_list for FIGURE_list later)
    fam <- regs2DAG_list
    fam_lavaan_ready <- LAV_list
    
  }
  
  ## MESSAGE: Early exit bc fam contains only one model; nothing to compare, nothing to compute
  if (length(fam_lavaan_ready) < 2) {
    message("An MVP Test could not be computed using the provided lavaan model stynax \nand path. No other models could be generated for comparison, thus halting \ncomputation of a test statistic.")
    return(invisible(NULL)) 
  }
  
  list(fam=fam, fam_lavaan_ready=fam_lavaan_ready)
}
