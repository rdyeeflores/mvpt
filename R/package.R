#' mvpt: Multiverse Path Test
#'
#' Researcher can often face uncertainty when multiple SEMs equally explain the same data. This can be problematic when an effect (path) between two variables can be modeled in different ways. To help manage path uncertainty due to having multiple, competing models, the SEM Multiverse Path (MVP) Test can auto-generate competing SEMs with the same specified path and test whether this one path changes across models. Significant path value changes will be traceable to alternate model specification changes.   
#'
#' @details
#' Main functions are \code{\link{mvpt}},
#' \code{\link{mvptZoom}}, and \code{\link{mvptRank}}.
#'
"_PACKAGE"

#' @importFrom lavaan sem lavInspect nobs vcov coef summary
#' @importFrom dagitty dagitty lavaanToGraph edges graphLayout equivalentDAGs convert
#' @importFrom ggplot2 ggplot aes
#' @importFrom stats na.omit pchisq
#' @importFrom grid unit
NULL
