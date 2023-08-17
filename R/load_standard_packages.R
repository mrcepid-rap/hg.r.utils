#' @import data.table
#' @import ggplot2
#' @import tidyr
#' @import stringr
#' @import patchwork
#' @import lemon
#' @import broom
load_standard_packages <- function() {

  requireNamespace(data.table)
  requireNamespace(ggplot2)
  requireNamespace(tidyr)
  requireNamespace(stringr)
  requireNamespace(patchwork)
  requireNamespace(lemon)
  requireNamespace(broom)

}
