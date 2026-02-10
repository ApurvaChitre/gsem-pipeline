#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# Example: Plot a 1-factor GenomicSEM usermodel path diagram
#
# This is NOT part of the pipeline. It's a standalone example that shows:
#   1) loading MPH_genomicSEM.RData (covstruc object: mph_out)
#   2) building a 1-factor model string
#   3) fitting GenomicSEM::usermodel()
#   4) plotting a path diagram via semPlot
#
# Requirements:
#   - conda env with GenomicSEM installed
#   - R packages: GenomicSEM, lavaan, semPlot
#
# Usage:
#   Rscript examples/path_diagram_1factor_example.R \
#     --mph_rdata /path/to/MPH_genomicSEM.RData \
#     --out_png /path/to/path_1factor.png
#
# Optional:
#   --traits trait1,trait2,trait3
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(GenomicSEM)
  library(lavaan)
  library(semPlot)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NA_character_) {
  hit <- grep(paste0('^', flag, '='), args, value = TRUE)
  if (length(hit) > 0) {
    return(sub(paste0('^', flag, '='), '', hit[1]))
  }
  idx <- which(args == flag)
  if (length(idx) > 0) {
    i <- idx[1]
    if (i == length(args)) stop('Missing value for ', flag)
    return(args[i + 1])
  }
  default
}

mph_rdata <- get_arg('--mph_rdata', NA_character_)
out_png   <- get_arg('--out_png',   'path_1factor.png')
traits_in <- get_arg('--traits',    '')

if (is.na(mph_rdata) || mph_rdata == '' || !file.exists(mph_rdata)) {
  stop('Missing/invalid --mph_rdata. Example: --mph_rdata /.../MPH_genomicSEM.RData')
}

# ---- Load covstruc ----
load(mph_rdata)  # expects an object named mph_out
if (!exists('mph_out')) {
  stop('Expected MPH RData to load an object named "mph_out".\n',
       'If your object name differs, edit this script accordingly.')
}

# ---- Determine traits ----
traits <- NULL
if (nzchar(traits_in)) {
  traits <- trimws(strsplit(traits_in, ',', fixed = TRUE)[[1]])
  traits <- traits[traits != '']
}

if (is.null(traits) || length(traits) == 0) {
  if (!is.null(mph_out$trait.names) && length(mph_out$trait.names) > 0) {
    traits <- mph_out$trait.names
  } else if (!is.null(mph_out$S)) {
    traits <- rownames(as.matrix(mph_out$S))
  }
}

if (is.null(traits) || length(traits) < 2) {
  stop('Could not infer trait names from mph_out.\n',
       'Pass them explicitly with: --traits trait1,trait2,...')
}

# ---- Build 1-factor model string ----
model_1f <- paste0('F1 =~ ', paste(traits, collapse = ' + '))

# ---- Fit model ----
fit <- usermodel(
  covstruc   = mph_out,
  model      = model_1f,
  estimation = 'DWLS',
  CFIcalc    = TRUE,
  std.lv     = TRUE,
  imp_cov    = FALSE
)

# -----------------------------------------------------------------------------
# Convert a GenomicSEM usermodel() fit into a semPlotModel
# (Lightly adapted from common GenomicSEM plotting snippets.)
# -----------------------------------------------------------------------------
semPlotModel_GSEM <- function(gsem_fit, est.label = 'STD_All') {
  object <- gsem_fit$results

  # semPlot expects lavaan-style @Pars and @Vars
  object$free <- 0
  numb <- 1:length(which(object$op != '~~'))
  object$free[which(object$op != '~~')] <- numb

  varNames  <- lavaanNames(object, type = 'ov')
  factNames <- lavaanNames(object, type = 'lv')
  factNames <- factNames[!factNames %in% varNames]

  if (is.null(object$label)) object$label <- rep('', nrow(object))
  if (is.null(object$group)) object$group <- ''

  semModel <- new('semPlotModel')
  object$est <- object[, est.label]

  semModel@Pars <- data.frame(
    label = object$label,
    lhs   = ifelse(object$op == '~' | object$op == '~1', object$rhs, object$lhs),
    edge  = '--',
    rhs   = ifelse(object$op == '~' | object$op == '~1', object$lhs, object$rhs),
    est   = object$est,
    std   = NA,
    group = object$group,
    fixed = object$free == 0,
    par   = object$free,
    stringsAsFactors = FALSE
  )

  semModel@Pars$edge[object$op == '~~']  <- '<->'
  semModel@Pars$edge[object$op == '~*~']<- '<->'
  semModel@Pars$edge[object$op == '~']  <- '~>'
  semModel@Pars$edge[object$op == '=~'] <- '->'
  semModel@Pars$edge[object$op == '~1'] <- 'int'
  semModel@Pars$edge[grepl('\\|', object$op)] <- '|'

  semModel@Thresholds <- semModel@Pars[grepl('\\|', semModel@Pars$edge), -(3:4)]
  semModel@Pars <- semModel@Pars[!object$op %in% c(':=', '<', '>', '==', '|', '<', '>'), ]

  semModel@Vars <- data.frame(
    name      = c(varNames, factNames),
    manifest  = c(varNames, factNames) %in% varNames,
    exogenous = NA,
    stringsAsFactors = FALSE
  )

  semModel@ObsCovs   <- list()
  semModel@ImpCovs   <- list()
  semModel@Computed  <- FALSE
  semModel@Original  <- list(object)

  semModel
}

pick_est_col <- function(fit) {
  # Different GenomicSEM versions / options may name the standardized column differently.
  if (!is.null(fit$results) && 'STD_All' %in% names(fit$results)) return('STD_All')
  if (!is.null(fit$results) && 'Standardized_Est' %in% names(fit$results)) return('Standardized_Est')
  stop('No standardized estimate column found in fit$results (expected STD_All or Standardized_Est).')
}

est_col <- pick_est_col(fit)
spm <- semPlotModel_GSEM(fit, est.label = est_col)

# ---- Plot ----
# If you want prettier labels, make a named vector like:
#   labels <- c(trait1 = 'Trait 1', trait2 = 'Trait 2')
# and substitute into nodeLabels.
node_names  <- spm@Vars$name
node_labels <- node_names

png(out_png, width = 1800, height = 1300, res = 300)
semPaths(
  spm,
  layout        = 'tree2',
  whatLabels    = 'est',
  nodeLabels    = node_labels,
  edge.width    = 2,
  node.width    = 2,
  edge.label.cex= 0.6,
  mar           = c(5, 4, 6, 2),
  nDigits       = 3
)
dev.off()

message('Wrote path diagram: ', out_png)
