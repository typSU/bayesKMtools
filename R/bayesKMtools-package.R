#' bayesKMtools: Bayesian Kaplan–Meier tools for discrete‐time survival analysis
#'
#' The **bayesKMtools** package implements a Bayesian analogue of the
#' Kaplan–Meier estimator for survival data measured in discrete days, as
#' described in our manuscript on Bayesian survival curves.  The key idea is
#' to treat each patient’s event time as a draw from a multinomial
#' distribution over a daily time grid.  A Dirichlet prior on these event
#' probabilities induces independent Beta priors on the day‐by‐day hazards,
#' allowing straightforward posterior sampling.
#'
#' Core functions include:
#' \itemize{
#'   \item \code{\link{compute_daily_hazards}}: draws Monte‑Carlo samples of the
#'     posterior hazard probability at each day, given right‑censored
#'     survival times and Beta hyper‑parameters.
#'   \item \code{\link{compute_hazard_ratio}}: computes the log hazard ratio
#'     between a treatment group and a control group using posterior hazard
#'     estimates under a discrete‑time proportional hazards model.
#'   \item \code{\link{plot_survival_with_km}}: plots the posterior mean
#'     survival curve and sampled posterior trajectories alongside the classical
#'     Kaplan–Meier estimate.
#' }
#'
#' These tools make it easy to reproduce the analyses in the paper, to explore
#' sensitivity to the prior, and to visualise uncertainty in the survival
#' curve.  The package assumes integer follow‑up times and focuses on a simple
#' Gamma/Beta prior structure, but the modular functions can be adapted for
#' more complex settings.
#'
#' @docType package
#' @name bayesKMtools
"_PACKAGE"
