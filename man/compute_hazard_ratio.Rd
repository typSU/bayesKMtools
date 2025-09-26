#' Estimate the log hazard ratio for two groups
#'
#' \code{compute_hazard_ratio()} estimates the log hazard ratio \eqn{\beta}
#' between a treatment and a control group when discrete hazards are known.
#' It solves the discrete-time Cox score equation by fixed-point iteration.
#'
#' \strong{Hazard ratio computation:}
#' Let \eqn{R^{(t)}(k)} and \eqn{R^{(c)}(k)} denote the numbers at risk on day
#' \eqn{k} in the treatment and control groups, and let \eqn{p^{(t)}_k} and
#' \eqn{p^{(c)}_k} be posterior draws of the discrete hazards.  The hazard ratio
#' \eqn{\exp(\beta)} satisfies
#' \deqn{N = \sum_{k=1}^K
#'       \frac{R^{(t)}(k)\, R^{(c)}(k)}{R^{(c)}(k) + \exp(\beta)\, R^{(t)}(k)}\,
#'       p^{(t)}_k,}
#' \deqn{D = \sum_{k=1}^K
#'       \frac{R^{(t)}(k)\, R^{(c)}(k)}{R^{(c)}(k) + \exp(\beta)\, R^{(t)}(k)}\,
#'       p^{(c)}_k,}
#' \deqn{\exp(\beta) = \frac{N}{D},}
#' where the sums run over the discrete time grid \eqn{k=1,\dots,K}.
#' For each posterior draw of \eqn{p^{(t)}_k} and \eqn{p^{(c)}_k}, the
#' algorithm iteratively updates \eqn{\beta} until convergence.  The resulting
#' set of \eqn{\beta}-values can then be summarised to obtain a posterior mean,
#' standard deviation, and credible intervals for the hazard ratio.
#'
#' @param trt Integer vector indicating group membership (1 = treatment,
#'   0 = control) for each subject.
#' @param time Integer vector of discrete follow‑up times, matching the length
#'   of \code{trt}.
#' @param hazard_trt Numeric vector of length \eqn{K} giving the posterior mean or
#'   sample hazards at each time point for the treatment group.
#' @param hazard_ctrl Numeric vector of length \eqn{K} giving the hazards for
#'   the control group.
#' @param beta_init Numeric; initial value of the log hazard ratio for the
#'   iteration (default \code{0}).
#' @param tol Numeric convergence tolerance (default \code{1e-6}).
#' @param max_iter Integer maximum number of iterations (default \code{10}).
#'
#' @return A numeric scalar giving the estimated log hazard ratio.  The hazard
#'   ratio itself is \code{exp(result)}.
#'
#' @examples
#' # Simulate hazards for demonstration
#' trt  <- c(1,0,1,0,1)
#' time <- c(2,3,3,5,5)
#' hazard_trt  <- c(0.05,0.04,0.03,0.02,0.01)
#' hazard_ctrl <- c(0.06,0.05,0.04,0.03,0.02)
#' log_hr <- compute_hazard_ratio(trt, time, hazard_trt, hazard_ctrl)
#' exp(log_hr)
