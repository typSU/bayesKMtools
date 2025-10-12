#' Compute posterior log hazard ratios for a two‑group comparison
#'
#' Given discrete hazard draws for a treatment and a control group,
#' \code{compute_hazard_ratio()} solves the discrete‑time Cox score
#' equation for the log hazard ratio \eqn{\beta} on a per‑draw basis.  For
#' each posterior draw \eqn{p^{(t)}_k} and \eqn{p^{(c)}_k}, the hazard ratio
#' estimate satisfies
#' \deqn{\exp(\beta) = \frac{N}{D},}
#' where
#' \deqn{N = \sum_{k=1}^K
#'       \frac{R^{(t)}(k)\, R^{(c)}(k)}{R^{(c)}(k) + \exp(\beta)\, R^{(t)}(k)}
#'       p^{(t)}_k,\quad
#' D = \sum_{k=1}^K
#'       \frac{R^{(t)}(k)\, R^{(c)}(k)}{R^{(c)}(k) + \exp(\beta)\, R^{(t)}(k)}
#'       p^{(c)}_k.}
#' Here \deqn{R^{(t)}(k)} and \deqn{R^{(c)}(k)} denote the numbers at risk in the
#' treatment and control groups at day \eqn{k}.
#' The function iterates to convergence for each draw and returns a vector of
#' log‑hazard ratios.  You can exponentiate these to obtain a posterior
#' distribution for the hazard ratio itself.
#'
#' @param trt Integer vector of group indicators (1 = treatment, 0 = control).
#' @param time Integer vector of observed times (days) for each subject; must
#'   be the same length as \code{trt}.
#' @param hazard_trt Numeric matrix of size \code{n_draws × K}, where each
#'   row contains a Monte‑Carlo draw of the discrete hazards for the treatment
#'   group.
#' @param hazard_ctrl Numeric matrix of size \code{n_draws × K}, containing
#'   hazard draws for the control group.
#' @param beta_init Numeric; initial value of the log hazard ratio for the
#'   iterative procedure (default \code{0}).
#' @param tol Numeric tolerance for convergence (default \code{1e-6}).
#' @param max_iter Integer; maximum number of iterations per draw (default \code{10}).
#'
#' @return A numeric vector of length \code{n_draws}; each element is the
#'   estimated log hazard ratio for one pair of hazard draws.  Exponentiate
#'   this vector to obtain the posterior distribution of the hazard ratio.
#'
#' @examples
#' # A simple illustration with two posterior draws of discrete hazards.
#' # Four subjects: two in the treatment group (trt = 1) and two in control (trt = 0).
#' trt  <- c(1, 0, 1, 0)
#' time <- c(2, 3, 3, 4)  # discrete follow‑up times
#'
#' # Each row of hazard_trt and hazard_ctrl is a draw of daily hazards.
#' # Here we have two posterior draws for a 4‑day grid.
#' hazard_trt  <- matrix(c(0.10, 0.05, 0.02, 0.02,
#'                         0.08, 0.04, 0.03, 0.01),
#'                       nrow = 2, byrow = TRUE)
#' hazard_ctrl <- matrix(c(0.12, 0.07, 0.04, 0.03,
#'                         0.10, 0.06, 0.02, 0.02),
#'                       nrow = 2, byrow = TRUE)
#'
#' # Compute log hazard ratios for each draw and exponentiate.
#' log_hr <- compute_hazard_ratio(trt, time, hazard_trt, hazard_ctrl)
#' exp(log_hr)
compute_hazard_ratio <- function(trt, time,
                                 hazard_trt, hazard_ctrl,
                                 beta_init = 0,
                                 tol = 1e-6,
                                 max_iter = 10) {
  if (!all(trt %in% c(0, 1))) {
    stop("`trt` must contain only 0 (control) and 1 (treatment) values")
  }
  if (length(trt) != length(time)) {
    stop("`trt` and `time` must have the same length")
  }
  if (!is.matrix(hazard_trt) || !is.matrix(hazard_ctrl)) {
    stop("`hazard_trt` and `hazard_ctrl` must be matrices of equal dimension")
  }
  if (nrow(hazard_trt) != nrow(hazard_ctrl) ||
      ncol(hazard_trt) != ncol(hazard_ctrl)) {
    stop("`hazard_trt` and `hazard_ctrl` must have the same dimensions")
  }
  K <- ncol(hazard_trt)
  # compute at-risk counts for each day
  at_risk_trt  <- sapply(seq_len(K), function(i) sum(trt[time >= i] == 1))
  at_risk_ctrl <- sapply(seq_len(K), function(i) sum(trt[time >= i] == 0))
  n_draws <- nrow(hazard_trt)
  beta_vec <- numeric(n_draws)
  # loop over hazard draws
  for (j in seq_len(n_draws)) {
    p_trt  <- hazard_trt[j, ]
    p_ctrl <- hazard_ctrl[j, ]
    beta_val <- beta_init
    for (iter in seq_len(max_iter)) {
      num <- sum(p_trt  * (at_risk_trt * at_risk_ctrl) /
                   (at_risk_ctrl + exp(beta_val) * at_risk_trt))
      den <- sum(p_ctrl * (at_risk_trt * at_risk_ctrl) /
                   (at_risk_ctrl + exp(beta_val) * at_risk_trt))
      beta_new <- log(num / den)
      if (abs(beta_new - beta_val) < tol) {
        beta_val <- beta_new
        break
      }
      beta_val <- beta_new
    }
    beta_vec[j] <- beta_val
  }
  beta_vec
}
