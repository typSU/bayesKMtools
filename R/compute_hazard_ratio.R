#' Estimate the log‑hazard ratio for two groups
#'
#' \code{compute_hazard_ratio()} estimates the log hazard ratio between a
#' treatment and a control group when the discrete hazard functions for
#' each group are known.  The function iteratively solves for the log hazard
#' ratio \eqn{\beta} by equating weighted sums of group‑specific hazards, where
#' weights reflect the number at risk in each group at each time point.  The
#' estimate returned is \eqn{\beta}; exponentiate it to obtain the hazard ratio.
#'
#' @param trt Integer vector of group indicators for each subject (1 = treatment,
#'   0 = control).
#' @param time Integer vector of discrete follow‑up times for each subject;
#'   must be the same length as \code{trt}.  Times should take values from 1 to
#'   the length of \code{hazard_trt}.
#' @param hazard_trt Numeric vector of length \eqn{K} giving the hazard
#'   probabilities at each discrete time point for the treatment group.
#' @param hazard_ctrl Numeric vector of length \eqn{K} giving the hazard
#'   probabilities at each discrete time point for the control group.
#' @param beta_init Numeric; initial value of the log hazard ratio for the
#'   iterative procedure (default \code{0}).
#' @param tol Numeric; convergence tolerance for the iterative update
#'   (default \code{1e-6}).  The algorithm stops when successive updates differ by less than
#'   \code{tol}.
#' @param max_iter Integer; maximum number of iterations to perform (default \code{10}).
#'
#' @return A numeric scalar giving the estimated log hazard ratio.  The hazard
#'   ratio itself is \code{exp( result )}.
#'
#' @examples
#' # Example with simulated discrete hazards
#' set.seed(42)
#' trt  <- rbinom(10, 1, 0.5)
#' time <- sample(1:5, 10, replace = TRUE)
#' haz_treat <- c(0.05, 0.04, 0.03, 0.02, 0.01)
#' haz_ctrl  <- c(0.06, 0.05, 0.04, 0.03, 0.02)
#' log_hr <- compute_hazard_ratio(trt, time, haz_treat, haz_ctrl)
#' hr <- exp(log_hr)
#' hr
compute_hazard_ratio <- function(trt, time,
                                 hazard_trt, hazard_ctrl,
                                 beta_init = 0,
                                 tol = 1e-6,
                                 max_iter = 10) {
  # input validation
  if (length(hazard_trt) != length(hazard_ctrl)) {
    stop("`hazard_trt` and `hazard_ctrl` must have the same length")
  }
  if (length(trt) != length(time)) {
    stop("`trt` and `time` must have the same length")
  }
  if (!all(trt %in% c(0, 1))) {
    stop("`trt` must contain only 0 (control) and 1 (treatment)")
  }
  K <- length(hazard_trt)

  # compute number at risk in each group at each time
  at_risk_trt  <- sapply(seq_len(K), function(i) sum(trt[time >= i] == 1))
  at_risk_ctrl <- sapply(seq_len(K), function(i) sum(trt[time >= i] == 0))

  beta_val <- beta_init
  for (iter in seq_len(max_iter)) {
    # numerator and denominator for the profile likelihood update
    num <- sum(hazard_trt  * (at_risk_trt * at_risk_ctrl) /
                 (at_risk_ctrl + exp(beta_val) * at_risk_trt))
    den <- sum(hazard_ctrl * (at_risk_trt * at_risk_ctrl) /
                 (at_risk_ctrl + exp(beta_val) * at_risk_trt))

    beta_new <- log(num / den)
    if (abs(beta_new - beta_val) < tol) {
      beta_val <- beta_new
      break
    }
    beta_val <- beta_new
  }
  beta_val
}
