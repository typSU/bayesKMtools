#' Sample posterior daily hazards on a discretised timeline
#'
#' For a discrete follow‑up measured in days, \code{compute_daily_hazards()}
#' draws Monte Carlo samples from the posterior distribution of the discrete hazard
#' at each day.  Let \eqn{n_i} be the number of observed events on day \eqn{i},
#' and \eqn{d_i} be the number of individuals still at risk at the start of day
#' \eqn{i}.  Under independent \eqn{\mathrm{Beta}(a_0, b_0)} priors on each
#' hazard, the posterior for the hazard on day \eqn{i} is
#' \eqn{\mathrm{Beta}(a_0 + n_i,\; b_0 + d_i - n_i)}.  This function samples
#' \code{n_samp} values from each of these Beta distributions and returns them
#' as a matrix.
#'
#' @param time Integer vector of follow‑up times (days) for each subject.  Times
#'   should be positive integers starting at 1.
#' @param status Integer vector of the same length as \code{time}, where
#'   \code{status = 1} indicates that the event occurred at that time and
#'   \code{status = 0} indicates a right‑censored observation.
#' @param a0 Positive numeric value giving the prior shape parameter for the number
#'   of events on each day (default is \code{0.001}).
#' @param b0 Positive numeric value giving the prior shape parameter for the number
#'   of non‑events on each day (default is \code{1}).
#' @param maxT A positive integer giving the upper bound of the discretized
#'   follow-up time. This value defines the largest day index in the time grid
#'   and therefore the maximum possible survival time considered in the analysis.
#' @param n_samp Integer; the number of Monte Carlo samples to draw from each
#'   Beta posterior distribution (default is \code{10000}).
#' @param progress Logical; if \code{TRUE} (the default) a progress bar is displayed.
#'
#' @return A numeric matrix with \code{n_samp} rows and \code{K} columns, where
#'   \code{K = max(time)}.  Column \code{i} contains draws from the posterior
#'   distribution of the hazard on day \code{i}.
#'
#' @examples
#' set.seed(1)
#' times  <- c(2, 4, 4, 5, 6, 7)
#' status <- c(1, 0, 1, 0, 1, 1)
#' hz_draws <- compute_daily_hazards(times, status, a0 = 0.01, b0 = 1, n_samp = 500)
#' # Posterior mean hazard for each day
#' colMeans(hz_draws)
compute_daily_hazards <- function(time, status,
                                  a0 = 0.001,
                                  b0 = 1,
                                  maxT = max(time),
                                  n_samp = 10000,
                                  progress = TRUE) {
  # basic input checks
  if (length(time) != length(status)) {
    stop("`time` and `status` must have the same length")
  }
  if (any(status != 0 & status != 1)) {
    stop("`status` must contain only 0 (censored) and 1 (event) values")
  }
  if (!is.numeric(a0) || a0 <= 0 || !is.numeric(b0) || b0 <= 0) {
    stop("`a0` and `b0` must be positive numeric values")
  }
  if (n_samp <= 0 || n_samp != as.integer(n_samp)) {
    stop("`n_samp` must be a positive integer")
  }

  # maximum day defines the number of columns (assumes integer days starting at 1)
  K <- maxT
  # number of events on each day
  num_events <- tabulate(time[status == 1], nbins = K)
  # number at risk at the start of each day: count subjects with time >= day
  den_atrisk <- sapply(seq_len(K), function(i) sum(time >= i))
  # preallocate result matrix
  hazard_samples <- matrix(0, nrow = n_samp, ncol = K)
  # optional progress bar
  if (progress) {
    pb <- utils::txtProgressBar(min = 0, max = K, style = 3)
  }
  for (i in seq_len(K)) {
    # posterior Beta parameters
    shape1 <- num_events[i] + a0
    shape2 <- den_atrisk[i] - num_events[i] + b0
    hazard_samples[, i] <- stats::rbeta(n_samp, shape1 = shape1, shape2 = shape2)
    if (progress) utils::setTxtProgressBar(pb, i)
  }
  if (progress) close(pb)
  colnames(hazard_samples) <- paste0("day_", seq_len(K))
  hazard_samples
}
