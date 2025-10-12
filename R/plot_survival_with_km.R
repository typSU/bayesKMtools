#' Plot Kaplan–Meier and posterior survival curves
#'
#' This function produces a survival plot that overlays the classical Kaplan–Meier
#' estimate, the posterior mean survival curve computed from a matrix of hazard
#' draws, and a set of posterior survival trajectories.  It is designed for
#' discrete‑time Bayesian survival analysis where hazards are sampled at each
#' integer time point.
#'
#' @param time Integer vector of observed times (days) for each subject.  Times
#'   should take values from 1 up to \code{max(time_grid)}.
#' @param status Integer vector of the same length as \code{time}, where
#'   \code{1} denotes an observed event and \code{0} denotes right‑censoring.
#' @param time_grid Integer vector defining the discretised time points.  Its
#'   length should match the number of columns in \code{hazard_draws}.
#' @param hazard_draws Numeric matrix of posterior hazard draws.  Each row is a
#'   Monte Carlo sample and each column corresponds to a time point in
#'   \code{time_grid}.
#' @param ylim Numeric vector of length two giving the y‑axis limits for survival
#'   probabilities.  Defaults to \code{c(0, 1)}.
#' @param burn Integer; number of initial rows in \code{hazard_draws} to discard
#'   as burn‑in (default \code{5000}).
#' @param n_sample Integer; number of posterior survival trajectories to display
#'   in the plot (default \code{100}).
#' @param seed Integer seed used when sampling the posterior trajectories
#'   (default \code{1}).
#' @param palette_option Character string specifying the viridis palette used by
#'   \code{ggplot2::scale_color_viridis_d()} (default \code{"viridis"}).
#' @param title Character string giving the plot title (default \code{"(a)"}).
#'
#' @return A \code{ggplot2} object.  You can further customise or save the plot
#'   using standard ggplot2 syntax.
#'
#' @importFrom survival Surv survfit
#' @importFrom tibble tibble
#' @importFrom ggplot2 ggplot geom_step scale_color_viridis_d coord_cartesian labs theme_minimal theme element_rect aes
#' @importFrom scales alpha
#' @export
#'
#' @examples
#' # Simulate a small data set and hazard draws for illustration
#' set.seed(2023)
#' times  <- c(3, 7, 7, 10, 12)
#' status <- c(1, 0, 1, 0, 1)
#' grid   <- seq_len(max(times))
#' # Fake posterior hazard draws: 200 samples by length(grid) columns
#' haz_draws <- matrix(runif(200 * length(grid), 0, 0.2), nrow = 200)
#' p <- plot_survival_with_km(times, status, grid, haz_draws,
#'                            burn = 50, n_sample = 5)
#' print(p)
plot_survival_with_km <- function(time, status,
                                  time_grid, hazard_draws,
                                  ylim = c(0, 1),
                                  burn = 0,
                                  n_sample = 100,
                                  seed = 1,
                                  palette_option = "viridis",
                                  title = "(a)") {
  if (length(time) != length(status)) {
    stop("`time` and `status` must be the same length")
  }
  if (ncol(hazard_draws) != length(time_grid)) {
    stop("Number of columns in `hazard_draws` must match length of `time_grid`")
  }
  # 1. Kaplan–Meier estimate
  km_fit <- survival::survfit(survival::Surv(time, status) ~ 1)
  km_df <- tibble::tibble(
    time     = c(0, km_fit$time),
    surv     = c(1, km_fit$surv),
    Estimate = "KM Estimate"
  )

  # 2. Posterior mean survival (with tail extension)
  K <- length(time_grid)
  # convert hazards to survival: survival = cumprod(1 - hazard)
  surv_mat <- 1 - hazard_draws[, seq_len(K), drop = FALSE]
  for (i in seq_len(nrow(surv_mat))) {
    surv_mat[i, ] <- cumprod(surv_mat[i, ])
  }
  max_t <- max(time_grid)
  post_mean_df <- tibble::tibble(
    time     = c(0, time_grid, max_t),
    surv     = c(1, colMeans(surv_mat), colMeans(surv_mat)[K]),
    Estimate = "Posterior Mean"
  )

  # 3. Sample a few posterior survival trajectories
  set.seed(seed)
  if (burn >= nrow(surv_mat)) {
    stop("`burn` must be less than the number of rows in `hazard_draws`")
  }
  available_idx <- seq(from = burn + 1, to = nrow(surv_mat))
  n_sample <- min(n_sample, length(available_idx))
  samp_idx <- sample(available_idx, n_sample, replace = FALSE)
  post_samp_df <- do.call(rbind, lapply(samp_idx, function(i) {
    surv_i <- surv_mat[i, ]
    tibble::tibble(
      time   = c(0, time_grid, max_t),
      surv   = c(1, surv_i, surv_i[K]),
      sample = factor(i)
    )
  }))

  # 4. Compose the plot
  ggplot2::ggplot() +
    # posterior sample trajectories
    ggplot2::geom_step(
      data  = post_samp_df,
      mapping = ggplot2::aes(x = time, y = surv, group = sample),
      color = "grey70", alpha = 0.3
    ) +
    # posterior mean and KM estimate
    ggplot2::geom_step(
      data = post_mean_df,
      mapping = ggplot2::aes(x = time, y = surv, color = Estimate),
      size = 1.2
    ) +
    ggplot2::geom_step(
      data = km_df,
      mapping = ggplot2::aes(x = time, y = surv, color = Estimate),
      size = 1.2
    ) +
    ggplot2::scale_color_viridis_d(option = palette_option) +
    ggplot2::coord_cartesian(ylim = ylim) +
    ggplot2::labs(
      title = title,
      x     = "Time (days)",
      y     = "Survival probability",
      color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      legend.position      = c(0.95, 0.95),
      legend.justification = c("right", "top"),
      legend.background    = ggplot2::element_rect(
        fill = scales::alpha("white", 0.5)
      )
    )
}
