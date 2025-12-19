#' Plot method for jmleIRT objects
#'
#' Quick diagnostic plots for Rasch JMLE results.
#'
#' @param x An object of class \code{"jmleIRT"}.
#' @param type Character string indicating the type of plot:
#'   \code{"hist"} for histograms of person and item parameters,
#'   \code{"icc"} for item characteristic curves of selected items.
#' @param items Integer vector of item indices to plot when \code{type = "icc"}.
#'   Defaults to the first 5 items (or all if fewer than 5).
#' @param ... Currently unused.
#' @method plot jmleIRT
#' @importFrom graphics hist legend lines par rug
#' @export
plot.jmleIRT <- function(x, type = c("hist", "icc"), items = NULL, ...) {
  type <- match.arg(type)
  theta <- x$theta
  beta <- x$beta

  if (type == "hist") {
    op <- par(mfrow = c(1, 2))
    on.exit(par(op), add = TRUE)

    hist(theta,
      main = "Person abilities (theta)",
      xlab = expression(theta), col = "grey80", border = "white"
    )
    rug(theta)

    hist(beta,
      main = "Item difficulties (beta)",
      xlab = expression(beta), col = "grey80", border = "white"
    )
    rug(beta)
  } else if (type == "icc") {
    if (is.null(items)) {
      items <- seq_len(min(5L, length(beta)))
    }
    items <- intersect(items, seq_along(beta))
    if (length(items) == 0L) {
      stop("No valid item indices in 'items'.")
    }

    # Grid of theta values
    th_grid <- seq(min(theta[is.finite(theta)], na.rm = TRUE) - 1,
      max(theta, na.rm = TRUE) + 1,
      length.out = 200
    )

    plot(NA,
      xlim = range(th_grid),
      ylim = c(0, 1),
      xlab = expression(theta),
      ylab = "P(X = 1)",
      main = "Item characteristic curves"
    )

    cols <- grDevices::rainbow(length(items))
    for (j in seq_along(items)) {
      i <- items[j]
      p <- 1 / (1 + exp(-(th_grid - beta[i])))
      lines(th_grid, p, col = cols[j], lwd = 2)
    }
    legend("bottomright",
      legend = paste0("Item ", items),
      col = cols, lwd = 2, bty = "n"
    )
  }

  invisible(x)
}
