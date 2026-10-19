# betabinom.R
#
# Beta-binomial building blocks on read level (concept doc, sec. 4.2):
# a read contributes k methylated out of n called CpGs in a block, and the
# beta-binomial over-dispersion rho is the correlation of CpG states on one
# molecule. rho -> 0 is the binomial (CpGs independent), rho -> 1 means every
# read is all-methylated or all-unmethylated (perfect epialleles).

# log-likelihood per read; k may be fractional (call_mode = "expected")
.bb_ll <- function(k, n, mu, rho) {
  mu  <- pmin(pmax(mu, 1e-6), 1 - 1e-6)
  rho <- pmin(pmax(rho, 1e-6), 1 - 1e-6)
  phi <- 1 / rho - 1
  a <- mu * phi
  b <- (1 - mu) * phi
  lchoose(n, k) + lbeta(k + a, n - k + b) - lbeta(a, b)
}

#' Fit an intercept-only beta-binomial to read counts
#'
#' @param k,n Methylated and called CpGs per read.
#' @param rho Optional fixed over-dispersion; NULL = estimate it.
#' @return list(mu, rho, loglik, converged)
#' @noRd
.bb_fit_null <- function(k, n, rho = NULL) {
  m0 <- sum(k) / sum(n)
  m0 <- min(max(m0, 1e-3), 1 - 1e-3)
  if (is.null(rho)) {
    nll <- function(p) -sum(.bb_ll(k, n, stats::plogis(p[1]), stats::plogis(p[2])))
    fit <- stats::optim(c(stats::qlogis(m0), stats::qlogis(0.1)), nll, method = "BFGS")
    list(mu = stats::plogis(fit$par[1]), rho = stats::plogis(fit$par[2]),
         loglik = -fit$value, converged = fit$convergence == 0L)
  } else {
    nll <- function(p) -sum(.bb_ll(k, n, stats::plogis(p), rho))
    fit <- stats::optimize(nll, c(-12, 12))
    list(mu = stats::plogis(fit$minimum), rho = rho, loglik = -fit$objective, converged = TRUE)
  }
}

#' Two-group beta-binomial fit with fixed over-dispersion
#'
#' logit(mu_r) = alpha + beta * g_r. Returns beta (log odds ratio of group 1
#' vs group 0), its Wald SE and the likelihood-ratio statistic against the
#' null beta = 0 at the same rho.
#' @noRd
.bb_fit_groups <- function(k, n, g, rho) {
  nll <- function(p) -sum(.bb_ll(k, n, stats::plogis(p[1] + p[2] * g), rho))
  m0 <- min(max(sum(k[g == 0]) / sum(n[g == 0]), 1e-3), 1 - 1e-3)
  m1 <- min(max(sum(k[g == 1]) / sum(n[g == 1]), 1e-3), 1 - 1e-3)
  start <- c(stats::qlogis(m0), stats::qlogis(m1) - stats::qlogis(m0))
  fit <- stats::optim(start, nll, method = "BFGS", hessian = TRUE)
  null <- .bb_fit_null(k, n, rho)
  se <- tryCatch(sqrt(diag(solve(fit$hessian)))[2], error = function(e) NA_real_)
  lrt <- max(0, 2 * (-fit$value - null$loglik))
  list(alpha = fit$par[1], beta = fit$par[2], se = se, lrt = lrt,
       p = stats::pchisq(lrt, df = 1, lower.tail = FALSE),
       converged = fit$convergence == 0L)
}

#' Two-component binomial mixture on reads (phasing-free bimodality)
#'
#' Each read is drawn from a "low" or a "high" epiallele with methylation
#' p_low / p_high; CpGs are independent given the component. EM with a fixed
#' number of iterations from a deterministic start. Compared with a single
#' beta-binomial (same number of parameters + 1) by BIC.
#' @return list(p_low, p_high, pi_high, separation, delta_bic, class)
#'   with `class` the posterior P(high) per read.
#' @noRd
.bimodality <- function(k, n, iter = 100L, tol = 1e-8) {
  m <- k / n
  pl <- max(0.01, stats::quantile(m, 0.25, names = FALSE))
  ph <- min(0.99, stats::quantile(m, 0.75, names = FALSE))
  if (ph - pl < 0.05) { pl <- max(0.01, pl - 0.1); ph <- min(0.99, ph + 0.1) }
  w <- 0.5
  ll_old <- -Inf
  for (it in seq_len(iter)) {
    l1 <- stats::dbinom(round(k), round(n), ph, log = TRUE) + log(w)
    l0 <- stats::dbinom(round(k), round(n), pl, log = TRUE) + log(1 - w)
    mx <- pmax(l1, l0)
    ll <- sum(mx + log(exp(l1 - mx) + exp(l0 - mx)))
    z <- exp(l1 - mx) / (exp(l1 - mx) + exp(l0 - mx))
    w <- min(max(mean(z), 1e-3), 1 - 1e-3)
    ph <- min(max(sum(z * k) / sum(z * n), 1e-4), 1 - 1e-4)
    pl <- min(max(sum((1 - z) * k) / sum((1 - z) * n), 1e-4), 1 - 1e-4)
    if (abs(ll - ll_old) < tol) break
    ll_old <- ll
  }
  if (pl > ph) { tmp <- pl; pl <- ph; ph <- tmp; w <- 1 - w; z <- 1 - z }
  bb <- .bb_fit_null(k, n)
  nr <- length(k)
  bic_mix <- -2 * ll + 3 * log(nr)
  bic_bb  <- -2 * bb$loglik + 2 * log(nr)
  list(p_low = pl, p_high = ph, pi_high = w, separation = ph - pl,
       delta_bic = bic_bb - bic_mix, class = z)
}
