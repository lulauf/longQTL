# calls.R
#
# Turning modkit-style per-call probabilities into the numbers the models
# count (concept doc, sec. 4.3). The probability column is P(modified); the
# confidence of a call is max(p, 1 - p).

#' Probability of "modified" from 5mC / 5hmC ML qualities
#'
#' ML integers q encode the probability bin [q/256, (q+1)/256); the bin
#' centre (q + 0.5) / 256 is used. Missing codes count as probability 0
#' when the other code is present.
#'
#' @param calls data.table with integer columns `q_m`, `q_h` (0-255, NA).
#' @param mod_combine `"m_plus_h"`: P(5mC) + P(5hmC) (what bisulfite
#'   measures); `"m_only"`: P(5mC), 5hmC counts as unmodified; `"drop_h"`:
#'   calls whose most likely state is 5hmC are removed, P(5mC) is
#'   renormalised over 5mC + unmodified.
#' @return `calls` with `p_m`, `p_h` and `prob_mod`.
#' @export
lq_mod_prob <- function(calls, mod_combine = c("m_plus_h", "m_only", "drop_h")) {
  mod_combine <- match.arg(mod_combine)
  dt <- data.table::copy(data.table::as.data.table(calls))
  dt <- dt[!(is.na(q_m) & is.na(q_h))]
  dt[, `:=`(p_m = data.table::fifelse(is.na(q_m), 0, (q_m + 0.5) / 256),
            p_h = data.table::fifelse(is.na(q_h), 0, (q_h + 0.5) / 256))]
  tot <- dt$p_m + dt$p_h
  over <- tot > 1
  if (any(over)) dt[over, `:=`(p_m = p_m / tot[over], p_h = p_h / tot[over])]
  if (mod_combine == "m_plus_h") {
    dt[, prob_mod := pmin(1, p_m + p_h)]
  } else if (mod_combine == "m_only") {
    dt[, prob_mod := p_m]
  } else {
    dt <- dt[!(p_h > p_m & p_h > 1 - p_m - p_h)]
    dt[, prob_mod := p_m / pmax(1e-12, 1 - p_h)]
  }
  dt[]
}

#' Apply the call mode to per-read CpG calls
#'
#' @param calls data.table with `prob_mod` (probability of the modified
#'   state) or ML qualities `q_m` / `q_h` (as from [lq_extract_reads()]),
#'   and, for `"hard_modkit"`, `sample_id`.
#' @param params An [lq_params()] object (uses `call_mode`, `min_call_prob`),
#'   or NULL with `mode` / `min_call_prob` given directly.
#' @param mode,min_call_prob Override without a params object.
#' @param mod_combine How 5hmC is counted when `calls` carries ML
#'   qualities (`q_m`, `q_h`) instead of `prob_mod`; see [lq_mod_prob()].
#'   Taken from `params` when given.
#' @param filter_percentile For `"hard_modkit"`: calls below this quantile
#'   of confidence within each sample are dropped (modkit's default
#'   behaviour uses the 10th percentile).
#'
#' @return `calls` with a numeric column `call` (0/1, or a probability for
#'   `"expected"`); dropped (ambiguous) calls are removed. Attribute
#'   `call_qc` gives the number and fraction dropped, per sample.
#' @export
lq_apply_call_mode <- function(calls, params = NULL,
                               mode = c("hard_fixed", "hard_modkit", "expected"),
                               min_call_prob = 0.8, filter_percentile = 0.1,
                               mod_combine = c("m_plus_h", "m_only", "drop_h")) {
  if (!is.null(params)) {
    stopifnot(inherits(params, "lq_params"))
    mode <- params$call_mode
    min_call_prob <- params$min_call_prob
    mod_combine <- params$mod_combine
  } else {
    mode <- match.arg(mode)
    mod_combine <- match.arg(mod_combine)
  }
  dt <- data.table::as.data.table(calls)
  if (!"prob_mod" %in% names(dt) && all(c("q_m", "q_h") %in% names(dt))) {
    dt <- lq_mod_prob(dt, mod_combine)
  }
  if (!"prob_mod" %in% names(dt)) stop("calls needs prob_mod, or q_m and q_h.")
  if (anyNA(dt$prob_mod) || any(dt$prob_mod < 0 | dt$prob_mod > 1)) {
    stop("prob_mod must be in [0, 1] without NA.")
  }
  dt <- data.table::copy(dt)
  dt[, conf_ := pmax(prob_mod, 1 - prob_mod)]

  keep <- switch(mode,
    hard_fixed  = dt$conf_ >= min_call_prob,
    hard_modkit = {
      if (!"sample_id" %in% names(dt)) stop("hard_modkit needs sample_id.")
      dt[, thr_ := stats::quantile(conf_, filter_percentile, names = FALSE), by = sample_id]
      dt$conf_ >= dt$thr_
    },
    expected = rep(TRUE, nrow(dt))
  )
  dt[, kept_ := keep]
  qc <- if ("sample_id" %in% names(dt)) {
    dt[, .(n_calls = .N, n_dropped = sum(!kept_), frac_dropped = mean(!kept_)), by = sample_id]
  } else {
    dt[, .(n_calls = .N, n_dropped = sum(!kept_), frac_dropped = mean(!kept_))]
  }

  out <- dt[kept_ == TRUE]
  out[, call := if (mode == "expected") prob_mod else as.numeric(prob_mod >= 0.5)]
  out[, intersect(c("conf_", "thr_", "kept_"), names(out)) := NULL]
  data.table::setattr(out, "call_qc", qc)
  data.table::setattr(out, "call_mode", mode)
  out[]
}
