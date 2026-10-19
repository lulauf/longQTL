# baselines.R
#
# Reference methods the longQTL model has to beat. Each represents an
# approach in current use:
#
#   fisher_pooled  - pool ALT vs REF CpG calls of all heterozygotes into one
#                    2x2 table (per-CpG Fisher logic as in many long-read ASM
#                    pipelines). Treats every CpG call as independent.
#   lm_mean        - one methylation value per sample, rank-normal
#                    transformed, regressed on dosage + stratum (fastQTL /
#                    MAGAR logic).
#   wilcoxon_delta - heterozygotes only, delta = meth(ALT reads) -
#                    meth(REF reads), Wilcoxon signed rank (current asmr
#                    call_haplotype_mqtl()).
#   glm_within     - heterozygotes only, counts per sample x allele,
#                    quasi-binomial GLM with a sample fixed effect.
#
# They differ in what they treat as the unit of replication (CpG call, read,
# person). The simulator shows where that choice breaks calibration.

#' Per-read summary of a block for one variant
#'
#' @param sim An `lq_sim` object.
#' @param variant_id Variant to test; default the causal variant, or the
#'   anchor variant when the scenario is null.
#' @param params [lq_params()] for the call mode.
#' @param evidence Allele evidence classes to accept (`"direct"`,
#'   `"phased"`).
#' @return data.table, one row per read: sample_id, read_id, k (methylated
#'   calls, fractional in `expected` mode), n (calls), allele (NA if no
#'   accepted evidence), dosage, stratum.
#' @export
lq_read_summary <- function(sim, variant_id = NULL, params = lq_params(),
                            evidence = "direct") {
  stopifnot(inherits(sim, "lq_sim"))
  vid <- .resolve_variant(sim, variant_id)
  cl <- lq_apply_call_mode(sim$calls, params)
  rd <- cl[, .(k = sum(call), n = .N), by = .(sample_id, read_id)]
  accepted <- match.arg(evidence, c("direct", "phased"), several.ok = TRUE)
  al <- sim$read_alleles[variant_id == vid & evidence %in% accepted, .(read_id, allele)]
  rd <- merge(rd, al, by = "read_id", all.x = TRUE)
  rd <- merge(rd, sim$genotypes[variant_id == vid, .(sample_id, dosage)], by = "sample_id")
  rd <- merge(rd, sim$samples[, .(sample_id, stratum)], by = "sample_id")
  data.table::setattr(rd, "variant_id", vid)
  rd[]
}

.resolve_variant <- function(sim, variant_id) {
  if (!is.null(variant_id)) {
    if (!variant_id %in% sim$variants$variant_id) stop("Unknown variant_id: ", variant_id)
    return(variant_id)
  }
  sc <- sim$scenario
  k <- if (sc$causal > 0) sc$causal else ceiling(length(sc$variant_offsets) / 2)
  sim$variants$variant_id[k]
}

#' Run the baseline methods on one simulated locus
#'
#' @inheritParams lq_read_summary
#' @param methods Subset of `"fisher_pooled"`, `"lm_mean"`,
#'   `"wilcoxon_delta"`, `"glm_within"`.
#' @param min_reads_per_allele For the heterozygote-based methods: minimum
#'   reads per allele for a sample to contribute.
#' @return data.table: method, estimate, p_value, n_units (what the method
#'   counts as independent: calls, samples).
#' @export
lq_baselines <- function(sim, variant_id = NULL, params = lq_params(),
                         evidence = "direct",
                         methods = c("fisher_pooled", "lm_mean", "wilcoxon_delta", "glm_within"),
                         min_reads_per_allele = 3L) {
  methods <- match.arg(methods, several.ok = TRUE)
  rd <- lq_read_summary(sim, variant_id, params, evidence)

  het <- rd[dosage == 1L & !is.na(allele)]
  per_allele <- het[, .(k = sum(k), n = sum(n), n_reads = .N), by = .(sample_id, allele)]
  ok <- per_allele[, .(ok = .N == 2L && all(n_reads >= min_reads_per_allele)), by = sample_id][ok == TRUE, sample_id]
  pa <- per_allele[sample_id %in% ok]

  res <- list()

  if ("fisher_pooled" %in% methods) {
    tot <- het[, .(k = sum(k), n = sum(n)), by = allele]
    if (nrow(tot) == 2L) {
      a <- tot[allele == 1L]; r <- tot[allele == 0L]
      m <- matrix(round(c(a$k, a$n - a$k, r$k, r$n - r$k)), 2L, byrow = TRUE)
      p <- tryCatch(stats::fisher.test(m)$p.value, error = function(e) NA_real_)
      res$fisher_pooled <- list(a$k / a$n - r$k / r$n, p, sum(tot$n))
    } else res$fisher_pooled <- list(NA_real_, NA_real_, 0L)
  }

  if ("lm_mean" %in% methods) {
    ps <- rd[, .(y = sum(k) / sum(n)), by = .(sample_id, dosage, stratum)]
    ps[, y_rn := stats::qnorm((rank(y) - 0.5) / .N)]
    f <- if (data.table::uniqueN(ps$stratum) > 1L) y_rn ~ dosage + stratum else y_rn ~ dosage
    if (data.table::uniqueN(ps$dosage) > 1L) {
      co <- summary(stats::lm(f, data = ps))$coefficients
      res$lm_mean <- list(co["dosage", 1L], co["dosage", 4L], nrow(ps))
    } else res$lm_mean <- list(NA_real_, NA_real_, nrow(ps))
  }

  if ("wilcoxon_delta" %in% methods) {
    if (length(ok) >= 2L) {
      w <- data.table::dcast(pa, sample_id ~ allele, value.var = c("k", "n"))
      d <- w$k_1 / w$n_1 - w$k_0 / w$n_0
      p <- suppressWarnings(stats::wilcox.test(d, mu = 0, exact = FALSE)$p.value)
      res$wilcoxon_delta <- list(mean(d), p, length(d))
    } else res$wilcoxon_delta <- list(NA_real_, NA_real_, length(ok))
  }

  if ("glm_within" %in% methods) {
    if (length(ok) >= 2L) {
      fit <- stats::glm(cbind(k, n - k) ~ factor(allele) + factor(sample_id),
                        family = stats::quasibinomial(), data = pa)
      co <- summary(fit)$coefficients
      res$glm_within <- list(co["factor(allele)1", 1L], co["factor(allele)1", 4L], length(ok))
    } else res$glm_within <- list(NA_real_, NA_real_, length(ok))
  }

  data.table::rbindlist(lapply(names(res), function(m) data.table::data.table(
    method = m, estimate = res[[m]][[1L]], p_value = res[[m]][[2L]],
    n_units = as.integer(res[[m]][[3L]]))))
}

#' Repeat a scenario and collect baseline p-values
#'
#' @param scenario An [lq_sim_scenario()] object.
#' @param n_rep Number of replicates; replicate `r` uses seed
#'   `scenario$seed + 7919 * r`, so results are reproducible and
#'   independent of the backend.
#' @param params,evidence,methods Passed to [lq_baselines()].
#' @param n_workers Forked workers (Unix); 1 = sequential.
#' @return data.table: rep, method, estimate, p_value, n_units.
#' @export
lq_sim_benchmark <- function(scenario, n_rep = 200L, params = lq_params(),
                             evidence = "direct",
                             methods = c("fisher_pooled", "lm_mean", "wilcoxon_delta", "glm_within"),
                             n_workers = 1L) {
  one <- function(r) {
    sim <- lq_simulate_cohort(scenario, seed = scenario$seed + 7919L * r)
    out <- lq_baselines(sim, params = params, evidence = evidence, methods = methods)
    out[, rep := r]
  }
  reps <- seq_len(n_rep)
  rows <- if (n_workers > 1L && .Platform$OS.type == "unix") {
    parallel::mclapply(reps, one, mc.cores = n_workers)
  } else lapply(reps, one)
  out <- data.table::rbindlist(rows)
  data.table::setcolorder(out, "rep")
  data.table::setattr(out, "scenario", scenario)
  out[]
}

#' Rejection rates and inflation per method
#'
#' Under a null scenario the rejection rate at level alpha should be close
#' to alpha; under an alternative it is the power.
#'
#' @param bench Output of [lq_sim_benchmark()].
#' @param alpha Levels.
#' @return data.table per method: n (non-missing p-values), reject_<alpha>
#'   columns, lambda_gc (median chi-square inflation), mean_estimate.
#' @export
lq_calibration <- function(bench, alpha = c(0.05, 0.01, 0.001)) {
  b <- bench[!is.na(p_value)]
  out <- b[, c(list(n = .N),
               stats::setNames(lapply(alpha, function(a) mean(p_value < a)),
                               paste0("reject_", alpha)),
               list(lambda_gc = stats::median(stats::qchisq(p_value, 1, lower.tail = FALSE)) /
                      stats::qchisq(0.5, 1),
                    mean_estimate = mean(estimate))),
           by = method]
  out[]
}
