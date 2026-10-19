# params.R
#
# One validated parameter object per run. Everything that changes a result
# lives here, so it can be written into the provenance of every output.

#' Analysis parameters for a longQTL run
#'
#' @section cis window and compute cost:
#' The genetic effect is split into two components (see concept doc, sec. 4):
#'
#' * **within** (allele-specific): ALT vs REF reads inside one heterozygous
#'   person. A read can only be assigned an allele at a variant it covers
#'   (`within_reach = "read"`) or, with phasing, at a variant in the same
#'   phase set (`within_reach = "phase_set"`). Its reach is therefore bounded
#'   by read length (~tens of kb) or phase-block length, never by
#'   `cis_window` alone. Cost per block ~ O(reads x variants they cover).
#' * **between** (dosage): needs only per-sample residual sums and genotype
#'   dosages, i.e. O(samples x variants) per block - the same cost class as
#'   fastQTL. This is the component that runs out to `cis_window`
#'   (default 1 Mb, as in QTLtools).
#'
#' Cost guards: variants must pass `min_maf` (between) or `min_het_carriers`
#' (within); at most `max_variants_per_block` variants closest to the block
#' are tested, and the number dropped is reported per block.
#'
#' @section Methylation calls:
#' * `"hard_fixed"` (default): a call counts only if its modkit probability
#'   for the called state is >= `min_call_prob`; the SAME threshold for
#'   every sample, chemistry and allele. Ambiguous calls are dropped and
#'   their rate is reported per block and allele.
#' * `"hard_modkit"`: use modkit's per-sample estimated pass threshold.
#'   Convenient, but the threshold differs between samples/chemistries,
#'   which biases the between component.
#' * `"expected"`: no threshold; each call contributes its probability of
#'   being modified (fractional k). Keeps all data, but variance is then
#'   only approximately beta-binomial; use as sensitivity analysis.
#'
#' `mod_combine` decides how 5hmC (`h`) is counted:
#' `"m_plus_h"` (default; equals what bisulfite arrays measure, needed for
#' replication against array mQTL catalogues), `"m_only"` (h counted as
#' unmodified), `"drop_h"` (h calls removed).
#'
#' @param cis_window Max distance in bp between a block and a tested variant
#'   (between component and phased-within). Default 1e6.
#' @param within_reach `"read"` (direct allele on the read only) or
#'   `"phase_set"` (also HP-tag propagation inside the same phase set).
#' @param min_maf Minimum cohort MAF for the between component.
#' @param min_het_carriers Minimum heterozygous samples for the within
#'   component in cohort mode.
#' @param max_variants_per_block Hard cap on variants per block.
#' @param call_mode,min_call_prob,mod_combine See section above.
#' @param block_mode `"comethyl"` or `"regions"`; see [lq_define_blocks()].
#' @param block_args Named list of arguments for the block mode; missing
#'   entries take the mode's defaults ([lq_block_modes()]).
#' @param seed Global seed. Per-block seeds are derived from it
#'   deterministically, so results do not depend on the number of workers.
#'
#' @return An object of class `lq_params`.
#' @export
lq_params <- function(cis_window             = 1e6,
                      within_reach           = c("read", "phase_set"),
                      min_maf                = 0.05,
                      min_het_carriers       = 3L,
                      max_variants_per_block = 10000L,
                      call_mode              = c("hard_fixed", "hard_modkit", "expected"),
                      min_call_prob          = 0.8,
                      mod_combine            = c("m_plus_h", "m_only", "drop_h"),
                      block_mode             = c("comethyl", "regions"),
                      block_args             = list(),
                      seed                   = 42L) {

  within_reach <- match.arg(within_reach)
  call_mode    <- match.arg(call_mode)
  mod_combine  <- match.arg(mod_combine)
  block_mode   <- match.arg(block_mode)

  .check_scalar_num(cis_window, "cis_window", lower = 0)
  if (cis_window > 5e6) {
    stop("cis_window = ", format(cis_window, big.mark = ","), " bp exceeds 5 Mb. ",
         "Beyond that the test is no longer 'cis' in any usual sense and the ",
         "multiple-testing burden grows without biological benefit.")
  }
  .check_scalar_num(min_maf, "min_maf", lower = 0, upper = 0.5)
  .check_scalar_num(min_het_carriers, "min_het_carriers", lower = 1)
  .check_scalar_num(max_variants_per_block, "max_variants_per_block", lower = 1)
  .check_scalar_num(min_call_prob, "min_call_prob", lower = 0.5, upper = 1)
  .check_scalar_num(seed, "seed")
  if (!is.list(block_args)) stop("block_args must be a list.")

  block_args <- .resolve_block_args(block_mode, block_args)

  structure(
    list(
      cis_window             = as.numeric(cis_window),
      within_reach           = within_reach,
      min_maf                = min_maf,
      min_het_carriers       = as.integer(min_het_carriers),
      max_variants_per_block = as.integer(max_variants_per_block),
      call_mode              = call_mode,
      min_call_prob          = min_call_prob,
      mod_combine            = mod_combine,
      block_mode             = block_mode,
      block_args             = block_args,
      seed                   = as.integer(seed)
    ),
    class = "lq_params"
  )
}

#' @export
print.lq_params <- function(x, ...) {
  cat("<lq_params>\n")
  cat("  cis window       :", format(x$cis_window, big.mark = ",", scientific = FALSE),
      "bp (between); within reach =", x$within_reach, "\n")
  cat("  variant filters  : MAF >=", x$min_maf, "| het carriers >=", x$min_het_carriers,
      "| max", x$max_variants_per_block, "per block\n")
  cat("  calls            :", x$call_mode,
      if (x$call_mode == "hard_fixed") paste0("(p >= ", x$min_call_prob, ")") else "",
      "| 5hmC:", x$mod_combine, "\n")
  cat("  blocks           :", x$block_mode, "\n")
  for (nm in names(x$block_args)) {
    cat("    ", format(nm, width = 16), ":", format(x$block_args[[nm]]), "\n")
  }
  cat("  seed             :", x$seed, "\n")
  invisible(x)
}

#' Genomic range in which variants are searched for one block
#'
#' @param block_start,block_end Block coordinates (1-based, inclusive).
#' @param params An [lq_params()] object.
#' @return Named numeric vector `c(start, end)`; start is clipped at 1.
#' @export
lq_variant_window <- function(block_start, block_end, params) {
  stopifnot(inherits(params, "lq_params"))
  c(start = pmax(1, block_start - params$cis_window),
    end   = block_end + params$cis_window)
}

# ---- internal validation helpers ---------------------------------------------

.check_scalar_num <- function(x, name, lower = -Inf, upper = Inf) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x)) {
    stop(name, " must be a single non-missing number.", call. = FALSE)
  }
  if (x < lower || x > upper) {
    stop(name, " must be in [", lower, ", ", upper, "], got ", x, ".", call. = FALSE)
  }
  invisible(TRUE)
}
