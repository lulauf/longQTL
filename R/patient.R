# patient.R
#
# M2: single-patient engine (concept doc, sec. 6, level 1).
#
# Question per methylation block: do the two haplotypes of THIS person differ,
# and which variants could be responsible? Steps (lq_patient_scan):
#
#   1. calls -> 0/1 (or probabilities) by the call mode
#   2. blocks: co-methylation blocks from the person's own reads, or regions
#   3. read x block table: k methylated of n called CpGs per read
#   4. over-dispersion: beta-binomial rho per block, a genome trend of rho
#      over the methylation level, empirical-Bayes shrinkage to that trend
#   5. test units per block
#        haplotype  reads split by HP inside one phase set (PS)
#        variant    reads split by the allele they carry directly, for
#                   variants the haplotype test cannot see (unphased SNVs,
#                   SVs - Sniffles SVs are unphased)
#   6. per unit: beta-binomial likelihood-ratio test with rho fixed at the
#      shrunken value; label permutation of reads for units that look
#      significant (exact under exchangeable reads, no rho assumption)
#   7. per block, phasing-free: two-epiallele mixture on reads (bimodality)
#   8. flags: CpG-destroying variant, variant next to a CpG, haplotype read
#      imbalance, overlap with known (e.g. imprinted) regions; candidates =
#      het variants on the same phase set near the block
#   9. multiple testing (BH) and calibration (lambda_GC)

#' Read x block methylation counts for one sample
#'
#' @param x An `lq_reads` object ([lq_extract_reads()]).
#' @param blocks An `lq_blocks` object.
#' @param params [lq_params()] (call mode, 5hmC handling).
#' @param min_cpgs_per_read Reads with fewer called CpGs in a block are dropped
#'   for that block.
#' @return data.table: read_id, block_id, k, n, hp, ps, strand.
#' @export
lq_block_reads <- function(x, blocks, params = lq_params(), min_cpgs_per_read = 3L) {
  stopifnot(inherits(x, "lq_reads"), inherits(blocks, "lq_blocks"))
  cl <- lq_apply_call_mode(x$calls, params)
  cl <- merge(cl[, .(read_id, chr, pos, call)], blocks$cpg_map, by = c("chr", "pos"))
  br <- cl[, .(k = sum(call), n = .N), by = .(read_id, block_id)]
  br <- br[n >= min_cpgs_per_read]
  br <- merge(br, x$reads[, .(read_id, hp, ps, strand)], by = "read_id")
  br[]
}

#' Over-dispersion per block with empirical-Bayes shrinkage
#'
#' Per block an intercept-only beta-binomial gives (mu, rho). Blocks are
#' binned by mu; the median logit(rho) per bin is the trend (CpG coupling on
#' a molecule depends strongly on the methylation level). A block's own
#' estimate is shrunk to the trend with weight n_reads / (n_reads + prior_reads).
#'
#' @param br Output of [lq_block_reads()].
#' @param min_reads Blocks with fewer reads only get the trend.
#' @param n_bins Bins of the methylation level for the trend.
#' @param prior_reads Prior strength in reads.
#' @return list(blocks = data.table(block_id, n_reads, mu, rho_hat),
#'   trend = function(mu) -> rho, shrink = function(rho_hat, mu, n_reads)).
#' @export
lq_dispersion <- function(br, min_reads = 10L, n_bins = 10L, prior_reads = 20) {
  per <- br[, {
    if (.N >= min_reads) {
      f <- .bb_fit_null(k, n)
      list(n_reads = .N, mu = f$mu, rho_hat = f$rho)
    } else list(n_reads = .N, mu = sum(k) / sum(n), rho_hat = NA_real_)
  }, by = block_id]

  ok <- per[!is.na(rho_hat)]
  if (nrow(ok) >= 3L) {
    br_q <- unique(stats::quantile(ok$mu, seq(0, 1, length.out = n_bins + 1L), names = FALSE))
    ok[, bin := cut(mu, br_q, include.lowest = TRUE)]
    tr <- ok[, .(mu_c = stats::median(mu), lr = stats::median(stats::qlogis(pmin(pmax(rho_hat, 1e-4), 1 - 1e-4)))),
             by = bin][order(mu_c)]
    trend <- function(mu) {
      if (nrow(tr) == 1L) return(rep(stats::plogis(tr$lr), length(mu)))
      stats::plogis(stats::approx(tr$mu_c, tr$lr, xout = mu, rule = 2)$y)
    }
  } else {
    glob <- if (nrow(ok)) stats::median(ok$rho_hat) else 0.1
    trend <- function(mu) rep(glob, length(mu))
  }
  shrink <- function(rho_hat, mu, n_reads) {
    t <- stats::qlogis(trend(mu))
    if (is.na(rho_hat)) return(stats::plogis(t))
    w <- n_reads / (n_reads + prior_reads)
    stats::plogis(w * stats::qlogis(min(max(rho_hat, 1e-4), 1 - 1e-4)) + (1 - w) * t)
  }
  # prior SD of logit(rho) around the trend: robust spread of the block
  # estimates (includes their sampling noise, so the prior is if anything
  # too wide, never too narrow); floor 0.25
  prior_sd <- if (nrow(ok) >= 10L) {
    max(0.25, stats::mad(stats::qlogis(pmin(pmax(ok$rho_hat, 1e-4), 1 - 1e-4)) -
                           stats::qlogis(trend(ok$mu))))
  } else 1
  list(blocks = per, trend = trend, shrink = shrink, prior_sd = prior_sd)
}

#' Scan one sample for allele- / haplotype-specific methylation
#'
#' @param x An `lq_reads` object ([lq_extract_reads()]).
#' @param params [lq_params()]: call mode, 5hmC handling, block mode and
#'   arguments, seed.
#' @param blocks Optional `lq_blocks`; default: defined from this sample's
#'   reads with `params$block_mode`.
#' @param min_block_cpgs Blocks with fewer CpGs are not tested.
#' @param min_cpgs_per_read See [lq_block_reads()].
#' @param min_reads_per_group Minimum reads on each side of a test.
#' @param candidate_window bp around a block in which heterozygous variants
#'   are listed as candidates (haplotype unit) or tested (variant unit).
#' @param known_regions Optional data.frame chr, start, end, name (1-based)
#'   of regions to flag, e.g. known imprinted DMRs.
#' @param n_perm Maximum label permutations per test unit (sequential:
#'   most null units stop after a few hundred). 0 = LRT only.
#' @param perm_hits Stop permuting after this many permuted statistics at
#'   least as extreme as the observed one (Besag & Clifford 1991).
#' @param prior_reads Shrinkage strength for rho (see [lq_dispersion()]).
#' @param mask_cpg_snp If `TRUE` (default), CpG sites overlapped by a
#'   heterozygous small variant of this sample (CpG-SNP, or an indel covering
#'   the CpG) are removed from every read before testing, so the variant
#'   itself can never create the methylation difference. The block is still
#'   tested on its remaining CpGs; `n_cpg_masked` reports how many were
#'   removed, `cpg_snp` whether the block had any. If `FALSE` (behaviour of
#'   longQTL <= M2), such blocks are tested unmasked and never called.
#' @param fdr,min_delta A unit is called (`asm_call = TRUE`) when q < `fdr`
#'   and |methylation difference| >= `min_delta` (and, with
#'   `mask_cpg_snp = FALSE`, no variant destroys a CpG of the block).
#'
#' @return Object of class `lq_patient`: list with `units` (one row per test),
#'   `blocks` (block-level bimodality and QC), `dispersion`, `qc`.
#' @export
lq_patient_scan <- function(x, params = lq_params(), blocks = NULL,
                            min_block_cpgs = 3L, min_cpgs_per_read = 3L,
                            min_reads_per_group = 5L, candidate_window = 5e4,
                            known_regions = NULL, n_perm = 10000L, perm_hits = 20L,
                            fdr = 0.05, min_delta = 0.1,
                            prior_reads = 20, mask_cpg_snp = TRUE) {
  stopifnot(inherits(x, "lq_reads"), inherits(params, "lq_params"))

  # het variants available in this sample
  v <- x$variants
  het <- v[!is.na(gt_hp1) & !is.na(gt_hp2) & gt_hp1 != gt_hp2]

  # 1-2: blocks from this sample's own reads (hard calls) unless given
  if (is.null(blocks)) {
    cl <- lq_apply_call_mode(x$calls, params)
    cl[, call := round(call)]
    blocks <- lq_define_blocks(cl[, .(chr, pos, read_id, call)], params = params)
  }
  bl <- blocks$blocks[n_cpg >= min_block_cpgs]

  # CpG sites destroyed (or covered) by a het small variant are removed from
  # every read after block definition, so blocks stay the same with and
  # without masking and only the test data change
  masked <- .cpg_snp_sites(x$calls, het)
  if (mask_cpg_snp && nrow(masked)) {
    x$calls <- x$calls[!masked, on = c("chr", "pos")]
  }

  # 3: read x block counts
  br <- lq_block_reads(x, blocks, params, min_cpgs_per_read)
  br <- br[block_id %in% bl$block_id]

  # 4: dispersion
  disp <- lq_dispersion(br, prior_reads = prior_reads)
  mu_blk <- stats::setNames(disp$blocks$mu, disp$blocks$block_id)

  direct <- x$read_alleles[evidence %in% c("direct", "direct_indel", "direct_sv") & allele %in% 0:1]

  # 5-6: test units
  units <- list()
  by_block <- split(br, by = "block_id")
  for (b in names(by_block)) {
    d <- by_block[[b]]
    info <- bl[block_id == b]
    seed <- lq_block_seed(b, params$seed)

    # haplotype unit: reads with HP inside the majority phase set
    h <- d[hp %in% 1:2 & !is.na(ps)]
    if (nrow(h)) {
      ps_major <- as.integer(names(which.max(table(h$ps))))
      h <- h[ps == ps_major]
      u <- .test_unit(h$k, h$n, as.integer(h$hp == 2L), min_reads_per_group, disp,
                      mu_blk[[b]], n_perm, perm_hits, seed)
      if (!is.null(u)) {
        cand <- het[phased == TRUE & ps == ps_major & chr == info$chr &
                      pos >= info$start - candidate_window & pos <= info$end + candidate_window]
        cand[, dist := pmax(0, pos - info$end, info$start - pos)]
        data.table::setorder(cand, dist)
        units[[length(units) + 1L]] <- data.table::data.table(
          block_id = b, unit_type = "haplotype", unit_id = paste0("PS:", ps_major),
          group_0 = "HP1", group_1 = "HP2", u,
          n_candidates = nrow(cand),
          nearest_candidate_bp = if (nrow(cand)) cand$dist[1] else NA_real_,
          candidates = paste(utils::head(cand$variant_id, 5L), collapse = ";"))
      }
    }

    # variant units: het variants not covered by the haplotype test
    vv <- het[chr == info$chr & pos >= info$start - candidate_window & pos <= info$end + candidate_window &
                (phased == FALSE | is.na(ps))]
    for (vid in vv$variant_id) {
      a <- direct[variant_id == vid, .(read_id, allele)]
      dv <- merge(d, a, by = "read_id")
      if (!nrow(dv)) next
      u <- .test_unit(dv$k, dv$n, dv$allele, min_reads_per_group, disp, mu_blk[[b]],
                      n_perm, perm_hits, seed)
      if (!is.null(u)) {
        units[[length(units) + 1L]] <- data.table::data.table(
          block_id = b, unit_type = "variant", unit_id = vid, group_0 = "REF", group_1 = "ALT", u,
          n_candidates = 1L,
          nearest_candidate_bp = { p <- vv[variant_id == vid, pos]; max(0, p - info$end, info$start - p) },
          candidates = vid)
      }
    }
  }
  ut <- data.table::rbindlist(units, use.names = TRUE, fill = TRUE)

  # 7-8: block-level bimodality and flags
  bq <- .block_qc(br, bl, blocks, het, known_regions)
  nm <- if (nrow(masked)) {
    blocks$cpg_map[masked, on = c("chr", "pos"), nomatch = NULL][, .(n_cpg_masked = .N), by = block_id]
  } else data.table::data.table(block_id = character(), n_cpg_masked = integer())
  bq <- merge(bq, nm, by = "block_id", all.x = TRUE)
  bq[is.na(n_cpg_masked), n_cpg_masked := 0L]
  if (!mask_cpg_snp) bq[, n_cpg_masked := 0L]

  if (nrow(ut)) {
    ut <- merge(bl[, .(block_id, chr, start, end, n_cpg)], ut, by = "block_id")
    ut <- merge(ut, bq, by = "block_id", all.x = TRUE)
    # 9: multiple testing over all units
    ut[, q := stats::p.adjust(p_final, "BH")]
    # a call needs significance and a methylation difference of at least
    # min_delta; CpG-SNP sites are masked, or (mask_cpg_snp = FALSE) the
    # whole block is excluded
    ut[, asm_call := q < fdr & abs(delta) >= min_delta]
    if (!mask_cpg_snp) ut[cpg_snp == TRUE, asm_call := FALSE]
    data.table::setorder(ut, p_final)
  }
  hap_p <- if (nrow(ut)) ut[unit_type == "haplotype", p_lrt] else numeric()
  qc <- list(n_blocks = nrow(bl), n_blocks_tested = data.table::uniqueN(ut$block_id),
             n_units = nrow(ut),
             lambda_gc_haplotype = .lambda(hap_p),
             call_mode = params$call_mode, block_mode = attr(blocks, "diagnostics")$mode,
             mask_cpg_snp = mask_cpg_snp, n_cpg_masked = if (mask_cpg_snp) nrow(masked) else 0L)
  structure(list(units = ut, blocks = bq, dispersion = disp$blocks, qc = qc),
            class = "lq_patient")
}

#' @export
print.lq_patient <- function(x, ...) {
  q <- x$qc
  cat("<lq_patient>", q$n_blocks_tested, "of", q$n_blocks, "blocks tested,", q$n_units, "units\n")
  cat("  lambda_GC (haplotype units):", sprintf("%.2f", q$lambda_gc_haplotype), "\n")
  if (q$n_units) {
    cat("  q < 0.05:", sum(x$units$q < 0.05), " | ASM calls:", sum(x$units$asm_call), " | top units:\n")
    print(utils::head(x$units[, .(block_id, unit_type, unit_id, n0 = n_0, n1 = n_1,
                                  m0 = round(meth_0, 2), m1 = round(meth_1, 2),
                                  p = signif(p_final, 2), q = signif(q, 2))], 10))
  }
  invisible(x)
}

# One test: k/n per read, g in {0, 1}
#   p_lrt  : beta-binomial LRT, rho fixed at its empirical-Bayes value
#            (power, effect size, ranking)
#   p_perm : sequential label permutation of reads (Besag & Clifford 1991):
#            exact under exchangeable reads, no assumption on rho; stops after
#            `perm_hits` exceedances or `n_perm` permutations
#   p_final: p_perm, unless the permutation ran out before reaching
#            `perm_hits` (p below its resolution) - then the smaller of
#            p_perm and the LRT p-value
.test_unit <- function(k, n, g, min_reads, disp, mu, n_perm, perm_hits, seed) {
  n0 <- sum(g == 0L); n1 <- sum(g == 1L)
  if (n0 < min_reads || n1 < min_reads) return(NULL)
  free <- .bb_fit_free(k, n, g)
  rho <- disp$shrink(free$rho, mu, length(k))
  f <- .bb_fit_groups(k, n, g, rho)
  m0 <- sum(k[g == 0L]) / sum(n[g == 0L]); m1 <- sum(k[g == 1L]) / sum(n[g == 1L])
  pp <- if (n_perm > 0L) .perm_seq(k, n, g, abs(m1 - m0), n_perm, perm_hits, seed) else
    list(p = NA_real_, B = 0L, exhausted = NA)
  p_final <- if (is.na(pp$p)) f$p else if (pp$exhausted) min(pp$p, f$p) else pp$p
  data.table::data.table(n_0 = n0, n_1 = n1, meth_0 = m0, meth_1 = m1, delta = m1 - m0,
                         beta = f$beta, se = f$se, rho = rho, p_lrt = f$p, p_perm = pp$p,
                         n_perm_used = pp$B, p_final = p_final)
}

# Sequential permutation of group labels; statistic |pooled meth 1 - 0|
.perm_seq <- function(k, n, g, obs, n_max, h, seed, batch = 200L) {
  set.seed(seed)
  nr <- length(k); n1 <- sum(g == 1L)
  K <- sum(k); N <- sum(n)
  hits <- 0L; B <- 0L
  while (B < n_max && hits < h) {
    b <- min(batch, n_max - B)
    r <- matrix(stats::runif(nr * b), nr)
    G <- apply(r, 2L, rank) <= n1               # random subsets of size n1
    k1 <- colSums(G * k); s1 <- colSums(G * n)
    d <- abs(k1 / s1 - (K - k1) / (N - s1))
    hits <- hits + sum(d >= obs - 1e-12)
    B <- B + b
  }
  list(p = (hits + 1) / (B + 1), B = B, exhausted = hits < h)
}

# Penalised-likelihood ratio test: logit(rho) ~ N(t, s^2) prior in both the
# null (mu common) and the alternative (mu by group) model.
.bb_test_penalized <- function(k, n, g, t, s) {
  pen <- function(lr) 0.5 * ((lr - t) / s)^2
  m0 <- min(max(sum(k[g == 0]) / sum(n[g == 0]), 1e-3), 1 - 1e-3)
  m1 <- min(max(sum(k[g == 1]) / sum(n[g == 1]), 1e-3), 1 - 1e-3)
  mp <- min(max(sum(k) / sum(n), 1e-3), 1 - 1e-3)
  nll_alt <- function(p) -sum(.bb_ll(k, n, stats::plogis(p[1] + p[2] * g), stats::plogis(p[3]))) + pen(p[3])
  nll_nul <- function(p) -sum(.bb_ll(k, n, stats::plogis(p[1]), stats::plogis(p[2]))) + pen(p[2])
  fa <- stats::optim(c(stats::qlogis(m0), stats::qlogis(m1) - stats::qlogis(m0), t), nll_alt,
                     method = "BFGS", hessian = TRUE)
  fn <- stats::optim(c(stats::qlogis(mp), t), nll_nul, method = "BFGS")
  stat <- max(0, 2 * (fn$value - fa$value))
  se <- tryCatch(sqrt(diag(solve(fa$hessian)))[2], error = function(e) NA_real_)
  list(beta = fa$par[2], se = se, rho = stats::plogis(fa$par[3]),
       p = stats::pchisq(stat, 1, lower.tail = FALSE))
}

.bb_fit_free <- function(k, n, g) {
  nll <- function(p) -sum(.bb_ll(k, n, stats::plogis(p[1] + p[2] * g), stats::plogis(p[3])))
  m0 <- min(max(sum(k[g == 0]) / sum(n[g == 0]), 1e-3), 1 - 1e-3)
  m1 <- min(max(sum(k[g == 1]) / sum(n[g == 1]), 1e-3), 1 - 1e-3)
  fit <- stats::optim(c(stats::qlogis(m0), stats::qlogis(m1) - stats::qlogis(m0), stats::qlogis(0.1)),
                      nll, method = "BFGS")
  list(rho = stats::plogis(fit$par[3]))
}

.lambda <- function(p) {
  p <- p[!is.na(p) & p > 0]
  if (length(p) < 20L) return(NA_real_)
  stats::median(stats::qchisq(p, 1, lower.tail = FALSE)) / stats::qchisq(0.5, 1)
}

# CpG sites (strand-combined position of the C) overlapped by a het small
# variant: the variant hits the C (pos) or the G (pos + 1), or an indel
# spans either base
.cpg_snp_sites <- function(calls, het) {
  hs <- het[class == "small"]
  if (!nrow(hs) || !nrow(calls)) return(data.table::data.table(chr = character(), pos = integer()))
  cp <- unique(calls[, .(chr, pos)])
  cp[, `:=`(c_start = pos, c_end = pos + 1L)]
  hv <- hs[, .(chr, s = as.integer(pos), e = as.integer(pmax(end, pos)))]
  data.table::setkey(cp, chr, c_start, c_end)
  data.table::setkey(hv, chr, s, e)
  ov <- data.table::foverlaps(hv, cp, by.x = c("chr", "s", "e"), by.y = c("chr", "c_start", "c_end"),
                              nomatch = NULL)
  unique(ov[, .(chr, pos)])
}

# Block-level: bimodality, haplotype read balance, variant flags, known regions
.block_qc <- function(br, bl, blocks, het, known_regions, ctx_bp = 5L) {
  bim <- br[n >= 5L, {
    if (.N >= 10L) {
      m <- .bimodality(k, n)
      hp_ok <- hp %in% 1:2
      conc <- if (sum(hp_ok) >= 10L) {
        a <- mean((m$class[hp_ok] > 0.5) == (hp[hp_ok] == 2L)); max(a, 1 - a)
      } else NA_real_
      list(bimodal_p_low = m$p_low, bimodal_p_high = m$p_high, bimodal_pi_high = m$pi_high,
           bimodal_separation = m$separation, bimodal_delta_bic = m$delta_bic,
           bimodal_hp_concordance = conc)
    } else list(bimodal_p_low = NA_real_, bimodal_p_high = NA_real_, bimodal_pi_high = NA_real_,
                bimodal_separation = NA_real_, bimodal_delta_bic = NA_real_,
                bimodal_hp_concordance = NA_real_)
  }, by = block_id]

  bal <- br[hp %in% 1:2, .(n_hp1 = sum(hp == 1L), n_hp2 = sum(hp == 2L)), by = block_id]
  bal[, hap_balance_p := mapply(function(a, b) stats::binom.test(a, a + b, 0.5)$p.value, n_hp1, n_hp2)]

  # variant inside a CpG of the block (CpG-SNP) or within ctx_bp of one
  cm <- blocks$cpg_map[block_id %in% bl$block_id]
  hs <- het[class == "small", .(chr, v_start = pos, v_end = end)]
  flags <- if (nrow(hs) && nrow(cm)) {
    cm2 <- cm[, .(chr, c_start = pos - ctx_bp, c_end = pos + 1L + ctx_bp, pos, block_id)]
    data.table::setkey(cm2, chr, c_start, c_end)
    hs[, `:=`(s = v_start, e = v_end)]
    data.table::setkey(hs, chr, s, e)
    ov <- data.table::foverlaps(hs, cm2, by.x = c("chr", "s", "e"), by.y = c("chr", "c_start", "c_end"),
                                nomatch = NULL)
    ov[, in_cpg := v_start <= pos + 1L & v_end >= pos]
    ov[, .(cpg_snp = any(in_cpg), context_snp = any(!in_cpg)), by = block_id]
  } else data.table::data.table(block_id = character(), cpg_snp = logical(), context_snp = logical())

  out <- merge(bl[, .(block_id)], bim, by = "block_id", all.x = TRUE)
  out <- merge(out, bal, by = "block_id", all.x = TRUE)
  out <- merge(out, flags, by = "block_id", all.x = TRUE)
  out[is.na(cpg_snp), cpg_snp := FALSE]
  out[is.na(context_snp), context_snp := FALSE]

  out[, known_region := NA_character_]
  if (!is.null(known_regions) && nrow(known_regions)) {
    kr <- data.table::as.data.table(known_regions)
    if (!"name" %in% names(kr)) kr[, name := paste0(chr, ":", start, "-", end)]
    b2 <- merge(out[, .(block_id)], bl[, .(block_id, chr, start, end)], by = "block_id")
    for (i in seq_len(nrow(kr))) {
      hit <- b2[chr == kr$chr[i] & end >= kr$start[i] & start <= kr$end[i], block_id]
      out[block_id %in% hit, known_region := kr$name[i]]
    }
  }
  out[]
}
