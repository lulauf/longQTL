# simulate_cohort.R
#
# Read-level cohort simulator: the ground truth against which every modelling
# choice is measured (concept doc, sec. 9). One locus = one methylation block
# plus a handful of variants around it. Generative model per read r of sample
# i, matching the analysis model of sec. 4:
#
#   logit mu_ir = alpha + beta_W (A_ir - Abar_i) + beta_B Abar_i
#                 + gamma[stratum_i] + u_i + v_i (A_ir - Abar_i)
#                 + imprint * [read comes from the maternal haplotype]
#   p_ir  ~ Beta(mean mu_ir, overdispersion rho)      (molecule-level state)
#   CpG state ~ Bernoulli(p_ir), independently per CpG given p_ir
#   modkit-like probability of "modified" drawn around the true state with a
#   stratum-specific sharpness (R9 calls are less confident than R10)
#
# A_ir is the read's allele at the causal variant, Abar_i = dosage / 2.
# Haplotype 1 of every sample is maternal, haplotype 2 paternal; alleles are
# drawn independently of origin, so imprinting has a random sign relative to
# ALT across heterozygous people - exactly the situation the tau term and the
# haplotype sign-flip null must handle.

#' Define a simulation scenario
#'
#' @param n_samples Number of individuals.
#' @param coverage Mean reads per sample overlapping the block (Poisson).
#' @param read_len_mean Mean read length in bp (gamma, shape 2).
#' @param n_cpg,cpg_spacing CpGs in the block and their spacing (bp).
#' @param block_start Position of the first CpG.
#' @param chr Chromosome name.
#' @param variant_offsets Variant positions relative to `block_start`.
#' @param causal Index of the causal variant in `variant_offsets`; 0 = no
#'   causal variant (global null for genetic effects).
#' @param maf Minor allele frequency of the anchor (causal) variant.
#' @param ld Probability that a non-anchor variant copies the anchor allele on
#'   a haplotype (crude LD control; 0 = independent, 1 = perfect LD).
#' @param base_meth Baseline methylation level (probability scale).
#' @param beta_W,beta_B Within (allelic) and between (dosage) effects, logit
#'   scale. Ignored when `causal = 0`.
#' @param sigma_u SD of the sample random intercept.
#' @param tau SD of the sample random allelic slope.
#' @param imprint Logit shift for reads from the maternal haplotype.
#' @param rho Read-level overdispersion (intra-molecule correlation of CpG
#'   states); 0 = independent CpGs.
#' @param strata Named probabilities of technical strata.
#' @param stratum_effect Named logit shifts per stratum.
#' @param call_sharpness Named list per stratum: Beta(a, b) parameters of the
#'   probability assigned to the TRUE state (higher a/b = sharper calls).
#' @param switch_error Probability that a phased (non-direct) allele
#'   assignment is flipped.
#' @param mapping_bias Probability that a read carrying the causal ALT allele
#'   and spanning it is lost (reference bias).
#' @param seed RNG seed.
#' @return Object of class `lq_scenario`.
#' @export
lq_sim_scenario <- function(n_samples       = 50L,
                            coverage        = 30,
                            read_len_mean   = 20000,
                            n_cpg           = 10L,
                            cpg_spacing     = 50L,
                            block_start     = 1e6,
                            chr             = "chrSim",
                            variant_offsets = c(-30000, -3000, 325, 4000, 40000),
                            causal          = 3L,
                            maf             = 0.3,
                            ld              = 0.8,
                            base_meth       = 0.5,
                            beta_W          = 0,
                            beta_B          = 0,
                            sigma_u         = 0.3,
                            tau             = 0,
                            imprint         = 0,
                            rho             = 0.1,
                            strata          = c(R9 = 0.5, R10 = 0.5),
                            stratum_effect  = c(R9 = 0, R10 = 0.2),
                            call_sharpness  = list(R9 = c(6, 1.5), R10 = c(12, 1.2)),
                            switch_error    = 0.01,
                            mapping_bias    = 0,
                            seed            = 1L) {
  .check_scalar_num(n_samples, "n_samples", lower = 1)
  .check_scalar_num(coverage, "coverage", lower = 0)
  .check_scalar_num(causal, "causal", lower = 0, upper = length(variant_offsets))
  .check_scalar_num(maf, "maf", lower = 0, upper = 1)
  .check_scalar_num(ld, "ld", lower = 0, upper = 1)
  .check_scalar_num(base_meth, "base_meth", lower = 1e-6, upper = 1 - 1e-6)
  .check_scalar_num(rho, "rho", lower = 0, upper = 0.999)
  .check_scalar_num(switch_error, "switch_error", lower = 0, upper = 0.5)
  .check_scalar_num(mapping_bias, "mapping_bias", lower = 0, upper = 1)
  if (!setequal(names(strata), names(stratum_effect)) ||
      !setequal(names(strata), names(call_sharpness))) {
    stop("strata, stratum_effect and call_sharpness must have the same names.")
  }
  structure(as.list(environment()), class = "lq_scenario")
}

#' @export
print.lq_scenario <- function(x, ...) {
  cat("<lq_scenario>", x$n_samples, "samples,", x$coverage, "reads/sample,",
      x$n_cpg, "CpGs\n")
  cat("  causal variant:", if (x$causal == 0) "none (null)" else
        paste0("#", x$causal, " at offset ", x$variant_offsets[x$causal]), "\n")
  cat("  beta_W =", x$beta_W, " beta_B =", x$beta_B, " tau =", x$tau,
      " imprint =", x$imprint, " rho =", x$rho, "\n")
  invisible(x)
}

#' Simulate a cohort at one locus, read by read
#'
#' @param scenario An [lq_sim_scenario()] object.
#' @param seed Optional seed overriding `scenario$seed`.
#' @param stratum Optional character vector (one per sample) that fixes the
#'   technical stratum of every sample, e.g. to keep strata constant across
#'   loci of one simulated cohort. The random draw still happens, so other
#'   random numbers do not change.
#' @return Object of class `lq_sim`, a list of data.tables:
#'   * `calls`: sample_id, read_id, chr, pos, strand, prob_mod, true_state
#'   * `reads`: sample_id, read_id, hap (1 maternal / 2 paternal), start, end,
#'     strand, allele_causal
#'   * `read_alleles`: sample_id, read_id, variant_id, allele (0/1),
#'     evidence (`direct` = read spans the variant, `phased` = via haplotype)
#'   * `variants`: variant_id, chr, pos, is_causal
#'   * `genotypes`: sample_id, variant_id, gt_hp1, gt_hp2, dosage
#'   * `samples`: sample_id, stratum, u, v
#'   and the `scenario`.
#' @export
lq_simulate_cohort <- function(scenario, seed = NULL, stratum = NULL) {
  stopifnot(inherits(scenario, "lq_scenario"))
  sc <- scenario
  set.seed(if (is.null(seed)) sc$seed else seed)

  S <- as.integer(sc$n_samples)
  V <- length(sc$variant_offsets)
  anchor <- if (sc$causal > 0) as.integer(sc$causal) else as.integer(ceiling(V / 2))
  sample_id <- sprintf("S%04d", seq_len(S))

  # ---- samples ----------------------------------------------------------------
  fixed_stratum <- stratum
  stratum <- sample(names(sc$strata), S, replace = TRUE, prob = sc$strata)
  if (!is.null(fixed_stratum)) {
    if (length(fixed_stratum) != S || !all(fixed_stratum %in% names(sc$strata))) {
      stop("stratum must have one entry per sample, each one of names(scenario$strata).")
    }
    stratum <- as.character(fixed_stratum)
  }
  samples <- data.table::data.table(sample_id = sample_id, stratum = stratum,
                                    u = stats::rnorm(S, 0, sc$sigma_u),
                                    v = stats::rnorm(S, 0, sc$tau))

  # ---- haplotypes (2S x V), row 2i-1 maternal, 2i paternal ---------------------
  anc <- stats::rbinom(2L * S, 1L, sc$maf)
  hap <- matrix(0L, 2L * S, V)
  for (v in seq_len(V)) {
    hap[, v] <- if (v == anchor) anc else
      ifelse(stats::runif(2L * S) < sc$ld, anc, stats::rbinom(2L * S, 1L, sc$maf))
  }
  var_pos <- as.integer(sc$block_start + sc$variant_offsets)
  variants <- data.table::data.table(variant_id = sprintf("%s:%d", sc$chr, var_pos),
                                     chr = sc$chr, pos = var_pos,
                                     is_causal = seq_len(V) == sc$causal)
  genotypes <- data.table::data.table(
    sample_id  = rep(sample_id, each = V),
    variant_id = rep(variants$variant_id, times = S),
    gt_hp1     = as.vector(t(hap[seq(1L, 2L * S, 2L), , drop = FALSE])),
    gt_hp2     = as.vector(t(hap[seq(2L, 2L * S, 2L), , drop = FALSE]))
  )
  genotypes[, dosage := gt_hp1 + gt_hp2]

  # ---- reads ---------------------------------------------------------------------
  cpg_pos <- as.integer(sc$block_start + (seq_len(sc$n_cpg) - 1L) * sc$cpg_spacing)
  block_end <- max(cpg_pos)
  n_reads <- stats::rpois(S, sc$coverage)
  R <- sum(n_reads)
  if (R == 0L) stop("Scenario produced no reads (coverage too low).")
  rs <- rep(seq_len(S), n_reads)
  h  <- sample.int(2L, R, replace = TRUE)
  len <- pmax(1000, stats::rgamma(R, shape = 2, scale = sc$read_len_mean / 2))
  start <- as.integer(stats::runif(R, sc$block_start - len, block_end + 1))
  end   <- as.integer(start + len)
  hap_row <- 2L * (rs - 1L) + h
  allele_c <- hap[cbind(hap_row, anchor)]

  # reference bias: ALT reads spanning the anchor variant are lost with prob
  spans_anchor <- start <= var_pos[anchor] & end >= var_pos[anchor]
  keep <- !(allele_c == 1L & spans_anchor & stats::runif(R) < sc$mapping_bias)
  rs <- rs[keep]; h <- h[keep]; start <- start[keep]; end <- end[keep]
  hap_row <- hap_row[keep]; allele_c <- allele_c[keep]
  R <- length(rs)

  # ---- read-level methylation probability ---------------------------------------
  dos <- hap[2L * rs - 1L, anchor] + hap[2L * rs, anchor]
  abar <- dos / 2
  dev  <- allele_c - abar
  eff_W <- if (sc$causal > 0) sc$beta_W else 0
  eff_B <- if (sc$causal > 0) sc$beta_B else 0
  eta <- stats::qlogis(sc$base_meth) + eff_W * dev + eff_B * abar +
    sc$stratum_effect[stratum[rs]] + samples$u[rs] + samples$v[rs] * dev +
    sc$imprint * (h == 1L)
  mu <- stats::plogis(eta)
  p_r <- if (sc$rho > 0) {
    phi <- 1 / sc$rho - 1
    stats::rbeta(R, mu * phi, (1 - mu) * phi)
  } else mu

  read_id <- sprintf("%s_r%05d", sample_id[rs], stats::ave(rs, rs, FUN = seq_along))
  strand <- sample(c("+", "-"), R, replace = TRUE)
  reads <- data.table::data.table(sample_id = sample_id[rs], read_id = read_id, hap = h,
                                  start = start, end = end, strand = strand,
                                  allele_causal = allele_c)

  # ---- CpG calls -------------------------------------------------------------------
  cover <- outer(start, cpg_pos, "<=") & outer(end, cpg_pos, ">=")
  idx <- which(cover, arr.ind = TRUE)
  r_i <- idx[, 1L]; c_i <- idx[, 2L]
  state <- stats::rbinom(length(r_i), 1L, p_r[r_i])
  st <- stratum[rs[r_i]]
  conf <- numeric(length(r_i))
  for (s in names(sc$call_sharpness)) {
    w <- which(st == s)
    ab <- sc$call_sharpness[[s]]
    conf[w] <- stats::rbeta(length(w), ab[1L], ab[2L])
  }
  prob_mod <- ifelse(state == 1L, conf, 1 - conf)
  calls <- data.table::data.table(sample_id = sample_id[rs[r_i]], read_id = read_id[r_i],
                                  chr = sc$chr, pos = cpg_pos[c_i], strand = strand[r_i],
                                  prob_mod = prob_mod, true_state = state)
  data.table::setorder(calls, sample_id, read_id, pos)

  # ---- read alleles at every variant ------------------------------------------------
  ra <- data.table::data.table(
    sample_id  = rep(sample_id[rs], each = V),
    read_id    = rep(read_id, each = V),
    variant_id = rep(variants$variant_id, times = R),
    allele     = as.vector(t(hap[hap_row, , drop = FALSE])),
    direct     = as.vector(t(outer(start, var_pos, "<=") & outer(end, var_pos, ">=")))
  )
  flip <- !ra$direct & stats::runif(nrow(ra)) < sc$switch_error
  ra[flip, allele := 1L - allele]
  ra[, evidence := data.table::fifelse(direct, "direct", "phased")]
  ra[, direct := NULL]

  structure(list(calls = calls, reads = reads, read_alleles = ra, variants = variants,
                 genotypes = genotypes, samples = samples, scenario = sc),
            class = "lq_sim")
}

#' One simulated sample as an `lq_reads` object
#'
#' Lets every function that works on extracted real data (e.g.
#' [lq_patient_scan()]) run on simulated data with known truth. All reads are
#' haplotagged (HP = simulated haplotype, one phase set); variants are phased.
#'
#' @param sim An `lq_sim` object.
#' @param sample Sample id (default: the first sample).
#' @return An `lq_reads` object.
#' @export
lq_sim_as_reads <- function(sim, sample = NULL) {
  stopifnot(inherits(sim, "lq_sim"))
  s <- if (is.null(sample)) sim$samples$sample_id[1] else sample
  reads <- sim$reads[sample_id == s, .(read_id, chr = sim$scenario$chr, start, end, strand,
                                       mapq = 60L, flag = 0L, hp = hap, ps = 1L,
                                       read_len = end - start + 1L)]
  calls <- sim$calls[sample_id == s, .(read_id, chr, pos, strand, prob_mod)]
  reads[, n_calls := calls[, .N, by = read_id][reads, on = "read_id", x.N]]
  reads[is.na(n_calls), n_calls := 0L]
  g <- sim$genotypes[sample_id == s]
  variants <- merge(sim$variants, g, by = "variant_id")
  variants <- variants[, .(variant_id, chr, pos, end = pos, ref = "A", alt = "G",
                           gt = paste0(gt_hp1, "|", gt_hp2), gt_hp1, gt_hp2,
                           phased = TRUE, ps = 1L, class = "small")]
  ra <- sim$read_alleles[sample_id == s, .(read_id, variant_id, allele, evidence, bq = NA_integer_)]
  structure(list(reads = reads, calls = calls, read_alleles = ra, variants = variants,
                 region = paste0(sim$scenario$chr, ":", min(calls$pos), "-", max(calls$pos)),
                 qc = list(n_reads = nrow(reads), n_calls = nrow(calls))),
            class = "lq_reads")
}

#' Simulate one person with many loci as an `lq_reads` object
#'
#' Every locus is an independent [lq_simulate_cohort()] draw for one sample,
#' placed 100 kb apart on one chromosome, so genome-wide steps (dispersion
#' trend, multiple testing) see many blocks as in real data.
#'
#' @param scenarios A list of [lq_sim_scenario()] objects (one per locus) or
#'   one scenario recycled `n_loci` times.
#' @param n_loci Number of loci if `scenarios` is a single scenario.
#' @param seed Base seed; locus `i` uses `seed + i`.
#' @return An `lq_reads` object with attribute `truth`: one row per locus
#'   (locus, block_start, block_end, imprint, beta_W, het_causal).
#' @export
lq_sim_person <- function(scenarios, n_loci = 100L, seed = 1L) {
  if (inherits(scenarios, "lq_scenario")) scenarios <- rep(list(scenarios), n_loci)
  parts <- lapply(seq_along(scenarios), function(i) {
    sc <- scenarios[[i]]
    sc$n_samples <- 1L
    sc$block_start <- 1e6 + (i - 1) * 1e5
    sc$chr <- "chrSim"
    s <- lq_simulate_cohort(sc, seed = seed + i)
    x <- lq_sim_as_reads(s)
    pre <- paste0("L", i, "_")
    x$reads[, read_id := paste0(pre, read_id)]
    x$calls[, read_id := paste0(pre, read_id)]
    x$read_alleles[, read_id := paste0(pre, read_id)]
    x$variants[, variant_id := paste0(pre, variant_id)]
    x$read_alleles[, variant_id := paste0(pre, variant_id)]
    cv <- s$variants[is_causal == TRUE, variant_id]
    het_c <- if (length(cv)) s$genotypes[variant_id == cv, gt_hp1 != gt_hp2] else NA
    x$truth <- data.table::data.table(locus = i, block_start = min(x$calls$pos), block_end = max(x$calls$pos),
                                      imprint = sc$imprint, beta_W = sc$beta_W, het_causal = het_c)
    x
  })
  out <- structure(list(
    reads = data.table::rbindlist(lapply(parts, `[[`, "reads")),
    calls = data.table::rbindlist(lapply(parts, `[[`, "calls")),
    read_alleles = data.table::rbindlist(lapply(parts, `[[`, "read_alleles")),
    variants = data.table::rbindlist(lapply(parts, `[[`, "variants")),
    region = "chrSim", qc = list()), class = "lq_reads")
  data.table::setattr(out, "truth", data.table::rbindlist(lapply(parts, `[[`, "truth")))
  out
}

#' @export
print.lq_sim <- function(x, ...) {
  cat("<lq_sim>", nrow(x$samples), "samples,", nrow(x$reads), "reads,",
      nrow(x$calls), "CpG calls,", nrow(x$variants), "variants\n")
  cv <- x$variants[is_causal == TRUE, variant_id]
  if (length(cv)) {
    g <- x$genotypes[variant_id == cv, table(dosage)]
    cat("  causal", cv, "genotypes 0/1/2:", paste(g, collapse = "/"), "\n")
  }
  invisible(x)
}


#' Simulate a cohort over many loci
#'
#' Every locus is an independent [lq_simulate_cohort()] draw for the same
#' `n_samples` people, placed `spacing` bp apart on one chromosome. Technical
#' strata are drawn once per person and kept across loci. The result plugs
#' into the cohort functions through [lq_io_sim()], so reference blocks and
#' the read-level store can be tested against known truth without BAM files.
#'
#' @param scenarios A list of [lq_sim_scenario()] objects (one per locus) or
#'   one scenario recycled `n_loci` times.
#' @param n_loci Number of loci if `scenarios` is a single scenario.
#' @param n_samples People in the cohort (overrides the scenarios).
#' @param seed Base seed; locus `i` uses `seed + i`.
#' @param spacing Distance between loci in bp.
#' @return Object of class `lq_sim_loci`: `samples` (sample_id, stratum),
#'   `reads` (named list of one `lq_reads` object per sample), `truth` (one
#'   row per locus).
#' @export
lq_sim_cohort_loci <- function(scenarios, n_loci = 20L, n_samples = 10L, seed = 1L,
                               spacing = 1e5) {
  if (inherits(scenarios, "lq_scenario")) scenarios <- rep(list(scenarios), n_loci)
  set.seed(seed)
  st <- scenarios[[1]]$strata
  sid <- sprintf("S%04d", seq_len(n_samples))
  strata <- sample(names(st), n_samples, replace = TRUE, prob = st)
  
  per_locus <- lapply(seq_along(scenarios), function(i) {
    sc <- scenarios[[i]]
    sc$n_samples <- as.integer(n_samples)
    sc$block_start <- 1e6 + (i - 1) * spacing
    sc$chr <- "chrSim"
    sim <- lq_simulate_cohort(sc, seed = seed + i, stratum = strata)
    cv <- sim$variants[is_causal == TRUE, variant_id]
    list(sim = sim,
         truth = data.table::data.table(
           locus = i, block_start = sc$block_start,
           block_end = sc$block_start + (sc$n_cpg - 1L) * sc$cpg_spacing,
           causal = if (length(cv)) cv else NA_character_, beta_W = sc$beta_W,
           beta_B = sc$beta_B, imprint = sc$imprint, base_meth = sc$base_meth))
  })
  
  reads <- lapply(sid, function(s) {
    parts <- lapply(seq_along(per_locus), function(i) {
      x <- lq_sim_as_reads(per_locus[[i]]$sim, s)
      pre <- paste0(s, ":L", i, "_")
      x$reads[, read_id := paste0(pre, read_id)]
      x$calls[, read_id := paste0(pre, read_id)]
      x$read_alleles[, read_id := paste0(pre, read_id)]
      x
    })
    out <- structure(list(
      reads = data.table::rbindlist(lapply(parts, `[[`, "reads")),
      calls = data.table::rbindlist(lapply(parts, `[[`, "calls")),
      read_alleles = data.table::rbindlist(lapply(parts, `[[`, "read_alleles")),
      variants = data.table::rbindlist(lapply(parts, `[[`, "variants")),
      region = "chrSim", qc = list()), class = "lq_reads")
    out$reads[, sample_id := s]
    out$calls[, sample_id := s]
    out
  })
  names(reads) <- sid
  structure(list(samples = data.table::data.table(sample_id = sid, stratum = strata),
                 reads = reads,
                 truth = data.table::rbindlist(lapply(per_locus, `[[`, "truth"))),
            class = "lq_sim_loci")
}