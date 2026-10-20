# cohort.R
#
# M3: everything the cohort scan (M4) needs, built once per cohort.
#
#   lq_cohort()            sample sheet + input/output hooks (BAM/VCF or
#                          simulator), technical stratum per sample
#   lq_reference_blocks()  pass A: per genome window and sample, CpG
#                          coverage, CpG-pair co-methylation and variants in
#                          CpGs; pooled over the cohort into ONE block set
#                          that every sample shares. CpGs destroyed by a
#                          variant in any sample are masked for everybody,
#                          so the phenotype is the same molecule set in
#                          carriers and non-carriers.
#   lq_sv_catalogue()      SVs merged across samples into shared ids
#   lq_cohort_store()      pass B: per chunk of reference blocks and sample,
#                          the sufficient statistics of the read model
#                          (k, n, strand, HP, PS per read and block), the
#                          read alleles with direct evidence, and sample x
#                          block summaries
#   lq_cohort_genotypes()  per-sample genotype tables for the between
#                          (dosage) part, SVs under catalogue ids
#   lq_stratum_check()     do technical strata (basecaller versions) differ
#                          more than people do?
#
# Block definition never sees alleles (see blocks.R); pooling over people
# does not change that.

# ---- cohort object -------------------------------------------------------------

#' Define a cohort
#'
#' @param samples data.frame with one row per sample: `sample_id` and, for
#'   the default BAM/VCF input, `bam` (BAM or CRAM), optionally `snv`
#'   (small-variant VCF), `sv` (Sniffles2 VCF), `reference` (FASTA, needed for
#'   CRAM), `vcf_sample` (sample column in multi-sample VCFs). `stratum`
#'   (technical stratum) is taken from the sheet if present, otherwise
#'   derived with [lq_bam_provenance()]. Further columns (tissue, sex, ...)
#'   are kept.
#' @param io Input hooks, default [lq_io_bam()]; [lq_io_sim()] for simulated
#'   cohorts.
#' @param check_files Check that all listed files exist.
#' @return Object of class `lq_cohort`.
#' @export
lq_cohort <- function(samples, io = lq_io_bam(), check_files = TRUE) {
  s <- data.table::as.data.table(samples)
  if (!"sample_id" %in% names(s)) stop("samples needs a sample_id column.")
  if (anyDuplicated(s$sample_id)) stop("Duplicated sample_id: ", s$sample_id[duplicated(s$sample_id)][1])
  s[, sample_id := as.character(sample_id)]
  if (isTRUE(io$needs_files)) {
    if (!"bam" %in% names(s)) stop("samples needs a bam column for BAM/VCF input.")
    if (check_files) {
      for (col in intersect(c("bam", "snv", "sv", "reference"), names(s))) {
        f <- s[[col]]
        miss <- !is.na(f) & nzchar(f) & !file.exists(f)
        if (any(miss)) stop(col, " file not found for ", s$sample_id[miss][1], ": ", f[miss][1])
      }
    }
  }
  if (!"stratum" %in% names(s)) {
    if ("bam" %in% names(s)) {
      prov <- lq_bam_provenance(s$bam, check_tags = FALSE,
                                reference = if ("reference" %in% names(s)) s$reference[1] else NULL)
      s[, stratum := prov$stratum]
    } else s[, stratum := "all"]
  }
  s[is.na(stratum) | !nzchar(stratum), stratum := "unknown"]
  structure(list(samples = s, io = io), class = "lq_cohort")
}

#' @export
print.lq_cohort <- function(x, ...) {
  cat("<lq_cohort>", nrow(x$samples), "samples, input:", x$io$name, "\n")
  print(x$samples[, .N, by = stratum])
  invisible(x)
}

#' Input hooks: haplotagged modBAM/modCRAM plus per-sample VCFs
#'
#' @param ... Passed to [lq_extract_reads()] (e.g. `min_mapq`,
#'   `min_base_qual`, `sv_window`).
#' @return A list of functions used by the cohort passes.
#' @export
lq_io_bam <- function(...) {
  ex_args <- list(...)
  col <- function(s, k) if (k %in% names(s) && !is.na(s[[k]]) && nzchar(s[[k]])) s[[k]] else NULL
  list(
    name = "bam",
    needs_files = TRUE,
    variants = function(s, region, type = c("small", "sv")) {
      vs <- if ("vcf_sample" %in% names(s) && !is.na(s$vcf_sample)) s$vcf_sample else 1L
      small <- if ("small" %in% type && !is.null(col(s, "snv"))) {
        lq_read_variants(col(s, "snv"), region, sample = vs, type = "small")
      } else .empty_variants()
      sv <- if ("sv" %in% type && !is.null(col(s, "sv"))) {
        lq_read_variants(col(s, "sv"), region, sample = vs, type = "sv")
      } else NULL
      list(small = small, sv = sv)
    },
    reads = function(s, region, small = NULL, sv = NULL, phased_max_dist = 0) {
      args <- c(list(bam = s$bam, region = region, small_variants = small, sv = sv,
                     sample_id = s$sample_id, reference = col(s, "reference"),
                     phased_max_dist = phased_max_dist), ex_args)
      suppressWarnings(do.call(lq_extract_reads, args))
    }
  )
}

#' Input hooks for a simulated cohort
#'
#' @param sim An `lq_sim_loci` object from [lq_sim_cohort_loci()].
#' @return A list of functions used by the cohort passes.
#' @export
lq_io_sim <- function(sim) {
  stopifnot(inherits(sim, "lq_sim_loci"))
  list(
    name = "simulation",
    needs_files = FALSE,
    variants = function(s, region, type = c("small", "sv")) {
      rg <- .parse_region(region)
      v <- sim$reads[[s$sample_id]]$variants[chr == rg$chr & pos >= rg$start & pos <= rg$end]
      v <- v[, .(variant_id, chr, pos, ref, alt, gt, gt_hp1, gt_hp2, phased, ps)]
      v[, multiallelic := rep(FALSE, .N)]
      list(small = v, sv = NULL)
    },
    reads = function(s, region, small = NULL, sv = NULL, phased_max_dist = 0) {
      rg <- .parse_region(region)
      x <- sim$reads[[s$sample_id]]
      r <- x$reads[chr == rg$chr & start <= rg$end & end >= rg$start]
      cl <- x$calls[read_id %in% r$read_id & pos >= rg$start & pos <= rg$end]
      vid <- if (is.null(small)) character() else small$variant_id
      ra <- x$read_alleles[read_id %in% r$read_id & variant_id %in% vid]
      ra <- ra[evidence == "direct" | phased_max_dist > 0]
      v <- if (is.null(small)) x$variants[0] else x$variants[variant_id %in% vid]
      structure(list(reads = r, calls = cl, read_alleles = ra, variants = v, region = region,
                     qc = list(n_reads = nrow(r), n_calls = nrow(cl))),
                class = "lq_reads")
    }
  )
}

# ---- genome windows ------------------------------------------------------------

#' Chromosome sizes from a BAM/CRAM header
#'
#' @param bam Path to a BAM or CRAM.
#' @param chroms Chromosomes to keep (default: chr1-22, chrX), in this order.
#' @return data.table chr, start = 1, end = length.
#' @export
lq_chrom_sizes <- function(bam, chroms = paste0("chr", c(1:22, "X"))) {
  h <- .read_header(bam)
  sq <- grep("^@SQ", h, value = TRUE)
  field <- function(key) vapply(strsplit(sq, "\t", fixed = TRUE), function(f) {
    v <- f[startsWith(f, key)]
    if (length(v)) substring(v[1L], nchar(key) + 1L) else NA_character_
  }, character(1))
  out <- data.table::data.table(chr = field("SN:"), end = as.integer(field("LN:")))
  out[, start := 1L]
  if (!is.null(chroms)) out <- out[match(chroms, chr)][!is.na(chr)]
  out[, .(chr, start, end)]
}

#' Cut regions into fixed windows for the reference pass
#'
#' @param regions Character vector "chr:start-end" or data.frame chr, start,
#'   end (1-based, inclusive), e.g. from [lq_chrom_sizes()].
#' @param window Window size in bp (default 5 Mb).
#' @return data.table in the chunk format of [lq_plan_chunks()].
#' @export
lq_cohort_windows <- function(regions, window = 5e6) {
  .check_scalar_num(window, "window", lower = 1e3)
  rg <- if (is.character(regions)) {
    data.table::rbindlist(lapply(regions, function(r) data.table::as.data.table(.parse_region(r))))
  } else data.table::as.data.table(regions)[, .(chr, start = as.integer(start), end = as.integer(end))]
  w <- rg[, {
    s <- seq(start, end, by = window)
    .(start = as.integer(s), end = as.integer(pmin(s + window - 1, end)))
  }, by = .(chr, r_start = start)]
  w[, r_start := NULL]
  w[, `:=`(var_start = start, var_end = end, n_blocks = NA_integer_, work = end - start + 1)]
  w[, chunk_id := sprintf("win_%05d", seq_len(.N))]
  data.table::setcolorder(w, c("chunk_id", "chr", "start", "end", "var_start", "var_end", "n_blocks", "work"))
  w[]
}

# ---- pass A: pooled CpG and pair statistics -------------------------------------

.passA_worker <- function(chunk, params, cohort, halo, min_cpg_reads) {
  core <- c(chunk$start, chunk$end)
  reg_core <- sprintf("%s:%d-%d", chunk$chr, core[1], core[2])
  reg_ext <- sprintf("%s:%d-%d", chunk$chr, core[1], core[2] + halo)
  comethyl <- params$block_mode == "comethyl"
  a <- params$block_args
  S <- cohort$samples

  cpg <- pairs <- sites <- info <- vector("list", nrow(S))
  for (i in seq_len(nrow(S))) {
    s <- S[i]
    v <- cohort$io$variants(s, reg_core, type = "small")$small
    x <- cohort$io$reads(s, reg_ext)
    info[[i]] <- data.table::data.table(sample_id = s$sample_id, n_reads = nrow(x$reads),
                                        n_calls = nrow(x$calls))
    v <- v[!(gt_hp1 %in% 0L & gt_hp2 %in% 0L)]
    if (nrow(v)) sites[[i]] <- v[, .(chr, s = pos, e = pos + pmax(nchar(ref), 1L) - 1L)]
    if (!nrow(x$calls)) next
    cl <- lq_apply_call_mode(x$calls, params)
    cpg[[i]] <- cl[pos >= core[1] & pos <= core[2],
                   .(n = .N, k = sum(call)), by = .(chr, pos)][, stratum := s$stratum]
    if (comethyl) {
      ps <- .pair_stats(cl[, .(chr, pos, read_id, call = round(call))], a$max_lag)
      ps <- ps[cpg_a >= core[1] & cpg_a <= core[2]]
      ps[, valid := n_pair >= a$min_pair_reads & !is.na(r2)]
      pairs[[i]] <- ps[, .(chr, cpg_a, cpg_b, n11, n10, n01, n00,
                           r2w = data.table::fifelse(valid, n_pair * r2, 0),
                           w = data.table::fifelse(valid, as.numeric(n_pair), 0),
                           ns = as.integer(valid))]
    }
  }

  cpg <- data.table::rbindlist(cpg, use.names = TRUE, fill = TRUE)
  cpg <- if (nrow(cpg)) {
    cpg[, .(n_samples = sum(n >= min_cpg_reads), n = sum(n), k = sum(k)), by = .(chr, pos, stratum)]
  } else data.table::data.table(chr = character(), pos = integer(), stratum = character(),
                                n_samples = integer(), n = integer(), k = numeric())
  pairs <- data.table::rbindlist(pairs, use.names = TRUE, fill = TRUE)
  if (nrow(pairs)) {
    pairs <- pairs[, .(n11 = sum(n11), n10 = sum(n10), n01 = sum(n01), n00 = sum(n00),
                       r2w = sum(r2w), w = sum(w), n_samples_r2 = sum(ns)),
                   by = .(chr, cpg_a, cpg_b)]
  }
  sites <- data.table::rbindlist(sites)
  if (nrow(sites)) sites <- sites[, .(n_carriers = .N), by = .(chr, s, e)]
  list(cpg = cpg, pairs = pairs, sites = sites, samples = data.table::rbindlist(info))
}

#' Reference blocks for a cohort
#'
#' Pass A of the cohort build. For every genome window and sample, reads are
#' extracted, calls are made with `params$call_mode`, and three things are
#' kept: per-CpG coverage and methylation, the 2x2 co-methylation counts and
#' r^2 of nearby CpG pairs (comethyl mode), and the variants that hit a CpG.
#' Pooled over the cohort, these give one block set for everybody:
#'
#' * CpGs covered by at least `min_cpg_reads` reads in at least
#'   `min_sample_frac` of the samples are kept.
#' * CpGs hit by a variant in at least `mask_min_carriers` samples (SNV on
#'   the C or G, indel covering either) are masked for all samples, so
#'   carriers and non-carriers are measured on the same CpGs. Rarer
#'   CpG-SNPs are masked in their carriers only ([lq_cohort_store()]). With
#'   the default of 1 every CpG-SNP seen in the cohort is masked for
#'   everybody; in large cohorts this can remove a noticeable share of
#'   CpGs, the print method reports it.
#' * comethyl: a pair is linked if the read-weighted mean of the per-sample
#'   r^2 reaches `r2_min` (per-sample r^2 avoids the spurious coupling that
#'   pooled counts get from methylation differences between people), or if
#'   it is near-constant and concordant in the pooled counts.
#'   regions: tiles or BED as in [lq_define_blocks()].
#'
#' @param cohort An [lq_cohort()].
#' @param params [lq_params()]: call mode, block mode and arguments.
#' @param regions Regions to process (see [lq_cohort_windows()]).
#' @param out_dir Directory for per-window results (resumable) and the
#'   reference (`reference.rds`, `reference_blocks.tsv`).
#' @param window Window size of pass A in bp.
#' @param min_cpg_reads,min_sample_frac CpG coverage filter (see above).
#' @param mask_min_carriers Carriers from which a CpG-SNP is masked
#'   cohort-wide (see above).
#' @param backend,n_workers See [lq_run_chunks()].
#' @return Object of class `lq_reference` (also `lq_blocks`): `blocks`
#'   (block_id, chr, start, end, n_cpg, n_cpg_masked, mean_meth,
#'   median_samples, mean_depth), `cpg_map`, `masked` (chr, pos,
#'   n_carriers, cohort_masked), `block_strata` (block_id, stratum, n_samples, meth),
#'   `samples` and `params`.
#' @export
lq_reference_blocks <- function(cohort, params, regions, out_dir, window = 5e6,
                                min_cpg_reads = 3L, min_sample_frac = 0.5, mask_min_carriers = 1L,
                                backend = c("sequential", "multicore", "future"), n_workers = 1L) {
  stopifnot(inherits(cohort, "lq_cohort"), inherits(params, "lq_params"))
  backend <- match.arg(backend)
  .check_scalar_num(min_sample_frac, "min_sample_frac", lower = 0, upper = 1)
  wins <- lq_cohort_windows(regions, window)
  halo <- if (params$block_mode == "comethyl") as.integer(params$block_args$max_bp) else 0L
  wdir <- file.path(out_dir, "passA")
  man <- lq_run_chunks(wins, .passA_worker, params, wdir, backend = backend, n_workers = n_workers,
                       cohort = cohort, halo = halo, min_cpg_reads = as.integer(min_cpg_reads))
  if (any(man$status == "failed")) {
    stop(sum(man$status == "failed"), " window(s) failed in pass A, e.g. ",
         man[status == "failed"]$error[1], ". Re-run to retry only those.")
  }
  res <- lapply(file.path(wdir, paste0(wins$chunk_id, ".rds")), readRDS)
  ref <- .build_reference(res, cohort, params, min_cpg_reads, min_sample_frac, mask_min_carriers)
  saveRDS(ref, file.path(out_dir, "reference.rds"))
  data.table::fwrite(ref$blocks, file.path(out_dir, "reference_blocks.tsv"), sep = "\t")
  ref
}

.build_reference <- function(res, cohort, params, min_cpg_reads, min_sample_frac, mask_min_carriers = 1L) {
  a <- params$block_args
  N <- nrow(cohort$samples)
  cs <- data.table::rbindlist(lapply(res, `[[`, "cpg"))
  tot <- cs[, .(n_samples = sum(n_samples), n = sum(n), k = sum(k)), by = .(chr, pos)]
  need <- max(1L, ceiling(min_sample_frac * N))
  keep <- tot[n_samples >= need]

  # population CpG-SNPs: a variant in any sample overlapping C or G
  sites <- data.table::rbindlist(lapply(res, `[[`, "sites"))
  masked <- data.table::data.table(chr = character(), pos = integer(), n_carriers = integer())
  if (nrow(sites) && nrow(keep)) {
    sites <- sites[, .(n_carriers = sum(n_carriers)), by = .(chr, s, e)]
    cp <- keep[, .(chr, c_s = pos, c_e = pos + 1L, pos)]
    data.table::setkey(cp, chr, c_s, c_e)
    data.table::setkey(sites, chr, s, e)
    ov <- data.table::foverlaps(sites, cp, by.x = c("chr", "s", "e"), by.y = c("chr", "c_s", "c_e"),
                                nomatch = NULL)
    if (nrow(ov)) masked <- ov[, .(n_carriers = max(n_carriers)), by = .(chr, pos)]
  }
  masked[, cohort_masked := n_carriers >= mask_min_carriers]
  keep <- keep[!masked[cohort_masked == TRUE], on = c("chr", "pos")]

  seg <- if (params$block_mode == "comethyl") {
    pr <- data.table::rbindlist(lapply(res, `[[`, "pairs"))
    if (nrow(pr)) {
      pr <- pr[keep[, .(chr, cpg_a = pos)], on = c("chr", "cpg_a"), nomatch = NULL]
      pr <- pr[keep[, .(chr, cpg_b = pos)], on = c("chr", "cpg_b"), nomatch = NULL]
      pr[, n_pair := n11 + n10 + n01 + n00]
      pr[, `:=`(p_a = (n11 + n10) / n_pair, p_b = (n11 + n01) / n_pair,
                r2 = data.table::fifelse(w > 0, r2w / w, NA_real_))]
      pr[, link := .link_pairs(n_pair, p_a, p_b, r2, a)]
    } else {
      pr <- data.table::data.table(chr = character(), cpg_a = integer(), cpg_b = integer(), link = integer())
    }
    .segment_cpgs(keep[, .(chr, pos)], pr, a)
  } else {
    .blocks_regions(keep[, .(chr, pos)], a)
  }
  blocks <- seg$blocks
  cm <- seg$cpg_map

  # block annotation from pooled CpG statistics
  bt <- merge(cm, tot, by = c("chr", "pos"))
  ann <- bt[, .(mean_meth = sum(k) / sum(n), median_samples = as.numeric(stats::median(n_samples)),
                mean_depth = sum(n) / (.N * N)), by = block_id]
  blocks <- merge(blocks, ann, by = "block_id", all.x = TRUE)
  blocks[, n_cpg_masked := 0L]
  if (nrow(masked) && nrow(blocks)) {
    mk <- masked[cohort_masked == TRUE, .(chr, s = pos, e = pos)]
    bb <- blocks[, .(chr, s = start, e = end, block_id)]
    data.table::setkey(bb, chr, s, e)
    ov <- data.table::foverlaps(mk, bb, nomatch = NULL)
    if (nrow(ov)) {
      nm <- ov[, .N, by = block_id]
      blocks[nm, n_cpg_masked := i.N, on = "block_id"]
    }
  }
  data.table::setorder(blocks, chr, start)
  blocks[, mode := params$block_mode]

  st <- merge(cm, cs, by = c("chr", "pos"))
  block_strata <- st[, .(n_samples = max(n_samples), meth = sum(k) / sum(n)), by = .(block_id, stratum)]

  diag <- list(mode = params$block_mode, n_samples = N, n_cpg_seen = nrow(tot),
               n_cpg_low_coverage = sum(tot$n_samples < need),
               n_cpg_masked = sum(masked$cohort_masked),
               n_cpg_masked_per_sample = sum(!masked$cohort_masked),
               mask_min_carriers = mask_min_carriers,
               n_cpg = nrow(cm), n_blocks = nrow(blocks),
               n_singletons = sum(blocks$n_cpg == 1L),
               median_cpgs = if (nrow(blocks)) stats::median(blocks$n_cpg) else NA_real_,
               median_bp = if (nrow(blocks)) stats::median(blocks$end - blocks$start + 1) else NA_real_,
               min_cpg_reads = min_cpg_reads, min_sample_frac = min_sample_frac)
  structure(list(blocks = blocks[], cpg_map = cm, masked = masked, block_strata = block_strata,
                 samples = cohort$samples, params = params),
            class = c("lq_reference", "lq_blocks"), diagnostics = diag, args = a)
}

#' @export
print.lq_reference <- function(x, ...) {
  d <- attr(x, "diagnostics")
  cat("<lq_reference> mode =", d$mode, "|", d$n_samples, "samples\n")
  cat("  CpGs seen:", d$n_cpg_seen, "| below coverage filter:", d$n_cpg_low_coverage,
      "| masked (variant in a CpG):", d$n_cpg_masked, "| in blocks:", d$n_cpg, "\n")
  if (d$n_cpg_masked_per_sample > 0L) {
    cat("  CpG-SNPs with <", d$mask_min_carriers, "carriers, masked in carriers only:",
        d$n_cpg_masked_per_sample, "\n")
  }
  cat("  blocks:", d$n_blocks, "| singletons:", d$n_singletons,
      "| median CpGs/block:", d$median_cpgs, "| median span:", d$median_bp, "bp\n")
  invisible(x)
}

# ---- SV catalogue ------------------------------------------------------------------

#' Merge structural variants across samples
#'
#' Per-sample Sniffles2 calls of the same event differ slightly in position
#' and length and carry sample-specific ids. Calls of the same type are
#' merged greedily along the chromosome when their start positions are
#' within `merge_bp` of the cluster's first call and their lengths agree
#' within `min_len_ratio`. The catalogue id is
#' `chr:pos:TYPE:len` of the cluster median, so it is reproducible.
#' (A Sniffles2 multi-sample run is the better input when available; then
#' each event already has one id and this step only renames.)
#'
#' @param cohort An [lq_cohort()].
#' @param regions Regions to read (see [lq_cohort_windows()]).
#' @param merge_bp,min_len_ratio Merging tolerances.
#' @return list(catalogue = data.table(variant_id, chr, pos, end, svtype,
#'   svlen, n_carriers), members = data.table(sample_id, orig_id,
#'   variant_id, gt_hp1, gt_hp2, phased, ps, support, rnames)).
#' @export
lq_sv_catalogue <- function(cohort, regions, merge_bp = 500L, min_len_ratio = 0.7) {
  stopifnot(inherits(cohort, "lq_cohort"))
  rg <- if (is.character(regions)) regions else sprintf("%s:%d-%d", regions$chr, regions$start, regions$end)
  calls <- list()
  for (i in seq_len(nrow(cohort$samples))) {
    s <- cohort$samples[i]
    for (r in rg) {
      sv <- cohort$io$variants(s, r, type = "sv")$sv
      if (is.null(sv) || !nrow(sv)) next
      sv <- sv[!(gt_hp1 %in% 0L & gt_hp2 %in% 0L) & !is.na(svtype) & svtype != "BND"]
      if (!nrow(sv)) next
      calls[[length(calls) + 1L]] <- sv[, .(sample_id = s$sample_id, orig_id = variant_id, chr, pos, end,
                                            svtype, len = abs(data.table::fifelse(is.na(svlen), end - pos, svlen)),
                                            gt_hp1, gt_hp2, phased, ps, support, rnames)]
    }
  }
  empty <- list(
    catalogue = data.table::data.table(variant_id = character(), chr = character(), pos = integer(),
                                       end = integer(), svtype = character(), svlen = integer(),
                                       n_carriers = integer()),
    members = data.table::data.table(sample_id = character(), orig_id = character(), variant_id = character(),
                                     gt_hp1 = integer(), gt_hp2 = integer(), phased = logical(),
                                     ps = integer(), support = integer(), rnames = list()))
  if (!length(calls)) return(empty)
  m <- unique(data.table::rbindlist(calls), by = c("sample_id", "orig_id"))
  data.table::setorder(m, chr, svtype, pos)
  m[, len := pmax(len, 1L)]
  m[, cluster := {
    id <- integer(.N); cur <- 0L; p0 <- -Inf; l0 <- NA_real_
    for (j in seq_len(.N)) {
      same <- pos[j] - p0 <= merge_bp && min(len[j], l0) / max(len[j], l0) >= min_len_ratio
      if (!isTRUE(same)) { cur <- cur + 1L; p0 <- pos[j]; l0 <- len[j] }
      id[j] <- cur
    }
    id
  }, by = .(chr, svtype)]
  cat_ <- m[, .(pos = as.integer(stats::median(pos)), end = as.integer(stats::median(end)),
                svlen = as.integer(stats::median(len)), n_carriers = data.table::uniqueN(sample_id)),
            by = .(chr, svtype, cluster)]
  cat_[, variant_id := sprintf("%s:%d:%s:%d", chr, pos, svtype, svlen)]
  # identical ids from different clusters (rare): make unique
  cat_[, variant_id := if (.N > 1L) paste0(variant_id, "_", seq_len(.N)) else variant_id, by = variant_id]
  m <- merge(m, cat_[, .(chr, svtype, cluster, variant_id)], by = c("chr", "svtype", "cluster"))
  # a sample with two calls in one cluster keeps the one with most support
  data.table::setorder(m, sample_id, variant_id, -support)
  m <- unique(m, by = c("sample_id", "variant_id"))
  list(catalogue = cat_[, .(variant_id, chr, pos, end, svtype, svlen, n_carriers)][order(chr, pos)],
       members = m[, .(sample_id, orig_id, variant_id, gt_hp1, gt_hp2, phased, ps, support, rnames)])
}

# Per-sample SV table in the lq_read_variants() schema: every catalogue SV of
# the region; the sample's own call contributes GT and RNAMES, otherwise the
# SV is genotyped from the alignment only (DEL/INS).
.sample_sv_table <- function(svcat, sid, chrom, from, to) {
  if (is.null(svcat) || !nrow(svcat$catalogue)) return(NULL)
  cc <- svcat$catalogue
  ct <- cc[cc$chr == chrom & cc$end >= from & cc$pos <= to]
  if (!nrow(ct)) return(NULL)
  mb <- svcat$members[svcat$members$sample_id == sid]
  out <- merge(ct, mb[, .(variant_id, gt_hp1, gt_hp2, phased, ps, support, rnames)],
               by = "variant_id", all.x = TRUE)
  out[is.na(phased), phased := FALSE]
  out[, rnames := lapply(rnames, function(x) if (is.null(x)) character() else x)]
  out[, `:=`(ref = "N", alt = paste0("<", svtype, ">"),
             gt = data.table::fifelse(is.na(gt_hp1), "./.", paste0(gt_hp1, "/", gt_hp2)),
             multiallelic = FALSE)]
  out[, .(variant_id, chr, pos, ref, alt, gt, gt_hp1, gt_hp2, phased, ps, multiallelic,
          svtype, end, svlen, support, rnames)]
}

# ---- pass B: read-level store ------------------------------------------------------

.store_format <- function(format) {
  if (format == "auto") format <- if (requireNamespace("arrow", quietly = TRUE)) "parquet" else "rds"
  format
}

.store_write <- function(dt, file, format) {
  if (format == "parquet") {
    arrow::write_parquet(dt, paste0(file, ".parquet"))
  } else saveRDS(dt, paste0(file, ".rds"))
}

.passB_worker <- function(chunk, params, cohort, reference, svcat, allele_flank,
                          min_cpgs_per_read, store_dir, format) {
  region <- sprintf("%s:%d-%d", chunk$chr, chunk$start, chunk$end)
  vs <- max(1L, as.integer(chunk$start - allele_flank))
  ve <- as.integer(chunk$end + allele_flank)
  var_region <- sprintf("%s:%d-%d", chunk$chr, vs, ve)
  bl <- reference$blocks[chr == chunk$chr & start >= chunk$start & end <= chunk$end]
  blk <- structure(list(blocks = bl, cpg_map = reference$cpg_map[block_id %in% bl$block_id]),
                   class = "lq_blocks")
  S <- cohort$samples
  br_all <- ra_all <- sb_all <- svg_all <- info <- vector("list", nrow(S))
  for (i in seq_len(nrow(S))) {
    s <- S[i]
    v <- cohort$io$variants(s, var_region, type = "small")$small
    sv <- .sample_sv_table(svcat, s$sample_id, chunk$chr, vs, ve)
    x <- cohort$io$reads(s, region, small = v, sv = sv, phased_max_dist = 0)
    # the sample's own CpG-destroying variants below the cohort masking
    # threshold: removed for this sample only
    own <- v[!(gt_hp1 %in% 0L & gt_hp2 %in% 0L)]
    if (nrow(own) && nrow(x$calls)) {
      own[, class := "small"]
      ms <- .cpg_snp_sites(x$calls, own)
      if (nrow(ms)) x$calls <- x$calls[!ms, on = c("chr", "pos")]
    }
    br <- if (nrow(x$calls)) lq_block_reads(x, blk, params, min_cpgs_per_read) else NULL
    if (!is.null(br) && nrow(br)) {
      br <- merge(br, x$reads[, .(read_id, read_len, mapq)], by = "read_id")
      br[, sample_id := s$sample_id]
      br_all[[i]] <- br
      sb_all[[i]] <- br[, .(n_reads = .N, k = sum(k), n = sum(n)), by = .(sample_id, block_id)]
      ra <- x$read_alleles[evidence != "phased" & read_id %in% br$read_id]
      if (nrow(ra)) ra_all[[i]] <- ra[, .(sample_id = s$sample_id, read_id, variant_id, allele, evidence)]
    }
    sva <- x$read_alleles[evidence == "direct_sv"]
    if (nrow(sva)) svg_all[[i]] <- sva[, .(sample_id = s$sample_id, read_id, variant_id, allele)]
    info[[i]] <- data.table::data.table(chunk_id = chunk$chunk_id, sample_id = s$sample_id,
                                        n_reads = nrow(x$reads),
                                        n_block_reads = if (is.null(br)) 0L else nrow(br),
                                        n_alleles = if (is.null(ra_all[[i]])) 0L else nrow(ra_all[[i]]))
  }
  d <- file.path(store_dir, chunk$chunk_id)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  .bind <- function(l, empty) { x <- data.table::rbindlist(l, use.names = TRUE, fill = TRUE); if (nrow(x)) x else empty }
  .store_write(.bind(br_all, data.table::data.table(read_id = character(), block_id = character(),
                                                    k = numeric(), n = integer(), hp = integer(), ps = integer(),
                                                    strand = character(), read_len = integer(), mapq = integer(),
                                                    sample_id = character())),
               file.path(d, "block_reads"), format)
  .store_write(.bind(ra_all, data.table::data.table(sample_id = character(), read_id = character(),
                                                    variant_id = character(), allele = integer(),
                                                    evidence = character())),
               file.path(d, "read_alleles"), format)
  .store_write(.bind(sb_all, data.table::data.table(sample_id = character(), block_id = character(),
                                                    n_reads = integer(), k = numeric(), n = integer())),
               file.path(d, "sample_blocks"), format)
  .store_write(.bind(svg_all, data.table::data.table(sample_id = character(), read_id = character(),
                                                     variant_id = character(), allele = integer())),
               file.path(d, "sv_reads"), format)
  data.table::rbindlist(info)
}

#' Build the read-level store of a cohort
#'
#' Pass B of the cohort build. Reference blocks are cut into load-balanced
#' chunks ([lq_plan_chunks()]); per chunk and sample the reads are extracted
#' once and reduced to the sufficient statistics of the read model:
#'
#' * `block_reads`: sample_id, read_id, block_id, k, n, strand, hp, ps,
#'   read_len, mapq (masked CpGs are not in the reference, so never counted)
#' * `read_alleles`: direct evidence only (SNV, indel, SV from RNAMES or the
#'   alignment) for variants within `allele_flank` of the chunk; phased
#'   alleles follow from hp/ps plus the phased genotypes and are derived in
#'   the scan, so they are not stored
#' * `sample_blocks`: n_reads, k, n per sample and block (between part,
#'   stratum check)
#' * `sv_reads`: SV allele of every read (also reads outside blocks) per
#'   catalogue SV and sample; [lq_sv_read_genotypes()] turns them into
#'   genotypes for samples without their own call
#'
#' CpGs hit by a variant are masked twice: cohort-wide in the reference
#' (see [lq_reference_blocks()]) and, for variants below that threshold,
#' in the carrier sample only.
#'
#' @param cohort An [lq_cohort()].
#' @param reference Output of [lq_reference_blocks()].
#' @param params [lq_params()] (call mode; cis_window for the chunk halo).
#' @param out_dir Store directory.
#' @param sv_catalogue Output of [lq_sv_catalogue()] or NULL (no SVs).
#' @param allele_flank bp beyond the chunk for which read alleles are stored.
#' @param target_work Chunk size, see [lq_plan_chunks()].
#' @param min_cpgs_per_read See [lq_block_reads()].
#' @param format `"rds"` (default), `"parquet"` (needs arrow; smaller and
#'   readable from Python, but not yet covered by the package tests) or
#'   `"auto"` (parquet if arrow is installed).
#' @param backend,n_workers See [lq_run_chunks()].
#' @return The manifest of [lq_run_chunks()] (invisibly). The chunk table is
#'   written to `out_dir/chunks.tsv`, per-sample counts to
#'   `out_dir/store_summary.tsv`.
#' @export
lq_cohort_store <- function(cohort, reference, params, out_dir, sv_catalogue = NULL,
                            allele_flank = 2e5, target_work = 2e5, min_cpgs_per_read = 3L,
                            format = c("rds", "parquet", "auto"),
                            backend = c("sequential", "multicore", "future"), n_workers = 1L) {
  stopifnot(inherits(cohort, "lq_cohort"), inherits(reference, "lq_reference"))
  backend <- match.arg(backend)
  format <- .store_format(match.arg(format))
  if (format == "parquet" && !requireNamespace("arrow", quietly = TRUE)) stop("format = 'parquet' needs arrow.")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  chunks <- lq_plan_chunks(reference, params, target_work = target_work)
  data.table::fwrite(chunks, file.path(out_dir, "chunks.tsv"), sep = "\t")
  saveRDS(list(format = format, samples = cohort$samples, params = params,
               allele_flank = allele_flank, reference_diag = attr(reference, "diagnostics")),
          file.path(out_dir, "store_info.rds"))
  if (!is.null(sv_catalogue)) saveRDS(sv_catalogue, file.path(out_dir, "sv_catalogue.rds"))
  man <- lq_run_chunks(chunks, .passB_worker, params, file.path(out_dir, "manifest"),
                       backend = backend, n_workers = n_workers,
                       cohort = cohort, reference = reference, svcat = sv_catalogue,
                       allele_flank = allele_flank, min_cpgs_per_read = as.integer(min_cpgs_per_read),
                       store_dir = file.path(out_dir, "chunks"), format = format)
  summ <- lq_collect(file.path(out_dir, "manifest"))
  if (nrow(summ)) data.table::fwrite(summ, file.path(out_dir, "store_summary.tsv"), sep = "\t")
  invisible(man)
}

#' Read tables from a cohort store
#'
#' @param store_dir Directory given to [lq_cohort_store()].
#' @param table One of `"block_reads"`, `"read_alleles"`, `"sample_blocks"`,
#'   `"sv_reads"`.
#' @param chunks Chunk ids to read (default all).
#' @return data.table.
#' @export
lq_store_read <- function(store_dir, table = c("block_reads", "read_alleles", "sample_blocks", "sv_reads"),
                          chunks = NULL) {
  table <- match.arg(table)
  info <- readRDS(file.path(store_dir, "store_info.rds"))
  ext <- if (info$format == "parquet") ".parquet" else ".rds"
  if (is.null(chunks)) chunks <- list.dirs(file.path(store_dir, "chunks"), full.names = FALSE, recursive = FALSE)
  files <- file.path(store_dir, "chunks", chunks, paste0(table, ext))
  files <- files[file.exists(files)]
  rd <- if (info$format == "parquet") function(f) data.table::as.data.table(arrow::read_parquet(f)) else readRDS
  out <- data.table::rbindlist(lapply(files, rd), use.names = TRUE, fill = TRUE)
  # reads crossing a chunk boundary are extracted in both chunks; their
  # allele records are identical, keep one
  if (table %in% c("read_alleles", "sv_reads") && nrow(out)) {
    out <- unique(out, by = c("sample_id", "read_id", "variant_id"))
  }
  out
}

#' SV genotypes from reads
#'
#' For catalogue SVs a sample did not call itself, the genotype is read off
#' its reads: ALT fraction below `max_alt_hom_ref` gives 0, above
#' `min_alt_hom_alt` gives 2, otherwise 1; fewer than `min_reads` informative
#' reads gives NA.
#'
#' @param store_dir Store directory.
#' @param min_reads,max_alt_hom_ref,min_alt_hom_alt Calling thresholds.
#' @return data.table sample_id, variant_id, n_ref, n_alt, dosage.
#' @export
lq_sv_read_genotypes <- function(store_dir, min_reads = 4L, max_alt_hom_ref = 0.15,
                                 min_alt_hom_alt = 0.85) {
  r <- lq_store_read(store_dir, "sv_reads")
  if (!nrow(r)) return(data.table::data.table(sample_id = character(), variant_id = character(),
                                              n_ref = integer(), n_alt = integer(), dosage = integer()))
  g <- r[, .(n_ref = sum(allele == 0L), n_alt = sum(allele == 1L)), by = .(sample_id, variant_id)]
  g[, f := n_alt / (n_ref + n_alt)]
  g[, dosage := data.table::fcase(n_ref + n_alt < min_reads, NA_integer_,
                                  f <= max_alt_hom_ref, 0L, f >= min_alt_hom_alt, 2L, default = 1L)]
  g[, f := NULL]
  g[]
}

# ---- genotypes -----------------------------------------------------------------------

#' Per-sample genotype tables for the cohort
#'
#' Non-reference small-variant calls of each sample's VCF plus its SV calls
#' under catalogue ids, one file per sample (`<out_dir>/genotypes/<id>.rds`).
#' A variant absent from a sample's VCF is homozygous reference in the
#' scan; this is only correct for single-sample VCFs from a caller that
#' reports every non-reference site (e.g. Clair3). A jointly genotyped
#' multi-sample VCF (missing = ./.) is the better input when available.
#' SVs not called in a sample are genotyped from its reads with
#' [lq_sv_read_genotypes()].
#'
#' @param cohort An [lq_cohort()].
#' @param regions Regions (see [lq_cohort_windows()]).
#' @param out_dir Store directory.
#' @param sv_catalogue Output of [lq_sv_catalogue()] or NULL.
#' @return data.table sample_id, n_small, n_sv (invisibly).
#' @export
lq_cohort_genotypes <- function(cohort, regions, out_dir, sv_catalogue = NULL) {
  stopifnot(inherits(cohort, "lq_cohort"))
  rg <- if (is.character(regions)) regions else sprintf("%s:%d-%d", regions$chr, regions$start, regions$end)
  d <- file.path(out_dir, "genotypes")
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  res <- lapply(seq_len(nrow(cohort$samples)), function(i) {
    s <- cohort$samples[i]
    sm <- data.table::rbindlist(lapply(rg, function(r) cohort$io$variants(s, r, type = "small")$small))
    sm <- if (nrow(sm)) sm[!(gt_hp1 %in% 0L & gt_hp2 %in% 0L) & !multiallelic,
                           .(variant_id, chr, pos, gt_hp1, gt_hp2, phased, ps, class = "small")] else NULL
    sv <- if (!is.null(sv_catalogue) && nrow(sv_catalogue$members)) {
      m <- sv_catalogue$members[sv_catalogue$members$sample_id == s$sample_id]
      merge(m[, .(variant_id, gt_hp1, gt_hp2, phased, ps)],
            sv_catalogue$catalogue[, .(variant_id, chr, pos)], by = "variant_id")[, class := "sv"]
    } else NULL
    g <- data.table::rbindlist(list(sm, sv), use.names = TRUE, fill = TRUE)
    if (!nrow(g)) g <- data.table::data.table(variant_id = character(), chr = character(), pos = integer(),
                                              gt_hp1 = integer(), gt_hp2 = integer(), phased = logical(),
                                              ps = integer(), class = character())
    g[, dosage := as.integer(gt_hp1 %in% 1L) + as.integer(gt_hp2 %in% 1L)]
    g[is.na(gt_hp1) | is.na(gt_hp2), dosage := NA_integer_]
    saveRDS(g, file.path(d, paste0(s$sample_id, ".rds")))
    data.table::data.table(sample_id = s$sample_id, n_small = sum(g$class == "small"), n_sv = sum(g$class == "sv"))
  })
  out <- data.table::rbindlist(res)
  data.table::fwrite(out, file.path(out_dir, "genotypes_summary.tsv"), sep = "\t")
  invisible(out)
}

# ---- stratum check -------------------------------------------------------------------

#' Do technical strata differ more than people do?
#'
#' For every pair of strata, the per-block difference of mean methylation
#' (sample-level means) is summarised as a systematic shift (mean over
#' blocks, overall and in low / intermediate / high methylation bins, where
#' call calibration differences show) and as the median absolute
#' difference. Both are compared with permutations of the stratum labels
#' among the samples of the two strata, which keeps the block structure and
#' answers whether the observed difference is larger than what the same
#' people would show in random groups of the same sizes.
#'
#' @param sample_blocks `lq_store_read(store, "sample_blocks")`.
#' @param samples data.frame sample_id, stratum.
#' @param min_reads Minimum reads for a sample's block value.
#' @param min_per_stratum Minimum samples with data per stratum and block.
#' @param max_blocks Blocks used (random subset, for speed).
#' @param n_perm Permutations.
#' @param seed RNG seed.
#' @return data.table, one row per stratum pair.
#' @export
lq_stratum_check <- function(sample_blocks, samples, min_reads = 5L, min_per_stratum = 2L,
                             max_blocks = 20000L, n_perm = 200L, seed = 1L) {
  sb <- data.table::as.data.table(sample_blocks)[n_reads >= min_reads & n > 0]
  st <- data.table::as.data.table(samples)[, .(sample_id, stratum)]
  sb <- merge(sb, st, by = "sample_id")
  if (!nrow(sb)) stop("No sample x block values with >= ", min_reads, " reads.")
  sb[, meth := k / n]
  M <- data.table::dcast(sb, block_id ~ sample_id, value.var = "meth")
  bid <- M$block_id
  M <- as.matrix(M[, -1L])
  lab <- st$stratum[match(colnames(M), st$sample_id)]
  strata <- names(which(table(lab) >= min_per_stratum))
  if (length(strata) < 2L) stop("Need at least two strata with >= ", min_per_stratum, " samples.")
  set.seed(seed)
  out <- list()
  for (ia in seq_along(strata)) for (ib in seq_along(strata)) if (ia < ib) {
    A <- strata[ia]; B <- strata[ib]
    cols <- which(lab %in% c(A, B))
    X <- M[, cols, drop = FALSE]
    g <- lab[cols] == A
    okb <- rowSums(!is.na(X[, g, drop = FALSE])) >= min_per_stratum &
           rowSums(!is.na(X[, !g, drop = FALSE])) >= min_per_stratum
    X <- X[okb, , drop = FALSE]
    if (nrow(X) > max_blocks) X <- X[sort(sample.int(nrow(X), max_blocks)), , drop = FALSE]
    if (!nrow(X)) next
    lvl <- rowMeans(X, na.rm = TRUE)
    bin <- cut(lvl, c(-Inf, 0.2, 0.8, Inf), labels = c("low", "mid", "high"))
    stat <- function(gg) {
      d <- rowMeans(X[, gg, drop = FALSE], na.rm = TRUE) - rowMeans(X[, !gg, drop = FALSE], na.rm = TRUE)
      c(shift = mean(d, na.rm = TRUE), mad = stats::median(abs(d), na.rm = TRUE),
        low = mean(d[bin == "low"], na.rm = TRUE), mid = mean(d[bin == "mid"], na.rm = TRUE),
        high = mean(d[bin == "high"], na.rm = TRUE))
    }
    obs <- stat(g)
    perm <- t(vapply(seq_len(n_perm), function(r) {
      gg <- sample(g)
      # permuted groups need data on both sides for every block
      stat(gg)
    }, numeric(5)))
    p_shift <- (1 + sum(abs(perm[, "shift"]) >= abs(obs["shift"]))) / (n_perm + 1)
    p_mid <- (1 + sum(abs(perm[, "mid"]) >= abs(obs["mid"]), na.rm = TRUE)) / (n_perm + 1)
    p_mad <- (1 + sum(perm[, "mad"] >= obs["mad"])) / (n_perm + 1)
    out[[length(out) + 1L]] <- data.table::data.table(
      stratum_a = A, stratum_b = B, n_a = sum(g), n_b = sum(!g), n_blocks = nrow(X),
      shift = obs[["shift"]], shift_low = obs[["low"]], shift_mid = obs[["mid"]], shift_high = obs[["high"]],
      p_shift = p_shift, p_shift_mid = p_mid,
      mad = obs[["mad"]], mad_perm_median = stats::median(perm[, "mad"]), p_mad = p_mad)
  }
  data.table::rbindlist(out)
}
