# plots.R
#
# Figures for longQTL results. Every function returns a ggplot (or a
# patchwork of ggplots), so it can be printed, saved with ggsave() or
# modified with further ggplot layers. ggplot2 and patchwork are suggested,
# not required: the analysis runs without them, the plot functions say what
# is missing.
#
#   single-sample ASM scan (lq_patient_scan()$units, also many windows bound
#   together):
#     lq_asm_regions()        neighbouring call blocks -> ASM regions
#     lq_plot_qq()            calibration of the haplotype test
#     lq_plot_manhattan()     -log10 p along the genome
#     lq_plot_chromosomes()   share of tested blocks with ASM per chromosome
#     lq_plot_effects()       HP1 vs HP2 methylation of every ASM region
#     lq_plot_known()         enrichment / recovery of known regions
#     lq_plot_locus_reads()   read-level view of one locus
#     lq_asm_report()         all of the above into one PDF + region table
#   cohort build:
#     lq_plot_blocks()        size of reference blocks
#     lq_plot_stratum_check() technical strata, from lq_stratum_check()

.need_plot_pkgs <- function(patchwork = FALSE) {
  need <- c("ggplot2", if (patchwork) "patchwork")
  miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) {
    stop("Plotting needs the package(s) ", paste(miss, collapse = ", "),
         ". Install with install.packages(c(", paste0('"', miss, '"', collapse = ", "), ")).", call. = FALSE)
  }
}

.lq_theme <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "bold"),
                   plot.subtitle = ggplot2::element_text(colour = "grey35"),
                   legend.position = "bottom")
}

.lq_cols <- c(call = "#2c7fb8", known = "#d95f0e", partial = "#7fbf7b",
              grey1 = "grey70", grey2 = "grey85")

.chr_order <- function(x) {
  x <- as.character(x)
  k <- sub("^chr", "", x)
  num <- suppressWarnings(as.integer(k))
  rank <- ifelse(!is.na(num), num, ifelse(k == "X", 23L, ifelse(k == "Y", 24L,
                 ifelse(k %in% c("M", "MT"), 25L, 26L))))
  factor(x, levels = unique(x[order(rank, x)]))
}

# known regions: BED path (0-based) or data.frame chr, start, end, name (1-based)
.read_known <- function(known) {
  if (is.null(known)) return(NULL)
  if (is.character(known)) {
    k <- data.table::fread(known, header = FALSE, select = 1:4, col.names = c("chr", "start", "end", "name"))
    k[, start := start + 1L]
  } else {
    k <- data.table::as.data.table(known)
    if (!"name" %in% names(k)) k[, name := sprintf("%s:%d-%d", chr, start, end)]
  }
  # methylated parental allele: a column `methylated_allele`, or the last
  # "|"-separated field of the name (as in "H19|somatic|P"); M = maternal,
  # P = paternal, anything else unknown
  al <- if ("methylated_allele" %in% names(k)) as.character(k$methylated_allele) else sub("^.*\\|", "", k$name)
  al[!al %in% c("M", "P")] <- NA_character_
  k[, .(chr = as.character(chr), start = as.integer(start), end = as.integer(end), name = as.character(name),
        methylated_allele = al)]
}

.overlap_name <- function(chr, start, end, k) {
  out <- rep(NA_character_, length(chr))
  if (is.null(k) || !nrow(k)) return(out)
  for (i in seq_len(nrow(k))) out[chr == k$chr[i] & end >= k$start[i] & start <= k$end[i]] <- k$name[i]
  out
}

# one row per block: its best haplotype test (variant test only where no
# haplotype test exists)
.best_per_block <- function(u) {
  u <- data.table::as.data.table(u)
  u[, pref_ := data.table::fifelse(unit_type == "haplotype", 0L, 1L)]
  b <- u[order(pref_, p_final)][, .SD[1L], by = block_id]
  b[, pref_ := NULL]
  b[]
}

# Points that land on the same spot of the figure are drawn once. `res` is
# the grid in axis units; with a few thousand grid cells per axis the figure
# is visually identical but a PDF of millions of tests stays small.
.thin_xy <- function(d, x, y, xres, yres) {
  key <- paste(round(d[[x]] / xres), round(d[[y]] / yres))
  d[!duplicated(key)]
}

.fmt_p <- function(p) if (p < 1e-300) "< 1e-300" else sprintf("%.1e", p)

.check_units <- function(u) {
  need <- c("block_id", "chr", "start", "end", "unit_type", "unit_id", "p_final", "q", "asm_call",
            "meth_0", "meth_1", "delta", "n_cpg")
  miss <- setdiff(need, names(u))
  if (length(miss)) stop("u is missing columns: ", paste(miss, collapse = ", "),
                         ". Pass the `units` table of lq_patient_scan().", call. = FALSE)
  data.table::as.data.table(u)
}

#' Merge neighbouring ASM calls into regions
#'
#' One biological ASM region (e.g. a 6 kb imprinting centre) is called in
#' many neighbouring blocks. Called haplotype tests in the same phase set
#' and with the same direction are merged when at most `max_gap` bp apart.
#' Regions, not blocks, are the unit to count and report.
#'
#' @param u `units` of [lq_patient_scan()] (several windows may be bound
#'   together; recompute `q` and `asm_call` over all of them first).
#' @param max_gap Maximum distance in bp between merged blocks.
#' @param known Optional known regions: BED file (0-based) or data.frame
#'   chr, start, end, name (1-based).
#' @return data.table, one row per region, sorted by `min_q`: chr, start,
#'   end, width, unit_id, n_blocks, n_cpg, meth_hp1, meth_hp2, delta, min_q,
#'   reads_hp1, reads_hp2, phasing_free_concordance, cpg_snp, context_snp,
#'   candidates, pattern (`"monoallelic"` when one haplotype is below 0.2
#'   and the other above 0.8, else `"partial"`), known_region.
#' @export
lq_asm_regions <- function(u, max_gap = 1000, known = NULL) {
  u <- .check_units(u)
  h <- u[asm_call == TRUE & unit_type == "haplotype"]
  if (!nrow(h)) return(data.table::data.table())
  h[, dir_ := sign(delta)]
  data.table::setorder(h, chr, unit_id, dir_, start)
  h[, region_ := cumsum(c(TRUE, start[-1L] - end[-.N] > max_gap)), by = .(chr, unit_id, dir_)]
  has <- function(col) if (col %in% names(h)) h[[col]] else rep(NA, nrow(h))
  h[, `:=`(conc_ = has("bimodal_hp_concordance"), cs_ = has("cpg_snp"), cx_ = has("context_snp"),
           cand_ = has("candidates"))]
  r <- h[, .(start = min(start), end = max(end), n_blocks = .N, n_cpg = sum(n_cpg),
             meth_hp1 = stats::weighted.mean(meth_0, n_cpg), meth_hp2 = stats::weighted.mean(meth_1, n_cpg),
             min_q = min(q), reads_hp1 = max(n_0), reads_hp2 = max(n_1),
             phasing_free_concordance = if (all(is.na(conc_))) NA_real_ else max(conc_, na.rm = TRUE),
             cpg_snp = any(cs_ %in% TRUE), context_snp = any(cx_ %in% TRUE),
             candidates = cand_[which.min(q)]),
         by = .(chr, unit_id, dir_, region_)]
  r[, `:=`(delta = meth_hp2 - meth_hp1, width = end - start + 1L)]
  r[, pattern := data.table::fifelse(pmin(meth_hp1, meth_hp2) < 0.2 & pmax(meth_hp1, meth_hp2) > 0.8,
                                     "monoallelic", "partial")]
  r[, known_region := .overlap_name(chr, start, end, .read_known(known))]
  r[, c("dir_", "region_") := NULL]
  data.table::setcolorder(r, c("chr", "start", "end", "width", "unit_id", "n_blocks", "n_cpg",
                               "meth_hp1", "meth_hp2", "delta", "min_q"))
  data.table::setorder(r, min_q)
  r[]
}

#' QQ plot of the haplotype test
#'
#' All haplotype tests, not only the calls. Points on the diagonal (inside
#' the 95 % band) mean no signal; a tail above it is real ASM. A bulk above
#' the diagonal (lambda well above 1) would point to a miscalibrated test.
#'
#' @inheritParams lq_asm_regions
#' @param cap Upper limit of the -log10 p axis.
#' @return A ggplot. Points that would overlap exactly in the figure are
#'   drawn once, so the plot stays small with millions of tests.
#' @export
lq_plot_qq <- function(u, cap = 50) {
  .need_plot_pkgs()
  h <- .check_units(u)[unit_type == "haplotype" & !is.na(p_final)]
  n <- nrow(h)
  if (!n) stop("No haplotype tests in u.")
  d <- data.table::data.table(obs = sort(-log10(pmax(h$p_final, 1e-300))), exp = -log10(rev(stats::ppoints(n))))
  i <- seq_len(n)
  band <- data.table::data.table(exp = -log10(stats::ppoints(n)),
                                 lo = -log10(stats::qbeta(0.975, i, n - i + 1)),
                                 hi = -log10(stats::qbeta(0.025, i, n - i + 1)))
  lam <- stats::median(stats::qchisq(h$p_final, 1, lower.tail = FALSE)) / stats::qchisq(0.5, 1)
  d[, obs := pmin(obs, cap)]
  rng <- max(d$exp, d$obs)
  d <- .thin_xy(d, "exp", "obs", rng / 2000, rng / 2000)
  band <- .thin_xy(band, "exp", "lo", rng / 2000, rng / 2000)
  ggplot2::ggplot(d, ggplot2::aes(.data$exp, .data$obs)) +
    ggplot2::geom_ribbon(data = band, ggplot2::aes(x = .data$exp, ymin = .data$lo, ymax = pmin(.data$hi, cap)),
                         inherit.aes = FALSE, fill = "grey88") +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey40", linetype = 2) +
    ggplot2::geom_point(size = 0.6, colour = .lq_cols[["call"]]) +
    ggplot2::labs(x = expression(expected ~ -log[10](p)), y = expression(observed ~ -log[10](p)),
                  title = "Calibration of the haplotype test",
                  subtitle = sprintf("%s tests, lambda_GC = %.2f", format(n, big.mark = ","), lam)) +
    .lq_theme()
}

#' Genome-wide view of an ASM scan
#'
#' Best test per block, -log10 p along the chromosomes. Calls in colour,
#' calls in known regions highlighted, the strongest regions labelled
#' (labels need ggrepel).
#'
#' @inheritParams lq_asm_regions
#' @param label_top Number of regions to label.
#' @param cap Upper limit of the -log10 p axis.
#' @return A ggplot.
#' @export
lq_plot_manhattan <- function(u, known = NULL, label_top = 10, cap = 50) {
  .need_plot_pkgs()
  b <- .best_per_block(.check_units(u))
  k <- .read_known(known)
  b[, known_ := .overlap_name(chr, start, end, k)]
  if ("known_region" %in% names(b)) b[is.na(known_) & !is.na(known_region) & nzchar(known_region), known_ := known_region]
  b[, chr := .chr_order(chr)]
  data.table::setorder(b, chr, start)
  off <- b[, .(len = as.numeric(max(end))), by = chr][, off := data.table::shift(cumsum(len), fill = 0)]
  b <- merge(b, off[, .(chr, off)], by = "chr")
  b[, `:=`(x = start + off, y = pmin(-log10(pmax(p_final, 1e-300)), cap))]
  b[, grp := data.table::fifelse(asm_call, data.table::fifelse(!is.na(known_), "call, known region", "call"),
                                 data.table::fifelse(as.integer(chr) %% 2L == 0L, "even", "odd"))]
  cols <- c(odd = .lq_cols[["grey1"]], even = .lq_cols[["grey2"]], call = .lq_cols[["call"]],
            `call, known region` = .lq_cols[["known"]])
  mids <- b[, .(mid = (min(x) + max(x)) / 2), by = chr]
  n_tested <- nrow(b); n_called <- sum(b$asm_call)
  # non-calls: one point per spot of the figure (calls are always all drawn)
  bg <- .thin_xy(b[asm_call == FALSE], "x", "y", max(b$x) / 1500, cap / 200)
  p <- ggplot2::ggplot(b, ggplot2::aes(.data$x, .data$y, colour = .data$grp)) +
    ggplot2::geom_point(data = bg, size = 0.4) +
    ggplot2::geom_point(data = b[asm_call == TRUE], size = 0.8) +
    ggplot2::scale_colour_manual(values = cols, breaks = c("call", "call, known region"), name = NULL) +
    ggplot2::scale_x_continuous(breaks = mids$mid, labels = sub("^chr", "", mids$chr), expand = c(0.01, 0)) +
    ggplot2::labs(x = NULL, y = expression(-log[10](p)), title = "Allele-specific methylation along the genome",
                  subtitle = sprintf("%s blocks tested, %s with ASM (q < 0.05, |delta| >= 0.1); y capped at %s",
                                     format(n_tested, big.mark = ","), format(n_called, big.mark = ","), cap)) +
    .lq_theme() + ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
  if (any(b$asm_call)) {
    p <- p + ggplot2::geom_hline(yintercept = -log10(max(b[asm_call == TRUE, p_final])),
                                 linetype = 2, colour = "grey40", linewidth = 0.3)
  }
  if (label_top > 0 && requireNamespace("ggrepel", quietly = TRUE)) {
    reg <- lq_asm_regions(u, known = known)
    if (nrow(reg)) {
      top <- utils::head(reg, label_top)
      top[, lab := data.table::fifelse(!is.na(known_region), known_region, sprintf("%s:%.2f Mb", chr, start / 1e6))]
      top <- merge(top[, .(chr = as.character(chr), start, lab, min_q)],
                   off[, .(chr = as.character(chr), off)], by = "chr")
      top[, `:=`(x = start + off, y = pmin(-log10(pmax(min_q, 1e-300)), cap))]
      p <- p + ggrepel::geom_text_repel(data = top, ggplot2::aes(.data$x, .data$y, label = .data$lab),
                                        inherit.aes = FALSE, size = 3, min.segment.length = 0, max.overlaps = 30)
    }
  }
  p
}

#' Share of tested blocks with ASM per chromosome
#'
#' In a female sample X inactivation makes chrX stand out; this is expected
#' and a useful second positive control, but chrX should be reported
#' separately from autosomal ASM.
#'
#' @inheritParams lq_asm_regions
#' @return A ggplot.
#' @export
lq_plot_chromosomes <- function(u) {
  .need_plot_pkgs()
  b <- .best_per_block(.check_units(u))
  s <- b[, .(tested = .N, called = sum(asm_call)), by = chr][, frac := called / tested]
  s[, chr := .chr_order(chr)]
  s[, sex := grepl("X|Y", chr)]
  med <- stats::median(s[sex == FALSE]$frac)
  ggplot2::ggplot(s, ggplot2::aes(.data$chr, 100 * .data$frac, fill = .data$sex)) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = c(`FALSE` = .lq_cols[["call"]], `TRUE` = .lq_cols[["known"]])) +
    ggplot2::geom_hline(yintercept = 100 * med, linetype = 2, colour = "grey40") +
    ggplot2::scale_x_discrete(labels = function(x) sub("^chr", "", x)) +
    ggplot2::labs(x = NULL, y = "% of tested blocks with ASM", title = "ASM per chromosome",
                  subtitle = "dashed: autosomal median") +
    .lq_theme()
}

#' Effect landscape of ASM regions
#'
#' Left: HP1 against HP2 methylation of every region. Corners (one allele
#' unmethylated, the other methylated) look like imprinting or X
#' inactivation; points between the corners are partial ASM, the usual
#' picture of sequence-dependent effects. Right: effect sizes.
#'
#' @param regions Output of [lq_asm_regions()].
#' @return A patchwork of two ggplots.
#' @export
lq_plot_effects <- function(regions) {
  .need_plot_pkgs(patchwork = TRUE)
  r <- data.table::as.data.table(regions)
  if (!nrow(r)) stop("No ASM regions.")
  r[, cls := data.table::fifelse(!is.na(known_region), "known region", pattern)]
  cols <- c(monoallelic = .lq_cols[["call"]], partial = .lq_cols[["partial"]], `known region` = .lq_cols[["known"]])
  p1 <- ggplot2::ggplot(r, ggplot2::aes(.data$meth_hp1, .data$meth_hp2, colour = .data$cls, size = .data$n_cpg)) +
    ggplot2::annotate("rect", xmin = 0, xmax = 0.2, ymin = 0.8, ymax = 1, alpha = 0.08) +
    ggplot2::annotate("rect", xmin = 0.8, xmax = 1, ymin = 0, ymax = 0.2, alpha = 0.08) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey60", linetype = 2) +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::scale_colour_manual(values = cols, name = NULL) +
    ggplot2::scale_size_area(max_size = 4, guide = "none") +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(x = "methylation HP1", y = "methylation HP2", title = "Effect landscape",
                  subtitle = "corners: one allele methylated, one not; between: partial ASM") +
    .lq_theme()
  p2 <- ggplot2::ggplot(r, ggplot2::aes(abs(.data$delta), fill = .data$cls)) +
    ggplot2::geom_histogram(binwidth = 0.05, boundary = 0, colour = "white") +
    ggplot2::scale_fill_manual(values = cols, name = NULL) +
    ggplot2::labs(x = "|methylation difference| between haplotypes", y = "ASM regions", title = "Effect sizes") +
    .lq_theme()
  patchwork::wrap_plots(p1, p2, nrow = 1)
}

#' Positive control: known regions
#'
#' Left: share of tested blocks with ASM inside the known regions against
#' the rest of the genome, with Fisher odds ratio. Right: which of the known
#' regions that were tested contain at least one ASM block.
#'
#' @inheritParams lq_asm_regions
#' @param known Known regions (required): BED file (0-based) or data.frame
#'   chr, start, end, name (1-based).
#' @param max_regions Up to this many tested regions are shown one by one;
#'   with more, recovery is summarised in one bar.
#' @return A patchwork of two ggplots.
#' @export
lq_plot_known <- function(u, known, max_regions = 40L) {
  .need_plot_pkgs(patchwork = TRUE)
  k <- .read_known(known)
  if (is.null(k) || !nrow(k)) stop("known regions are required.")
  b <- .best_per_block(.check_units(u))
  b[, known_ := .overlap_name(chr, start, end, k)]
  b[, inside := !is.na(known_)]
  tab <- b[, .(tested = .N, called = sum(asm_call)), by = inside]
  if (nrow(tab) < 2L) stop("Need tested blocks both inside and outside the known regions.")
  data.table::setorder(tab, -inside)
  tab[, `:=`(frac = called / tested, lab = data.table::fifelse(inside, "known regions", "rest of genome"))]
  ft <- stats::fisher.test(matrix(c(tab$called, tab$tested - tab$called), 2))
  rec <- b[inside == TRUE, .(hit = any(asm_call)), by = known_]
  n_untested <- nrow(k) - nrow(rec)
  p1 <- ggplot2::ggplot(tab, ggplot2::aes(.data$lab, 100 * .data$frac, fill = .data$inside)) +
    ggplot2::geom_col(width = 0.6, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = c(`TRUE` = .lq_cols[["known"]], `FALSE` = .lq_cols[["call"]])) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%d / %d", .data$called, .data$tested)), vjust = -0.4, size = 3.2) +
    ggplot2::labs(x = NULL, y = "% of tested blocks with ASM", title = "Positive control",
                  subtitle = sprintf("odds ratio %.1f (95%% CI %.1f-%.1f), Fisher p %s",
                                     ft$estimate, ft$conf.int[1], ft$conf.int[2],
                                     sub("^([0-9])", "= \\1", .fmt_p(ft$p.value)))) +
    .lq_theme()
  if (nrow(rec) <= max_regions) {
    # few regions: one tile per region
    rec[, known_ := factor(known_, levels = known_[order(hit, known_)])]
    p2 <- ggplot2::ggplot(rec, ggplot2::aes(.data$known_, 1, fill = .data$hit)) +
      ggplot2::geom_tile(colour = "white") +
      ggplot2::scale_fill_manual(values = c(`TRUE` = .lq_cols[["known"]], `FALSE` = "grey85"),
                                 labels = c(`TRUE` = "ASM found", `FALSE` = "not found"), name = NULL) +
      ggplot2::coord_flip() +
      ggplot2::labs(x = NULL, y = NULL, title = sprintf("Recovered: %d of %d tested regions", sum(rec$hit), nrow(rec)),
                    subtitle = if (n_untested > 0L) sprintf("%d further regions had no testable block", n_untested)) +
      .lq_theme() + ggplot2::theme(axis.text.x = ggplot2::element_blank(), panel.grid = ggplot2::element_blank())
  } else {
    # many regions: one bar, found / not found / not testable
    st <- data.table::data.table(
      status = factor(c("ASM found", "tested, not found", "no testable block"),
                      levels = c("no testable block", "tested, not found", "ASM found")),
      n = c(sum(rec$hit), sum(!rec$hit), n_untested))
    p2 <- ggplot2::ggplot(st, ggplot2::aes(x = "known regions", y = .data$n, fill = .data$status)) +
      ggplot2::geom_col(width = 0.5) +
      ggplot2::geom_text(ggplot2::aes(label = .data$n), position = ggplot2::position_stack(vjust = 0.5), size = 3.2) +
      ggplot2::scale_fill_manual(values = c(`ASM found` = .lq_cols[["known"]], `tested, not found` = "grey80",
                                            `no testable block` = "grey93"), name = NULL) +
      ggplot2::labs(x = NULL, y = "regions", title = sprintf("Recovered: %d of %d tested regions", sum(rec$hit), nrow(rec)),
                    subtitle = sprintf("per-region view: use a list with <= %d regions", max_regions)) +
      .lq_theme()
  }
  patchwork::wrap_plots(p1, p2, nrow = 1)
}

#' Read-level view of one locus
#'
#' Every read a row, every CpG a dot (dark = methylated), grouped by
#' haplotype. Heterozygous variants as orange lines; with `u`, blocks with
#' an ASM call are shaded. This is the picture only long reads give: the
#' methylation state of whole molecules, per allele.
#'
#' @param bam BAM or CRAM.
#' @param region "chr:start-end".
#' @param reference FASTA for CRAM.
#' @param small_variants,sv VCF paths or tables, see [lq_extract_reads()].
#' @param u Optional `units` table; its calls in the region are shaded.
#' @param params [lq_params()] (call mode).
#' @param max_reads Reads shown per haplotype (random subset).
#' @param x Optional `lq_reads` object instead of reading `bam` again.
#' @param title Plot title.
#' @return A ggplot.
#' @export
lq_plot_locus_reads <- function(bam = NULL, region, reference = NULL, small_variants = NULL, sv = NULL,
                                u = NULL, params = lq_params(), max_reads = 60L, x = NULL, title = region) {
  .need_plot_pkgs()
  if (is.null(x)) x <- lq_extract_reads(bam, region, small_variants = small_variants, sv = sv, reference = reference)
  rg <- .parse_region(region)
  cl <- lq_apply_call_mode(x$calls[chr == rg$chr & pos >= rg$start & pos <= rg$end], params)
  if (!nrow(cl)) stop("No CpG calls in ", region)
  cl <- merge(cl[, .(read_id, pos, call)], x$reads[, .(read_id, hp, start, end)], by = "read_id")
  cl[, hp := data.table::fifelse(is.na(hp), "untagged", paste0("HP", hp))]
  rd <- unique(cl[, .(read_id, hp, start, end)])
  set.seed(1L)
  rd <- rd[, .SD[sample(.N, min(.N, max_reads))], by = hp]
  data.table::setorder(rd, hp, start)
  rd[, row := seq_len(.N), by = hp]
  cl <- merge(cl, rd[, .(read_id, row)], by = "read_id")
  hm <- cl[, .(m = mean(call)), by = hp]
  rd[, hp_lab := sprintf("%s (%d reads, mean %d%%)", hp, .N, round(100 * hm$m[match(hp, hm$hp)])), by = hp]
  cl <- merge(cl, unique(rd[, .(hp, hp_lab)]), by = "hp")
  cl[, state := factor(data.table::fifelse(call >= 0.5, "methylated", "unmethylated"),
                       levels = c("unmethylated", "methylated"))]
  p <- ggplot2::ggplot()
  if (!is.null(u)) {
    cb <- data.table::as.data.table(u)[asm_call == TRUE & chr == rg$chr & end >= rg$start & start <= rg$end,
                                       .(s = max(min(start), rg$start), e = min(max(end), rg$end)), by = block_id]
    if (nrow(cb)) p <- p + ggplot2::geom_rect(data = cb, ggplot2::aes(xmin = .data$s, xmax = .data$e, ymin = -Inf, ymax = Inf),
                                              fill = "#fdae6b", alpha = 0.25)
  }
  het <- x$variants[!is.na(gt_hp1) & !is.na(gt_hp2) & gt_hp1 != gt_hp2 & chr == rg$chr &
                      pos >= rg$start & pos <= rg$end]
  if (nrow(het)) p <- p + ggplot2::geom_vline(data = het, ggplot2::aes(xintercept = .data$pos),
                                              colour = .lq_cols[["known"]], linewidth = 0.25, alpha = 0.6)
  p +
    ggplot2::geom_segment(data = rd, ggplot2::aes(x = pmax(.data$start, rg$start), xend = pmin(.data$end, rg$end),
                                                  y = .data$row, yend = .data$row), colour = "grey80", linewidth = 0.3) +
    ggplot2::geom_point(data = cl, ggplot2::aes(.data$pos, .data$row, colour = .data$state), size = 0.7, shape = 16) +
    ggplot2::scale_colour_manual(values = c(unmethylated = "#c6dbef", methylated = "#08306b"), name = NULL, drop = FALSE) +
    ggplot2::facet_grid(hp_lab ~ ., scales = "free_y", space = "free_y") +
    ggplot2::scale_x_continuous(labels = function(v) sprintf("%.3f", v / 1e6)) +
    ggplot2::coord_cartesian(xlim = c(rg$start, rg$end)) +
    ggplot2::labs(x = sprintf("%s (Mb)", rg$chr), y = NULL, title = title,
                  subtitle = "each row one molecule; shaded: blocks with ASM; orange lines: heterozygous variants") +
    .lq_theme() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_blank(),
                   strip.text.y = ggplot2::element_text(angle = 0, hjust = 0))
}

#' All ASM figures in one PDF
#'
#' Writes the QQ plot, the genome-wide view, ASM per chromosome, the effect
#' landscape and (with `known`) the positive control to `pdf`, the region
#' table to `<pdf without .pdf>_regions.tsv`, and with `bam` read-level
#' views of the `n_loci` strongest regions to `<...>_loci.pdf`. If the known
#' regions carry their methylated parental allele (BED names ending in `|M`
#' / `|P`), the parental origin of haplotypes is assigned with
#' [lq_parent_of_origin()]: one more page plus `<...>_phase_sets.tsv` and
#' `<...>_parent_of_origin.tsv`.
#'
#' @inheritParams lq_asm_regions
#' @param pdf Output PDF path.
#' @param bam,reference,small_variants Optional; with `bam` the strongest
#'   regions are drawn read by read.
#' @param n_loci Number of read-level locus views.
#' @param flank bp added on both sides of a region in the locus views.
#' @param width,height PDF page size in inches.
#' @return The region table (invisibly).
#' @export
lq_asm_report <- function(u, pdf, known = NULL, bam = NULL, reference = NULL, small_variants = NULL,
                          n_loci = 6L, flank = 2000L, width = 11, height = 6.5) {
  .need_plot_pkgs(patchwork = TRUE)
  u <- .check_units(u)
  reg <- lq_asm_regions(u, known = known)
  base <- sub("\\.pdf$", "", pdf)
  data.table::fwrite(reg, paste0(base, "_regions.tsv"), sep = "\t")
  grDevices::pdf(pdf, width = width, height = height)
  tryCatch({
    print(lq_plot_qq(u))
    print(lq_plot_manhattan(u, known = known))
    print(lq_plot_chromosomes(u))
    if (nrow(reg)) print(lq_plot_effects(reg))
    if (!is.null(known)) print(tryCatch(lq_plot_known(u, known), error = function(e) {
      message("positive control skipped: ", conditionMessage(e)); NULL
    }))
    # parent of origin, when the known regions carry their methylated allele
    kk <- .read_known(known)
    if (!is.null(kk) && any(!is.na(kk$methylated_allele))) {
      poo <- tryCatch(lq_parent_of_origin(u, kk), error = function(e) {
        message("parent of origin skipped: ", conditionMessage(e)); NULL })
      if (!is.null(poo)) {
        data.table::fwrite(poo$phase_sets, paste0(base, "_phase_sets.tsv"), sep = "\t")
        data.table::fwrite(poo$regions[!is.na(methylated_parent)], paste0(base, "_parent_of_origin.tsv"), sep = "\t")
        if (any(!is.na(poo$regions$methylated_parent))) print(lq_plot_parent_of_origin(poo))
      }
    }
  }, finally = grDevices::dev.off())
  if (!is.null(bam) && nrow(reg) && n_loci > 0L) {
    grDevices::pdf(paste0(base, "_loci.pdf"), width = width, height = height)
    tryCatch(for (i in seq_len(min(n_loci, nrow(reg)))) {
      r <- reg[i]
      region <- sprintf("%s:%d-%d", r$chr, max(1L, r$start - flank), r$end + flank)
      ttl <- sprintf("%s  |  HP1 %.0f%% vs HP2 %.0f%%, q = %.1e%s", region, 100 * r$meth_hp1, 100 * r$meth_hp2,
                     r$min_q, if (!is.na(r$known_region)) paste0("  |  ", r$known_region) else "")
      print(lq_plot_locus_reads(bam, region, reference = reference, small_variants = small_variants,
                                u = u, title = ttl))
    }, finally = grDevices::dev.off())
  }
  invisible(reg)
}

#' Size of reference blocks
#'
#' @param reference Output of [lq_reference_blocks()] (or any `lq_blocks`).
#' @return A patchwork of two ggplots: CpGs per block and block span.
#' @export
lq_plot_blocks <- function(reference) {
  .need_plot_pkgs(patchwork = TRUE)
  b <- data.table::as.data.table(reference$blocks)
  b[, span := end - start + 1L]
  d <- attr(reference, "diagnostics")
  sub <- if (!is.null(d$n_cpg_masked)) sprintf("%s blocks, %s CpGs; %s CpGs masked (CpG-SNPs)",
                                               format(nrow(b), big.mark = ","), format(sum(b$n_cpg), big.mark = ","),
                                               format(d$n_cpg_masked, big.mark = ",")) else NULL
  p1 <- ggplot2::ggplot(b, ggplot2::aes(.data$n_cpg)) +
    ggplot2::geom_histogram(binwidth = 1, fill = .lq_cols[["call"]], colour = "white") +
    ggplot2::labs(x = "CpGs per block", y = "blocks", title = "Reference blocks", subtitle = sub) + .lq_theme()
  p2 <- ggplot2::ggplot(b[n_cpg > 1L], ggplot2::aes(.data$span)) +
    ggplot2::geom_histogram(bins = 50, fill = .lq_cols[["call"]], colour = "white") +
    ggplot2::scale_x_log10() +
    ggplot2::labs(x = "block span (bp, blocks with > 1 CpG)", y = "blocks", title = " ") + .lq_theme()
  patchwork::wrap_plots(p1, p2, nrow = 1)
}

#' Technical strata compared
#'
#' One bar per pair of strata and methylation level: mean methylation
#' difference (stratum a - stratum b) over blocks. Labels give the
#' permutation p-value of the overall and the intermediate-level shift.
#'
#' @param check Output of [lq_stratum_check()].
#' @return A ggplot.
#' @export
lq_plot_stratum_check <- function(check) {
  .need_plot_pkgs()
  s <- data.table::as.data.table(check)
  if (!nrow(s)) stop("Empty stratum check.")
  s[, pair := sprintf("%s - %s\n(n = %d vs %d)", stratum_a, stratum_b, n_a, n_b)]
  long <- data.table::melt(s, id.vars = c("pair", "p_shift", "p_shift_mid"),
                           measure.vars = c("shift_low", "shift_mid", "shift_high", "shift"),
                           variable.name = "level", value.name = "diff")
  long[, level := factor(level, levels = c("shift_low", "shift_mid", "shift_high", "shift"),
                         labels = c("low (< 0.2)", "intermediate", "high (> 0.8)", "all blocks"))]
  lab <- s[, .(pair, txt = sprintf("p = %.2g (all), %.2g (intermediate)", p_shift, p_shift_mid))]
  ggplot2::ggplot(long, ggplot2::aes(.data$level, 100 * .data$diff, fill = .data$level)) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40") +
    ggplot2::facet_wrap(~ pair) +
    ggplot2::geom_text(data = lab, ggplot2::aes(x = 2.5, y = Inf, label = .data$txt), inherit.aes = FALSE,
                       vjust = 1.5, size = 3, colour = "grey30") +
    ggplot2::scale_fill_manual(values = c("#9ecae1", "#3182bd", "#08519c", "grey50")) +
    ggplot2::labs(x = "methylation level of the block", y = "methylation difference a - b (percentage points)",
                  title = "Technical strata",
                  subtitle = "p from permuting people between strata; a calibration difference shows at intermediate levels") +
    .lq_theme()
}
