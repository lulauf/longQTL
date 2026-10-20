# parent_of_origin.R
#
# Parent of origin of haplotypes in one person, without parental data.
#
# Known imprinted DMRs are methylated on a fixed parental allele (e.g. H19:
# paternal, PEG3: maternal). Where such a DMR shows ASM, the methylated
# haplotype therefore reveals which haplotype (HP1 or HP2) of that phase set
# is maternal. HP labels are consistent within a phase set, so the assignment
# carries over to every other ASM region in the same phase set
# (Akbari et al. 2022, eLife 11:e77898, used the same idea with nanopore
# data).
#
# What it gives: for ASM regions in an anchored phase set, whether the
# maternal or the paternal allele is methylated. A new region that is
# methylated on the same parental allele in every person is an imprinting
# candidate; one whose parental side differs between people (as chr11:1.921 Mb
# between b241 and b262) behaves like a sequence-dependent effect (mQTL).
#
# Limits: only phase sets that contain an anchor are assigned; a switch error
# inside a phase set mislabels everything beyond it, which shows up as
# disagreeing anchors (status "conflict").

#' Assign the parental origin of haplotypes from known imprinted DMRs
#'
#' @param u `units` of [lq_patient_scan()] (several windows bound together
#'   are fine; recompute `q` and `asm_call` over all of them first).
#' @param known Known imprinted DMRs with their methylated parental allele:
#'   BED file whose name ends in `|M` or `|P` (as in
#'   `imprinted_dmrs_hg38_somatic.bed`), or a data.frame chr, start, end
#'   (1-based) and `methylated_allele` (`"M"`/`"P"`).
#' @param min_anchor_delta Minimum |methylation difference| for a known DMR to
#'   serve as anchor (weak, partially imprinted DMRs are not used).
#' @param min_agreement Share of anchor weight that has to agree for a phase
#'   set to be assigned; below it the phase set is a `"conflict"`.
#' @param max_gap Passed to [lq_asm_regions()].
#' @return A list:
#'   * `phase_sets`: chr, unit_id, start, end (span of its tests), n_anchors,
#'     n_agree, agreement, maternal_hp (1/2), status (`"assigned"`,
#'     `"conflict"`), anchors (names);
#'   * `regions`: [lq_asm_regions()] plus maternal_hp, is_anchor,
#'     methylated_parent (`"maternal"`/`"paternal"`), meth_maternal,
#'     meth_paternal; regions in phase sets without assignment have NA.
#' @export
lq_parent_of_origin <- function(u, known, min_anchor_delta = 0.3, min_agreement = 0.8, max_gap = 1000) {
  u <- .check_units(u)
  k <- .read_known(known)
  if (is.null(k) || !any(!is.na(k$methylated_allele))) {
    stop("known needs the methylated parental allele (name ending in |M or |P, or a methylated_allele column).")
  }
  reg <- lq_asm_regions(u, max_gap = max_gap)
  empty_ps <- data.table::data.table(chr = character(), unit_id = character(), start = integer(), end = integer(),
                                     n_anchors = integer(), n_agree = integer(), agreement = numeric(),
                                     maternal_hp = integer(), status = character(), anchors = character())
  if (!nrow(reg)) return(list(phase_sets = empty_ps, regions = reg))

  # anchors: ASM regions overlapping a known DMR with known allele
  reg[, `:=`(anchor_name = NA_character_, anchor_allele = NA_character_)]
  ka <- k[!is.na(methylated_allele)]
  for (i in seq_len(nrow(ka))) {
    reg[chr == ka$chr[i] & end >= ka$start[i] & start <= ka$end[i] & is.na(anchor_name),
        `:=`(anchor_name = ka$name[i], anchor_allele = ka$methylated_allele[i])]
  }
  reg[, is_anchor := !is.na(anchor_allele) & abs(delta) >= min_anchor_delta]
  reg[, meth_hp_ := data.table::fifelse(meth_hp2 > meth_hp1, 2L, 1L)]          # methylated haplotype
  # maternal haplotype implied by an anchor
  reg[is_anchor == TRUE, implied_ := data.table::fifelse(anchor_allele == "M", meth_hp_, 3L - meth_hp_)]
  reg[, w_ := -log10(pmax(min_q, 1e-300))]

  an <- reg[is_anchor == TRUE]
  ps <- if (nrow(an)) {
    an[, {
      w1 <- sum(w_[implied_ == 1L]); w2 <- sum(w_[implied_ == 2L])
      mhp <- if (w1 >= w2) 1L else 2L
      .(n_anchors = .N, n_agree = sum(implied_ == mhp), agreement = max(w1, w2) / (w1 + w2),
        maternal_hp = mhp, anchors = paste(unique(sub("\\|.*", "", anchor_name)), collapse = ","))
    }, by = .(chr, unit_id)]
  } else empty_ps[, .(chr, unit_id, n_anchors, n_agree, agreement, maternal_hp, anchors)]
  ps[, status := data.table::fifelse(agreement >= min_agreement, "assigned", "conflict")]
  span <- u[unit_type == "haplotype", .(start = min(start), end = max(end)), by = .(chr, unit_id)]
  ps <- merge(ps, span, by = c("chr", "unit_id"), all.x = TRUE)
  data.table::setcolorder(ps, c("chr", "unit_id", "start", "end", "n_anchors", "n_agree", "agreement",
                                "maternal_hp", "status", "anchors"))

  ok <- ps[status == "assigned", .(chr, unit_id, maternal_hp)]
  reg <- merge(reg, ok, by = c("chr", "unit_id"), all.x = TRUE)
  reg[, methylated_parent := data.table::fifelse(is.na(maternal_hp), NA_character_,
                               data.table::fifelse(meth_hp_ == maternal_hp, "maternal", "paternal"))]
  reg[, `:=`(meth_maternal = data.table::fifelse(maternal_hp == 1L, meth_hp1, meth_hp2),
             meth_paternal = data.table::fifelse(maternal_hp == 1L, meth_hp2, meth_hp1))]
  reg[, known_region := anchor_name]
  reg[, c("anchor_name", "anchor_allele", "meth_hp_", "implied_", "w_") := NULL]
  data.table::setorder(reg, min_q)
  data.table::setorder(ps, chr, start)
  list(phase_sets = ps[], regions = reg[])
}

#' Parental origin of ASM regions
#'
#' Left: maternal against paternal methylation of every ASM region in a phase
#' set whose parental origin is known (anchors marked). Right: how many
#' non-anchor regions are methylated on the maternal and on the paternal
#' allele, split into monoallelic-like and partial ASM.
#'
#' @param poo Output of [lq_parent_of_origin()].
#' @return A patchwork of two ggplots.
#' @export
lq_plot_parent_of_origin <- function(poo) {
  .need_plot_pkgs(patchwork = TRUE)
  r <- data.table::as.data.table(poo$regions)[!is.na(methylated_parent)]
  if (!nrow(r)) stop("No ASM region lies in a phase set with known parental origin.")
  ps <- data.table::as.data.table(poo$phase_sets)
  r[, kind := data.table::fifelse(is_anchor, "anchor (known DMR)", "other ASM region")]
  p1 <- ggplot2::ggplot(r, ggplot2::aes(.data$meth_maternal, .data$meth_paternal, colour = .data$kind)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey60", linetype = 2) +
    ggplot2::geom_point(ggplot2::aes(size = .data$n_cpg), alpha = 0.6) +
    ggplot2::scale_colour_manual(values = c(`anchor (known DMR)` = .lq_cols[["known"]],
                                            `other ASM region` = .lq_cols[["call"]]), name = NULL) +
    ggplot2::scale_size_area(max_size = 4, guide = "none") +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(x = "methylation, maternal allele", y = "methylation, paternal allele",
                  title = "Parent of origin",
                  subtitle = sprintf("%d phase sets assigned (%d with conflicting anchors), %d ASM regions placed",
                                     sum(ps$status == "assigned"), sum(ps$status == "conflict"), nrow(r))) +
    .lq_theme()
  cnt <- r[is_anchor == FALSE, .N, by = .(methylated_parent, pattern)]
  p2 <- ggplot2::ggplot(cnt, ggplot2::aes(.data$methylated_parent, .data$N, fill = .data$pattern)) +
    ggplot2::geom_col(width = 0.6) +
    ggplot2::scale_fill_manual(values = c(monoallelic = .lq_cols[["call"]], partial = .lq_cols[["partial"]]),
                               name = NULL) +
    ggplot2::labs(x = "methylated allele", y = "ASM regions (without anchors)",
                  title = "Which parent is methylated?",
                  subtitle = "monoallelic regions on one parental side: imprinting candidates") +
    .lq_theme()
  patchwork::wrap_plots(p1, p2, nrow = 1, widths = c(1.2, 1))
}
