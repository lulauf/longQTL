# extract.R
#
# M1: one call per region and sample turns a haplotagged modBAM/modCRAM
# plus its VCFs into the three read-level tables every model downstream
# uses (concept doc, sec. 4.1 and 4.4):
#
#   reads        one row per primary alignment
#   calls        one row per read x CpG: 5mC and 5hmC probability
#   read_alleles one row per read x variant: allele and evidence class
#                  direct     SNV/indel base read off the alignment
#                  direct_sv  read is in Sniffles RNAMES (ALT) or spans the
#                             SV without supporting it (REF)
#                  phased     HP tag of the read + phased GT, same phase set

#' Extract read-level methylation and alleles for one region
#'
#' @param bam Path to an indexed, haplotagged BAM or CRAM with MM/ML tags.
#' @param region "chr:start-end", 1-based inclusive.
#' @param small_variants Output of [lq_read_variants()] for SNVs/indels, or
#'   a path to the VCF (read for this region). NULL = no small variants.
#' @param sv Same for structural variants (Sniffles2 with RNAMES). NULL = none.
#' @param sample_id Prefix for read ids (keeps reads unique across samples).
#' @param reference FASTA for CRAM input.
#' @param min_mapq Minimum mapping quality (default 10).
#' @param skip_supplementary Skip supplementary alignments (default TRUE;
#'   their SEQ is hard-clipped while MM/ML refer to the full read).
#' @param indel_window bp around an indel in which an alignment gap of the
#'   exact length counts as the ALT allele (ONT indels are often shifted).
#' @param sv_flank bp a non-supporting read must extend beyond an SV on both
#'   sides to count as REF.
#' @param sv_window bp around an SV breakpoint in which an alignment event
#'   (deletion / insertion of matching length, +-30 %) or a clip counts.
#' @param min_clip Minimum clip length (bp) at a breakpoint that makes a read
#'   unusable as REF for that SV.
#' @param min_base_qual Minimum base quality for a direct SNV allele.
#' @param phased_max_dist Variants farther than this from a read are not
#'   assigned through phasing.
#'
#' @return Object of class `lq_reads`: list of data.tables `reads`, `calls`
#'   (read_id, chr, pos, strand, q_m, q_h as ML integers 0-255),
#'   `read_alleles` (read_id, variant_id, allele, evidence, bq),
#'   `variants` (all variants used, with `class` small/sv) and `qc`.
#' @export
lq_extract_reads <- function(bam, region, small_variants = NULL, sv = NULL,
                             sample_id = NULL, reference = NULL, min_mapq = 10L,
                             skip_supplementary = TRUE, indel_window = 5L,
                             sv_flank = 500L, sv_window = 500L, min_clip = 100L,
                             min_base_qual = 10L,
                             phased_max_dist = 1e6) {
  if (!file.exists(bam)) stop("BAM/CRAM not found: ", bam)
  rg <- .parse_region(region)
  if (grepl("\\.cram$", bam, ignore.case = TRUE) && is.null(reference)) {
    stop("CRAM input needs `reference`.")
  }
  if (is.character(small_variants)) small_variants <- lq_read_variants(small_variants, region, type = "small")
  if (is.character(sv)) sv <- lq_read_variants(sv, region, type = "sv")

  sv_ok <- if (is.null(sv)) .empty_variants() else sv[chr == rg$chr]
  sm <- if (is.null(small_variants)) .empty_variants() else small_variants[chr == rg$chr & multiallelic == FALSE]
  n_multi <- if (is.null(small_variants)) 0L else sum(small_variants$multiallelic)
  data.table::setorder(sm, pos)

  sv_code <- if (nrow(sv_ok)) data.table::fcase(sv_ok$svtype == "DEL", 1L, sv_ok$svtype == "INS", 2L,
                                                default = 0L) else integer()
  raw <- extract_region_cpp(normalizePath(bam), region, if (is.null(reference)) "" else normalizePath(reference),
                            as.integer(min_mapq), isTRUE(skip_supplementary),
                            sm$pos, sm$ref, sm$alt, as.integer(indel_window),
                            if (nrow(sv_ok)) as.integer(sv_ok$pos) else integer(),
                            if (nrow(sv_ok)) as.integer(sv_ok$end) else integer(),
                            sv_code,
                            if (nrow(sv_ok)) as.integer(abs(ifelse(is.na(sv_ok$svlen), sv_ok$end - sv_ok$pos, sv_ok$svlen))) else integer(),
                            as.integer(sv_window), as.integer(sv_flank), as.integer(min_clip))

  prefix <- if (is.null(sample_id)) "" else paste0(sample_id, ":")
  reads <- data.table::as.data.table(raw$reads)
  reads[, read_id := paste0(prefix, read_name)]
  reads[, strand := ifelse(bitwAnd(flag, 16L) > 0L, "-", "+")]
  reads[, chr := rg$chr]
  if (!is.null(sample_id)) reads[, sample_id := sample_id]

  calls <- data.table::as.data.table(raw$calls)
  calls[, `:=`(read_id = reads$read_id[read], chr = rg$chr)]
  calls <- calls[pos >= rg$start & pos <= rg$end]
  calls[, read := NULL]
  data.table::setcolorder(calls, c("read_id", "chr", "pos", "strand", "q_m", "q_h"))
  if (!is.null(sample_id)) calls[, sample_id := sample_id]

  # direct alleles of small variants
  al <- data.table::as.data.table(raw$alleles)
  direct <- if (nrow(al)) {
    al[, .(read_id = reads$read_id[read], variant_id = sm$variant_id[var],
           allele = allele, evidence = "direct", bq = bq)]
  } else NULL
  if (!is.null(direct)) {
    is_snv <- nchar(sm$ref[match(direct$variant_id, sm$variant_id)]) == 1L &
              nchar(sm$alt[match(direct$variant_id, sm$variant_id)]) == 1L
    direct[is_snv & !is.na(bq) & bq < min_base_qual, allele := NA_integer_]
    direct[!is_snv, evidence := "direct_indel"]
    direct <- direct[!is.na(allele)]
  }

  sv_res <- .sv_alleles(reads, sv_ok, data.table::as.data.table(raw$sv))
  sv_al <- sv_res$alleles
  ph <- .phased_alleles(reads, rbind(sm[, .(variant_id, pos, end = pos, gt_hp1, gt_hp2, phased, ps)],
                                     if (nrow(sv_ok)) sv_ok[, .(variant_id, pos, end, gt_hp1, gt_hp2, phased, ps)]),
                        phased_max_dist)
  read_alleles <- data.table::rbindlist(list(direct, sv_al, ph), use.names = TRUE, fill = TRUE)
  if (!nrow(read_alleles)) {
    read_alleles <- data.table::data.table(read_id = character(), variant_id = character(),
                                           allele = integer(), evidence = character(), bq = integer())
  }

  variants <- data.table::rbindlist(list(
    if (nrow(sm)) sm[, .(variant_id, chr, pos, end = pos + nchar(ref) - 1L, ref, alt, gt,
                         gt_hp1, gt_hp2, phased, ps, class = "small")],
    if (nrow(sv_ok)) sv_ok[, .(variant_id, chr, pos, end, ref, alt, gt, gt_hp1, gt_hp2,
                               phased, ps, class = "sv", svtype, svlen, support)]
  ), use.names = TRUE, fill = TRUE)
    if (!nrow(variants)) {
      variants <- data.table::data.table(variant_id = character(), chr = character(), pos = integer(),
                                       end = integer(), ref = character(), alt = character(),
                                       gt = character(), gt_hp1 = integer(), gt_hp2 = integer(),
                                       phased = logical(), ps = integer(), class = character())
    }

  qc <- c(raw$qc, list(
    n_reads = nrow(reads), n_calls = nrow(calls),
    frac_reads_with_mm = if (nrow(reads)) mean(reads$has_mm == 1L) else NA_real_,
    frac_reads_hp = if (nrow(reads)) mean(!is.na(reads$hp)) else NA_real_,
    n_multiallelic_skipped = n_multi,
    sv_without_rnames = sv_res$no_rnames,
    sv_rnames_unmatched = sv_res$unmatched
  ))
  if (length(sv_res$unmatched)) {
    warning(length(sv_res$unmatched), " SV(s) list supporting reads (RNAMES) but none of them is ",
            "among the reads of this BAM (e.g. ", sv_res$unmatched[1], "). Was the SV VCF called ",
            "on a different BAM? ALT alleles for these SVs rely on the alignment only.", call. = FALSE)
  }

  structure(list(reads = reads[, .(read_id, chr, start, end, strand, mapq, flag, hp, ps,
                                   read_len, n_calls)],
                 calls = calls, read_alleles = read_alleles, variants = variants,
                 region = region, qc = qc),
            class = "lq_reads")
}

#' @export
print.lq_reads <- function(x, ...) {
  q <- x$qc
  cat("<lq_reads>", x$region, "\n")
  cat("  reads:", q$n_reads, " (skipped by flag:", q$skipped_flag, ", by MAPQ:", q$skipped_mapq, ")\n")
  cat("  CpG calls:", q$n_calls, " | reads with MM/ML:", sprintf("%.0f%%", 100 * q$frac_reads_with_mm),
      " | haplotagged:", sprintf("%.0f%%", 100 * q$frac_reads_hp), "\n")
  if (nrow(x$read_alleles)) {
    cat("  read-variant alleles by evidence:\n")
    print(x$read_alleles[, .N, by = evidence])
  }
  invisible(x)
}

# SV alleles from two independent sources:
#   RNAMES  - Sniffles' list of supporting reads
#   CIGAR   - the alignment itself (deletion / insertion of matching length,
#             clips at the breakpoints), computed in extract_region_cpp()
# ALT = in RNAMES or carries the event. REF = spans the SV with flank,
# carries no event, is not clipped at a breakpoint and is not in RNAMES.
# DEL / INS can be genotyped without RNAMES; INV / DUP only with RNAMES
# (an inversion leaves no event in a linear alignment), otherwise skipped.
.sv_alleles <- function(reads, sv, ev) {
  out <- list(alleles = NULL, no_rnames = character(), unmatched = character())
  if (!nrow(sv) || !nrow(reads)) return(out)
  sv <- data.table::copy(sv)
  sv[, idx := .I]
  sv[, has_rn := lengths(rnames) > 0L]
  out$no_rnames <- sv[has_rn == FALSE & svtype != "BND", variant_id]

  if (!nrow(ev)) ev <- data.table::data.table(read = integer(), sv = integer(), event = integer(),
                                               clipped = integer(), spans = integer())
  ev[, read_name := reads$read_name[read]]
  ev[, read_id := reads$read_id[read]]
  ev <- merge(ev, sv[, .(sv = idx, variant_id, svtype, has_rn)], by = "sv")

  # RNAMES reads that are in this BAM region but not in ev (e.g. split reads
  # whose primary alignment does not reach the breakpoint window)
  rn <- if (any(sv$has_rn)) {
    sv[has_rn == TRUE, .(read_name = as.character(unlist(rnames))), by = .(sv = idx, variant_id, svtype, has_rn)]
  } else {
    data.table::data.table(sv = integer(), variant_id = character(), svtype = character(),
                           has_rn = logical(), read_name = character())
  }
  rn <- rn[read_name %in% reads$read_name]
  rn[, read_id := reads$read_id[match(read_name, reads$read_name)]]
  matched <- unique(rn$variant_id)
  out$unmatched <- sv[has_rn == TRUE & !variant_id %in% matched & svtype != "BND" &
                        (is.na(support) | support > 0), variant_id]

  ev[, in_rn := paste(read_name, variant_id) %chin% paste(rn$read_name, rn$variant_id)]
  ev[, allele := data.table::fcase(
    svtype == "BND", NA_integer_,
    in_rn | event == 1L, 1L,
    spans == 1L & clipped == 0L & (svtype %in% c("DEL", "INS") | has_rn), 0L,
    default = NA_integer_)]
  extra_alt <- rn[!paste(read_name, variant_id) %chin% paste(ev$read_name, ev$variant_id) & svtype != "BND",
                  .(read_id, variant_id, allele = 1L)]
  al <- data.table::rbindlist(list(ev[!is.na(allele), .(read_id, variant_id, allele)], extra_alt))
  if (nrow(al)) {
    al[, `:=`(evidence = "direct_sv", bq = NA_integer_)]
    out$alleles <- unique(al, by = c("read_id", "variant_id"))
  }
  out
}

# Phased alleles: read HP (1/2) and PS must match a phased heterozygous
# variant's phase set; allele = GT of that haplotype.
.phased_alleles <- function(reads, v, max_dist) {
  v <- v[phased == TRUE & !is.na(ps) & gt_hp1 != gt_hp2,
         .(variant_id, v_pos = pos, v_end = end, gt_hp1, gt_hp2, ps)]
  r <- reads[!is.na(hp) & hp %in% 1:2 & !is.na(ps), .(read_id, r_start = start, r_end = end, hp, ps)]
  if (!nrow(v) || !nrow(r)) return(NULL)
  j <- merge(r, v, by = "ps", allow.cartesian = TRUE)
  j <- j[pmax(0, v_pos - r_end, r_start - v_end) <= max_dist]
  if (!nrow(j)) return(NULL)
  j[, .(read_id, variant_id, allele = ifelse(hp == 1L, gt_hp1, gt_hp2),
        evidence = "phased", bq = NA_integer_)]
}
