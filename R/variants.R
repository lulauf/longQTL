# variants.R
#
# VCF input for small variants (Clair3, phased) and structural variants
# (Sniffles2, with RNAMES). Plain text parsing with data.table, no
# Bioconductor dependency; region access through tabix when available.

#' Read variants from a VCF
#'
#' @param vcf Path to a (bgzipped) VCF.
#' @param region Optional region "chr:start-end" (1-based, inclusive).
#'   Uses tabix if installed and indexed, otherwise filters after reading.
#' @param sample Sample column (name or 1-based index among samples).
#' @param type `"auto"` (SV if an SVTYPE INFO field is present), `"small"`
#'   or `"sv"`.
#' @param pass_only Keep only FILTER == PASS (or ".").
#'
#' @return data.table: variant_id, chr, pos, ref, alt, gt, gt_hp1, gt_hp2,
#'   phased, ps, and for SVs additionally svtype, end, svlen, support,
#'   rnames (list column of supporting read names). Multi-allelic records
#'   are kept but flagged (`multiallelic`).
#' @export
lq_read_variants <- function(vcf, region = NULL, sample = 1L,
                             type = c("auto", "small", "sv"), pass_only = TRUE) {
  type <- match.arg(type)
  if (!file.exists(vcf)) stop("VCF not found: ", vcf)

  use_tabix <- !is.null(region) && nzchar(Sys.which("tabix")) &&
    (file.exists(paste0(vcf, ".tbi")) || file.exists(paste0(vcf, ".csi")))
  cmd <- if (use_tabix) {
    paste("tabix -h", shQuote(vcf), shQuote(region), "| grep -v '^##'")
  } else {
    paste(if (grepl("\\.gz$", vcf)) "gzip -dc" else "cat", shQuote(vcf), "| grep -v '^##'")
  }
  dt <- data.table::fread(cmd = cmd, sep = "\t", header = TRUE, colClasses = "character",
                          showProgress = FALSE)
  if (!nrow(dt)) return(.empty_variants())
  data.table::setnames(dt, 1L, "CHROM")
  samples <- setdiff(names(dt), c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"))
  scol <- if (is.numeric(sample)) samples[sample] else sample
  if (is.na(scol) || !scol %in% samples) stop("Sample not found in VCF: ", sample)

  if (!is.null(region) && !use_tabix) {
    rg <- .parse_region(region)
    dt <- dt[CHROM == rg$chr & as.integer(POS) >= rg$start & as.integer(POS) <= rg$end]
  }
  if (pass_only) dt <- dt[FILTER %in% c("PASS", ".")]
  if (!nrow(dt)) return(.empty_variants())

  if (type == "auto") type <- if (any(grepl("(^|;)SVTYPE=", dt$INFO))) "sv" else "small"

  gt <- .format_field(dt$FORMAT, dt[[scol]], "GT")
  ps <- .format_field(dt$FORMAT, dt[[scol]], "PS")
  phased <- grepl("|", gt, fixed = TRUE)
  al <- strsplit(gt, "[|/]")
  a1 <- suppressWarnings(as.integer(vapply(al, `[`, "", 1L)))
  a2 <- suppressWarnings(as.integer(vapply(al, function(x) if (length(x) > 1) x[2L] else NA_character_, "")))

  out <- data.table::data.table(
    chr = dt$CHROM, pos = as.integer(dt$POS), ref = toupper(dt$REF), alt = toupper(dt$ALT),
    gt = gt, gt_hp1 = a1, gt_hp2 = a2, phased = phased,
    ps = suppressWarnings(as.integer(ps)),
    multiallelic = grepl(",", dt$ALT, fixed = TRUE)
  )
  out[ps %in% 0L, ps := NA_integer_]

  if (type == "sv") {
    info <- dt$INFO
    out[, `:=`(
      svtype  = .info_field(info, "SVTYPE"),
      end     = suppressWarnings(as.integer(.info_field(info, "END"))),
      svlen   = suppressWarnings(as.integer(.info_field(info, "SVLEN"))),
      support = suppressWarnings(as.integer(.info_field(info, "SUPPORT")))
    )]
    out[is.na(end), end := pos]
    rn <- .info_field(info, "RNAMES")
    out[, rnames := lapply(rn, function(x) if (is.na(x)) character() else strsplit(x, ",", fixed = TRUE)[[1]])]
    out[, variant_id := ifelse(dt$ID != ".", dt$ID, sprintf("%s:%d:%s", chr, pos, svtype))]
    out[, alt := ifelse(nchar(alt) > 50, paste0("<", svtype, ">"), alt)]
    out[, ref := ifelse(nchar(ref) > 50, substr(ref, 1, 1), ref)]
  } else {
    out[, variant_id := sprintf("%s:%d:%s:%s", chr, pos, ref, alt)]
  }
  data.table::setcolorder(out, "variant_id")
  data.table::setorder(out, chr, pos)
  data.table::setattr(out, "variant_type", type)
  out[]
}

.empty_variants <- function() {
  data.table::data.table(variant_id = character(), chr = character(), pos = integer(),
                         ref = character(), alt = character(), gt = character(),
                         gt_hp1 = integer(), gt_hp2 = integer(), phased = logical(),
                         ps = integer(), multiallelic = logical())
}

.format_field <- function(format, values, key) {
  f <- strsplit(format, ":", fixed = TRUE)
  v <- strsplit(values, ":", fixed = TRUE)
  vapply(seq_along(f), function(i) {
    j <- match(key, f[[i]])
    if (is.na(j) || j > length(v[[i]])) NA_character_ else v[[i]][j]
  }, character(1))
}

.info_field <- function(info, key) {
  m <- regmatches(info, regexec(paste0("(?:^|;)", key, "=([^;]*)"), info, perl = TRUE))
  vapply(m, function(x) if (length(x) == 2L) x[2L] else NA_character_, character(1))
}

.parse_region <- function(region) {
  m <- regmatches(region, regexec("^([^:]+):([0-9,]+)-([0-9,]+)$", region))[[1]]
  if (length(m) != 4L) stop("Region must look like chr:start-end, got ", region)
  list(chr = m[2L], start = as.integer(gsub(",", "", m[3L])),
       end = as.integer(gsub(",", "", m[4L])))
}
