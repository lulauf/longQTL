# blocks.R
#
# Methylation phenotype units ("blocks"). Two interchangeable modes behind one
# interface, so both can be run on the same data and compared:
#
#   "comethyl" - blocks from CpG-CpG coupling on the SAME reads (r^2 of the
#                2x2 call table over reads covering both CpGs), in the spirit
#                of methylation haplotype blocks (Guo et al. 2017, Nat Genet
#                49:635). Needs no cohort: works for one sample.
#   "regions"  - fixed tiles of `tile_width` bp, or a supplied BED (CpG
#                islands, cCREs, promoters, ...). Reproducible and comparable
#                across studies, blind to the data.
#
# Both return the same schema, so everything downstream is mode-agnostic.
#
# Genotype-blindness: block definition never sees alleles. If the coupling
# used to build "comethyl" blocks is partly created by an mQTL (ALT reads
# unmethylated, REF reads methylated), blocks may enrich around mQTLs, but
# the null distribution stays valid because permutations of allele labels
# leave the blocks unchanged.

# ---- mode registry -------------------------------------------------------------

.BLOCK_MODES <- list(
  comethyl = list(
    defaults = list(
      max_lag          = 5L,     # compare each CpG with its next 5 CpGs
      r2_min           = 0.5,    # linkage threshold (Guo et al. 2017 use 0.5)
      min_pair_reads   = 6L,     # reads covering both CpGs for a pair to count
      merge_invariant  = TRUE,   # link pairs that are near-constant and concordant
      invariant_maf    = 0.1,    # near-constant: min(p, 1 - p) below this
      invariant_tol    = 0.1,    # |p_a - p_b| tolerance for invariant merging
      max_gap          = 500L,   # bp between adjacent CpGs inside a block
      max_cpgs         = 50L,    # cap CpGs per block
      max_bp           = 5000L,  # cap block span (must be coverable by reads)
      min_cpgs         = 2L,     # smaller blocks become singletons
      keep_singletons  = TRUE    # FALSE drops CpGs outside blocks
    ),
    fun = function(calls, args) .blocks_comethyl(calls, args)
  ),
  regions = list(
    defaults = list(
      bed        = NULL,         # data.frame chr/start/end (0-based start); NULL = tiles
      tile_width = 1000L,
      min_cpgs   = 1L
    ),
    fun = function(calls, args) .blocks_regions(calls, args)
  )
)

#' Available block-definition modes and their default arguments
#' @return Named list: mode -> list of default arguments.
#' @export
lq_block_modes <- function() lapply(.BLOCK_MODES, `[[`, "defaults")

.resolve_block_args <- function(mode, args) {
  defaults <- .BLOCK_MODES[[mode]]$defaults
  unknown <- setdiff(names(args), names(defaults))
  if (length(unknown) > 0L) {
    stop("Unknown block_args for mode '", mode, "': ", paste(unknown, collapse = ", "),
         ". Valid: ", paste(names(defaults), collapse = ", "), call. = FALSE)
  }
  utils::modifyList(defaults, args, keep.null = TRUE)
}

# ---- public entry point ------------------------------------------------------------

#' Define methylation blocks
#'
#' @param calls data.frame/data.table of read-level CpG calls with columns
#'   `chr`, `pos` (1-based CpG position), `read_id`, `call` (0/1). An
#'   optional `sample_id` column is used to make read ids unique across
#'   samples. `"regions"` mode only needs `chr` and `pos`.
#' @param params An [lq_params()] object (uses `block_mode`, `block_args`),
#'   or `NULL` together with `mode`/`args`.
#' @param mode,args Override mode and arguments without a params object.
#'
#' @return An object of class `lq_blocks`: a list with
#'   `blocks` (block_id, chr, start, end, n_cpg, mode) and
#'   `cpg_map` (chr, pos, block_id). Attribute `diagnostics` holds summary
#'   counts.
#' @export
lq_define_blocks <- function(calls, params = NULL, mode = NULL, args = list()) {
  if (!is.null(params)) {
    stopifnot(inherits(params, "lq_params"))
    mode <- params$block_mode
    args <- params$block_args
  } else {
    if (is.null(mode)) stop("Give either params or mode.")
    mode <- match.arg(mode, names(.BLOCK_MODES))
    args <- .resolve_block_args(mode, args)
  }

  calls <- data.table::as.data.table(calls)
  need <- if (mode == "regions") c("chr", "pos") else c("chr", "pos", "read_id", "call")
  miss <- setdiff(need, names(calls))
  if (length(miss)) stop("calls is missing columns: ", paste(miss, collapse = ", "))

  out <- .BLOCK_MODES[[mode]]$fun(calls, args)
  out$blocks[, mode := mode]

  diag <- list(
    mode         = mode,
    n_cpg        = nrow(out$cpg_map),
    n_blocks     = nrow(out$blocks),
    n_singletons = sum(out$blocks$n_cpg == 1L),
    median_cpgs  = stats::median(out$blocks$n_cpg),
    median_bp    = stats::median(out$blocks$end - out$blocks$start + 1)
  )
  structure(out, class = "lq_blocks", diagnostics = diag, args = args)
}

#' @export
print.lq_blocks <- function(x, ...) {
  d <- attr(x, "diagnostics")
  cat("<lq_blocks> mode =", d$mode, "\n")
  cat("  CpGs:", d$n_cpg, "| blocks:", d$n_blocks,
      "| singletons:", d$n_singletons, "\n")
  cat("  median CpGs/block:", d$median_cpgs, "| median span:", d$median_bp, "bp\n")
  invisible(x)
}

# ---- mode: comethyl --------------------------------------------------------------

#' Pairwise co-methylation statistics between nearby CpGs on the same reads
#'
#' @inheritParams lq_define_blocks
#' @param max_lag Compare each CpG with its next `max_lag` CpGs.
#' @return data.table: chr, cpg_a, cpg_b (positions), lag, n11, n10, n01, n00,
#'   n_pair, p_a, p_b, r2 (NA when a margin is empty).
#' @export
lq_pair_comethylation <- function(calls, max_lag = 5L) {
  calls <- data.table::as.data.table(calls)
  .pair_stats(calls, max_lag)
}

.pair_stats <- function(calls, max_lag) {
  dt <- calls[call %in% c(0L, 1L), .(chr, pos = as.integer(pos),
                                    read_id = if ("sample_id" %in% names(calls))
                                      paste(sample_id, read_id, sep = ":") else
                                      as.character(read_id),
                                    call = as.integer(call))]
  dt <- unique(dt, by = c("chr", "pos", "read_id"))

  res <- lapply(split(dt, by = "chr", keep.by = TRUE), function(d) {
    cpgs <- sort(unique(d$pos))
    d[, cpg_idx := match(pos, cpgs) - 1L]
    d[, read_idx := match(read_id, unique(read_id))]
    data.table::setorder(d, read_idx, cpg_idx)

    pc <- pair_counts_cpp(d$read_idx, d$cpg_idx, d$call, length(cpgs), as.integer(max_lag))
    a   <- rep(seq_along(cpgs), each = max_lag)
    lg  <- rep(seq_len(max_lag), times = length(cpgs))
    b   <- a + lg
    keep <- b <= length(cpgs)
    out <- data.table::data.table(
      chr = d$chr[1L], cpg_a = cpgs[a[keep]], cpg_b = cpgs[b[keep]], lag = lg[keep],
      n11 = pc$n11[keep], n10 = pc$n10[keep], n01 = pc$n01[keep], n00 = pc$n00[keep]
    )
    out
  })
  ps <- data.table::rbindlist(res)
  ps[, n_pair := n11 + n10 + n01 + n00]
  ps[, `:=`(p_a = (n11 + n10) / n_pair, p_b = (n11 + n01) / n_pair)]
  # r^2 = D^2 / (pA (1-pA) pB (1-pB)), D = p11 - pA pB
  ps[, r2 := {
    den <- p_a * (1 - p_a) * p_b * (1 - p_b)
    D <- n11 / n_pair - p_a * p_b
    data.table::fifelse(n_pair > 0 & den > 0, D^2 / den, NA_real_)
  }]
  ps[]
}

.blocks_comethyl <- function(calls, a) {
  ps <- .pair_stats(calls, a$max_lag)

  # pair linkage: 1 linked, 0 broken, -1 unobserved
  ps[, link := {
    observed  <- n_pair >= a$min_pair_reads
    # r^2 carries no information when both CpGs are (nearly) constant: one
    # miscalled read then drives it to ~0. Such pairs are linked when their
    # levels agree.
    invariant <- observed & pmin(p_a, 1 - p_a) < a$invariant_maf &
      pmin(p_b, 1 - p_b) < a$invariant_maf
    lk <- data.table::fifelse(!observed, -1L,
            data.table::fifelse(!is.na(r2) & r2 >= a$r2_min, 1L, 0L))
    if (isTRUE(a$merge_invariant)) {
      lk[invariant & abs(p_a - p_b) <= a$invariant_tol] <- 1L
    }
    lk
  }]

  cpg_all <- unique(calls[, .(chr, pos = as.integer(pos))])
  data.table::setorder(cpg_all, chr, pos)

  maps <- lapply(split(cpg_all, by = "chr", keep.by = TRUE), function(cc) {
    n <- nrow(cc)
    link <- rep(-1L, n * a$max_lag)
    pp <- ps[chr == cc$chr[1L]]
    if (nrow(pp)) {
      ia <- match(pp$cpg_a, cc$pos)
      cell <- (ia - 1L) * a$max_lag + pp$lag
      link[cell] <- pp$link
    }
    bs <- segment_blocks_cpp(cc$pos, link, as.integer(a$max_lag), as.integer(a$max_gap),
                             as.integer(a$max_cpgs), as.integer(a$max_bp))
    cc[, block_start_idx := bs]
    cc
  })
  cpg_map <- data.table::rbindlist(maps)

  cpg_map[, n_cpg := .N, by = .(chr, block_start_idx)]
  small <- cpg_map$n_cpg < a$min_cpgs
  if (any(small)) {
    if (isTRUE(a$keep_singletons)) {
      # every CpG of an undersized block becomes its own block
      cpg_map[small, block_start_idx := -seq_len(sum(small))]
      cpg_map[small, n_cpg := 1L]
    } else {
      cpg_map <- cpg_map[!small]
    }
  }
  cpg_map[, block_id := sprintf("%s:%d-%d", chr, min(pos), max(pos)),
          by = .(chr, block_start_idx)]
  blocks <- cpg_map[, .(start = min(pos), end = max(pos), n_cpg = .N),
                    by = .(block_id, chr)]
  data.table::setorder(blocks, chr, start)
  list(blocks = blocks[], cpg_map = cpg_map[, .(chr, pos, block_id)])
}

# ---- mode: regions ------------------------------------------------------------------

.blocks_regions <- function(calls, a) {
  cpg <- unique(calls[, .(chr, pos = as.integer(pos))])
  data.table::setorder(cpg, chr, pos)

  if (is.null(a$bed)) {
    w <- as.integer(a$tile_width)
    if (w < 1L) stop("tile_width must be >= 1")
    cpg[, tile := (pos - 1L) %/% w]
    cpg[, block_id := sprintf("%s:%d-%d", chr, tile * w + 1L, (tile + 1L) * w)]
    cpg[, tile := NULL]
  } else {
    bed <- data.table::as.data.table(a$bed)[, 1:3]
    data.table::setnames(bed, c("chr", "start", "end"))
    bed[, `:=`(start = as.integer(start), end = as.integer(end))]
    # BED is 0-based half-open: CpG at 1-based pos is inside if start < pos <= end
    bed[, block_id := sprintf("%s:%d-%d", chr, start + 1L, end)]
    q <- cpg[, .(chr, s = pos - 1L, e = pos, pos)]
    b <- bed[, .(chr, s = start, e = end - 1L, block_id)]
    q[, e := s]
    data.table::setkey(b, chr, s, e)
    hit <- data.table::foverlaps(q, b, nomatch = NULL, mult = "first")
    cpg <- hit[, .(chr, pos, block_id)]
  }

  cpg[, n_cpg := .N, by = block_id]
  cpg <- cpg[n_cpg >= a$min_cpgs]
  blocks <- cpg[, .(start = min(pos), end = max(pos), n_cpg = .N), by = .(block_id, chr)]
  data.table::setorder(blocks, chr, start)
  list(blocks = blocks[], cpg_map = cpg[, .(chr, pos, block_id)])
}
