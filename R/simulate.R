# simulate.R
#
# Read-level simulator. This first version produces CpG call tables with a
# known block structure (ground truth for the block-definition modes). The
# full M0 simulator adds samples, alleles, allelic effects (beta_W, beta_B),
# imprinting (tau^2), cell-type mixtures, phasing and call errors on top of
# the same read generator.

#' Simulate read-level CpG calls with known co-methylation segments
#'
#' CpGs sit on a regular grid. The region is split into segments of two
#' kinds: in a **coupled** segment every read draws one epiallele state
#' (methylated with probability `p_meth`) and all its CpGs follow that state
#' up to call noise; in an **independent** segment every CpG on every read is
#' an independent Bernoulli(`p_meth`) draw.
#'
#' @param segments data.frame with columns `n_cpg` and `type`
#'   (`"coupled"` / `"independent"`), in genomic order.
#' @param spacing Distance in bp between neighbouring CpGs.
#' @param n_reads Number of reads.
#' @param read_bp Read length in bp.
#' @param p_meth Methylation probability (state probability for coupled
#'   segments, per-CpG probability for independent ones).
#' @param noise Call error rate.
#' @param chr Chromosome name.
#' @param seed RNG seed.
#'
#' @return list: `calls` (chr, pos, read_id, call) and `truth`
#'   (chr, pos, segment, type).
#' @export
lq_simulate_reads <- function(segments = data.frame(n_cpg = c(10, 10, 10),
                                                    type = c("coupled", "independent", "coupled")),
                              spacing = 50L, n_reads = 200L, read_bp = 3000L,
                              p_meth = 0.5, noise = 0.02, chr = "chrSim", seed = 1L) {
  segments <- data.table::as.data.table(segments)
  stopifnot(all(segments$type %in% c("coupled", "independent")))
  set.seed(seed)

  seg_id <- rep(seq_len(nrow(segments)), segments$n_cpg)
  type   <- rep(segments$type, segments$n_cpg)
  n <- length(seg_id)
  pos <- as.integer(seq_len(n) * spacing)
  span <- max(pos) + spacing

  rows <- vector("list", n_reads)
  for (r in seq_len(n_reads)) {
    s <- sample.int(span + read_bp, 1L) - read_bp
    covered <- which(pos >= s & pos < s + read_bp)
    if (!length(covered)) next
    state <- stats::rbinom(nrow(segments), 1L, p_meth)          # one draw per segment
    truth_call <- ifelse(type[covered] == "coupled",
                         state[seg_id[covered]],
                         stats::rbinom(length(covered), 1L, p_meth))
    flip <- stats::rbinom(length(covered), 1L, noise) == 1L
    rows[[r]] <- data.table::data.table(chr = chr, pos = pos[covered],
                                        read_id = sprintf("read%05d", r),
                                        call = as.integer(xor(truth_call == 1L, flip)))
  }
  list(
    calls = data.table::rbindlist(rows),
    truth = data.table::data.table(chr = chr, pos = pos, segment = seg_id, type = type)
  )
}
