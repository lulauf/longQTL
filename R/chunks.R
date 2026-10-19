# chunks.R
#
# Parallelisation. A cis analysis is local: the result for one block only
# needs the reads of that block and the variants within cis_window of it.
# Blocks are therefore independent units of work. Instead of one job per
# chromosome (MAGAR's unit; chr1 is ~10x chr21 and dominates wall time),
# the genome is cut into chunks of consecutive blocks with roughly equal
# work. Each chunk records its own variant range (block range +/- cis
# window, the "halo"), so a worker loads only what it needs.
#
# Reproducibility does not depend on the number of workers: every block
# gets a seed derived from (global seed, block_id) by a fixed hash.

#' Split blocks into load-balanced chunks with a cis halo
#'
#' @param blocks `lq_blocks` object or data.table with block_id, chr, start,
#'   end, n_cpg.
#' @param params [lq_params()] object (uses `cis_window`).
#' @param target_work Approximate work per chunk, in units of `n_cpg`
#'   (default 2e5 CpGs). Smaller chunks = better load balance and smaller
#'   memory per worker, more scheduling overhead.
#' @param max_span_bp Hard cap on the block span of one chunk (default 20 Mb),
#'   which bounds the variant range a worker has to load.
#'
#' @return data.table, one row per chunk: chunk_id, chr, start, end (span of
#'   its blocks), var_start, var_end (variant range incl. halo), n_blocks,
#'   work, first_block, last_block.
#' @export
lq_plan_chunks <- function(blocks, params, target_work = 2e5, max_span_bp = 2e7) {
  stopifnot(inherits(params, "lq_params"))
  b <- if (inherits(blocks, "lq_blocks")) blocks$blocks else data.table::as.data.table(blocks)
  miss <- setdiff(c("block_id", "chr", "start", "end", "n_cpg"), names(b))
  if (length(miss)) stop("blocks is missing: ", paste(miss, collapse = ", "))
  .check_scalar_num(target_work, "target_work", lower = 1)
  .check_scalar_num(max_span_bp, "max_span_bp", lower = 1)

  b <- data.table::copy(b)
  data.table::setorder(b, chr, start)

  b[, chunk_in_chr := {
    # greedy left-to-right: open a new chunk when work or span would overflow
    id <- integer(.N); cur <- 1L; w <- 0; s0 <- start[1L]
    for (k in seq_len(.N)) {
      if (k > 1L && (w + n_cpg[k] > target_work || end[k] - s0 > max_span_bp)) {
        cur <- cur + 1L; w <- 0; s0 <- start[k]
      }
      w <- w + n_cpg[k]
      id[k] <- cur
    }
    id
  }, by = chr]

  ch <- b[, .(start = min(start), end = max(end), n_blocks = .N, work = sum(n_cpg),
              first_block = block_id[1L], last_block = block_id[.N]),
          by = .(chr, chunk_in_chr)]
  ch[, `:=`(var_start = pmax(1, start - params$cis_window),
            var_end   = end + params$cis_window)]
  data.table::setorder(ch, -work)                 # largest first: better packing
  ch[, chunk_id := sprintf("chunk_%05d", seq_len(.N))]
  ch[, chunk_in_chr := NULL]
  data.table::setcolorder(ch, c("chunk_id", "chr", "start", "end", "var_start",
                                "var_end", "n_blocks", "work"))
  ch[]
}

#' Deterministic per-block seed
#'
#' 32-bit FNV-1a hash of `block_id`, mixed with the global seed. Identical on
#' every platform and independent of worker count or processing order.
#'
#' @param block_id Character vector.
#' @param seed Global seed (integer).
#' @return Integer vector of seeds.
#' @export
lq_block_seed <- function(block_id, seed) {
  vapply(as.character(block_id), function(s) {
    h <- 2166136261                                     # FNV offset basis
    for (byte in as.integer(charToRaw(enc2utf8(s)))) {
      h <- .xor32(h, byte)
      # h * 16777619 mod 2^32, exactly, in doubles: 16777619 = 2^24 + 403
      h <- ((h %% 256) * 16777216 + h * 403) %% 4294967296
    }
    h <- .xor32(h, as.numeric(seed) %% 4294967296)
    as.integer(h - 2147483648)                          # map to signed int range
  }, integer(1), USE.NAMES = FALSE)
}

.xor32 <- function(a, b) {
  # bitwise xor of two numbers in [0, 2^32) without 32-bit signed overflow
  hi <- bitwXor(as.integer(a %/% 65536), as.integer(b %/% 65536))
  lo <- bitwXor(as.integer(a %% 65536), as.integer(b %% 65536))
  hi * 65536 + lo
}

#' Run a worker over chunks, in parallel and resumable
#'
#' @param chunks Output of [lq_plan_chunks()].
#' @param worker Function `function(chunk, params, ...)` returning a
#'   data.frame, NULL or any other R object (stored as is). It receives one chunk row (a one-row data.table).
#' @param params [lq_params()] object, passed to the worker.
#' @param out_dir Directory for per-chunk results (`<chunk_id>.rds`) and the
#'   manifest. Existing results are reused, so an interrupted run resumes.
#' @param backend `"sequential"`, `"multicore"` (forked, Unix) or `"future"`
#'   (any `future::plan()`, e.g. a SLURM cluster via future.batchtools).
#' @param n_workers Workers for `"multicore"`.
#' @param ... Passed to `worker`.
#'
#' @return The manifest: chunk_id, status (done / failed / cached), seconds,
#'   n_rows, error. Results stay on disk; read them with [lq_collect()].
#' @export
lq_run_chunks <- function(chunks, worker, params, out_dir,
                          backend = c("sequential", "multicore", "future"),
                          n_workers = 1L, ...) {
  backend <- match.arg(backend)
  stopifnot(is.function(worker), inherits(params, "lq_params"))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  run_one <- function(k) {
    ch <- chunks[k]
    f <- file.path(out_dir, paste0(ch$chunk_id, ".rds"))
    if (file.exists(f)) {
      return(data.table::data.table(chunk_id = ch$chunk_id, status = "cached",
                                    seconds = NA_real_, n_rows = NA_integer_,
                                    error = NA_character_))
    }
    t0 <- proc.time()[["elapsed"]]
    res <- tryCatch(worker(ch, params, ...), error = function(e) e)
    dt <- proc.time()[["elapsed"]] - t0
    if (inherits(res, "error")) {
      return(data.table::data.table(chunk_id = ch$chunk_id, status = "failed",
                                    seconds = dt, n_rows = NA_integer_,
                                    error = conditionMessage(res)))
    }
    tmp <- paste0(f, ".tmp")                     # atomic write: no half files
    saveRDS(res, tmp)
    file.rename(tmp, f)
    data.table::data.table(chunk_id = ch$chunk_id, status = "done", seconds = dt,
                           n_rows = if (is.null(res)) 0L else if (is.data.frame(res)) nrow(res) else NA_integer_,
                           error = NA_character_)
  }

  idx <- seq_len(nrow(chunks))
  rows <- switch(backend,
    sequential = lapply(idx, run_one),
    multicore  = {
      if (.Platform$OS.type != "unix") stop("multicore backend needs a Unix OS.")
      parallel::mclapply(idx, run_one, mc.cores = n_workers, mc.preschedule = FALSE)
    },
    future = {
      if (!requireNamespace("future.apply", quietly = TRUE)) {
        stop("backend = 'future' needs the future.apply package.")
      }
      future.apply::future_lapply(idx, run_one, future.seed = TRUE)
    }
  )
  manifest <- data.table::rbindlist(rows)
  data.table::fwrite(manifest, file.path(out_dir, "manifest.tsv"), sep = "\t")

  n_fail <- sum(manifest$status == "failed")
  if (n_fail > 0L) {
    warning(n_fail, " chunk(s) failed; see manifest.tsv. Re-run to retry only those.")
  }
  manifest[]
}

#' Collect per-chunk results from disk
#' @param out_dir Directory used in [lq_run_chunks()].
#' @return data.table of all chunk results row-bound.
#' @export
lq_collect <- function(out_dir) {
  files <- list.files(out_dir, pattern = "^chunk_.*\\.rds$", full.names = TRUE)
  data.table::rbindlist(lapply(files, readRDS), use.names = TRUE, fill = TRUE)
}
