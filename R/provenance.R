# provenance.R
#
# Sequencing provenance per BAM: pore chemistry, basecaller and its models.
# These define the technical strata that must enter the model as covariates
# (between component) and the strata for replication. Within-person
# comparisons share chemistry by construction; between-person comparisons
# do not.

#' Read pore, basecaller and modification models from BAM headers
#'
#' Parses `@RG` and `@PG` header lines (Dorado writes `basecall_model=` and
#' `modbase_models=` into `@RG DS`; command lines in `@PG CL` are also
#' searched). Optionally peeks at the first reads to confirm that MM/ML
#' (modification) and HP (haplotag) tags are present.
#'
#' @param bams Character vector of BAM/CRAM paths. If named, the names are
#'   used as `sample_id`.
#' @param check_tags Inspect the first `n_reads` alignments for MM/ML and
#'   HP tags (needs samtools). Only those reads are streamed; the file is
#'   never read in full.
#' @param n_reads Number of alignments to inspect.
#' @param reference FASTA of the reference genome. Needed to read records
#'   from CRAM files (headers work without it). Ignored for BAM.
#' @param header_text Optional list of character vectors (one per BAM) with
#'   header lines; if given, samtools is not called (useful for tests or
#'   when headers were exported elsewhere).
#'
#' @return data.table, one row per BAM: sample_id, file, pore, kit,
#'   basecaller, basecaller_version, model_tier (sup/hac/fast),
#'   model_version (e.g. v5.0.0), mod_model (e.g. 5mCG_5hmCG),
#'   basecall_model, modbase_models, info_source (`"header"` or `"path"`:
#'   where the pore was found), aligner, aligner_version, haplotagger,
#'   has_mm_ml, has_hp, mod_codes (MM codes seen in the reads, e.g.
#'   `"C+h?;C+m?"`), stratum (pore/tier/model version).
#'
#' @details Dorado writes model names into `@RG` and `@PG`. When reads are
#'   converted to FASTQ before alignment those lines are lost; the
#'   information then survives only in file names, which are parsed as a
#'   fallback (`info_source = "path"`). `mod_codes` is read from the reads
#'   themselves and is the authoritative answer to "was 5hmC called
#'   separately?" (`h` present).
#' @export
lq_bam_provenance <- function(bams, check_tags = TRUE, n_reads = 1000L,
                              reference = NULL, header_text = NULL) {
  sample_ids <- if (!is.null(names(bams))) names(bams) else rep(NA_character_, length(bams))
  bams <- unname(bams)
  if (!is.null(header_text)) check_tags <- FALSE

  # validate everything before touching any file
  is_cram <- grepl("\\.cram$", bams, ignore.case = TRUE)
  if (check_tags && any(is_cram)) {
    if (is.null(reference)) {
      stop("CRAM input needs `reference` (the FASTA used for alignment) to ",
           "read records; or set check_tags = FALSE to parse headers only.")
    }
    if (!file.exists(reference)) stop("reference not found: ", reference)
  }

  if (is.null(header_text)) {
    if (!nzchar(Sys.which("samtools"))) {
      stop("samtools not found on PATH; install it or pass header_text.")
    }
    missing_files <- bams[!file.exists(bams)]
    if (length(missing_files)) stop("Files not found: ", paste(missing_files, collapse = ", "))
    header_text <- lapply(bams, .read_header)
  }
  if (length(header_text) != length(bams)) stop("header_text must match bams.")

  rows <- Map(function(f, h, sid, cram) {
    p <- .parse_bam_header(h, path = f)
    p$sample_id <- sid
    p$file <- f
    if (check_tags) {
      recs <- .head_records(f, n_reads, if (cram) reference else NULL)
      p$has_mm_ml <- any(grepl("\tM[Mm]:Z:", recs)) && any(grepl("\tM[Ll]:B:C", recs))
      p$has_hp    <- any(grepl("\tHP:i:", recs))
      p$mod_codes <- .mod_codes(recs)
    } else {
      p$has_mm_ml <- NA
      p$has_hp    <- NA
      p$mod_codes <- NA_character_
    }
    data.table::as.data.table(p)
  }, bams, header_text, sample_ids, is_cram)

  out <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  out[, stratum := paste(pore, model_tier, model_version, sep = "/")]
  n_unknown <- sum(is.na(out$pore))
  if (n_unknown > 0L) {
    warning(n_unknown, " file(s) without any pore information in header or path; ",
            "fill `pore` / `model_version` by hand before using `stratum`.")
  }
  data.table::setcolorder(out, c("sample_id", "file"))
  n_strata <- data.table::uniqueN(out$stratum)
  if (n_strata > 1L) {
    message(n_strata, " technical strata found: ",
            paste(unique(out$stratum), collapse = ", "),
            ". Use 'stratum' as a covariate and replicate within strata.")
  }
  out[]
}

# Header lines of one file; a file samtools cannot read is an error that
# names the file, not a warning from system2.
.read_header <- function(f) {
  h <- suppressWarnings(system2("samtools", c("view", "-H", shQuote(f)),
                                stdout = TRUE, stderr = FALSE))
  status <- attr(h, "status")
  if (!is.null(status) && status != 0L) {
    stop("samtools could not read the header of ", f,
         " (corrupt, truncated or not a BAM/CRAM file).", call. = FALSE)
  }
  h
}

# Stream only the first n alignment records; `head` closes the pipe, so
# samtools stops after n lines instead of decoding the whole file.
.head_records <- function(f, n, reference = NULL) {
  ref <- if (is.null(reference)) "" else paste("-T", shQuote(reference))
  cmd <- sprintf("samtools view %s %s 2>/dev/null | head -n %d", ref, shQuote(f), as.integer(n))
  suppressWarnings(system2("sh", c("-c", shQuote(cmd)), stdout = TRUE))
}

#' Find alignment files in a per-sample folder layout
#'
#' Expands a path pattern containing `*` (e.g.
#' `"/data/1kg/*/alignment/*.bam"`) and takes the sample id from one path
#' component, by default the one where the first `*` stands.
#'
#' @param pattern Glob pattern; may contain several `*`.
#' @param id_level Which path component (counted from the start, after
#'   splitting on "/") holds the sample id. Default: the component that
#'   contains the first `*` in `pattern`.
#' @param exclude Regex of files to drop (default drops index files).
#' @return Named character vector: names = sample ids, values = paths. Warns
#'   if a sample has more than one file.
#' @export
lq_find_alignments <- function(pattern, id_level = NULL,
                               exclude = "\\.(bai|crai|csi)$") {
  files <- Sys.glob(pattern)
  files <- files[!grepl(exclude, files)]
  if (!length(files)) stop("No files match: ", pattern)
  parts <- strsplit(pattern, "/", fixed = TRUE)[[1]]
  if (is.null(id_level)) id_level <- which(grepl("*", parts, fixed = TRUE))[1L]
  ids <- vapply(strsplit(files, "/", fixed = TRUE), `[`, character(1), id_level)
  dup <- unique(ids[duplicated(ids)])
  if (length(dup)) {
    warning(length(dup), " sample(s) have several files (e.g. ", dup[1L],
            "); all are kept, check before running.")
  }
  stats::setNames(files, ids)
}

.first_match <- function(x, pattern, group = 1L) {
  m <- regmatches(x, regexec(pattern, x, perl = TRUE))
  hits <- vapply(m, function(v) if (length(v) > group) v[group + 1L] else NA_character_,
                 character(1))
  hits <- hits[!is.na(hits)]
  if (length(hits)) hits[1L] else NA_character_
}

.pg_field <- function(pg_lines, program_regex, field) {
  hit <- pg_lines[grepl(paste0("\tPN:", program_regex), pg_lines, ignore.case = TRUE)]
  if (!length(hit)) return(NA_character_)
  .first_match(hit, paste0("\t", field, ":([^\t]+)"))
}

# Provenance is read in three layers, most reliable first:
#   1. structured header fields: @RG DS basecall_model= / modbase_models=,
#      @PG PN:dorado VN:
#   2. the model name itself (dna_r10.4.1_e8.2_400bps_sup@v5.0.0 encodes
#      pore, tier and version)
#   3. naming tokens anywhere in the header text or the file path, e.g.
#      "R10-LSK114-dorado090_sup_5mCG_5hmCG_v500" - the only trace left when
#      reads went through FASTQ before alignment and @RG/@PG were lost.
# `info_source` records which layer the pore came from.
.parse_bam_header <- function(lines, path = "") {
  lines <- as.character(lines)
  rg <- lines[startsWith(lines, "@RG")]
  pg <- lines[startsWith(lines, "@PG")]
  all <- paste(c(lines, path), collapse = "\n")

  # layer 1
  basecall_model <- .first_match(c(rg, pg), "basecall_model=([^\\s\t]+)")
  if (is.na(basecall_model)) {
    basecall_model <- .first_match(c(pg, rg), "(dna_r[0-9.]+_[A-Za-z0-9_.]+@v[0-9.]+)")
  }
  modbase <- .first_match(c(rg, pg), "modbase_models=([^\\s\t]+)")
  if (is.na(modbase)) modbase <- .first_match(pg, "--modified-bases(?:-models)?[ =]([^\t ]+)")
  # placeholders such as "modbase_model_version_id" are not model names
  if (!is.na(modbase) && !grepl("5m|6mA|4mC|@v", modbase)) modbase <- NA_character_

  dorado_pg <- any(grepl("\tPN:dorado", pg, ignore.case = TRUE))
  bc_version <- if (dorado_pg) .pg_field(pg, "dorado", "VN") else NA_character_

  # layer 2: decode the model name
  pore <- tier <- version <- NA_character_
  if (!is.na(basecall_model)) {
    pore    <- .first_match(basecall_model, "dna_r([0-9.]+)_")
    if (!is.na(pore)) pore <- paste0("R", pore)
    tier    <- .first_match(basecall_model, "_(sup|hac|fast)@")
    version <- .first_match(basecall_model, "@(v[0-9.]+)")
  }
  source <- if (!is.na(pore)) "header" else NA_character_

  # layer 3: naming tokens
  tok <- .parse_tokens(all)
  if (is.na(pore) && !is.na(tok$pore)) { pore <- tok$pore; source <- "path" }
  if (is.na(tier))       tier <- tok$tier
  if (is.na(version))    version <- tok$version
  if (is.na(bc_version)) bc_version <- tok$dorado_version
  mod_model <- if (!is.na(modbase)) .first_match(modbase, "(5mCG_5hmCG|5mC_5hmC|5mCG|5mC|6mA)") else tok$mod_model

  basecaller <- if (dorado_pg || !is.na(tok$dorado_version) ||
                    grepl("@v[0-9]", basecall_model %||% "")) "dorado" else
                if (grepl("guppy", all, ignore.case = TRUE)) "guppy" else NA_character_
  if (identical(basecaller, "guppy")) bc_version <- .pg_field(pg, "guppy[^\t]*", "VN")

  list(
    pore               = pore,
    kit                = tok$kit,
    basecaller         = basecaller,
    basecaller_version = bc_version,
    model_tier         = tier,
    model_version      = version,
    mod_model          = mod_model,
    basecall_model     = basecall_model,
    modbase_models     = modbase,
    info_source        = source,
    aligner            = if (grepl("PN:minimap2", all)) "minimap2" else NA_character_,
    aligner_version    = .pg_field(pg, "minimap2", "VN"),
    haplotagger        = if (grepl("PN:whatshap", all, ignore.case = TRUE)) "whatshap" else
                         if (grepl("PN:longphase", all, ignore.case = TRUE)) "longphase" else NA_character_
  )
}

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1L && is.na(a))) b else a

.parse_tokens <- function(text) {
  pore <- .first_match(text, "[Rr](9\\.4\\.1|10\\.4\\.1)")
  if (!is.na(pore)) {
    pore <- paste0("R", pore)
  } else if (grepl("(^|[^A-Za-z0-9])R10([^0-9.]|$)", text)) {
    # LSK114 / e8.2 / 400bps only exist for R10.4.1
    pore <- if (grepl("LSK114|e8\\.2|400bps", text)) "R10.4.1" else "R10"
  } else if (grepl("(^|[^A-Za-z0-9])R9([^0-9.]|$)", text)) {
    pore <- if (grepl("LSK109|LSK110|450bps", text)) "R9.4.1" else "R9"
  }
  dv <- regmatches(text, regexec("dorado[_-]?v?([0-9])\\.?([0-9])\\.?([0-9])", text,
                                 ignore.case = TRUE))[[1]]
  list(
    pore           = pore,
    kit            = .first_match(text, "(?:SQK-)?((?:LSK|ULK|RBK|NBD|LRK)[0-9]{3})"),
    dorado_version = if (length(dv) == 4L) paste(dv[2:4], collapse = ".") else NA_character_,
    tier           = tolower(.first_match(text, "(?i)(?:^|[_@./-])(sup|hac|fast)(?=[_@./-]|$)")),
    version        = {
      v <- regmatches(text, regexec("_v([0-9])([0-9])([0-9])(?:[^0-9]|$)", text))[[1]]
      if (length(v) == 4L) paste0("v", paste(v[2:4], collapse = ".")) else NA_character_
    },
    mod_model      = .first_match(text, "(5mCG_5hmCG|5mC_5hmC|5mCG|5mC|6mA)")
  )
}

# Modification codes present in MM tags, e.g. "C+h?;C+m?"
.mod_codes <- function(recs) {
  mm <- regmatches(recs, regexpr("\tM[Mm]:Z:[^\t]*", recs))
  if (!length(mm)) return(NA_character_)
  parts <- unlist(strsplit(sub("^\tM[Mm]:Z:", "", mm), ";", fixed = TRUE))
  parts <- parts[nzchar(parts)]
  codes <- sort(unique(sub(",.*$", "", parts)))
  paste(codes, collapse = ";")
}
