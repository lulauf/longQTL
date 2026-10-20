#' longQTL: read-level methylation QTL mapping from long reads
#'
#' Module map (see the concept document for the scientific rationale):
#'
#' * `lq_params()`            - one validated parameter object for a run
#' * `lq_define_blocks()`     - methylation phenotype units, two modes:
#'                              `"comethyl"` (read co-methylation) and
#'                              `"regions"` (fixed tiles or a supplied BED)
#' * `lq_plan_chunks()`       - genome -> load-balanced chunks with a cis halo
#' * `lq_run_chunks()`        - parallel, resumable execution of any chunk worker
#' * `lq_bam_provenance()`    - pore / basecaller / modification model per BAM
#' * `lq_simulate_reads()`    - read-level simulator (ground truth for tests)
#'
#' @keywords internal
#' @useDynLib longQTL, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @import data.table
"_PACKAGE"

# data.table non-standard evaluation column names used across the package
utils::globalVariables(c(
  ".", ".N", ".I", ".SD", "chr", "pos", "read_id", "call", "start", "end",
  "block_id", "n_cpg", "cpg_idx", "read_idx", "chunk_id", "i.block_id",
  "n11", "n10", "n01", "n00", "n_pair", "r2", "lag", "cpg_a", "cpg_b",
  "cum", "var_start", "var_end", "status", "mode", "tile", "region_id",
  "first_pos", "last_pos", "cpg_in_block", "work", "link", "block_start_idx",
  "p_a", "p_b", "s", "e", "sample_id", "chunk_in_chr", "stratum", "pore",
  "basecaller", "dosage", "gt_hp1", "gt_hp2", "allele", "direct", "evidence",
  "is_causal", "conf_", "thr_", "kept_", "prob_mod", "k", "n", "n_reads", "y",
  "y_rn", "p_value", "estimate", "method", "variant_id", "ok", "allele_causal",
  "true_state", "rep", "u", "v", "hap", "strand", "model_tier", "model_version",
  "CHROM", "POS", "FILTER", "multiallelic", "ps", "phased", "svtype", "rnames",
  "read_name", "flag", "read", "bq", "mapq", "hp", "read_len", "n_calls", "var",
  "q_m", "q_h", "support", "svlen", "ref", "alt", "gt", "has_mm", "p_m", "p_h",
  "v_pos", "v_end", "r_start", "r_end", "idx", "has_rn", "in_rn", "event", "clipped",
  "spans", "bin", "mu", "rho_hat", "mu_c", "lr", "dist", "n_0", "n_1", "meth_0",
  "meth_1", "p_final", "p_lrt", "p_perm", "q", "perm_confirmed", "unit_type",
  "n_hp1", "n_hp2", "hap_balance_p", "c_start", "c_end", "in_cpg", "cpg_snp",
  "context_snp", "known_region", "name", "x.N", "unit_id", "asm_call", "delta",
  "n_cpg_masked", "valid", "r2w", "w", "ns", "n_samples", "n_samples_r2", "c_s", "c_e",
  "i.N", "mean_meth", "median_samples", "mean_depth", "cluster", "len", "orig_id",
  "meth", "n_carriers", "class", "r_start", "multiallelic", "block_start_idx",
  "cohort_masked", "n_ref", "n_alt", "f", "v_start",
  ".data", "pref_", "dir_", "region_", "conc_", "cs_", "cx_", "cand_", "known_", "grp", "x",
  "lab", "min_q", "meth_hp1", "meth_hp2", "pattern", "cls", "inside", "frac", "called",
  "tested", "sex", "hp_lab", "state", "row", "span", "pair", "stratum_a", "stratum_b",
  "n_a", "n_b", "p_shift", "p_shift_mid", "level", "diff", "off", "len", "mid", "y", "hit", "txt",
  "is_anchor", "meth_hp_", "implied_", "w_", "anchor_name", "anchor_allele", "maternal_hp",
  "methylated_parent", "meth_maternal", "meth_paternal", "kind", "N", "agreement", "methylated_allele",
  "n_anchors", "n_agree", "anchors", "obs"
))
