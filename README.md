# longQTL

Read-level methylation QTL mapping from Oxford Nanopore long reads.

Every ONT read carries both its allele and its CpG methylation calls. longQTL
models that coupling directly (beta-binomial GLMM over reads) and separates the
**within-person** allele-specific effect from the **between-person** dosage
effect.

## Status (v0.0.1)

| Module | Function | State |
| --- | --- | --- |
| Parameters | `lq_params()` | done: cis window (default 1 Mb, cap 5 Mb), call modes, 5hmC handling |
| Blocks | `lq_define_blocks()` | done: modes `comethyl` and `regions` |
| Parallelisation | `lq_plan_chunks()`, `lq_run_chunks()`, `lq_block_seed()` | done |
| Provenance | `lq_bam_provenance()` | done: pore, basecaller, models, MM/ML + HP tags |
| Simulator | `lq_simulate_reads()`, `lq_sim_scenario()`, `lq_simulate_cohort()` | done (M0): alleles, LD, beta_W/beta_B, tau, imprinting, rho, strata, call probabilities, phasing errors, mapping bias |
| Call modes | `lq_apply_call_mode()` | done: hard_fixed, hard_modkit, expected |
| Baselines | `lq_baselines()`, `lq_sim_benchmark()`, `lq_calibration()` | done: fisher_pooled, lm_mean, wilcoxon_delta, glm_within; gate script in `inst/benchmarks/m0_gate.R` |
| Read extraction (htslib) | `lq_extract_reads()`, `lq_read_variants()`, `lq_mod_prob()` | done (M1): own MM/ML parser (validated call-by-call against pysam), alleles direct / direct_indel / direct_sv / phased; check script `inst/scripts/m1_check_sample.R` |
| Single-patient engine | `lq_patient_scan()`, `lq_block_reads()`, `lq_dispersion()` | done (M2): beta-binomial LRT with empirical-Bayes rho, sequential permutation p-values, phasing-free bimodality, CpG-SNP / context flags; simulators `lq_sim_as_reads()`, `lq_sim_person()`; script `inst/scripts/m2_patient_scan.R` |
| Cohort scan | — | M4 |

## Quick start

```r
library(longQTL)

p <- lq_params(cis_window = 1e6, block_mode = "comethyl")

sim <- lq_simulate_reads()
blocks <- lq_define_blocks(sim$calls, params = p)

chunks <- lq_plan_chunks(blocks, p, target_work = 2e5)
worker <- function(chunk, params) {
  # load reads for chunk$chr:chunk$start-chunk$end,
  # variants for chunk$chr:chunk$var_start-chunk$var_end, then test
  NULL
}
manifest <- lq_run_chunks(chunks, worker, p, out_dir = "results",
                          backend = "multicore", n_workers = 16)
```

## Technical strata

```r
prov <- lq_bam_provenance(list.files("bams", "\\.bam$", full.names = TRUE))
table(prov$stratum)   # e.g. R9.4.1/guppy vs R10.4.1/dorado
```

## Development

```sh
Rscript -e 'Rcpp::compileAttributes(); roxygen2::roxygenise()'
R CMD INSTALL .
Rscript -e 'testthat::test_local()'
```
