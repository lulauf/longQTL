// blocks.cpp
//
// C++ kernels for read co-methylation block definition.
//
// pair_counts_cpp: for every pair of CpGs (a, a + lag), lag = 1..max_lag,
//   counts the reads that carry calls at BOTH CpGs, split into the 2x2 table
//   (1,1) (1,0) (0,1) (0,0). Input rows must be sorted by read, then CpG
//   index, with at most one call per (read, CpG). Cost: O(calls x max_lag).
//
// segment_blocks_cpp: left-to-right segmentation of CpGs into blocks given a
//   pair-linkage matrix (1 linked, 0 not linked, -1 not observed).

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List pair_counts_cpp(IntegerVector read_idx, IntegerVector cpg_idx,
                     IntegerVector call, int n_cpg, int max_lag) {
  const R_xlen_t n = read_idx.size();
  if (cpg_idx.size() != n || call.size() != n) {
    stop("read_idx, cpg_idx and call must have the same length");
  }
  if (max_lag < 1) stop("max_lag must be >= 1");

  const R_xlen_t cells = static_cast<R_xlen_t>(n_cpg) * max_lag;
  IntegerVector n11(cells), n10(cells), n01(cells), n00(cells);

  R_xlen_t i = 0;
  while (i < n) {
    // [i, j) is one read
    R_xlen_t j = i + 1;
    while (j < n && read_idx[j] == read_idx[i]) {
      if (cpg_idx[j] <= cpg_idx[j - 1]) {
        stop("rows must be sorted by read, then strictly increasing CpG index "
             "(duplicate (read, CpG) at row %d)", (int)(j + 1));
      }
      ++j;
    }
    for (R_xlen_t a = i; a < j; ++a) {
      const int ca = cpg_idx[a];
      if (ca < 0 || ca >= n_cpg) stop("cpg_idx out of range");
      for (R_xlen_t b = a + 1; b < j; ++b) {
        const int lag = cpg_idx[b] - ca;
        if (lag > max_lag) break;
        const R_xlen_t cell = static_cast<R_xlen_t>(ca) * max_lag + (lag - 1);
        const int x = call[a], y = call[b];
        if (x == 1 && y == 1)      n11[cell]++;
        else if (x == 1 && y == 0) n10[cell]++;
        else if (x == 0 && y == 1) n01[cell]++;
        else if (x == 0 && y == 0) n00[cell]++;
        else stop("call must be 0 or 1");
      }
    }
    i = j;
  }
  return List::create(_["n11"] = n11, _["n10"] = n10,
                      _["n01"] = n01, _["n00"] = n00);
}

// link: integer vector of length n_cpg * max_lag, cell (a, lag) at
//   a * max_lag + lag - 1. 1 = linked, 0 = not linked, -1 = unobserved.
// Rule for adding CpG j to the open block:
//   - gap pos[j] - pos[j-1] <= max_gap
//   - adjacent pair (j-1, j) observed AND linked
//   - every OBSERVED pair (k, j) with k in the block and j - k <= max_lag
//     is linked (unobserved non-adjacent pairs do not break the block)
//   - block stays within max_cpgs CpGs and max_bp base pairs
// Returns the 0-based index of each CpG's block start.
// [[Rcpp::export]]
IntegerVector segment_blocks_cpp(IntegerVector pos, IntegerVector link,
                                 int max_lag, int max_gap, int max_cpgs,
                                 int max_bp) {
  const int n = pos.size();
  IntegerVector block_start(n);
  if (n == 0) return block_start;
  if (link.size() != static_cast<R_xlen_t>(n) * max_lag) {
    stop("link has wrong length");
  }
  int bs = 0;
  block_start[0] = 0;
  for (int j = 1; j < n; ++j) {
    bool ok = (pos[j] - pos[j - 1]) <= max_gap &&
              (j - bs) < max_cpgs &&
              (pos[j] - pos[bs]) <= max_bp;
    if (ok) {
      const int adj = link[static_cast<R_xlen_t>(j - 1) * max_lag];
      ok = (adj == 1);
    }
    if (ok) {
      const int k0 = std::max(bs, j - max_lag);
      for (int k = k0; k < j - 1; ++k) {
        const int v = link[static_cast<R_xlen_t>(k) * max_lag + (j - k - 1)];
        if (v == 0) { ok = false; break; }
      }
    }
    if (!ok) bs = j;
    block_start[j] = bs;
  }
  return block_start;
}
