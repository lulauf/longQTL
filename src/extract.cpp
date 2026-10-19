// extract.cpp
//
// Read-level extraction from a modBAM / modCRAM region with htslib:
//   * per read: name, flag, mapq, reference span, strand, HP / PS tags
//   * per read x CpG: modification qualities for 5mC (m) and 5hmC (h)
//   * per read x small variant: the allele the read carries (direct evidence)
//
// Coordinates: CpG positions are reported strand-combined as the 1-based
// position of the C of the CpG on the forward reference strand (as modkit
// with --combine-strands). A reverse-strand read measures the C of the
// complementary strand, which sits on the G of the forward CpG, so its call
// is shifted one base to the left.
//
// MM modes: for '?' (unknown) only listed positions are calls; for '.' or no
// mode flag (implicit), unlisted C's in CpG context are calls of
// "unmodified" (quality 0) - see SAMtags.pdf, section on MM/ML. MM/ML are
// parsed by parse_mm_ml() below, not by htslib's basemod API.
//
// Supplementary alignments are skipped by default: their SEQ is
// hard-clipped while MM/ML describe the full read, so positions do not match.

#include <Rcpp.h>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <string>
#include <vector>
#include <cstring>
#include <algorithm>
#include <cctype>
#include <cstdlib>

using namespace Rcpp;

namespace {

inline char seq_base(const uint8_t* seq, int i) {
  return seq_nt16_str[bam_seqi(seq, i)];
}

// true if no query position before q is aligned yet
inline bool q2r_empty(const std::vector<int>& q2r, int q) {
  for (int i = 0; i < q && i < (int)q2r.size(); ++i) if (q2r[i] >= 0) return false;
  return true;
}

struct ModCalls {
  bool has_m = false, has_h = false;
  bool m_implicit = false, h_implicit = false;
  bool ok = true;                 // false = malformed tags, read gets no calls
  std::vector<int> qm, qh;        // ML quality per SEQ index, -1 = not listed
};

// Own MM/ML parser following the SAM tags specification (section 1.7):
//   MM:Z:C+h?,d1,d2,...;C+m?,d1,...;   ML:B:C,q1,q2,...
// * positions count occurrences of the canonical base in the ORIGINAL read
//   orientation; for reverse-mapped reads that is the reverse complement of
//   SEQ, i.e. G's in SEQ counted from the end
// * an entry with k codes and n positions consumes n*k ML values,
//   interleaved by position; entries consume ML in the order they appear
// * mode '?' = unlisted bases are unknown; '.' or no flag = unmodified
// Written in-house because htslib < 1.16 reads the '?' / '.' flag as an
// extra modification code and rejects every Dorado record.
ModCalls parse_mm_ml(const bam1_t* b) {
  ModCalls mc;
  const int qlen = b->core.l_qseq;
  uint8_t* mm = bam_aux_get(b, "MM");
  if (!mm) mm = bam_aux_get(b, "Mm");
  if (!mm) return mc;
  const char* s = bam_aux2Z(mm);
  if (!s) return mc;
  uint8_t* mn = bam_aux_get(b, "MN");
  if (mn && bam_aux2i(mn) != qlen) { mc.ok = false; return mc; }   // SEQ trimmed
  uint8_t* ml = bam_aux_get(b, "ML");
  if (!ml) ml = bam_aux_get(b, "Ml");
  const uint32_t n_ml = ml ? bam_auxB_len(ml) : 0;

  mc.qm.assign(qlen, -1);
  mc.qh.assign(qlen, -1);
  const bool rev = bam_is_rev(b);
  const uint8_t* seq = bam_get_seq(b);

  // SEQ indices of C (forward) / G read from the end (reverse), lazily built
  std::vector<int> c_index;
  bool c_built = false;
  auto build_c = [&]() {
    c_index.clear();
    if (!rev) {
      for (int i = 0; i < qlen; ++i) if (seq_base(seq, i) == 'C') c_index.push_back(i);
    } else {
      for (int i = qlen - 1; i >= 0; --i) if (seq_base(seq, i) == 'G') c_index.push_back(i);
    }
    c_built = true;
  };

  std::string str(s);
  size_t start = 0;
  uint32_t ml_cursor = 0;
  while (start < str.size()) {
    size_t end = str.find(';', start);
    if (end == std::string::npos) end = str.size();
    const std::string item = str.substr(start, end - start);
    start = end + 1;
    if (item.size() < 3) continue;

    const size_t comma = item.find(',');
    const std::string head = item.substr(0, comma);       // e.g. C+h?  C+mh  C+m.
    const char base = head[0];
    const char strand = head[1];
    std::string codes = head.substr(2);
    bool implicit = true;
    if (!codes.empty() && (codes.back() == '?' || codes.back() == '.')) {
      implicit = codes.back() == '.';
      codes.pop_back();
    }
    // ChEBI numeric codes count as one code
    int k = 0;
    if (!codes.empty() && isdigit((unsigned char)codes[0])) k = 1; else k = (int)codes.size();
    std::vector<int> deltas;
    if (comma != std::string::npos) {
      size_t p = comma + 1;
      while (p < item.size()) {
        size_t q = item.find(',', p);
        if (q == std::string::npos) q = item.size();
        if (q > p) deltas.push_back(std::atoi(item.substr(p, q - p).c_str()));
        p = q + 1;
      }
    }
    const uint32_t need = (uint32_t)deltas.size() * (uint32_t)k;
    if (ml_cursor + need > n_ml) { mc.ok = false; return mc; }

    const bool relevant = (base == 'C' && strand == '+');
    if (relevant) {
      for (int j = 0; j < k; ++j) {
        if (codes[j] == 'm') { mc.has_m = true; mc.m_implicit = implicit; }
        if (codes[j] == 'h') { mc.has_h = true; mc.h_implicit = implicit; }
      }
      if (!c_built) build_c();
      long idx = -1;
      for (size_t i = 0; i < deltas.size(); ++i) {
        idx += (long)deltas[i] + 1;
        if (idx >= (long)c_index.size()) { mc.ok = false; return mc; }
        const int qpos = c_index[idx];
        for (int j = 0; j < k; ++j) {
          const int qv = (int)bam_auxB2i(ml, ml_cursor + (uint32_t)(i * k + j));
          if (codes[j] == 'm') mc.qm[qpos] = qv;
          else if (codes[j] == 'h') mc.qh[qpos] = qv;
        }
      }
    }
    ml_cursor += need;
  }
  return mc;
}

}  // namespace

// [[Rcpp::export]]
List extract_region_cpp(std::string path, std::string region, std::string reference,
                        int min_mapq, bool skip_supplementary,
                        IntegerVector var_pos, CharacterVector var_ref,
                        CharacterVector var_alt, int indel_window,
                        IntegerVector sv_pos, IntegerVector sv_end,
                        IntegerVector sv_type, IntegerVector sv_len,
                        int sv_window, int sv_flank, int min_clip) {
  // sv_type: 1 = DEL, 2 = INS, 0 = other (INV, DUP: breakpoints only)
  const int ns = sv_pos.size();
  samFile* fp = sam_open(path.c_str(), "r");
  if (!fp) stop("cannot open " + path);
  if (!reference.empty()) {
    if (hts_set_fai_filename(fp, reference.c_str()) != 0) {
      sam_close(fp);
      stop("cannot use reference " + reference);
    }
  }
  sam_hdr_t* hdr = sam_hdr_read(fp);
  if (!hdr) { sam_close(fp); stop("cannot read header of " + path); }
  hts_idx_t* idx = sam_index_load(fp, path.c_str());
  if (!idx) { sam_hdr_destroy(hdr); sam_close(fp); stop("no index (.bai/.crai) for " + path); }
  hts_itr_t* itr = sam_itr_querys(idx, hdr, region.c_str());
  if (!itr) {
    hts_idx_destroy(idx); sam_hdr_destroy(hdr); sam_close(fp);
    stop("invalid region or contig not in header: " + region);
  }

  // variants: sorted by position (checked)
  const int nv = var_pos.size();
  std::vector<int> vpos(nv);
  std::vector<std::string> vref(nv), valt(nv);
  for (int i = 0; i < nv; ++i) {
    vpos[i] = var_pos[i];
    vref[i] = as<std::string>(var_ref[i]);
    valt[i] = as<std::string>(var_alt[i]);
    if (i > 0 && vpos[i] < vpos[i - 1]) {
      hts_itr_destroy(itr); hts_idx_destroy(idx); sam_hdr_destroy(hdr); sam_close(fp);
      stop("variants must be sorted by position");
    }
  }

  // read table
  std::vector<std::string> r_name;
  std::vector<int> r_flag, r_mapq, r_start, r_end, r_hp, r_ps, r_len, r_ncall;
  std::vector<int> r_has_mm;
  // call table
  std::vector<int> c_read, c_pos, c_qm, c_qh;
  std::vector<char> c_strand;
  // allele table
  std::vector<int> a_read, a_var, a_allele, a_bq;
  // SV evidence table: per read x SV, alignment-level evidence
  std::vector<int> s_read, s_sv, s_event, s_clipped, s_spans;

  bam1_t* b = bam_init1();
  int n_skipped_flag = 0, n_skipped_mapq = 0, n_bad_mm = 0;

  int ret;
  while ((ret = sam_itr_next(fp, itr, b)) >= 0) {
    const uint16_t flag = b->core.flag;
    if (flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP)) { ++n_skipped_flag; continue; }
    if (skip_supplementary && (flag & BAM_FSUPPLEMENTARY)) { ++n_skipped_flag; continue; }
    if (b->core.qual < min_mapq) { ++n_skipped_mapq; continue; }

    const int ridx = (int)r_name.size();   // 0-based read index
    const int qlen = b->core.l_qseq;
    const hts_pos_t rstart = b->core.pos;  // 0-based
    const hts_pos_t rend = bam_endpos(b);  // 0-based exclusive
    const bool rev = bam_is_rev(b);
    const uint8_t* seq = bam_get_seq(b);
    const uint8_t* qual = bam_get_qual(b);

    // query position -> reference position (0-based), -1 if not aligned
    std::vector<int> q2r(qlen, -1);
    // reference offset -> query position, -1 deleted; and per reference
    // offset the length of an insertion following it / deletion starting at it
    const int span = (int)(rend - rstart);
    std::vector<int> r2q(span > 0 ? span : 0, -1);
    std::vector<int> ins_after(span > 0 ? span : 0, 0), del_at(span > 0 ? span : 0, 0);
    // alignment gaps (>= 10 bp) and clips, for SV evidence
    std::vector<int> big_del_pos, big_del_len, big_ins_pos, big_ins_len;
    int clip_left = 0, clip_right = 0;
    {
      const uint32_t* cig = bam_get_cigar(b);
      int q = 0, r = (int)rstart;
      for (uint32_t k = 0; k < b->core.n_cigar; ++k) {
        const int op = bam_cigar_op(cig[k]);
        const int len = bam_cigar_oplen(cig[k]);
        switch (op) {
          case BAM_CMATCH: case BAM_CEQUAL: case BAM_CDIFF:
            for (int j = 0; j < len; ++j) {
              q2r[q + j] = r + j;
              r2q[r + j - rstart] = q + j;
            }
            q += len; r += len; break;
          case BAM_CINS:
            if (r - 1 - rstart >= 0 && r - 1 - rstart < span) ins_after[r - 1 - rstart] += len;
            if (len >= 10) { big_ins_pos.push_back(r); big_ins_len.push_back(len); }
            q += len; break;
          case BAM_CSOFT_CLIP: case BAM_CHARD_CLIP:
            // clips before the first aligned base are left clips (H and S
            // can both occur at one end)
            if (r == (int)rstart && q2r_empty(q2r, q)) clip_left += len; else clip_right += len;
            if (op == BAM_CSOFT_CLIP) q += len;
            break;
          case BAM_CDEL:
            if (r - rstart >= 0 && r - rstart < span) del_at[r - rstart] = len;
            if (len >= 10) { big_del_pos.push_back(r); big_del_len.push_back(len); }
            r += len; break;
          case BAM_CREF_SKIP:
            r += len; break;
          default: break;  // pad
        }
      }
    }

    // --- structural variants: evidence from this alignment ----------------------
    // event   : deletions / insertions (>= 10 bp each) inside the breakpoint
    //           window add up to the SV length (+-30 %, at least +-20 bp)
    //           -> carries ALT; a partial sum (>= 30 % of the length) counts
    //           like a clip (unclear, never REF)
    // clipped : a clip >= min_clip ending within sv_window of a breakpoint ->
    //           the read may carry the SV as a split alignment; never REF
    // spans   : the alignment covers the SV with sv_flank on both sides
    for (int si = 0; si < ns; ++si) {
      const int p = sv_pos[si];                       // 1-based anchor
      const int e = (sv_type[si] == 2) ? p : sv_end[si];
      if (rend < p - sv_window - sv_flank || rstart > e + sv_window + sv_flank) continue;
      const int L = std::abs(sv_len[si]);
      const int tol_len = std::max(20, (int)(0.3 * L));
      // sum gaps of the SV's type inside the breakpoint window: ONT often
      // splits one large insertion / deletion into several smaller gaps
      int gap_sum = 0;
      if (sv_type[si] == 1) {
        for (size_t k = 0; k < big_del_pos.size(); ++k)
          if (big_del_pos[k] >= p - sv_window && big_del_pos[k] <= e + sv_window) gap_sum += big_del_len[k];
      } else if (sv_type[si] == 2) {
        for (size_t k = 0; k < big_ins_pos.size(); ++k)
          if (std::abs(big_ins_pos[k] - p) <= sv_window) gap_sum += big_ins_len[k];
      }
      int event = 0;
      int partial = 0;
      if ((sv_type[si] == 1 || sv_type[si] == 2) && L > 0) {
        if (std::abs(gap_sum - L) <= tol_len) event = 1;
        else if (gap_sum >= 0.3 * L) partial = 1;   // some evidence: never REF
      }
      int clipped = partial;
      for (int bp : {p, e}) {
        if (clip_left >= min_clip && std::abs((int)rstart + 1 - bp) <= sv_window) clipped = 1;
        if (clip_right >= min_clip && std::abs((int)rend - bp) <= sv_window) clipped = 1;
      }
      const int spans = (rstart + 1 <= p - sv_flank && rend >= e + sv_flank) ? 1 : 0;
      if (event || clipped || spans) {
        s_read.push_back(ridx + 1); s_sv.push_back(si + 1);
        s_event.push_back(event); s_clipped.push_back(clipped); s_spans.push_back(spans);
      }
    }

    // --- modification calls ---------------------------------------------------
    int ncall = 0;
    const ModCalls mi = parse_mm_ml(b);
    if (!mi.ok) ++n_bad_mm;
    const bool has_mm = mi.ok && (mi.has_m || mi.has_h);
    if (has_mm) {
      const std::vector<int>& qm = mi.qm;
      const std::vector<int>& qh = mi.qh;
      // walk CpG sites of the read in alignment space
      for (int q = 0; q < qlen; ++q) {
        int cpg_pos0;  // 0-based forward-strand position of the CpG C
        if (!rev) {
          if (q + 1 >= qlen) break;
          if (seq_base(seq, q) != 'C' || seq_base(seq, q + 1) != 'G') continue;
          if (q2r[q] < 0 || q2r[q + 1] != q2r[q] + 1) continue;
          cpg_pos0 = q2r[q];
        } else {
          if (q < 1) continue;
          if (seq_base(seq, q) != 'G' || seq_base(seq, q - 1) != 'C') continue;
          if (q2r[q] < 0 || q2r[q - 1] != q2r[q] - 1) continue;
          cpg_pos0 = q2r[q] - 1;
        }
        int m = qm[q], h = qh[q];
        if (m < 0 && mi.has_m && mi.m_implicit) m = 0;
        if (h < 0 && mi.has_h && mi.h_implicit) h = 0;
        if (m < 0 && h < 0) continue;              // not called ('?' mode)
        c_read.push_back(ridx + 1);
        c_pos.push_back(cpg_pos0 + 1);
        c_qm.push_back(m < 0 ? NA_INTEGER : m);
        c_qh.push_back(h < 0 ? NA_INTEGER : h);
        c_strand.push_back(rev ? '-' : '+');
        ++ncall;
      }
    }

    // --- alleles at small variants --------------------------------------------
    {
      const int lo = (int)rstart + 1;          // 1-based inclusive
      const int hi = (int)rend;                // 1-based inclusive
      auto it = std::lower_bound(vpos.begin(), vpos.end(), lo);
      for (int vi = (int)(it - vpos.begin()); vi < nv && vpos[vi] <= hi; ++vi) {
        const std::string& R = vref[vi];
        const std::string& A = valt[vi];
        const int p0 = vpos[vi] - 1;             // 0-based position of REF[0]
        int allele = NA_INTEGER, bq = NA_INTEGER;
        if (R.size() == 1 && A.size() >= 1 && A.find(',') == std::string::npos && A.size() == 1) {
          // SNV
          const int off = p0 - (int)rstart;
          if (off >= 0 && off < span && r2q[off] >= 0) {
            const char base = seq_base(seq, r2q[off]);
            bq = qual[r2q[off]];
            allele = (base == R[0]) ? 0 : (base == A[0]) ? 1 : 2;
          }
        } else if (A.find(',') == std::string::npos && R.size() != A.size() &&
                   (R.size() == 1 || A.size() == 1)) {
          // simple indel, left-anchored: look for an event of the exact
          // length within +-indel_window bp of the anchor
          const bool is_ins = A.size() > R.size();
          const int len = (int)(is_ins ? A.size() - R.size() : R.size() - A.size());
          const int w0 = p0 - indel_window - (int)rstart;
          const int w1 = p0 + (is_ins ? 0 : len) + indel_window - (int)rstart;
          if (w0 >= 0 && w1 < span) {             // read spans the window
            bool found = false;
            for (int off = w0; off <= w1 && !found; ++off) {
              if (is_ins ? ins_after[off] == len : del_at[off] == len) found = true;
            }
            allele = found ? 1 : 0;
          }
        }
        if (allele != NA_INTEGER) {
          a_read.push_back(ridx + 1);
          a_var.push_back(vi + 1);
          a_allele.push_back(allele);
          a_bq.push_back(bq);
        }
      }
    }

    // --- read record -------------------------------------------------------------
    uint8_t* hp = bam_aux_get(b, "HP");
    uint8_t* ps = bam_aux_get(b, "PS");
    r_name.emplace_back(bam_get_qname(b));
    r_flag.push_back(flag);
    r_mapq.push_back(b->core.qual);
    r_start.push_back((int)rstart + 1);
    r_end.push_back((int)rend);
    r_hp.push_back(hp ? (int)bam_aux2i(hp) : NA_INTEGER);
    r_ps.push_back(ps ? (int)bam_aux2i(ps) : NA_INTEGER);
    r_len.push_back(qlen);
    r_ncall.push_back(ncall);
    r_has_mm.push_back(has_mm ? 1 : 0);
  }

  bam_destroy1(b);
  hts_itr_destroy(itr);
  hts_idx_destroy(idx);
  sam_hdr_destroy(hdr);
  sam_close(fp);
  if (ret < -1) stop("error while reading " + path + " (truncated file?)");

  CharacterVector strand(c_strand.size());
  for (size_t i = 0; i < c_strand.size(); ++i) strand[i] = std::string(1, c_strand[i]);

  return List::create(
    _["reads"] = List::create(
      _["read_name"] = wrap(r_name), _["flag"] = wrap(r_flag), _["mapq"] = wrap(r_mapq),
      _["start"] = wrap(r_start), _["end"] = wrap(r_end), _["hp"] = wrap(r_hp),
      _["ps"] = wrap(r_ps), _["read_len"] = wrap(r_len), _["n_calls"] = wrap(r_ncall),
      _["has_mm"] = wrap(r_has_mm)),
    _["calls"] = List::create(
      _["read"] = wrap(c_read), _["pos"] = wrap(c_pos), _["strand"] = strand,
      _["q_m"] = wrap(c_qm), _["q_h"] = wrap(c_qh)),
    _["alleles"] = List::create(
      _["read"] = wrap(a_read), _["var"] = wrap(a_var), _["allele"] = wrap(a_allele),
      _["bq"] = wrap(a_bq)),
    _["sv"] = List::create(
      _["read"] = wrap(s_read), _["sv"] = wrap(s_sv), _["event"] = wrap(s_event),
      _["clipped"] = wrap(s_clipped), _["spans"] = wrap(s_spans)),
    _["qc"] = List::create(
      _["skipped_flag"] = n_skipped_flag, _["skipped_mapq"] = n_skipped_mapq,
      _["bad_mm"] = n_bad_mm));
}
