/*------------------------------------------------------------------------*/
/*! \file fsa_search.cpp
    \brief contains a search-based detection of the final-stage adder

  Part of AMulet2 : AIG Multiplier Verification Tool.
  Copyright(C) 2020, 2021 Daniela Kaufmann, Johannes Kepler University Linz
  Copyright(C) 2026 Ben Jensen (final-stage-adder search)
  Released under the MIT License, see LICENSE.
*/
/*------------------------------------------------------------------------*/
#include <vector>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <stdint.h>
#include <array>
#include "fsa_search.h"
#include "signal_statistics.h"

namespace {

const int NONE = -1;
struct Row { int x, y; };           // literals or NONE
// kind: 0 double, 1 single, 2 double-nocarry, 3 pass, 4 carry-only
struct Opt { int kind; Row r; };
const char * kind_name[] = {"double", "single", "double-nocarry", "pass",
                            "carry-only"};

unsigned n_out;
std::vector<Row> rows;
std::vector<int> kinds;
std::vector<unsigned> region;       // AND variables between rows and outputs
unsigned evals, eval_budget;
uint64_t rng = 0x9e3779b97f4a7c15ull;
uint64_t next_rand() {               // splitmix64
  uint64_t z = (rng += 0x9e3779b97f4a7c15ull);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
  return z ^ (z >> 31);
}
inline unsigned var_of(unsigned lit) { return lit >> 1; }

// lit == x XOR y for an AIG XOR pattern v = AND(!A, !B), A = AND(a0,a1),
// B = AND(!a0,!a1)
bool xor_ins(unsigned lit, unsigned & x, unsigned & y) {
  if (lit < 2) return false;
  aiger_and * v = is_model_and(lit & ~1u);
  if (!v) return false;
  if (!(v->rhs0 & 1) || !(v->rhs1 & 1)) return false;
  aiger_and * A = is_model_and(v->rhs0 & ~1u);
  aiger_and * B = is_model_and(v->rhs1 & ~1u);
  if (!A || !B) return false;
  unsigned a0 = A->rhs0, a1 = A->rhs1, b0 = B->rhs0, b1 = B->rhs1;
  bool ok = (b0 == (a0 ^ 1) && b1 == (a1 ^ 1)) ||
            (b0 == (a1 ^ 1) && b1 == (a0 ^ 1));
  if (!ok) return false;
  x = a0 ^ (lit & 1);
  y = a1;
  return true;
}

// check columns 0..i; on success region holds the cone (sorted)
// marks/values are flat arrays indexed by variable (stamped, so no clearing
// between evaluations)
std::vector<unsigned> row_mark, reg_mark;
std::vector<uint64_t> val;
unsigned stamp = 0;
bool check(unsigned i, const char ** why) {
  evals++;
  unsigned nv = get_model_maxvar() + 1;
  if (row_mark.size() < nv) {
    row_mark.assign(nv, 0);
    reg_mark.assign(nv, 0);
    val.assign(nv, 0);
  }
  stamp++;
  std::vector<unsigned> rowvars;
  for (unsigned k = 0; k <= i; k++) {
    int r2[2] = {rows[k].x, rows[k].y};
    for (int j = 0; j < 2; j++)
      if (r2[j] != NONE && r2[j] >= 2 && row_mark[var_of(r2[j])] != stamp) {
        row_mark[var_of(r2[j])] = stamp; rowvars.push_back(var_of(r2[j]));
      }
  }
  std::vector<unsigned> gs, st;
  for (unsigned k = 0; k <= i; k++) {
    unsigned o = slit(k);
    if (o < 2) { *why = "constant output"; return false; }
    if (row_mark[var_of(o)] != stamp) st.push_back(var_of(o));
  }
  while (!st.empty()) {
    unsigned v = st.back(); st.pop_back();
    if (reg_mark[v] == stamp) continue;
    if (is_model_input(2*v)) { *why = "path to primary input"; return false; }
    aiger_and * a = is_model_and(2*v);
    if (!a) { *why = "unknown signal"; return false; }
    reg_mark[v] = stamp; gs.push_back(v);
    unsigned c[2] = {a->rhs0, a->rhs1};
    for (int j = 0; j < 2; j++) {
      unsigned w = var_of(c[j]);
      if (w == 0 || row_mark[w] == stamp || reg_mark[w] == stamp) continue;
      st.push_back(w);
    }
  }
  // AIGER: lhs > rhs, so ascending = topological
  std::sort(gs.begin(), gs.end());
  val[0] = 0;
  auto get = [&](unsigned lit) -> uint64_t {
    uint64_t r = val[var_of(lit)];
    return (lit & 1) ? ~r : r;
  };
  for (int rep = 0; rep < 4; rep++) {
    for (unsigned v : rowvars) val[v] = next_rand();
    for (unsigned v : gs) {
      aiger_and * a = is_model_and(2*v);
      val[v] = get(a->rhs0) & get(a->rhs1);
    }
    uint64_t c = 0;
    for (unsigned k = 0; k <= i; k++) {
      uint64_t x = rows[k].x != NONE ? get(rows[k].x) : 0;
      uint64_t y = rows[k].y != NONE ? get(rows[k].y) : 0;
      uint64_t s = x ^ y ^ c;
      c = (x & y) | (x & c) | (y & c);
      if (get(slit(k)) != s) { *why = "sum mismatch"; return false; }
    }
  }
  region = gs;
  return true;
}

// An AIG XOR pattern only fixes its inputs up to a joint negation
// (x^y == !x^!y), but the carry of an adder column depends on the true
// polarity, so both polarities are offered.
std::vector<Opt> options(unsigned i, bool carry_possible) {
  std::vector<Opt> o;
  unsigned s = slit(i), u, v;
  if (xor_ins(s, u, v)) {
    unsigned pc[2][2] = {{u, v}, {v, u}};
    for (int j = 0; j < 2; j++) {
      unsigned p0, p1;
      // the propagate gate's own polarity is ambiguous too: 4 sign
      // combinations
      if (xor_ins(pc[j][0], p0, p1))
        for (int sg = 0; sg < 4; sg++)
          o.push_back({0, {static_cast<int>(p0 ^ (sg & 1)),
                           static_cast<int>(p1 ^ (sg >> 1))}});
    }
    if (carry_possible) {
      o.push_back({1, {static_cast<int>(u), NONE}});
      o.push_back({1, {static_cast<int>(u ^ 1), NONE}});
      o.push_back({1, {static_cast<int>(v), NONE}});
      o.push_back({1, {static_cast<int>(v ^ 1), NONE}});
    }
    o.push_back({2, {static_cast<int>(u), static_cast<int>(v)}});
    o.push_back({2, {static_cast<int>(u ^ 1), static_cast<int>(v ^ 1)}});
  }
  o.push_back({3, {static_cast<int>(s), NONE}});
  if (carry_possible) o.push_back({4, {NONE, NONE}});
  return o;
}

bool dfs(unsigned i, bool carry_possible) {
  if (evals > eval_budget) return false;
  if (i == n_out) return true;
  const char * why = 0;
  for (const Opt & op : options(i, carry_possible)) {
    rows[i] = op.r;
    bool okc = check(i, &why);
    if (verbose >= 3)
      msg("fsa-search: column %u try %s rows %d %d -> %s", i,
          kind_name[op.kind], op.r.x, op.r.y, okc ? "ok" : why);
    if (!okc) continue;
    kinds[i] = op.kind;
    bool cp = carry_possible || op.kind == 0 || op.kind == 2;
    if (dfs(i+1, cp)) return true;
    if (evals > eval_budget) return false;
  }
  rows[i] = {NONE, NONE};
  return false;
}

// --- construction of the ripple-carry adder over the rows ---------------
unsigned next_var;
// lhs, rhs0, rhs1 (shared by miter and rewritten)
std::vector<std::array<unsigned, 3>> new_ands;
// structural hashing: a ripple/miter AND that already exists (in the adder
// region or among the new gates) is reused, so shared propagate/generate
// gates are merged before the SAT call
std::unordered_map<uint64_t, unsigned> strash;
bool use_strash = true;
inline uint64_t skey(unsigned a, unsigned b) {
  if (a > b) std::swap(a, b);
  return (static_cast<uint64_t>(a) << 32) | b;
}
unsigned mk_and(unsigned a, unsigned b) {
  if (a == 0 || b == 0 || a == (b ^ 1)) return 0;
  if (a == 1 || a == b) return b;
  if (b == 1) return a;
  if (use_strash) {
    auto it = strash.find(skey(a, b));
    if (it != strash.end()) return it->second;
  }
  unsigned l = 2*(next_var++);
  new_ands.push_back({l, a, b});
  if (use_strash) strash[skey(a, b)] = l;
  return l;
}
unsigned mk_or(unsigned a, unsigned b) { return mk_and(a ^ 1, b ^ 1) ^ 1; }
// 3 ANDs, same shape as btor_ha
unsigned mk_xor(unsigned a, unsigned b) {
  if (a == 0) return b;
  if (b == 0) return a;
  if (a == 1) return b ^ 1;
  if (b == 1) return a ^ 1;
  unsigned one = mk_and(a ^ 1, b ^ 1), two = mk_and(a, b);
  return mk_and(one ^ 1, two ^ 1);
}
std::vector<unsigned> ripple_out;

void build_ripple() {
  new_ands.clear(); ripple_out.clear();
  next_var = get_model_maxvar() + 1;
  strash.clear();
  use_strash = !getenv("AMULET_NO_STRASH");
  if (use_strash)
    for (unsigned v : region) {
      aiger_and * a = is_model_and(2*v);
      strash[skey(a->rhs0, a->rhs1)] = a->lhs;
    }
  unsigned c = 0;
  for (unsigned k = 0; k < n_out; k++) {
    unsigned x = rows[k].x != NONE ? static_cast<unsigned>(rows[k].x) : 0;
    unsigned y = rows[k].y != NONE ? static_cast<unsigned>(rows[k].y) : 0;
    unsigned t = mk_xor(x, y);
    ripple_out.push_back(mk_xor(t, c));
    if (k + 1 < n_out) c = mk_or(mk_and(x, y), mk_and(t, c));
  }
}

std::vector<unsigned> row_vars_sorted() {
  std::unordered_set<unsigned> rv;
  for (const Row & r : rows) {
    if (r.x != NONE && r.x >= 2) rv.insert(var_of(r.x));
    if (r.y != NONE && r.y >= 2) rv.insert(var_of(r.y));
  }
  std::vector<unsigned> v(rv.begin(), rv.end());
  std::sort(v.begin(), v.end());
  return v;
}

std::vector<unsigned> miter_vars_in;   // row variables (miter inputs)
std::vector<std::array<unsigned, 3>> miter_ands;
unsigned miter_out;

}  // namespace

bool fsa_search_substitute() {
  n_out = get_model_num_outputs();
  rows.assign(n_out, {NONE, NONE});
  kinds.assign(n_out, -1);
  evals = 0; eval_budget = 200 * n_out;
  msg("fsa-search: searching final stage adder rows (%u columns)", n_out);
  bool ok = dfs(0, false);
  if (!ok) {
    msg("substitution not possible - fsa-search: %s (%u evaluations)",
        evals > eval_budget ? "evaluation budget exhausted"
                            : "no consistent adder hypothesis", evals);
    return 0;
  }
  const char * why = 0;
  if (!check(n_out - 1, &why)) {
    msg("substitution not possible - fsa-search: "
        "hypothesis failed re-simulation (%s)", why);
    return 0;
  }
  unsigned doubles = 0;
  std::string shape;
  for (unsigned k = 0; k < n_out; k++) {
    if (kinds[k] == 0 || kinds[k] == 2) doubles++;
    shape += static_cast<char>('0' + (rows[k].x != NONE) + (rows[k].y != NONE));
  }
  if (!doubles) {
    msg("substitution not possible - fsa-search: no two-operand column");
    return 0;
  }
  msg("fsa-search: found adder, rows per column (LSB first) %s, "
      "%u region gates, %u evaluations",
      shape.c_str(), static_cast<unsigned>(region.size()), evals);
  if (verbose >= 2)
    for (unsigned k = 0; k < n_out; k++)
      msg("fsa-search: column %u kind %s rows %d %d", k, kind_name[kinds[k]],
          rows[k].x, rows[k].y);

  build_ripple();

  // miter: inputs = row variables; region gates; ripple gates;
  // out = OR_k (s_k XOR r_k)
  miter_vars_in = row_vars_sorted();
  unsigned acc = 0;   // OR of differences
  for (unsigned k = 0; k < n_out; k++) {
    unsigned d = mk_xor(slit(k), ripple_out[k]);
    acc = mk_or(acc, d);
  }
  // new_ands now also holds the difference/OR gates built just above
  miter_ands.clear();
  for (unsigned v : region) {
    aiger_and * a = is_model_and(2*v);
    miter_ands.push_back({a->lhs, a->rhs0, a->rhs1});
  }
  for (auto & g : new_ands) miter_ands.push_back(g);
  miter_out = acc;
  return 1;
}

bool fsa_search_miter_to_file(FILE * file) {
  if (!file) return 0;
  // map AIG variables to DIMACS variables; constant FALSE gets its own variable
  std::unordered_map<unsigned, int> dv;
  int m = 0;
  auto mapv = [&](unsigned v) -> int {
    auto it = dv.find(v);
    if (it != dv.end()) return it->second;
    dv[v] = ++m;
    return m;
  };
  auto L = [&](unsigned lit) -> int {
    int x = mapv(var_of(lit));
    return (lit & 1) ? -x : x;
  };
  std::vector<std::vector<int>> cls;
  cls.push_back({-mapv(0)});                       // constant false
  for (unsigned v : miter_vars_in) mapv(v);
  // XOR-aware Tseitin: v = AND(!A,!B), A = AND(x,y), B = AND(!x,!y), A and B
  // used only by v -> v <-> x XOR y as 4 clauses, A and B not encoded
  // (AMULET_NO_XORCNF=1 disables)
  std::unordered_map<unsigned, const std::array<unsigned, 3> *> gl;
  std::unordered_map<unsigned, unsigned> fo;
  for (auto & g : miter_ands) {
    gl[g[0]] = &g;
    fo[g[1] & ~1u]++;
    fo[g[2] & ~1u]++;
  }
  fo[miter_out & ~1u]++;
  std::unordered_set<unsigned> skip;
  std::vector<std::array<unsigned, 3>> xors;   // v, x, y
  if (!getenv("AMULET_NO_XORCNF"))
    for (auto & g : miter_ands) {
      if (!(g[1] & 1) || !(g[2] & 1)) continue;
      auto ia = gl.find(g[1] & ~1u), ib = gl.find(g[2] & ~1u);
      if (ia == gl.end() || ib == gl.end()) continue;
      const auto & A = *ia->second; const auto & B = *ib->second;
      bool ok = (B[1] == (A[1] ^ 1) && B[2] == (A[2] ^ 1)) ||
                (B[1] == (A[2] ^ 1) && B[2] == (A[1] ^ 1));
      if (!ok || fo[A[0]] != 1 || fo[B[0]] != 1 || A[0] == B[0]) continue;
      if (skip.count(A[0]) || skip.count(B[0])) continue;
      skip.insert(A[0]); skip.insert(B[0]);
      xors.push_back({g[0], A[1], A[2]});
    }
  std::unordered_set<unsigned> isx;
  for (auto & x : xors) isx.insert(x[0]);
  for (auto & g : miter_ands) {
    if (skip.count(g[0]) || isx.count(g[0])) continue;
    int l = L(g[0]), a = L(g[1]), b = L(g[2]);
    cls.push_back({-l, a});
    cls.push_back({-l, b});
    cls.push_back({l, -a, -b});
  }
  for (auto & x : xors) {
    int v = L(x[0]), a = L(x[1]), b = L(x[2]);
    cls.push_back({-v, a, b});
    cls.push_back({-v, -a, -b});
    cls.push_back({v, -a, b});
    cls.push_back({v, a, -b});
  }
  cls.push_back({L(miter_out)});
  // the header counts what is written
  fprintf(file, "p cnf %d %u\n", m, static_cast<unsigned>(cls.size()));
  for (auto & c : cls) {
    for (int x : c) fprintf(file, "%d ", x);
    fprintf(file, "0\n");
  }
  msg("fsa-search: adder miter with %u row inputs, %u gates, %d vars, "
      "%u clauses", static_cast<unsigned>(miter_vars_in.size()),
      static_cast<unsigned>(miter_ands.size()), m,
      static_cast<unsigned>(cls.size()));
  return 1;
}

bool fsa_search_write_rewritten(FILE * file) {
  // rewritten AIG = original inputs + original ANDs in the cones of the row
  // variables + ripple adder
  for (unsigned i = 0; i < get_model_num_inputs(); i++)
    aiger_add_input(rewritten, get_model_inputs_lit(i),
                    get_model_inputs_name(i));
  std::unordered_set<unsigned> keep;
  std::vector<unsigned> st = row_vars_sorted();
  while (!st.empty()) {
    unsigned v = st.back(); st.pop_back();
    if (keep.count(v)) continue;
    aiger_and * a = is_model_and(2*v);
    if (!a) continue;                 // primary input
    keep.insert(v);
    if (var_of(a->rhs0)) st.push_back(var_of(a->rhs0));
    if (var_of(a->rhs1)) st.push_back(var_of(a->rhs1));
  }
  std::vector<unsigned> kv(keep.begin(), keep.end());
  std::sort(kv.begin(), kv.end());
  for (unsigned v : kv) {
    aiger_and * a = is_model_and(2*v);
    aiger_add_and(rewritten, a->lhs, a->rhs0, a->rhs1);
  }
  // ripple gates only (the diff gates of the miter are not needed): they
  // precede the diff gates
  std::unordered_set<unsigned> need;
  std::vector<unsigned> st2(ripple_out.begin(), ripple_out.end());
  std::unordered_map<unsigned, std::array<unsigned, 3>> byl;
  for (auto & g : new_ands) byl[g[0]] = g;
  // region gates reused by the ripple adder through structural hashing
  std::vector<unsigned> regneed;
  while (!st2.empty()) {
    unsigned l = st2.back() & ~1u; st2.pop_back();
    if (l < 2 || need.count(l) || keep.count(var_of(l))) continue;
    if (!byl.count(l)) {
      aiger_and * a = is_model_and(l);
      if (!a) continue;               // primary input
      need.insert(l); regneed.push_back(var_of(l));
      st2.push_back(a->rhs0); st2.push_back(a->rhs1);
      continue;
    }
    need.insert(l);
    st2.push_back(byl[l][1]); st2.push_back(byl[l][2]);
  }
  std::sort(regneed.begin(), regneed.end());
  for (unsigned v : regneed) {
    aiger_and * a = is_model_and(2*v);
    aiger_add_and(rewritten, a->lhs, a->rhs0, a->rhs1);
  }
  for (auto & g : new_ands)
    if (need.count(g[0])) aiger_add_and(rewritten, g[0], g[1], g[2]);
  char buf[32];
  for (unsigned k = 0; k < n_out; k++) {
    snprintf(buf, sizeof buf, "O%u", k);
    aiger_add_output(rewritten, ripple_out[k], buf);
  }
  const char * err = aiger_check(rewritten);
  if (err) { msg("fsa-search: rewritten aig check failed: %s", err); return 0; }
  aiger_reencode(rewritten);
  return aiger_write_to_file(rewritten, aiger_binary_mode, file);
}
