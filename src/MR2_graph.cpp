// MR2_graph.cpp
// Rcpp/RcppArmadillo translation of MR2's chordal-graph MCMC engine.
//
// Profiling MR2() showed >50% of total runtime inside check_chordal(), the ripclique/
// junction-tree construction, and ln_h_ratio_zo() -- all called graphParNiter times PER MR2
// MCMC iteration. These are small-vector, many-tiny-function-call operations, i.e. exactly
// the pattern where R's per-call interpreter overhead dominates over actual FLOPs, and where
// Rcpp gives the largest win: not just faster inner loops, but eliminating the R<->C++
// boundary crossing for the whole graphParNiter loop (generate_graph_cpp runs it entirely in
// C++ and returns only the final graph + acceptance rate).
//
// All vertex sets are represented as 0/1 (n_responses x 1) arma::vec ("zero-one" vectors, zo),
// mirroring MR2's original R representation, to keep the translation mechanical and checkable
// line-by-line against the R source. Vertex *labels* are 0-indexed internally (R is 1-indexed);
// this file performs the conversion at the R/C++ boundary only.
//
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// ---------------------------------------------------------------------------------------
// Small zero-one vector helpers (union_zo, intersect_zo)
// ---------------------------------------------------------------------------------------

static inline vec intersect_zo_cpp(const vec &a, const vec &b) { return a % b; }

static inline vec union_zo_cpp(const vec &a, const vec &b)
{
  vec s = a + b;
  s.for_each([](vec::elem_type &x) { if (x > 1.0) x = 1.0; });
  return s;
}

// ---------------------------------------------------------------------------------------
// check_chordal: maximum cardinality search + clique test (Tarjan-Yannakakis)
// ---------------------------------------------------------------------------------------

struct ChordalResult { int chordal; ivec order; };

static ChordalResult check_chordal_cpp(const mat &g_in)
{
  mat g = g_in;
  g.diag().ones();
  int p = g.n_rows;

  ivec order(p, fill::zeros);
  int chordal = 1;

  std::vector<int> numbered;
  numbered.reserve(p);
  std::vector<bool> is_numbered(p, false);

  numbered.push_back(0);
  is_numbered[0] = true;
  order(0) = 0;

  for (int i = 1; i < p; i++)
  {
    int best_u = -1, best_score = -1;
    for (int v = 0; v < p; v++)
    {
      if (is_numbered[v]) continue;
      int score = 0;
      for (int nb : numbered) if (g(v, nb) != 0.0) score++;
      if (score > best_score) { best_score = score; best_u = v; }
    }
    int u = best_u;
    numbered.push_back(u);
    is_numbered[u] = true;
    order(i) = u;

    std::vector<int> pa;
    for (int t = 0; t < i; t++) if (g(u, order(t)) != 0.0) pa.push_back(order(t));

    bool all_ones = true;
    for (size_t a1 = 0; a1 < pa.size() && all_ones; a1++)
      for (size_t a2 = 0; a2 < pa.size() && all_ones; a2++)
        if (g(pa[a1], pa[a2]) != 1.0) all_ones = false;

    if (!all_ones) { chordal = 0; break; }
  }
  return ChordalResult{chordal, order};
}

// ---------------------------------------------------------------------------------------
// chordal_to_ripcliques: perfect-elimination-ordering to running-intersection cliques
// ---------------------------------------------------------------------------------------

static mat chordal_to_ripcliques_cpp(const mat &G, const ivec &order)
{
  int p = G.n_rows;
  mat pa(p, p, fill::zeros);

  for (int i = 1; i < p; i++)
  {
    int v = order(i);
    vec pre_v(p, fill::zeros);
    for (int t = 0; t < i; t++) pre_v(order(t)) = 1.0;
    vec pa_i = intersect_zo_cpp(G.col(v), pre_v);
    pa.col(v) = pa_i;
  }

  ivec ladder(p); ladder.fill(-1);
  for (int i = 0; i < p; i++)
  {
    bool take = (i == p - 1) ? true : (accu(pa.col(order(i))) >= accu(pa.col(order(i + 1))));
    if (take) ladder(i) = order(i);
  }

  mat cliques(p, p, fill::zeros);
  int idx = -1;
  for (int i = 0; i < p; i++)
  {
    if (ladder(i) != -1)
    {
      idx++;
      vec v_col(p, fill::zeros);
      v_col(ladder(i)) = 1.0;
      cliques.col(idx) = union_zo_cpp(v_col, pa.col(ladder(i)));
    }
  }
  return cliques;
}

static inline int num_cliques_of(const mat &cliques)
{
  rowvec sizes = sum(cliques, 0);
  int n = 0;
  for (uword c = 0; c < sizes.n_elem; c++) if (sizes(c) != 0.0) n++;
  return n;
}

// ---------------------------------------------------------------------------------------
// ripcliques_to_jtree / sepsize / reachability
// ---------------------------------------------------------------------------------------

static mat ripcliques_to_jtree_cpp(const mat &cliques)
{
  int p = cliques.n_rows;
  int num_cl = num_cliques_of(cliques);
  mat jtree(p, p, fill::zeros);

  if (num_cl >= 2)
  {
    vec score(p, fill::zeros);
    for (int i = 1; i < num_cl; i++)
    {
      for (int k = 0; k < i; k++) score(k) = accu(intersect_zo_cpp(cliques.col(i), cliques.col(k)));
      double mx = 0.0;
      for (int k = 0; k < i; k++) if (score(k) > mx) mx = score(k);
      if (mx != 0.0)
      {
        int j = -1;
        for (int k = 0; k < i; k++) if (score(k) == mx) { j = k; break; }
        jtree(i, j) = 1.0; jtree(j, i) = 1.0;
      }
    }
  }
  return jtree;
}

static mat sepsize_cpp(const mat &cliques, const mat &jtree)
{
  int p = cliques.n_rows;
  int num_cl = num_cliques_of(cliques);
  mat sepsize(p, p, fill::zeros);

  for (int i = 0; i < num_cl; i++)
    for (int j = i + 1; j < num_cl; j++)
      if (jtree(i, j) == 1.0)
      {
        double s = accu(intersect_zo_cpp(cliques.col(i), cliques.col(j)));
        sepsize(i, j) = s; sepsize(j, i) = s;
      }
  return sepsize;
}

static mat reachability_graph_cpp(const mat &g)
{
  int p = g.n_rows;
  mat A = g;
  mat reach(p, p, fill::zeros);
  for (int i = 0; i < p - 1; i++) { reach += A; A = A * g; }
  reach.for_each([](mat::elem_type &x) { x = (x > 0.0) ? 1.0 : 0.0; });
  return reach;
}

// ---------------------------------------------------------------------------------------
// find_clique_containing / parents_vertex
// ---------------------------------------------------------------------------------------

static inline int find_clique_containing_cpp(int a, const mat &cliques)
{
  for (uword c = 0; c < cliques.n_cols; c++) if (cliques(a, c) != 0.0) return (int)c;
  return -1;
}

static inline vec parents_vertex_cpp(const mat &g, int i)
{
  vec ps = g.col(i);
  for (uword t = i; t < ps.n_elem; t++) ps(t) = 0.0;
  return ps;
}

static inline int which_nonzero_single(const vec &v)
{
  // faithful to R's `which(x != 0)` when the caller guarantees at most one nonzero entry
  for (uword t = 0; t < v.n_elem; t++) if (v(t) != 0.0) return (int)t;
  return -1;   // "integer(0)" in R -> sum()==0, handled by callers via -1 sentinel
}

// ---------------------------------------------------------------------------------------
// check_edge_delete_zo
// ---------------------------------------------------------------------------------------

struct EdgeDeleteResult { int delete_ok; vec C; };

static EdgeDeleteResult check_edge_delete_cpp(int i, int j, const mat &cliques)
{
  int p = cliques.n_rows;
  int delete_ok = 1;
  vec C(p, fill::zeros);
  int count_cliques = 0;
  int num_cl = num_cliques_of(cliques);

  vec edge_col_vec(p, fill::zeros);
  edge_col_vec(i) = 1.0; edge_col_vec(j) = 1.0;

  for (int k = 0; k < num_cl; k++)
  {
    vec clique_k = cliques.col(k);
    if (accu(intersect_zo_cpp(edge_col_vec, clique_k)) == 2.0)
    {
      C = clique_k;
      count_cliques++;
      if (count_cliques > 1) { delete_ok = 0; C = vec(3, fill::value(9.0)); break; }   // matches R's rep(9,3) sentinel exactly (unused downstream when delete_ok==0)
    }
  }
  return EdgeDeleteResult{delete_ok, C};
}

// ---------------------------------------------------------------------------------------
// check_edge_add_same_component_zo -- direct, careful 0-indexed port of the R original
// ---------------------------------------------------------------------------------------

struct EdgeAddResult { int edge_ok; vec C; bool has_C; };

static EdgeAddResult check_edge_add_same_component_cpp(int a, int b, const mat &jtree,
                                                        const mat &sepsize, const mat &cliques)
{
  int p = cliques.n_rows;
  int edge_ok = 0, quit = 0;
  int CASE_same_branch = 0, CASE_sats_on_b2_branch = 0, locate_a2 = 0;

  vec ab_edge_vec(p, fill::zeros);
  ab_edge_vec(a) = 1.0; ab_edge_vec(b) = 1.0;

  int index_b = find_clique_containing_cpp(b, cliques);
  int index_a = find_clique_containing_cpp(a, cliques);
  int index_b2 = std::max(index_a, index_b);
  int a2 = (index_b2 == index_b) ? a : b;

  int index_a2 = -1;
  int index = index_b2;
  int next_index = which_nonzero_single(parents_vertex_cpp(jtree, index));

  if (index_b2 - 1 >= 1)
  {
    for (int dummy = 1; dummy <= index_b2 - 1; dummy++)
    {
      if (next_index == -1 && locate_a2 == 0) { CASE_same_branch = 0; break; }
      if (next_index == -1) break;
      if (cliques(a2, next_index) == 1.0)
      { index_a2 = next_index; locate_a2 = 1; CASE_same_branch = 1; break; }
      index = next_index;
      next_index = which_nonzero_single(parents_vertex_cpp(jtree, index));
    }
  }
  if (CASE_same_branch == 0) index_a2 = find_clique_containing_cpp(a2, cliques);

  vec clique_a2_int_clique_b2 = intersect_zo_cpp(cliques.col(index_a2), cliques.col(index_b2));
  double s = accu(clique_a2_int_clique_b2);

  vec C(p, fill::zeros);
  bool has_C = false;

  if (s == 0.0) { edge_ok = 0; quit = 1; }
  if (jtree(index_a2, index_b2) == 1.0)
  { edge_ok = 1; C = union_zo_cpp(ab_edge_vec, clique_a2_int_clique_b2); has_C = true; quit = 1; }

  if (quit == 0)
  {
    if (CASE_same_branch == 1)
    {
      int bottom = index_b2;
      int next_parent_b2 = which_nonzero_single(parents_vertex_cpp(jtree, index_b2));
      while (next_parent_b2 != -1 && next_parent_b2 >= index_a2)
      {
        if (sepsize(next_parent_b2, bottom) == s)
        { edge_ok = 1; C = union_zo_cpp(ab_edge_vec, clique_a2_int_clique_b2); has_C = true; quit = 1; break; }
        bottom = next_parent_b2;
        next_parent_b2 = which_nonzero_single(parents_vertex_cpp(jtree, bottom));
      }
    }
    else
    {
      std::vector<int> ancestors_a2, ancestors_b2;
      vec ancestors_a2_col_vec(p, fill::zeros), ancestors_b2_col_vec(p, fill::zeros);

      int next_ancestor_a2 = index_a2;
      for (int count = 0; count < index_a2 + 1 && next_ancestor_a2 != -1; count++)
      {
        ancestors_a2.push_back(next_ancestor_a2);
        ancestors_a2_col_vec(next_ancestor_a2) = 1.0;
        next_ancestor_a2 = which_nonzero_single(parents_vertex_cpp(jtree, next_ancestor_a2));
      }
      int next_ancestor_b2 = index_b2;
      for (int count = 0; count < index_b2 + 1 && next_ancestor_b2 != -1; count++)
      {
        ancestors_b2.push_back(next_ancestor_b2);
        ancestors_b2_col_vec(next_ancestor_b2) = 1.0;
        next_ancestor_b2 = which_nonzero_single(parents_vertex_cpp(jtree, next_ancestor_b2));
      }

      vec fork_set = intersect_zo_cpp(ancestors_a2_col_vec, ancestors_b2_col_vec);
      int fork = -1;
      for (int t = p - 1; t >= 0; t--) if (fork_set(t) != 0.0) { fork = t; break; }

      std::vector<int> branch_b2_gt_fork;
      for (int v : ancestors_b2) if (v > fork) branch_b2_gt_fork.push_back(v);

      if (!branch_b2_gt_fork.empty())
      {
        for (size_t count = 0; count + 1 < ancestors_b2.size() && count < branch_b2_gt_fork.size(); count++)
        {
          int index_row = ancestors_b2[count], index_col = ancestors_b2[count + 1];
          if (sepsize(index_row, index_col) == s)
          {
            edge_ok = 1; C = union_zo_cpp(ab_edge_vec, clique_a2_int_clique_b2); has_C = true;
            quit = 1; CASE_sats_on_b2_branch = 1; break;
          }
        }
      }
      if (quit == 0 && CASE_sats_on_b2_branch == 0)
      {
        std::vector<int> branch_a2_gt_fork;
        for (int v : ancestors_a2) if (v > fork) branch_a2_gt_fork.push_back(v);
        if (!branch_a2_gt_fork.empty())
        {
          for (size_t count = 0; count + 1 < ancestors_a2.size() && count < branch_a2_gt_fork.size(); count++)
          {
            int index_row = ancestors_a2[count + 1], index_col = ancestors_a2[count];
            if (sepsize(index_row, index_col) == s)
            {
              edge_ok = 1; C = union_zo_cpp(ab_edge_vec, clique_a2_int_clique_b2); has_C = true;
              quit = 1; break;
            }
          }
        }
      }
    }
  }
  if (quit == 0) { edge_ok = 0; has_C = false; }
  return EdgeAddResult{edge_ok, C, has_C};
}

// ---------------------------------------------------------------------------------------
// next_edge_candidate / next_graph_candidate_zo
// ---------------------------------------------------------------------------------------

static inline std::pair<int,int> next_edge_candidate_cpp(int p)
{
  IntegerVector idx = Rcpp::sample(p, 2, false);   // 1-indexed, distinct
  return {idx[0] - 1, idx[1] - 1};
}

struct GraphCandidateResult { mat g_proposal; int a, b; vec C; bool has_C; int CASE_add, CASE_delete; };

static GraphCandidateResult next_graph_candidate_cpp(const mat &g_current, const mat &jtree,
                                                     const mat &sepsize, const mat &cliques,
                                                     const mat &reach_graph)
{
  int p = g_current.n_rows;
  mat g_proposal = g_current;
  int CASE_add = 0, CASE_delete = 0;
  int a = -1, b = -1;
  vec C; bool has_C = false;

  while (approx_equal(g_proposal, g_current, "absdiff", 0.0))
  {
    auto edge = next_edge_candidate_cpp(p);
    int i = edge.first, j = edge.second;

    if (g_current(i, j) == 1.0)
    {
      EdgeDeleteResult ss = check_edge_delete_cpp(i, j, cliques);
      CASE_delete = ss.delete_ok;
      if (CASE_delete == 1)
      { g_proposal(i, j) = 0.0; g_proposal(j, i) = 0.0; a = i; b = j; C = ss.C; has_C = true; }
    }
    else
    {
      if (reach_graph(i, j) == 0.0)
      {
        CASE_add = 1; a = i; b = j;
        C = vec(p, fill::zeros); C(a) = 1.0; C(b) = 1.0; has_C = true;
        g_proposal(i, j) = 1.0; g_proposal(j, i) = 1.0;
      }
      else
      {
        EdgeAddResult ss = check_edge_add_same_component_cpp(i, j, jtree, sepsize, cliques);
        CASE_add = ss.edge_ok;
        if (CASE_add == 1)
        { a = i; b = j; C = ss.C; has_C = ss.has_C; g_proposal(i, j) = 1.0; g_proposal(j, i) = 1.0; }
      }
    }
  }
  return GraphCandidateResult{g_proposal, a, b, C, has_C, CASE_add, CASE_delete};
}

// ---------------------------------------------------------------------------------------
// ln_h_ratio_zo
// ---------------------------------------------------------------------------------------

static double ln_h_ratio_cpp(const vec &C, int a, int b, double delta, const mat &Phi)
{
  int p = C.n_elem;
  std::vector<int> C_nodes;
  for (int t = 0; t < p; t++) if (C(t) != 0.0) C_nodes.push_back(t);

  std::vector<int> Sq2;
  for (int v : C_nodes) if (v != a && v != b) Sq2.push_back(v);
  int numSq2 = (int)Sq2.size();
  double delta_star = delta - 1.0;

  std::vector<int> idxAll = Sq2;
  idxAll.push_back(a); idxAll.push_back(b);
  int total = (int)idxAll.size();

  mat Phi_CC(total, total);
  for (int r = 0; r < total; r++)
    for (int c = 0; c < total; c++)
      Phi_CC(r, c) = Phi(idxAll[r], idxAll[c]);

  mat U = chol(Phi_CC);          // R's chol() is upper-triangular by default
  mat L = trans(U);              // t(chol(Phi_CC)) in R

  double l_aa = L(numSq2, numSq2);
  double l_bb = L(numSq2 + 1, numSq2 + 1);
  double l_ab = L(numSq2 + 1, numSq2);

  double ln_top = (delta_star + numSq2 + 2) * (std::log(l_aa) + std::log(l_bb));
  double ln_bottom = (delta_star + numSq2 + 1) / 2.0 *
                      std::log(l_aa * l_aa * l_ab * l_ab + l_aa * l_aa * l_bb * l_bb);
  double ln_gam_bit = -std::log(2.0 * std::sqrt(M_PI)) +
                       R::lgammafn((delta_star + numSq2 + 1) / 2.0) -
                       R::lgammafn((delta_star + numSq2 + 2) / 2.0);

  return ln_top - ln_bottom + ln_gam_bit;
}

// ---------------------------------------------------------------------------------------
// generate_graph_cpp: the whole graphParNiter loop, run entirely in C++
// ---------------------------------------------------------------------------------------

// [[Rcpp::export]]
List generate_graph_cpp(NumericMatrix g_current_r, double delta, NumericMatrix Phi_r,
                        double N, NumericMatrix S_N_r, int warmup, int iters_to_keep)
{
  mat g_next = as<mat>(g_current_r);
  mat Phi    = as<mat>(Phi_r);
  mat S_N    = as<mat>(S_N_r);

  if (iters_to_keep <= warmup) stop("not enough iterations specified");

  int accept_count_sample = 0;

  for (int it = 1; it <= iters_to_keep; it++)
  {
    mat g_curr = g_next;

    ChordalResult cc = check_chordal_cpp(g_curr);
    if (cc.chordal == 0) stop("graph not decomposable (should not happen during MCMC exploration)");

    mat cliques = chordal_to_ripcliques_cpp(g_curr, cc.order);
    mat jtree   = ripcliques_to_jtree_cpp(cliques);
    mat sepsz   = sepsize_cpp(cliques, jtree);
    mat reach   = reachability_graph_cpp(g_curr);

    double u = R::unif_rand();

    GraphCandidateResult out = next_graph_candidate_cpp(g_curr, jtree, sepsz, cliques, reach);

    double likelihood_ratio = 0.0;
    if (out.CASE_add == 1)
    {
      double ln_lr = ln_h_ratio_cpp(out.C, out.a, out.b, delta, Phi) -
                     ln_h_ratio_cpp(out.C, out.a, out.b, delta + N, Phi + S_N);
      likelihood_ratio = std::exp(ln_lr);
    }
    else if (out.CASE_delete == 1)
    {
      double ln_lr = ln_h_ratio_cpp(out.C, out.a, out.b, delta + N, Phi + S_N) -
                     ln_h_ratio_cpp(out.C, out.a, out.b, delta, Phi);
      likelihood_ratio = std::exp(ln_lr);
    }

    double accept_prob = std::min(1.0, likelihood_ratio);
    if (u <= accept_prob)
    {
      g_next = out.g_proposal;
      if (it > warmup) accept_count_sample++;
    }
    else
    {
      g_next = g_curr;
    }
  }

  return List::create(Named("g_next") = wrap(g_next),
                      Named("accept_rate") = (double)accept_count_sample / iters_to_keep);
}

// ---------------------------------------------------------------------------------------
// Exported wrappers for the individual pieces, used only for R-vs-C++ correctness testing
// ---------------------------------------------------------------------------------------

// [[Rcpp::export]]
List check_chordal_cpp_export(NumericMatrix g_r)
{
  ChordalResult cc = check_chordal_cpp(as<mat>(g_r));
  IntegerVector order_r(cc.order.n_elem);
  for (uword t = 0; t < cc.order.n_elem; t++) order_r[t] = cc.order(t) + 1;   // back to 1-indexed
  return List::create(Named("chordal") = cc.chordal, Named("order") = order_r);
}

// [[Rcpp::export]]
NumericMatrix chordal_to_ripcliques_cpp_export(NumericMatrix G_r, IntegerVector order_r)
{
  ivec order(order_r.size());
  for (int t = 0; t < order_r.size(); t++) order(t) = order_r[t] - 1;
  mat cliques = chordal_to_ripcliques_cpp(as<mat>(G_r), order);
  return wrap(cliques);
}

// [[Rcpp::export]]
NumericMatrix ripcliques_to_jtree_cpp_export(NumericMatrix cliques_r)
{
  return wrap(ripcliques_to_jtree_cpp(as<mat>(cliques_r)));
}

// [[Rcpp::export]]
NumericMatrix sepsize_cpp_export(NumericMatrix cliques_r, NumericMatrix jtree_r)
{
  return wrap(sepsize_cpp(as<mat>(cliques_r), as<mat>(jtree_r)));
}

// [[Rcpp::export]]
NumericMatrix reachability_graph_cpp_export(NumericMatrix g_r)
{
  return wrap(reachability_graph_cpp(as<mat>(g_r)));
}

// [[Rcpp::export]]
double ln_h_ratio_cpp_export(NumericVector C_r, int a, int b, double delta, NumericMatrix Phi_r)
{
  vec C = as<vec>(C_r);
  return ln_h_ratio_cpp(C, a - 1, b - 1, delta, as<mat>(Phi_r));
}

// [[Rcpp::export]]
List check_edge_delete_cpp_export(int i, int j, NumericMatrix cliques_r)
{
  EdgeDeleteResult r = check_edge_delete_cpp(i - 1, j - 1, as<mat>(cliques_r));
  return List::create(Named("delete_ok") = r.delete_ok, Named("C") = wrap(r.C));
}

// [[Rcpp::export]]
List check_edge_add_same_component_cpp_export(int a, int b, NumericMatrix jtree_r,
                                              NumericMatrix sepsize_r, NumericMatrix cliques_r)
{
  EdgeAddResult r = check_edge_add_same_component_cpp(a - 1, b - 1, as<mat>(jtree_r),
                                                      as<mat>(sepsize_r), as<mat>(cliques_r));
  if (r.has_C) return List::create(Named("edge_ok") = r.edge_ok, Named("C") = wrap(r.C));
  return List::create(Named("edge_ok") = r.edge_ok, Named("C") = R_NilValue);
}
