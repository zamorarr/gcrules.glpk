#include <Rcpp.h>
#include "glpk.h"

void add_max_support_vars(glp_prob* mip, int n, int d, int R_max, double ws) {
  int num_vars = n + d + 2;
  glp_add_cols(mip, num_vars);

  // alpha binary variables
  std::string varname;
  for (int i = 1; i <= d; i++) {
    varname = "alpha" + std::to_string(i);
    glp_set_col_name(mip, i, varname.c_str());
    glp_set_col_kind(mip, i, GLP_BV);
  }

  for (int i = 1; i <= n; i++) {
    varname = "r" + std::to_string(i);
    glp_set_col_name(mip, d + i, varname.c_str());
    glp_set_col_kind(mip, d + i, GLP_BV);
  }

  // R (total columns selected)
  glp_set_col_name(mip, d + n + 1, "R");
  glp_set_col_kind(mip, d + n + 1, GLP_IV);
  glp_set_col_bnds(mip, d + n + 1, GLP_DB, 0, R_max);

  // S (support - number of rows matching rule)
  glp_set_col_name(mip, d + n + 2, "S");
  glp_set_col_kind(mip, d + n + 2, GLP_IV);
  glp_set_col_bnds(mip, d + n + 2, GLP_DB, 0, n);

  // objective
  glp_set_obj_coef(mip, d + n + 1, -1);
  glp_set_obj_coef(mip, d + n + 2, ws);
}

void add_max_support_constraints(glp_prob* mip, int n, int d, int s, int R_max) {
  //int num_constraints = 2 + s;
  glp_add_rows(mip, n + s + 4);

  // consistency constraint
  for (int i = 1; i <= s; i++) {
    glp_set_row_bnds(mip, i, GLP_LO, 1.0, 0.0);
  }

  // rule constraint
  int M = d + 1;
  for (int i = 1; i <= n; i++) {
    glp_set_row_bnds(mip, i + s, GLP_UP, 0.0, M);
  }

  // relevance constraint: sum(alpha*alpha_e) - R == 0
  glp_set_row_bnds(mip, n + s + 1, GLP_FX, 0.0, 0.0);

  // auxiliary constraint (sum(alpha) == R)
  glp_set_row_bnds(mip, n + s + 2, GLP_FX, 0.0, 0.0);

  // max sparsity constraint
  glp_set_row_bnds(mip, n + s + 3, GLP_UP, 0.0, R_max);

  // auxillary constraint (sum(r) == S)
  glp_set_row_bnds(mip, n + s + 4, GLP_FX, 0.0, 0.0);
}

void add_max_support_matrix(glp_prob* mip, Rcpp::NumericVector mi, Rcpp::NumericVector mj, Rcpp::NumericVector mv) {
  int num_entries = mi.size();
  int m_row[1 + num_entries];
  int m_col[1 + num_entries];
  double m_val[1 + num_entries];

  // set contraints
  for (int k = 1; k <= num_entries; k++) {
    m_row[k] = mi(k - 1);
    m_col[k] = mj(k - 1);
    m_val[k] = mv(k - 1);
  }

  // load constraint matrix
  glp_load_matrix(mip, num_entries, m_row, m_col, m_val);
}

// [[Rcpp::export]]
Rcpp::List max_support_cpp(Rcpp::NumericVector mi, Rcpp::NumericVector mj, Rcpp::NumericVector mv, int n, int d, int s, int R_max, double ws) {
  // create problem
  glp_prob *mip;
  mip = glp_create_prob();
  glp_set_prob_name(mip, "max_support");
  glp_set_obj_dir(mip, GLP_MAX);

  // lp parameters
  glp_smcp lp_parm;
  glp_init_smcp(&lp_parm);
  lp_parm.meth = GLP_DUALP;

  // ip parameters
  glp_iocp ip_parm;
  glp_init_iocp(&ip_parm);
  ip_parm.tm_lim = 30 * 1000; // milliseconds

  // setup variables and constraints
  add_max_support_vars(mip, n, d, R_max, ws);
  add_max_support_constraints(mip, n, d, s, R_max);
  add_max_support_matrix(mip, mi, mj, mv);

  // write problem
  glp_write_lp(mip, NULL, "problem.lp");

  // solve problem
  int ret;
  ret = glp_simplex(mip, &lp_parm);
  if (ret == 0) {
    glp_intopt(mip, &ip_parm);
  } else {
    return Rcpp::List::create(
      Rcpp::Named("error_code") = ret
    );
  }

  // extract solution
  int status = glp_mip_status(mip);
  double obj_val = glp_mip_obj_val(mip);

  std::vector<int> alpha;
  for (int i = 1; i <= d; i++) {
    if (glp_mip_col_val(mip, i) > 0) {
      alpha.push_back(i);
    }
  }

  std::vector<int> r;
  for (int i = 1; i <= n; i++) {
    if (glp_mip_col_val(mip, d + i) > 0) {
      r.push_back(i);
    }
  }

  int R = glp_mip_col_val(mip, 1 + d + n);
  int S = glp_mip_col_val(mip, 2 + d + n);

  // cleanup problem
  glp_delete_prob(mip);

  // return value
  return Rcpp::List::create(
    Rcpp::Named("is_optimal") = status == GLP_OPT,
    Rcpp::Named("obj") = glp_mip_obj_val(mip),
    Rcpp::Named("sparsity") = R,
    Rcpp::Named("support") = S,
    Rcpp::Named("alpha") = alpha,
    Rcpp::Named("r") = r
  );
}
