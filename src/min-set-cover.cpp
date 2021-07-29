#include <Rcpp.h>
#include "glpk.h"

void add_vars(glp_prob* mip, int n, int d) {
  int num_vars = d + 1;
  glp_add_cols(mip, num_vars);

  // alpha binary variables
  std::string varname;
  for (int i = 1; i <= d; i++) {
    varname = "alpha" + std::to_string(i);
    glp_set_col_name(mip, i, varname.c_str());
    glp_set_col_kind(mip, i, GLP_BV);
  }

  // R
  glp_set_col_name(mip, d + 1, "R");
  glp_set_col_kind(mip, d + 1, GLP_IV);
  glp_set_col_bnds(mip, d + 1, GLP_DB, 1, d);
  glp_set_obj_coef(mip, d + 1, 1.0);
}

void add_constraints(glp_prob* mip, int s) {
  //int num_constraints = 2 + s;
  glp_add_rows(mip, s + 2);

  // consistency constraint
  for (int i = 1; i <= s; i++) {
    glp_set_row_bnds(mip, i, GLP_LO, 1.0, 0.0);
  }

  // relevance constraint: sum(alpha*alpha_e) - R == 0
  glp_set_row_bnds(mip, s + 1, GLP_FX, 0.0, 0.0);

  // auxiliary constraint
  glp_set_row_bnds(mip, s + 2, GLP_FX, 0.0, 0.0);
}

void add_matrix(glp_prob* mip, Rcpp::NumericVector mi, Rcpp::NumericVector mj, Rcpp::NumericVector mv) {
  int num_entries = mi.size();
  int m_row[1 + num_entries];
  int m_col[1 + num_entries];
  double m_val[1 + num_entries];

  // constraint 1: sum(alpha) - R == 0 -> 1*a1 + 1*a2 + .. 1*ad - 1*R == 0
  for (int k = 1; k <= num_entries; k++) {
    m_row[k] = mi(k - 1);
    m_col[k] = mj(k - 1);
    m_val[k] = mv(k - 1);
  }

  // load constraint matrix
  glp_load_matrix(mip, num_entries, m_row, m_col, m_val);
}

// [[Rcpp::export]]
Rcpp::List min_set_cover_cpp(Rcpp::NumericVector mi, Rcpp::NumericVector mj, Rcpp::NumericVector mv, int n, int d, int s) {
  // create problem
  glp_prob *mip;
  mip = glp_create_prob();
  glp_set_prob_name(mip, "min_set_cover");
  glp_set_obj_dir(mip, GLP_MIN);

  // lp parameters
  glp_smcp lp_parm;
  glp_init_smcp(&lp_parm);
  lp_parm.meth = GLP_DUALP;

  // ip parameters
  glp_iocp ip_parm;
  glp_init_iocp(&ip_parm);
  ip_parm.tm_lim = 30 * 1000; // milliseconds

  // setup variables and constraints
  add_vars(mip, n, d);
  add_constraints(mip, s);
  add_matrix(mip, mi, mj, mv);

  // write problem
  //glp_write_lp(mip, NULL, "problem.lp");

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
  int R = glp_mip_col_val(mip, d + 1);

  // cleanup problem
  glp_delete_prob(mip);

  // return value
  return Rcpp::List::create(
    Rcpp::Named("is_optimal") = status == GLP_OPT,
    Rcpp::Named("obj") = glp_mip_obj_val(mip),
    Rcpp::Named("sparsity") = R,
    Rcpp::Named("alpha") = alpha
  );
}
