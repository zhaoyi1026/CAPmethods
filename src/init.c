// Auto-generated native-routine registration for the CAPmethods package.
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <R_ext/Rdynload.h>

extern SEXP hd_covls_cpp(SEXP);
extern SEXP hd_score_cpp(SEXP, SEXP);
extern SEXP hd_capbeta_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP hd_capd1_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP hd_objfunc_cpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP lcap_covls_cpp(SEXP);
extern SEXP lcap_gammasolve_cpp(SEXP, SEXP);
extern SEXP lcap_score_cpp(SEXP, SEXP);
extern SEXP lcap_recar_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP lcap_capd1_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP lcap_objfunc_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP mlcap_compute_scores_cpp(SEXP, SEXP, SEXP);
extern SEXP mlcap_obj_func_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP coc_cov_ls_cpp(SEXP);
extern SEXP coc_cov_sk_x_cpp(SEXP);
extern SEXP coc_cov_sk_y_cpp(SEXP, SEXP, SEXP);
extern SEXP coc_score_cpp(SEXP, SEXP);
extern SEXP coc_accum_cpp(SEXP, SEXP);
extern SEXP coc_eigen_solve_cpp(SEXP, SEXP);
extern SEXP med_score_cpp(SEXP, SEXP);
extern SEXP med_accum_cpp(SEXP, SEXP);
extern SEXP med_eigen_solve_cpp(SEXP, SEXP);
extern SEXP med_cov_cpp(SEXP);
extern SEXP hdcov_precompute_cpp(SEXP);
extern SEXP hdcov_score_cpp(SEXP, SEXP);
extern SEXP hdcov_accum_cpp(SEXP, SEXP);
extern SEXP cluster_smat_cpp(SEXP);
extern SEXP cluster_score_cpp(SEXP, SEXP);
extern SEXP cluster_accum_cpp(SEXP, SEXP);
extern SEXP cluster_eigen_solve_cpp(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"hd_covls_cpp", (DL_FUNC) &hd_covls_cpp, 1},
    {"hd_score_cpp", (DL_FUNC) &hd_score_cpp, 2},
    {"hd_capbeta_cpp", (DL_FUNC) &hd_capbeta_cpp, 8},
    {"hd_capd1_cpp", (DL_FUNC) &hd_capd1_cpp, 9},
    {"hd_objfunc_cpp", (DL_FUNC) &hd_objfunc_cpp, 5},
    {"lcap_covls_cpp", (DL_FUNC) &lcap_covls_cpp, 1},
    {"lcap_gammasolve_cpp", (DL_FUNC) &lcap_gammasolve_cpp, 2},
    {"lcap_score_cpp", (DL_FUNC) &lcap_score_cpp, 2},
    {"lcap_recar_cpp", (DL_FUNC) &lcap_recar_cpp, 16},
    {"lcap_capd1_cpp", (DL_FUNC) &lcap_capd1_cpp, 15},
    {"lcap_objfunc_cpp", (DL_FUNC) &lcap_objfunc_cpp, 12},
    {"mlcap_compute_scores_cpp", (DL_FUNC) &mlcap_compute_scores_cpp, 3},
    {"mlcap_obj_func_cpp", (DL_FUNC) &mlcap_obj_func_cpp, 15},
    {"coc_cov_ls_cpp", (DL_FUNC) &coc_cov_ls_cpp, 1},
    {"coc_cov_sk_x_cpp", (DL_FUNC) &coc_cov_sk_x_cpp, 1},
    {"coc_cov_sk_y_cpp", (DL_FUNC) &coc_cov_sk_y_cpp, 3},
    {"coc_score_cpp", (DL_FUNC) &coc_score_cpp, 2},
    {"coc_accum_cpp", (DL_FUNC) &coc_accum_cpp, 2},
    {"coc_eigen_solve_cpp", (DL_FUNC) &coc_eigen_solve_cpp, 2},
    {"med_score_cpp", (DL_FUNC) &med_score_cpp, 2},
    {"med_accum_cpp", (DL_FUNC) &med_accum_cpp, 2},
    {"med_eigen_solve_cpp", (DL_FUNC) &med_eigen_solve_cpp, 2},
    {"med_cov_cpp", (DL_FUNC) &med_cov_cpp, 1},
    {"hdcov_precompute_cpp", (DL_FUNC) &hdcov_precompute_cpp, 1},
    {"hdcov_score_cpp", (DL_FUNC) &hdcov_score_cpp, 2},
    {"hdcov_accum_cpp", (DL_FUNC) &hdcov_accum_cpp, 2},
    {"cluster_smat_cpp", (DL_FUNC) &cluster_smat_cpp, 1},
    {"cluster_score_cpp", (DL_FUNC) &cluster_score_cpp, 2},
    {"cluster_accum_cpp", (DL_FUNC) &cluster_accum_cpp, 2},
    {"cluster_eigen_solve_cpp", (DL_FUNC) &cluster_eigen_solve_cpp, 2},
    {NULL, NULL, 0}
};

void R_init_CAPmethods(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, FALSE);
}
