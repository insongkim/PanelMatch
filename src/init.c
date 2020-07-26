#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP _PanelMatch_check_control_units_for_treatment_restriction(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _PanelMatch_check_missing_data_control_units(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _PanelMatch_check_missing_data_treated_units(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _PanelMatch_check_treated_units_for_treatment_reversion(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _PanelMatch_do_exact_matching_refinement(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _PanelMatch_equality_four_cpp(SEXP, SEXP, SEXP);
extern SEXP _PanelMatch_get_comparison_histories(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _PanelMatch_get_dits(SEXP, SEXP);
extern SEXP _PanelMatch_get_msets_helper(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _PanelMatch_get_treated_indices(SEXP, SEXP, SEXP, SEXP);
extern SEXP _PanelMatch_get_vit_index(SEXP, SEXP, SEXP);
extern SEXP _PanelMatch_get_vit_index_unsigned(SEXP, SEXP, SEXP);
extern SEXP _PanelMatch_get_yearly_dmats(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _PanelMatch_handle_vits(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _PanelMatch_multiply_weights_msm(SEXP, SEXP);
extern SEXP _PanelMatch_needs_renormalization(SEXP);
extern SEXP _PanelMatch_non_matching_matcher(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _PanelMatch_prep_lead_years(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_PanelMatch_check_control_units_for_treatment_restriction", (DL_FUNC) &_PanelMatch_check_control_units_for_treatment_restriction, 6},
  {"_PanelMatch_check_missing_data_control_units",              (DL_FUNC) &_PanelMatch_check_missing_data_control_units,              5},
  {"_PanelMatch_check_missing_data_treated_units",              (DL_FUNC) &_PanelMatch_check_missing_data_treated_units,              6},
  {"_PanelMatch_check_treated_units_for_treatment_reversion",   (DL_FUNC) &_PanelMatch_check_treated_units_for_treatment_reversion,   6},
  {"_PanelMatch_do_exact_matching_refinement",                  (DL_FUNC) &_PanelMatch_do_exact_matching_refinement,                  6},
  {"_PanelMatch_equality_four_cpp",                             (DL_FUNC) &_PanelMatch_equality_four_cpp,                             3},
  {"_PanelMatch_get_comparison_histories",                      (DL_FUNC) &_PanelMatch_get_comparison_histories,                      7},
  {"_PanelMatch_get_dits",                                      (DL_FUNC) &_PanelMatch_get_dits,                                      2},
  {"_PanelMatch_get_msets_helper",                              (DL_FUNC) &_PanelMatch_get_msets_helper,                              5},
  {"_PanelMatch_get_treated_indices",                           (DL_FUNC) &_PanelMatch_get_treated_indices,                           4},
  {"_PanelMatch_get_vit_index",                                 (DL_FUNC) &_PanelMatch_get_vit_index,                                 3},
  {"_PanelMatch_get_vit_index_unsigned",                        (DL_FUNC) &_PanelMatch_get_vit_index_unsigned,                        3},
  {"_PanelMatch_get_yearly_dmats",                              (DL_FUNC) &_PanelMatch_get_yearly_dmats,                              6},
  {"_PanelMatch_handle_vits",                                   (DL_FUNC) &_PanelMatch_handle_vits,                                   7},
  {"_PanelMatch_multiply_weights_msm",                          (DL_FUNC) &_PanelMatch_multiply_weights_msm,                          2},
  {"_PanelMatch_needs_renormalization",                         (DL_FUNC) &_PanelMatch_needs_renormalization,                         1},
  {"_PanelMatch_non_matching_matcher",                          (DL_FUNC) &_PanelMatch_non_matching_matcher,                          6},
  {"_PanelMatch_prep_lead_years",                               (DL_FUNC) &_PanelMatch_prep_lead_years,                               2},
  {NULL, NULL, 0}
};

void R_init_PanelMatch(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}