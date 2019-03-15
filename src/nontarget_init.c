#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void mass(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP adduct_search(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP combine_tuple(SEXP, SEXP, SEXP, SEXP);
extern SEXP homol_triplet(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP kdtree(SEXP);
extern SEXP kdtree4(SEXP, SEXP, SEXP);
extern SEXP metagroup(SEXP, SEXP);
extern SEXP node_delete(SEXP, SEXP, SEXP);
extern SEXP peak_search(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP search_boxtree(SEXP, SEXP, SEXP, SEXP);
extern SEXP search_kdtree(SEXP, SEXP, SEXP);
extern SEXP search_kdtree_homol(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP search_kdtree3(SEXP, SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"mass", (DL_FUNC) &mass, 19},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"adduct_search",       (DL_FUNC) &adduct_search,        8},
    {"combine_tuple",       (DL_FUNC) &combine_tuple,        4},
    {"homol_triplet",       (DL_FUNC) &homol_triplet,        9},
    {"kdtree",              (DL_FUNC) &kdtree,               1},
    {"kdtree4",             (DL_FUNC) &kdtree4,              3},
    {"metagroup",           (DL_FUNC) &metagroup,            2},
    {"node_delete",         (DL_FUNC) &node_delete,          3},
    {"peak_search",         (DL_FUNC) &peak_search,         10},
    {"search_boxtree",      (DL_FUNC) &search_boxtree,       4},
    {"search_kdtree",       (DL_FUNC) &search_kdtree,        3},
    {"search_kdtree_homol", (DL_FUNC) &search_kdtree_homol,  7},
    {"search_kdtree3",      (DL_FUNC) &search_kdtree3,       4},
    {NULL, NULL, 0}
};

void R_init_nontarget(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
