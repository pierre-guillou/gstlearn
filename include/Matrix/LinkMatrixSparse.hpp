/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Basic/VectorNumT.hpp"

#ifndef SWIG

namespace gstlrn
{
class MatrixSparse;
class cs;
class css;
class csn;
class csd;
class NF_Triplet;

class GSTLEARN_EXPORT QChol
{
public:
  MatrixSparse* Q;
  css* S;
  csn* N;
};

class GSTLEARN_EXPORT cs_MG
{
public:
  int nh;
  int nH;
  double* sumrow;
  MatrixSparse* IhH;
  QChol* A;
};

class GSTLEARN_EXPORT cs_MGS
{
public:
  int flag_cg;     /* Apply Conjugate-gradient */
  int nlevels;     /* Number of multigrid levels */
  int npath;       /* Number of paths */
  int type_coarse; /* Type of coarsening algorithm */
  int ngc;         /* Maximum count of Conjugate-Gradient iters */
  int nmg;         /* Maximum count of Multigrid operations*/
  int ngs;         /* Number of Gauss_Siedel relaxation cycles */
  int ncur;        /* Number of mesh vertices */
  int* path;       /* Path description */
  double tolnmg;   /* Tolerance for Multigrid */
  double tolcg;    /* Tolerance for Conjugate-Gradient */
  double* diag;    /* Normation diagonal */
  cs_MG** mg;      /* Array of cs_MG structures */
};

/// TODO : cs_*2 functions to be removed (encapsulation)
GSTLEARN_EXPORT cs* cs_spfree2(cs* A);
GSTLEARN_EXPORT css* cs_sfree2(css* S);
GSTLEARN_EXPORT csn* cs_nfree2(csn* N);

GSTLEARN_EXPORT void cs_vector_xM(const cs* A, int nout, const double* x, double* y);
GSTLEARN_EXPORT void cs_mulvec_uptri(const cs* A, int nout, const double* x, double* y, int flag_diag);
GSTLEARN_EXPORT void cs_mulvec_lowtri(const cs* A, int nout, const double* x, double* y, int flag_diag);
GSTLEARN_EXPORT cs* cs_matvecnorm(const cs* A, const double* x, int oper);
GSTLEARN_EXPORT void cs_matvecnorm_inplace(cs* A, const double* x, int oper);
GSTLEARN_EXPORT double* cs_col_sumrow(const cs* A, int* ncol, int* nrow);

GSTLEARN_EXPORT void cs_print_nice(const char* title, const cs* A, int maxrow = -1, int maxcol = -1);
GSTLEARN_EXPORT void cs_print_dim(const char* title, const cs* A);
GSTLEARN_EXPORT void cs_print_file(const char* radix, int rank, const cs* A);

GSTLEARN_EXPORT void cs_force_dimension(cs* T, int nrow, int ncol);
GSTLEARN_EXPORT cs* cs_diag(VectorDouble diag, double tol = EPSILON10);

GSTLEARN_EXPORT int cs_lsolve_lowtri(const cs* L, const double* x, double* y);
GSTLEARN_EXPORT int cs_lsolve_uptri(const cs* L, const double* x, double* y);

// Qchol operations
GSTLEARN_EXPORT int qchol_cholesky(int verbose, QChol* QC);
GSTLEARN_EXPORT void cs_chol_invert(QChol* qctt, double* xcr, double* rhs, double* work);
GSTLEARN_EXPORT void cs_chol_simulate(QChol* qctt, double* simu, double* work);

// Multigrid operations
GSTLEARN_EXPORT cs_MGS* cs_multigrid_manage(cs_MGS* mgs, int mode, int nlevels, int path_type);
GSTLEARN_EXPORT void cs_multigrid_params(cs_MGS* mgs, int flag_cg, int type_coarse, int ngc, int nmg, int ngs, double tolcg, double tolnmg);
GSTLEARN_EXPORT void cs_multigrid_print(cs_MGS* mgs);
GSTLEARN_EXPORT int cs_multigrid_setup(cs_MGS* mgs, QChol* Qctt, int flag_sel, int verbose, double** sel);
GSTLEARN_EXPORT int cs_multigrid_process(cs_MGS* mgs, QChol* qctt, int verbose, double* x, double* b, double* work);

GSTLEARN_EXPORT NF_Triplet csToTriplet(const cs* A, int shiftRow = 0, int shiftCol = 0, double tol = EPSILON10);

GSTLEARN_EXPORT void cs_print_range(const char* title, const cs* C);
GSTLEARN_EXPORT cs* cs_eye(int number, double value);
GSTLEARN_EXPORT cs* cs_extract_diag(const cs* C, int oper_choice = 1);
GSTLEARN_EXPORT double* csd_extract_diag(const cs* C, int oper_choice = 1);

GSTLEARN_EXPORT void cs_rowcol(const cs* A, int* nrows, int* ncols, int* count, double* percent);
GSTLEARN_EXPORT cs* cs_duplicate(const cs* b1);
GSTLEARN_EXPORT cs* cs_normalize_by_diag_and_release(cs* Q, int flag_rel);
GSTLEARN_EXPORT cs* cs_prod_norm_and_release(cs* b1, cs* lambda, int flag_rel);
GSTLEARN_EXPORT int cs_coarsening(const cs* Q, int type, int** indCo, cs** L);
GSTLEARN_EXPORT cs* cs_interpolate(const cs* AA, const cs* Lt, const int* Co);
GSTLEARN_EXPORT int cs_get_nrow(const cs* A);
GSTLEARN_EXPORT int cs_get_ncol(const cs* A);

GSTLEARN_EXPORT void cs_set_status_update_nonzero_value(int status = 2);
GSTLEARN_EXPORT int cs_get_status_update_nonzero_value();

#endif // SWIG
}