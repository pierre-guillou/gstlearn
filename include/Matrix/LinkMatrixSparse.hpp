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

#include "LinearOp/CholeskySparse.hpp"
#include "gstlearn_export.hpp"

#ifndef SWIG

namespace gstlrn
{
class MatrixSparse;
class cs;

class GSTLEARN_EXPORT QChol
{
public:
  MatrixSparse* Q;
  CholeskySparse* chol;
};

GSTLEARN_EXPORT void cs_print_dim(const char* title, const cs* A);
GSTLEARN_EXPORT cs* cs_duplicate(const cs* b1);

// Qchol operations
GSTLEARN_EXPORT bool is_chol_ready(QChol* QC);
GSTLEARN_EXPORT QChol* qchol_free(QChol* QC);
GSTLEARN_EXPORT int qchol_getNCols(QChol* QC);
GSTLEARN_EXPORT int qchol_getNRows(QChol* QC);
GSTLEARN_EXPORT int qchol_cholesky(int verbose, QChol* QC);
GSTLEARN_EXPORT void cs_chol_invert(QChol* qctt, double* xcr, const double* rhs, const double* work);
GSTLEARN_EXPORT void cs_chol_simulate(QChol* qctt, double* simu, const double* work);

#endif // SWIG
}