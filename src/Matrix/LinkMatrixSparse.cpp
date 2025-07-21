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
#include "Matrix/LinkMatrixSparse.hpp"

#include "Basic/AStringable.hpp"
#include "Basic/OptDbg.hpp"

#include "geoslib_old_f.h"

// External library
#include "csparse_f.h"

#define MAX_NEIGH      100
#define XCR(ilevel, i) (xcr[(ilevel) * ncur + (i)])
#define RHS(ilevel, i) (rhs[(ilevel) * ncur + (i)])
#define MAT(i, j)      (mat[(i) * n + (j)])
#define DEBUG          0
#define MAX_NEIGH      100
#define XCR(ilevel, i) (xcr[(ilevel) * ncur + (i)])
#define RHS(ilevel, i) (rhs[(ilevel) * ncur + (i)])
#define MAT(i, j)      (mat[(i) * n + (j)])
#define DEBUG          0

namespace gstlrn
{
/****************************************************************************/
/*!
 **  Finalize the construction of the QChol structure.
 **  Perform the Cholesky decomposition
 **
 ** \return  Error return code
 **
 ** \param[in]  verbose   Verbose flag
 ** \param[in]  QC   QChol structure to be finalized
 **
 ** \remarks In case of problem the message is issued in this function
 ** \remarks If the decomposition is already performed, nothing is done
 **
 *****************************************************************************/
int qchol_cholesky(int verbose, QChol* QC)
{
  /* Check that the Q matrix has already been defined */

  if (QC->Q == nullptr) return (1);

  /* Cholesky decomposition is only valid for square matric */

  if (qchol_getNRows(QC) != qchol_getNCols(QC))
  {
    messerr("You wish to perform a Cholesky Decomposition of a Matrix");
    messerr("which is not square: %d x %d", qchol_getNRows(QC), qchol_getNCols(QC));
    messerr("This must be an error");
    return (1);
  }

  if (verbose) message("  Cholesky Decomposition... ");

  // Perform the Cholesky decomposition (new style)

  if (QC->chol == nullptr)
    // Perform the Cholesky decomposition (new style)

    if (QC->chol == nullptr)
    {
      QC->chol = new CholeskySparse(*QC->Q);
      if (QC->chol == nullptr)
      {
        messerr("Error in Cholesky decompostion (new version)");
        messerr("Error in Cholesky decompostion (new version)");
        goto label_err;
      }
    }

  if (verbose) message("Finished\n");

  /* Optional printout */

  if (OptDbg::query(EDbg::KRIGING) || OptDbg::force())
  {
    message("Q Sparse Matrix\n");
    print_matrix("", 1, *QC->Q);
    QC->Q->dumpRange("Q");
  }
  return (0);

label_err:
  delete QC->chol;
  return (1);
}

bool is_chol_ready(QChol* QC)
{
  return QC->chol != nullptr;
}

QChol* qchol_free(QChol* QC)
{
  if (QC == nullptr) return QC;
  delete QC->chol;
  delete QC;
  return QC;
}
int qchol_getNCols(QChol* QC)
{
  if (QC == nullptr) return 0;
  if (QC->Q == nullptr) return 0;
  return QC->Q->getNCols();
}
int qchol_getNRows(QChol* QC)
{
  if (QC == nullptr) return 0;
  if (QC->Q == nullptr) return 0;
  return QC->Q->getNRows();
}

/****************************************************************************/
/*!
 **  Inversion using Cholesky
 **
 ** \param[in]  qctt     Qchol structure
 ** \param[in,out]  xcr  Current vector
 ** \param[in]  rhs      Current R.H.S. vector
 **
 ** \param[out] work     Working array
 **
 *****************************************************************************/
void cs_chol_invert(QChol* qctt, double* xcr, const double* rhs, const double* work)
{
  DECLARE_UNUSED(work);
  DECLARE_UNUSED(work);
  if (DEBUG) message("Cholesky Inversion\n");
  int n = qchol_getNCols(qctt);

  VectorDouble rhsVD(n);
  for (int i = 0; i < n; i++) rhsVD[i] = rhs[i];
  VectorDouble xcrVD = qctt->chol->solveX(rhsVD);
  for (int i = 0; i < n; i++) xcr[i] = xcrVD[i];
}

/****************************************************************************/
/*!
 **  Simulate using Cholesky
 **
 ** \param[in]  qctt     Qchol structure
 **
 ** \param[out] simu     Simulated array
 ** \param[out] work     Working array
 **
 *****************************************************************************/
void cs_chol_simulate(QChol* qctt, double* simu, const double* work)
{
  if (DEBUG) message("Cholesky Simulation\n");
  int n = qchol_getNCols(qctt);

  VectorDouble simuVD(n);
  VectorDouble workVD(n);
  for (int i = 0; i < n; i++) workVD[i] = work[i];
  (void)qctt->chol->addSimulateToDest(workVD, simuVD);
  for (int i = 0; i < n; i++) simu[i] = simuVD[i];
}


} // namespace gstlrn