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
#include "Matrix/NF_Triplet.hpp"

#include "Basic/AStringable.hpp"
#include "Basic/File.hpp"
#include "Basic/Memory.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/String.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/VectorHelper.hpp"
#include "Core/Keypair.hpp"

#include "geoslib_old_f.h"

#include <set>

// External library
#include "csparse_f.h"

#define MAX_NEIGH      100
#define XCR(ilevel, i) (xcr[(ilevel) * ncur + (i)])
#define RHS(ilevel, i) (rhs[(ilevel) * ncur + (i)])
#define MAT(i, j)      (mat[(i) * n + (j)])
#define DEBUG          0

namespace gstlrn
{
static int flagUpdateNonzero = 1;

/**
 * Define the status when modifying the value of a nonzero element of the sparse matrix
 * @param status 0 (no test); 1 (warning issued); 2 (throw issued)
 */
void cs_set_status_update_nonzero_value(int status)
{
  flagUpdateNonzero = status;
}

int cs_get_status_update_nonzero_value()
{
  return flagUpdateNonzero;
}

cs* cs_spfree2(cs* A)
{
  return cs_spfree(A);
}
css* cs_sfree2(css* S)
{
  return cs_sfree(S);
}
csn* cs_nfree2(csn* N)
{
  return cs_nfree(N);
}
void cs_force_dimension(cs* T, int nrow, int ncol)
{
  if (cs_get_nrow(T) > nrow)
    messageAbort("Forcing CS dimension: NRows current(%d) is larger than forecast(%d)",
                 cs_get_nrow(T), nrow);
  if (cs_getncol(T) > ncol)
    messageAbort("Forcing CS dimension: NCols current(%d) is larger than forecast(%d)",
                 cs_getncol(T), ncol);
  T->m = nrow;
  T->n = ncol;
}

/* Transform VectorDouble to cs diagonal */
cs* cs_diag(VectorDouble diag, double tol)
{
  int number = (int)diag.size();

  cs* Striplet = cs_spalloc(0, 0, 1, 1, 1);
  for (int i = 0; i < number; i++)
  {
    if (ABS(diag[i]) < tol) continue;
    if (!cs_entry(Striplet, i, i, diag[i]))
    {
      return nullptr;
    }
  }

  cs_force_dimension(Striplet, number, number);
  cs* Q = cs_triplet(Striplet);
  cs_spfree(Striplet);
  return Q;
}

/****************************************************************************/
/*!
 **  Update the selection vector
 **
 ** \param[in] ncur   Dimension of array 'sel'
 ** \param[in] sel   Array containing the current selection
 ** \param[in] indCo   Array of selected samples
 **
 ** \remarks The array 'indCo' is relative to the input selection
 **
 *****************************************************************************/
static void st_selection_update(int ncur, double* sel, const int* indCo)
{
  int ecr, nval;

  ecr = nval = 0;
  for (int iech = 0; iech < ncur; iech++)
  {
    if (sel[iech] == 0) continue;
    sel[iech] = indCo[ecr++];
    nval += static_cast<int>(sel[iech]);
  }
}

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
  int nmax = 8;

  /* Check that the Q matrix has already been defined */

  if (QC->Q == nullptr) return (1);

  /* Cholesky decomposition is only valid for square matric */

  if (QC->Q->getNRows() != QC->Q->getNCols())
  {
    messerr("You wish to perform a Cholesky Decomposition of a Matrix");
    messerr("which is not square: %d x %d", QC->Q->getNRows(), QC->Q->getNCols());
    messerr("This must be an error");
    return (1);
  }

  /* Perform the Cholesky decomposition */

  if (verbose) message("  Cholesky Decomposition... ");

  if (QC->S == nullptr)
  {
    if (verbose) message("Ordering... ");
    QC->S = cs_schol(QC->Q->getCS(), 0);
    if (QC->S == nullptr)
    {
      messerr("Error in cs_schol function");
      goto label_err;
    }
  }

  if (QC->N == nullptr)
  {
    if (verbose) message("Factorization... ");
    QC->N = cs_chol(QC->Q->getCS(), QC->S);
    if (QC->N == nullptr)
    {
      messerr("Error in cs_chol function");
      goto label_err;
    }
  }
  if (verbose) message("Finished\n");

  /* Optional printout */

  if (OptDbg::query(EDbg::KRIGING) || OptDbg::force())
  {
    message("Q Sparse Matrix\n");
    cs_print(QC->Q->getCS(), 1);
    cs_print_range("Q", QC->Q->getCS());
  }
  return (0);

label_err:
  if (verbose)
    cs_print_nice("Cholesky Decomposition of QC", QC->Q->getCS(), nmax, nmax);
  QC->N = cs_nfree(QC->N);
  QC->S = cs_sfree(QC->S);
  return (1);
}

/****************************************************************************/
/*!
 **  Define the path according to the number of levels
 **
 ** \param[in]  mgs    cs_MGS structure
 ** \param[in]  nlevels    Number of coarsening levels
 ** \param[in]  path_type Type of the path
 **
 *****************************************************************************/
static void st_path_define(cs_MGS* mgs, int nlevels, int path_type)
{
  int ecr, n;
  static int niw  = 6;
  static int iw[] = {1, 1, -1, 1, -1, -1};

  /* Dispatch */

  switch (path_type)
  {
    case 1: // V type
      mgs->npath = 2 * nlevels;
      mgs->path  = (int*)mem_alloc(sizeof(int) * mgs->npath, 1);

      ecr = 0;
      for (int i = 0; i < nlevels; i++)
        mgs->path[ecr++] = 1;
      for (int i = 0; i < nlevels; i++)
        mgs->path[ecr++] = -1;
      break;

    case 2: // W type
      if (nlevels == 1)
      {
        mgs->npath   = 2 * nlevels;
        mgs->path    = (int*)mem_alloc(sizeof(int) * mgs->npath, 1);
        mgs->path[0] = 1;
        mgs->path[1] = -1;
      }
      else
      {
        mgs->npath = niw;
        mgs->path  = (int*)mem_alloc(sizeof(int) * mgs->npath, 1);
        for (int i = 0; i < niw; i++)
          mgs->path[i] = iw[i];

        for (int ilevel = 2; ilevel < nlevels; ilevel++)
        {
          n          = mgs->npath;
          mgs->npath = 2 + 2 * n;
          mgs->path  = (int*)mem_realloc((char*)mgs->path,
                                         sizeof(int) * mgs->npath, 1);
          for (int i = n - 1; i >= 0; i--)
          {
            mgs->path[i + 1]     = mgs->path[i];
            mgs->path[i + 1 + n] = mgs->path[i];
          }
          mgs->path[0]         = 1;
          mgs->path[2 * n + 1] = -1;
        }
      }
      break;

    case 3: // F type
      mgs->npath = nlevels;
      mgs->path  = (int*)mem_alloc(sizeof(int) * mgs->npath, 1);
      for (int i = 0; i < nlevels; i++)
        mgs->path[i] = 1;

      for (int ilevel = 1; ilevel < nlevels; ilevel++)
      {
        n = mgs->npath;
        mgs->npath += 2 * ilevel;
        mgs->path = (int*)mem_realloc((char*)mgs->path,
                                      sizeof(int) * mgs->npath, 1);
        for (int i = 0; i < ilevel; i++)
        {
          mgs->path[n + i]          = -1;
          mgs->path[n + i + ilevel] = 1;
        }
      }

      n = mgs->npath;
      mgs->npath += nlevels;
      mgs->path = (int*)mem_realloc((char*)mgs->path,
                                    sizeof(int) * mgs->npath, 1);
      for (int i = 0; i < nlevels; i++)
        mgs->path[i + n] = -1;
      break;

    default:
      break;
  }
}

/****************************************************************************/
/*!
 **  Define the default parameters for the Multigrid option
 **
 ** \param[in]  mgs    cs_MGS structure
 **
 *****************************************************************************/
static void st_multigrid_set_default_params(cs_MGS* mgs)
{
  mgs->flag_cg     = 1;
  mgs->type_coarse = 1;
  mgs->ngc         = 50;
  mgs->nmg         = 4;
  mgs->ngs         = 1;
  mgs->tolcg       = 1.e-07;
  mgs->tolnmg      = 1.e-07;
}

/****************************************************************************/
/*!
 **  Define the parameters for the multigrid manipulation
 **
 ** \param[in]  mgs     cs_MGS structure
 ** \param[in]  flag_cg     1 for Conjugate-Gradient use; 0 otherwise
 ** \param[in]  type_coarse Type of coarsening algorithm
 ** \param[in]  ngc     Maximum number of Conjugate-Gradient iterations
 ** \param[in]  nmg     Maximum number of mutligrid iterations
 ** \param[in]  ngs     Number of Gauss-Siedel relaxation cycles
 ** \param[in]  tolcg     Tolerance for the Conjugate-Gradient algorithm
 ** \param[in]  tolnmg     Tolerance for the Multigrid algorithm
 **
 *****************************************************************************/
void cs_multigrid_params(cs_MGS* mgs,
                         int flag_cg,
                         int type_coarse,
                         int ngc,
                         int nmg,
                         int ngs,
                         double tolcg,
                         double tolnmg)
{
  if (mgs == nullptr) return;
  mgs->flag_cg = flag_cg;
  if (mgs->nlevels == 0) mgs->flag_cg = 0;
  mgs->type_coarse = type_coarse;
  mgs->ngc         = ngc;
  mgs->nmg         = nmg;
  mgs->ngs         = ngs;
  mgs->tolcg       = tolcg;
  mgs->tolnmg      = tolnmg;
}

/****************************************************************************/
/*!
 **  Allocate one sub-structure for the multigrid manipulation
 **
 ** \param[in]  mode     1 for allocation; -1 for deallocation
 ** \param[in]  mg     cs_MG structure (used for deallocation)
 **
 *****************************************************************************/
static cs_MG* st_monogrid_manage(int mode, cs_MG* mg)
{
  /* Dispatch */

  if (mode > 0)
  {
    mg         = (cs_MG*)mem_alloc(sizeof(cs_MG), 1);
    mg->nh     = 0;
    mg->nH     = 0;
    mg->sumrow = nullptr;
    mg->IhH    = nullptr;
    mg->A      = qchol_manage(1, NULL);
  }
  else
  {
    if (mg != nullptr)
    {
      mg->sumrow = (double*)mem_free((char*)mg->sumrow);
      delete mg->IhH;
      mg->IhH = nullptr;
      mg->A   = qchol_manage(-1, mg->A);
      mg      = (cs_MG*)mem_free((char*)mg);
    }
  }
  return (mg);
}

/****************************************************************************/
/*!
 **  Print the path
 **
 *****************************************************************************/
static void st_path_print(int nlevels, int npath, int* path)
{
  if (npath == 0) return;
  message("MultiGrid Path =");
  for (int i = 0; i < npath; i++)
    message(" %d", path[i]);
  message(" -> Number of levels = %d\n", nlevels);
}

/****************************************************************************/
/*!
 **  Allocate the structure for the multigrid manipulation
 **
 ** \param[in]  mgs     cs_MGS structure
 ** \param[in]  mode     1 for allocation; -1 for deallocation
 ** \param[in]  nlevels     Number of levels of the multigrid (only for mode > 0)
 ** \param[in]  path_type   Type of the Path (1:V; 2:W, 3:F) (only for mode > 0)
 **
 *****************************************************************************/
cs_MGS* cs_multigrid_manage(cs_MGS* mgs, int mode, int nlevels, int path_type)
{
  /* Dispatch */

  if (mode > 0)
  {
    mgs          = (cs_MGS*)mem_alloc(sizeof(cs_MGS), 1);
    mgs->nlevels = nlevels;
    mgs->diag    = nullptr;
    st_path_define(mgs, nlevels, path_type);
    st_multigrid_set_default_params(mgs);

    mgs->mg = (cs_MG**)mem_alloc(sizeof(cs_MG*) * (1 + nlevels), 1);
    for (int i = 0; i <= nlevels; i++)
      mgs->mg[i] = st_monogrid_manage(1, NULL);
  }
  else
  {
    if (mgs == nullptr) return (mgs);
    mgs->path = (int*)mem_free((char*)mgs->path);
    mgs->diag = (double*)mem_free((char*)mgs->diag);
    for (int i = 0; i <= nlevels; i++)
      mgs->mg[i] = st_monogrid_manage(-1, mgs->mg[i]);
    mgs->mg = (cs_MG**)mem_free((char*)mgs->mg);
    mgs     = (cs_MGS*)mem_free((char*)mgs);
  }
  return (mgs);
}

/****************************************************************************/
/*!
 **  Print one cs_MG structure
 **
 *****************************************************************************/
static void st_mg_print(cs_MGS* mgs, int rank)
{
  cs_MG* mg;

  mg = mgs->mg[rank];
  mestitle(2, "Contents of the MG structure for level %d", rank);
  if (mg->nh > 0 && mg->nH > 0)
    message("Transition between %d and %d vertices\n", mg->nh, mg->nH);
  if (mg->IhH != NULL) cs_print_range("Range of IhH", mg->IhH->getCS());
  if (mg->A != NULL) cs_print_range("Range of A", mg->A->Q->getCS());
  if (mg->sumrow != NULL)
    print_range("Range of sumrow", mg->nh, mg->sumrow, NULL);
}

/****************************************************************************/
/*!
 **  Print the cs_MGS structure
 **
 ** \param[in] mgs     cs_MGS structure
 **
 *****************************************************************************/
void cs_multigrid_print(cs_MGS* mgs)
{
  mestitle(1, "Multigrid Levels");
  st_path_print(mgs->nlevels, mgs->npath, mgs->path);
  print_range("Range of diag", mgs->ncur, mgs->diag, NULL);
  for (int i = 0; i <= mgs->nlevels; i++)
    st_mg_print(mgs, i);
}

/****************************************************************************/
/*!
 **  Normalize / Denormalize the initial solution and RHS
 **  for the multigrid kriging
 **
 ** \param[in] mgs    cs_MGS structure
 ** \param[in] mode    1 : Normalize; -1 : Denormalize
 ** \param[in] z      Solution vector
 ** \param[in] b      Right-hand side
 **
 ** \remarks: When mode>0, z and b are divided by diagonal
 ** \remarks: When mode<0, z is multiplied by diagonal
 **
 *****************************************************************************/
static void st_multigrid_scale(cs_MGS* mgs, int mode, double* z, double* b)
{

  /* Dispatch */

  if (mgs->diag == nullptr) return;
  if (mode > 0)
  {
    for (int icur = 0; icur < mgs->ncur; icur++)
    {
      b[icur] *= mgs->diag[icur];
      z[icur] /= mgs->diag[icur];
    }
  }
  else
  {
    for (int icur = 0; icur < mgs->ncur; icur++)
      z[icur] *= mgs->diag[icur];
  }
}

/****************************************************************************/
/*!
 **  Ascent step
 **
 ** \param[in]  mgs     cs_MGS structure
 ** \param[in]  level     Current level
 ** \param[in]  flag_init   1 if the output vector must be initialized
 ** \param[in]  flag_scale  1 if the output vector must be scaled
 ** \param[in]  zin     Input vector
 **
 ** \param[out] zout     Output vector
 ** \param[out] work     Working array
 **
 ** \remark Arguments 'zin' and 'zout' may coincide (if correctly dimensionned)
 **
 *****************************************************************************/
static void st_multigrid_ascent(cs_MGS* mgs,
                                int level,
                                int flag_init,
                                int flag_scale,
                                double* zin,
                                double* zout,
                                double* work)
{
  cs_MG* mg;

  if (DEBUG)
    message("Ascending from %d to %d (init=%d scale=%d)\n", level + 1, level,
            flag_init, flag_scale);
  mg = mgs->mg[level];
  mg->IhH->prodMatVecInPlacePtr(zin, work, true);
  if (flag_init)
    for (int icur = 0; icur < mg->nh; icur++)
      zout[icur] = work[icur];
  else
    for (int icur = 0; icur < mg->nh; icur++)
      zout[icur] += work[icur];

  if (flag_scale)
    for (int icur = 0; icur < mg->nh; icur++)
      if (mg->sumrow[icur] != 0.) zout[icur] /= mg->sumrow[icur];
}

/****************************************************************************/
/*!
 **  Descent step
 **
 ** \param[in]  mgs     cs_MGS structure
 ** \param[in]  level     Current level
 ** \param[in]  zin     Input vector
 ** \param[in]  rhsin     Input R.H.S. vector
 **
 ** \param[out] rhsout     Output R.H.S. vector
 ** \param[out] work     Working array
 **
 *****************************************************************************/
static void st_multigrid_descent(cs_MGS* mgs,
                                 int level,
                                 double* zin,
                                 const double* rhsin,
                                 double* rhsout,
                                 double* work)
{
  cs_MG* mg;

  if (DEBUG) message("Descending from %d to %d\n", level - 1, level);
  mg = mgs->mg[level - 1];
  mg->A->Q->prodVecMatInPlacePtr(zin, work, false);
  for (int icur = 0; icur < mg->nh; icur++)
    work[icur] = rhsin[icur] - work[icur];
  mg->IhH->prodVecMatInPlacePtr(work, rhsout, false);
}

/****************************************************************************/
/*!
 **  Setup the Multigrid system
 **
 ** \return  Error return code
 **
 ** \param[in] mgs    cs_MGS structure
 ** \param[in] qctt    QChol matrix of the upper level
 ** \param[in] flag_sel    If 1, an output selection is created
 ** \param[in] verbose    Verbose flag
 **
 ** \param[out] sel_arg    Vector of selection (in double as used in db)
 **
 *****************************************************************************/
int cs_multigrid_setup(cs_MGS* mgs,
                       QChol* qctt,
                       int flag_sel,
                       int verbose,
                       double** sel_arg)
{
  int *indCo, error, flag_print;
  cs* L;
  double* sel;
  cs_MG *mg, *mgold;

  // Initializations

  error      = 1;
  flag_print = (int)get_keypone("MG_Flag_Print", 0.);
  indCo      = nullptr;
  L          = nullptr;
  sel        = nullptr;
  if (verbose) mestitle(1, "Coarsening %d levels", mgs->nlevels);
  if (flag_print) cs_print_file("QTT_avant", ITEST, qctt->Q->getCS());

  // Define the size of the system

  mgs->ncur = qctt->Q->getNCols();

  // Initialize the Selection vector

  if (flag_sel)
  {
    sel = (double*)mem_alloc(sizeof(double) * mgs->ncur, 0);
    if (sel == nullptr) goto label_end;
    for (int i = 0; i < mgs->ncur; i++)
      sel[i] = 1;
  }

  // Normalize the initial Qctt matrix

  if (mgs->nlevels > 0)
  {
    mgs->diag = csd_extract_diag(qctt->Q->getCS(), -3);
    if (mgs->diag == nullptr) return (1);
    cs* local = cs_normalize_by_diag_and_release(qctt->Q->getCSUnprotected(), 1);
    qctt->Q->freeCS();
    qctt->Q->setCS(local);
    local = cs_spfree(local);
  }
  if (flag_print) cs_print_file("QTT_apres", ITEST, qctt->Q->getCS());

  // Loop on the levels

  for (int ilevel = 0; ilevel <= mgs->nlevels; ilevel++)
  {
    mg = mgs->mg[ilevel];

    // Store A

    if (ilevel == 0)
    {
      mg->A->Q = qctt->Q;
      mg->nh   = mg->A->Q->getNCols();
    }
    else
    {
      mgold    = mgs->mg[ilevel - 1];
      mg->A->Q = prodNormMatMat(mgold->A->Q, mgold->IhH, false);
      if (mg->A->Q == nullptr) goto label_end;
    }
    if (flag_print) cs_print_file("A", ilevel, mg->A->Q->getCS());

    if (ilevel < mgs->nlevels)
    {
      // Extract the coarse grid

      if (cs_coarsening(mg->A->Q->getCS(), mgs->type_coarse, &indCo, &L)) goto label_end;
      if (flag_print) cs_print_file("L", ilevel, L);
      if (flag_print)
        for (int ik = 0; ik < mg->A->Q->getNCols(); ik++)
          message("indco[%d] = %d\n", ik, indCo[ik]);

      // Interpolation

      mg->IhH->setCS(cs_interpolate(mg->A->Q->getCS(), L, indCo));
      if (flag_print) cs_print_file("IhH", ilevel, mg->IhH->getCS());
      delete L;
      L = nullptr;

      // Calculate the sum of rows per column

      mg->sumrow = cs_col_sumrow(mg->IhH->getCS(), &mg->nh, &mg->nH);
      if (mg->sumrow == nullptr) goto label_end;

      // Update the selection

      if (flag_sel) st_selection_update(mgs->ncur, sel, indCo);
      indCo = (int*)mem_free((char*)indCo);
    }

    // Cholesky of the last level

    if (ilevel == mgs->nlevels)
    {
      if (qchol_cholesky(verbose, mg->A)) goto label_end;
    }
  }

  // Optional printout

  if (verbose) cs_multigrid_print(mgs);

  // Set the error return code

  if (flag_sel) *sel_arg = sel;
  error = 0;

label_end:
  if (error) mem_free((char*)sel);
  mem_free((char*)indCo);
  delete L;
  return (error);
}

/****************************************************************************/
/*!
 **  Inversion using Cholesky
 **
 ** \param[in]  qctt     Qchol structure
 ** \param[in,out]  xcr Current vector
 ** \param[in]  rhs     Current R.H.S. vector
 **
 ** \param[out] work     Working array
 **
 *****************************************************************************/
void cs_chol_invert(QChol* qctt, double* xcr, double* rhs, double* work)
{
  if (DEBUG) message("Cholesky Inversion\n");
  int n = qctt->Q->getNCols();
  cs_ipvec(n, qctt->S->Pinv, rhs, work);
  cs_lsolve(qctt->N->L, work);
  cs_ltsolve(qctt->N->L, work);
  cs_pvec(n, qctt->S->Pinv, work, xcr);
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
void cs_chol_simulate(QChol* qctt, double* simu, double* work)
{
  if (DEBUG) message("Cholesky Simulation\n");
  int n = qctt->Q->getNCols();
  cs_ltsolve(qctt->N->L, work);
  cs_pvec(n, qctt->S->Pinv, work, simu);
}

/****************************************************************************/
/*!
 **  Relaxation step
 **
 ** \param[in]  mgs     cs_MGS structure
 ** \param[in]  level     Current level
 ** \param[in]  mode     1 for Descending and -1 for Ascending
 ** \param[in,out]  xcr     Current vector
 ** \param[in]  rhs     Current R.H.S. vector
 **
 ** \param[out] work     Working array
 **
 *****************************************************************************/
static void st_relaxation(cs_MGS* mgs,
                          int level,
                          int mode,
                          double* xcr,
                          const double* rhs,
                          double* work)
{
  cs_MG* mg;

  mg = mgs->mg[level];

  if (mode > 0)
  {

    /* Descending case */

    if (DEBUG) message("Relaxation Descending\n");

    if (level != 0)
      for (int icur = 0; icur < mg->nh; icur++)
        xcr[icur] = 0.;

    for (int i = 0; i < mgs->ngs; i++)
    {
      cs_mulvec_uptri(mg->A->Q->getCS(), mg->nh, xcr, work, 0);
      for (int icur = 0; icur < mg->nh; icur++)
        xcr[icur] = rhs[icur] - work[icur];
      cs_lsolve_lowtri(mg->A->Q->getCS(), xcr, work);
      for (int icur = 0; icur < mg->nh; icur++)
        xcr[icur] = work[icur];
    }
  }
  else
  {

    /* Ascending case */

    if (DEBUG) message("Relaxation Ascending\n");
    for (int i = 0; i < mgs->ngs; i++)
    {
      cs_mulvec_lowtri(mg->A->Q->getCS(), mg->nh, xcr, work, 0);
      for (int icur = 0; icur < mg->nh; icur++)
        xcr[icur] = rhs[icur] - work[icur];
      cs_lsolve_uptri(mg->A->Q->getCS(), xcr, work);
      for (int icur = 0; icur < mg->nh; icur++)
        xcr[icur] = work[icur];
    }
  }
}

/****************************************************************************/
/*!
 **  Iterative phases for the multigrid kriging
 **
 ** \return  Error return code
 **
 ** \param[in]    mgs     cs_MGS structure
 ** \param[in]    verbose  Verbose flag
 ** \param[in,out] x     Input vector
 ** \param[in]    b     R.H.S. vector
 **
 ** \param[out] work     Working array
 **
 ** \remark 'verbose' is passed as argument as it may be different from one
 ** \remark case to the other when calling this function
 **
 *****************************************************************************/
static int st_multigrid_kriging_prec(cs_MGS* mgs,
                                     int verbose,
                                     double* x,
                                     double* b,
                                     double* work)
{
  double *xcr, *rhs, *scores, norm, delta, score;
  int error, nlevels, level, ncur, niter, mode, flag_sym, nfois;

  /* Initializations */

  error = 1;
  xcr = rhs = scores = nullptr;
  ncur               = mgs->ncur;
  nlevels            = mgs->nlevels;
  norm               = VH::innerProduct(b, b, ncur);
  flag_sym           = (int)get_keypone("MG_Flag_Symmetric", 1.);
  if (verbose)
    message("Pre-conditioning phase (Niter=%d Tol=%15.10lf)\n", mgs->nmg,
            mgs->tolnmg);

  /* Core allocation */

  xcr = (double*)mem_alloc(sizeof(double) * ncur * (1 + nlevels), 0);
  if (xcr == nullptr) goto label_end;
  rhs = (double*)mem_alloc(sizeof(double) * ncur * (1 + nlevels), 0);
  if (rhs == nullptr) goto label_end;
  scores = (double*)mem_alloc(sizeof(double) * mgs->nmg, 0);
  if (scores == nullptr) goto label_end;

  /* Initialize the internal arrays */

  for (level = 0; level <= nlevels; level++)
    for (int icur = 0; icur < ncur; icur++)
    {
      XCR(level, icur) = (level == 0) ? x[icur] : 0.;
      RHS(level, icur) = (level == 0) ? b[icur] : 0.;
    }

  /* Loop on the multigrid iterations */

  niter = 0;
  nfois = (flag_sym) ? 2 : 1;

  for (int iter = 0; iter < mgs->nmg; iter++)
  {
    level = 0;

    if (mgs->nlevels > 0)
    {
      /* Loop for symmetrization */

      for (int ifois = 0; ifois < nfois; ifois++)
      {

        /* Loop on the path levels */

        for (int k = 0; k < mgs->npath; k++)
        {
          if (level != nlevels)
          {
            mode = (k > 0) ? mgs->path[k - 1] : 1;
            if (flag_sym) mode = 2 * ifois - 1;
            st_relaxation(mgs, level, mode, &XCR(level, 0), &RHS(level, 0), work);
          }
          else
            cs_chol_invert(mgs->mg[level]->A, &XCR(level, 0), &RHS(level, 0), work);

          level += mgs->path[k];

          if (mgs->path[k] > 0)
            st_multigrid_descent(mgs, level, &XCR(level - 1, 0),
                                 &RHS(level - 1, 0), &RHS(level, 0), work);
          else
            st_multigrid_ascent(mgs, level, 0, 0, &XCR(level + 1, 0),
                                &XCR(level, 0), work);
        }
        mode = (!flag_sym) ? -1 : 1;
        st_relaxation(mgs, 0, mode, &XCR(level, 0), &RHS(level, 0), work);
      }
    }
    else
    {
      cs_chol_invert(mgs->mg[level]->A, &XCR(level, 0), &RHS(level, 0), work);
    }

    // Calculate the score

    score = 0.;
    cs_vector_xM(mgs->mg[0]->A->Q->getCS(), mgs->mg[0]->A->Q->getNCols(), &XCR(0, 0), work);
    for (int icur = 0; icur < ncur; icur++)
    {
      delta = (b[icur] - work[icur]);
      score += delta * delta;
    }
    score           = sqrt(score / norm);
    scores[niter++] = score;
    if (verbose)
      message("Iteration %3d -> Score = %15.10lf\n", iter + 1, score);
    if (score < mgs->tolnmg) break;
  }

  // Save the result

  for (int icur = 0; icur < ncur; icur++)
    x[icur] = XCR(0, icur);

  // Store the scores

  set_keypair("Multigrid_Scores", 1, 1, niter, scores);

  // Set the error return code

  error = 0;

label_end:
  mem_free((char*)xcr);
  mem_free((char*)rhs);
  mem_free((char*)scores);
  return (error);
}

/****************************************************************************/
/*!
 **  Conjugate Gradient algorithm for the multigrid kriging
 **
 ** \return  Error return code
 **
 ** \param[in]    mgs    cs_MGS structure
 ** \param[in]    verbose Verbose flag
 ** \param[in,out] x    Input vector
 ** \param[in]    b    R.H.S. vector
 **
 ** \param[out] work    Working array
 **
 *****************************************************************************/
static int st_multigrid_kriging_cg(cs_MGS* mgs,
                                   int verbose,
                                   double* x,
                                   double* b,
                                   double* work)
{
  cs_MG* mg;
  double *resid, *p, *z, *scores, *temp, sn, s, alpha, beta, norm, score;
  int ncur, error, niter;

  // Initializations

  error = 1;
  ncur  = mgs->ncur;
  norm  = VH::innerProduct(b, b, ncur);
  resid = p = z = scores = temp = nullptr;
  if (verbose)
    message("Conjugate-Gradient Phase (Nmax=%d Tol=%15.10lf)\n", mgs->ngc,
            mgs->tolcg);

  // Core allocation

  p = (double*)mem_alloc(sizeof(double) * ncur, 0);
  if (p == nullptr) goto label_end;
  z = (double*)mem_alloc(sizeof(double) * ncur, 0);
  if (z == nullptr) goto label_end;
  resid = (double*)mem_alloc(sizeof(double) * ncur, 0);
  if (resid == nullptr) goto label_end;
  temp = (double*)mem_alloc(sizeof(double) * ncur, 0);
  if (temp == nullptr) goto label_end;
  scores = (double*)mem_alloc(sizeof(double) * mgs->ngc, 0);
  if (scores == nullptr) goto label_end;

  // Calculate initial residual

  mg = mgs->mg[0];
  cs_vector_xM(mg->A->Q->getCS(), mg->nh, x, work);
  for (int icur = 0; icur < mg->nh; icur++)
    resid[icur] = b[icur] - work[icur];
  for (int icur = 0; icur < ncur; icur++)
    z[icur] = 0.;
  st_multigrid_kriging_prec(mgs, 0, z, resid, work);
  matrix_product_safe(1, ncur, 1, resid, z, &sn);
  for (int icur = 0; icur < ncur; icur++)
    p[icur] = z[icur];

  // Loop

  niter = 0;
  for (int iter = 0; iter < mgs->ngc; iter++)
  {
    s = sn;
    cs_vector_xM(mg->A->Q->getCS(), mg->nh, p, temp);
    matrix_product_safe(1, ncur, 1, p, temp, &alpha);
    alpha = s / alpha;

    for (int icur = 0; icur < ncur; icur++)
    {
      x[icur] += alpha * p[icur];
      resid[icur] -= alpha * temp[icur];
    }

    score           = sqrt(VH::innerProduct(resid, resid, ncur) / norm);
    scores[niter++] = score;
    if (verbose)
      message("Iteration Gradient %3d -> Score = %15.10lf\n", iter + 1, score);
    if (score < mgs->tolcg) break;

    for (int icur = 0; icur < ncur; icur++)
      z[icur] = 0.;
    st_multigrid_kriging_prec(mgs, 0, z, resid, work);

    matrix_product_safe(1, ncur, 1, resid, z, &sn);
    matrix_product_safe(1, ncur, 1, temp, z, &beta);
    beta *= -alpha / s;
    for (int icur = 0; icur < ncur; icur++)
      p[icur] = z[icur] + beta * p[icur];
  }

  // Store the scores

  set_keypair("Multigrid_Gradient_Scores", 1, 1, niter, scores);

  // Set the error return code

  error = 0;

label_end:
  mem_free((char*)p);
  mem_free((char*)z);
  mem_free((char*)resid);
  mem_free((char*)temp);
  mem_free((char*)scores);
  return (error);
}

/****************************************************************************/
/*!
 **  Perform the multigrid kriging
 **
 ** \return  Error return code
 **
 ** \param[in]  mgs     cs_MGS structure
 ** \param[in]  qctt     QChol matrix of the upper level
 ** \param[in]  verbose     Verbose flag
 ** \param[in]  x0     Input vector
 ** \param[in]  b     R.H.S. vector
 **
 ** \param[out] work     Working array
 **
 *****************************************************************************/
int cs_multigrid_process(cs_MGS* mgs,
                         QChol* qctt,
                         int verbose,
                         double* x0,
                         double* b,
                         double* work)
{
  // Perform the setup (if not already done)

  if (mgs->diag == nullptr)
  {
    if (cs_multigrid_setup(mgs, qctt, 0, verbose, NULL)) return (1);
  }
  else
  {
    if (mgs->ncur != qctt->Q->getNCols())
      messageAbort("Check that multigrid has been setup up correctly");
  }

  // Initial normalization

  st_multigrid_scale(mgs, 1, x0, b);

  // Dispatch

  if (mgs->flag_cg)
  {

    /* Multigrid using Conjugate-gradient */

    if (st_multigrid_kriging_cg(mgs, verbose, x0, b, work)) return (1);
  }
  else
  {

    /* Iterative multigrid processing */

    if (st_multigrid_kriging_prec(mgs, verbose, x0, b, work)) return (1);
  }

  // Denormalize the results

  st_multigrid_scale(mgs, -1, x0, NULL);

  return (0);
}

NF_Triplet csToTriplet(const cs* A, int shiftRow, int shiftCol, double tol)
{
  NF_Triplet NF_T;

  if (A == nullptr) return NF_T;
  int ncols = cs_getncol(A);

  cs* AT = cs_transpose(A, 1);
  if (AT != nullptr)
  {
    cs_spfree(AT);
  }

  int* Ap    = A->p;
  int* Ai    = A->i;
  double* Ax = A->x;
  int nz     = A->nz;
  if (nz >= 0) return NF_T;

  for (int j = 0; j < ncols; j++)
    for (int p = Ap[j]; p < Ap[j + 1]; p++)
    {
      double value = (Ax) ? Ax[p] : 1.;
      if (ABS(value) <= tol) continue;
      NF_T.add(Ai[p] + shiftRow, j + shiftCol, value);
    }
  return NF_T;
}

int cs_get_nrow(const cs* A)
{
  if (A == nullptr) return 0;
  int nrow = cs_getnrow(A);
  return (nrow);
}

int cs_get_ncol(const cs* A)
{
  if (A == nullptr) return 0;
  int ncol = cs_getncol(A);
  return (ncol);
}

void cs_print_dim(const char* title, const cs* A)
{
  if (A == nullptr) return;
  message("%s: Nrow=%d Ncol=%d NNZ=%d\n", title, cs_getnrow(A), cs_getncol(A), A->nzmax);
}

void cs_print_range(const char* title, const cs* C)
{
  if (C == nullptr) return;

  /* Convert the contents of the sparse matrix into triplets (eigen formet) */

  NF_Triplet NF_T = csToTriplet(C);

  /* Calculate the extreme values */

  StatResults stats = ut_statistics(NF_T.getNElements(), NF_T.getValues().data());

  /* Printout */

  if (title != NULL)
    message("%s\n", title);
  else
    message("Sparse matrix\n");
  message(" Descr: m=%d n=%d nnz=%d\n", cs_getnrow(C), cs_getncol(C), C->nzmax);
  if (NF_T.getNElements() > 0)
    message(" Range: [%lf ; %lf] (%d/%d)\n", stats.mini, stats.maxi, stats.nvalid, NF_T.getNElements());
  else
    message(" All terms are set to zero\n");
}

/* Construct the sparse diagonal matrix */
cs* cs_eye(int number, double value)
{
  cs *Atriplet, *A;

  /* Initializations */

  A = nullptr;

  /* Fill the new sparse triplet */

  Atriplet = cs_spalloc(0, 0, 1, 1, 1);
  if (Atriplet == nullptr) goto label_end;
  for (int i = 0; i < number; i++)
  {
    if (!cs_entry(Atriplet, i, i, value)) goto label_end;
  }

  A = cs_triplet(Atriplet);

label_end:
  cs_spfree(Atriplet);
  return (A);
}

// Build a sparse matrix containing the diagonal of a sparse matrix
// oper_choice: Operation on the diagonal term (see Utilities::operate_XXX)
cs* cs_extract_diag(const cs* C, int oper_choice)
{
  operate_function oper_func = operate_Identify(oper_choice);
  cs* A                      = nullptr;
  cs* Atriplet               = cs_spalloc(0, 0, 1, 1, 1);
  int* Cp                    = C->p;
  int* Ci                    = C->i;
  double* Cx                 = C->x;
  if (Atriplet == nullptr) goto label_end;

  /* Loop on the rows */

  for (int j = 0; j < cs_getncol(C); j++)
  {
    for (int p = Cp[j]; p < Cp[j + 1]; p++)
    {
      if (Ci[p] != j) continue;
      if (!cs_entry(Atriplet, j, j, oper_func(Cx[p]))) goto label_end;
    }
  }
  A = cs_triplet(Atriplet);

label_end:
  cs_spfree(Atriplet);
  return (A);
}

// Extract the (transformed) diagonal of a sparse matrix
// oper_choice: Operation on the diagonal term (see Utilities::operate_XXX)
double* csd_extract_diag(const cs* C, int oper_choice)
{
  operate_function oper_func = operate_Identify(oper_choice);
  double* diag               = nullptr;
  int size                   = cs_getncol(C);
  int* Cp                    = C->p;
  int* Ci                    = C->i;
  double* Cx                 = C->x;

  /* Core allocation */

  diag = (double*)mem_alloc(sizeof(double) * size, 0);
  if (diag == nullptr) goto label_end;
  for (int i = 0; i < size; i++)
    diag[i] = 0.;

  /* Loop on the rows */

  for (int j = 0; j < cs_getncol(C); j++)
  {
    for (int p = Cp[j]; p < Cp[j + 1]; p++)
    {
      if (Ci[p] != j) continue;
      diag[j] = oper_func(Cx[p]);
    }
  }

label_end:
  return (diag);
}

/* Return the number of rows and columns */
/* as well as the percentage of filled terms */
void cs_rowcol(const cs* A, int* nrows, int* ncols, int* count, double* percent)
{
  cs* AT;
  int* Ap;
  double* Ax;

  (*nrows) = (*ncols) = (*count) = 0;
  (*percent)                     = 0.;
  if (!A) return;

  Ap = A->p;
  Ax = A->x;
  if (A->nz >= 0) return;

  for (int j = 0; j < cs_getncol(A); j++)
    for (int p = Ap[j]; p < Ap[j + 1]; p++)
    {
      if (ABS(Ax[p]) > 0) (*count)++;
    }

  *ncols = cs_getncol(A);
  AT     = cs_transpose(A, 1);
  *nrows = cs_getncol(AT);
  cs_spfree(AT);

  if ((*nrows) > 0 && (*ncols) > 0)
    (*percent) = ((100. * (double)(*count)) / ((double)(*nrows) * (double)(*ncols)));
}

/* Print a nice sparse matrix */
/* The format is copied from Matrix package */
void cs_print_nice(const char* title, const cs* A, int maxrow, int maxcol)
{
  int p, j, m, n, *Ap, *Ai, npass, jdeb, jfin, found, num_line;
  double* Ax;
  int nbypass = 7;

  if (!A)
  {
    message("(null)\n");
    return;
  }
  m  = cs_getnrow(A);
  n  = cs_getncol(A);
  Ap = A->p;
  Ai = A->i;
  Ax = A->x;
  if (A->nz >= 0) return;
  if (maxcol >= 0) n = maxcol;
  if (maxrow >= 0) m = maxrow;

  npass = (int)ceil((double)n / (double)nbypass);

  /* Print the title (optional) */

  if (title != NULL)
    message("%s", title);
  else
    message("Print Sparse Matrix");
  if (maxrow >= 0) message(" nrows<=%d", maxrow);
  if (maxcol >= 0) message(" ncols<=%d", maxcol);
  message("\n");

  /* Loop on the passes */

  for (int ipass = 0; ipass < npass; ipass++)
  {
    jdeb = ipass * nbypass;
    jfin = MIN(jdeb + nbypass, n);

    /* Title of the columns */

    message("      ");
    for (j = jdeb; j < jfin; j++)
      message("    [,%3d]", j + 1);
    message("\n");

    /* Loop on the lines */

    for (int i = 0; i < m; i++)
    {
      message("[%3d,] ", i + 1);

      /* Loop on the columns */

      for (j = jdeb; j < jfin; j++)
      {

        /* Search for the correct line number */

        found = -1;
        for (p = Ap[j]; p < Ap[j + 1] && found < 0; p++)
        {
          num_line = Ai[p];
          if (num_line == i) found = p;
        }

        if (found < 0)
          message(" .        ");
        else
          message("%9.4lf ", Ax[found]);
      }
      message("\n");
    }
    message("\n");
  }
}

cs* cs_duplicate(const cs* b1)
{
  cs* bres;
  bres = cs_add(b1, b1, 1., 0.);
  return (bres);
}

cs* cs_prod_norm_and_release(cs* b1, cs* lambda, int flag_release)
{
  cs *bres1, *bres2;

  bres1 = bres2 = nullptr;

  bres1 = cs_multiply(lambda, b1);
  if (bres1 == nullptr) goto label_end;
  bres2 = cs_multiply(bres1, lambda);
  if (bres2 == nullptr) goto label_end;
  if (flag_release) cs_spfree(b1);

label_end:
  cs_spfree(bres1);
  return (bres2);
}

/* Return per column, the sum along the row */
/* Note: Check existence of returned argument */
/* Note: The returned array must be freed by calling function */
double* cs_col_sumrow(const cs* A, int* ncol, int* nrow)
{
  int *Ap, error;
  cs* AT; // Will store t(A)
  double *Ax, *vect, sumval;

  error = 1;
  *ncol = *nrow = 0;
  AT            = nullptr;
  vect          = nullptr;
  if (!A) goto label_end;
  if (A->nz >= 0) goto label_end;
  Ap = A->p;
  Ax = A->x;

  /* Transpose A into AT to follow ordering per column */

  AT = cs_transpose(A, 1);
  if (AT == nullptr) goto label_end;

  /* Core allocation */

  vect = (double*)mem_alloc(sizeof(double) * cs_getncol(A), 0);
  if (vect == nullptr) goto label_end;

  /* Loop on the rows */

  for (int j = 0; j < cs_getncol(A); j++)
  {
    sumval = 0.;
    for (int p = Ap[j]; p < Ap[j + 1]; p++)
      sumval += Ax[p];
    vect[j] = sumval;
  }

  *ncol = cs_getncol(A);
  *nrow = cs_getncol(AT);
  error = 0;

label_end:
  cs_spfree(AT);
  if (error) vect = (double*)mem_free((char*)vect);
  return (vect);
}

/* Operate the product of a vector by a sparse matrix */
/* y = x %*% A */
void cs_vector_xM(const cs* A, int nout, const double* x, double* y)
{
  int *Ap, *Ai, n;
  double* Ax;

  n  = cs_getncol(A);
  Ap = A->p;
  Ai = A->i;
  Ax = A->x;

  for (int j = 0; j < nout; j++)
    y[j] = 0.;

  for (int j = 0; j < n; j++)
    for (int p = Ap[j]; p < Ap[j + 1]; p++)
      y[j] += x[Ai[p]] * Ax[p];
}

cs* cs_normalize_by_diag_and_release(cs* Q, int flag_release)
{
  cs *diag, *Qp;

  Qp = diag = nullptr;

  diag = cs_extract_diag(Q, -3);
  if (diag == nullptr) goto label_end;
  Qp = cs_prod_norm_and_release(Q, diag, flag_release);
  if (Qp == nullptr) goto label_end;

label_end:
  cs_spfree(diag);
  return (Qp);
}

/* Operate norm-product of sparse matrix A by diagonal matrix X */
/* (entered as a vector): B = oper(X) %*% A %*% oper(X) */
/* oper_choice: Operation on the diagonal term (see Utilities::operate_XXX) */
cs* cs_matvecnorm(const cs* A, const double* x, int oper_choice)
{
  operate_function oper_func = operate_Identify(oper_choice);
  cs* B                      = cs_duplicate(A);
  if (B == nullptr) return (B);
  int n      = cs_getncol(A);
  int* Ap    = A->p;
  int* Ai    = A->i;
  double* Ax = A->x;
  double* Bx = B->x;

  for (int j = 0; j < n; j++)
    for (int p = Ap[j]; p < Ap[j + 1]; p++)
      Bx[p] = Ax[p] * oper_func(x[j]) * oper_func(x[Ai[p]]);
  return (B);
}

/* Same of cs_matvecnorm ... but in place                                   */
/* oper_choice: Operation on the diagonal term (see Utilities::operate_XXX) */
void cs_matvecnorm_inplace(cs* A, const double* x, int oper_choice)
{
  operate_function oper_func = operate_Identify(oper_choice);
  int n                      = cs_getncol(A);
  int* Ap                    = A->p;
  int* Ai                    = A->i;
  double* Ax                 = A->x;

  for (int j = 0; j < n; j++)
    for (int p = Ap[j]; p < Ap[j + 1]; p++)
      Ax[p] = Ax[p] * oper_func(x[j]) * oper_func(x[Ai[p]]);
}

static void st_get_FiCo(cs* L,
                        cs* Lt,
                        int* lambda,
                        int* indUd,
                        int* indFi,
                        int* indCo)
{
  int n, nU, icol, newCo, *Ltp, *Lti, *LLp, *LLi, nFi, newFi;
  double* Ltx;

  std::map<int, std::set<int>> tab;
  std::map<int, std::set<int>>::reverse_iterator rit;
  std::map<int, std::set<int>>::iterator it;
  std::set<int>* vec;

  /* Initializations */

  Ltp = L->p;
  Ltx = L->x;
  Lti = L->i;
  LLp = Lt->p;
  LLi = Lt->i;

  n  = cs_getncol(L);
  nU = 0;
  for (int i = 0; i < n; i++) nU += indUd[i];

  for (int i = 0; i < n; i++)
  {
    it = tab.find(lambda[i]);
    if (it == tab.end()) tab[lambda[i]] = std::set<int>();
    tab[lambda[i]].insert(i);
  }

  // Update indices

  while (nU > 0)
  {
    // Find the sample mostly connected

    rit          = tab.rbegin();
    newCo        = *rit->second.begin();
    indCo[newCo] = 1;
    indUd[newCo] = 0;
    nU--;
    vec = &rit->second;
    vec->erase(newCo);
    if (vec->empty()) tab.erase(lambda[newCo]);

    // Add the samples to the fine set

    nFi = 0;
    for (int p = Ltp[newCo]; p < Ltp[newCo + 1]; p++)
    {
      icol = Lti[p];
      if (Ltx[p] == 1 && indUd[icol] == 1)
      {
        indFi[nFi++] = icol;
        indCo[icol]  = 0;
        indUd[icol]  = 0;
        nU--;
        it = tab.find(lambda[icol]);
        if (it == tab.end())
        {
          messageAbort("%d not found for %d (%d)", lambda[icol], icol, nU);
        }
        vec = &it->second;
        vec->erase(icol);
        if (vec->empty()) tab.erase(lambda[icol]);
      }
    }

    // Update Lambda (due to newCo)

    for (int p = LLp[newCo]; p < LLp[newCo + 1]; p++)
    {
      if (indUd[LLi[p]])
      {
        it = tab.find(lambda[LLi[p]]);
        if (it == tab.end())
        {
          messageAbort("Fi : %lf not found for %d", lambda[LLi[p]], LLi[p]);
        }
        vec = &it->second;
        vec->erase(LLi[p]);
        if (vec->empty()) tab.erase(lambda[LLi[p]]);
        lambda[LLi[p]]--;
        it = tab.find(lambda[LLi[p]]);
        if (it == tab.end()) // means lambda[LLi[p]] is not in tab as a key.
          tab[lambda[LLi[p]]] = std::set<int>();
        tab[lambda[LLi[p]]].insert(LLi[p]);
      }
    }

    // Update Lambda (due to newFi)

    if (nFi > 0)
    {
      for (int k = 0; k < nFi; k++)
      {
        newFi = indFi[k];
        for (int p = LLp[newFi]; p < LLp[newFi + 1]; p++)
        {
          if (indUd[LLi[p]])
          {
            it  = tab.find(lambda[LLi[p]]);
            vec = &it->second;
            vec->erase(LLi[p]);
            if (vec->empty()) tab.erase(lambda[LLi[p]]);
            lambda[LLi[p]]++;
            it = tab.find(lambda[LLi[p]]);
            if (it == tab.end()) tab[lambda[LLi[p]]] = std::set<int>();
            tab[lambda[LLi[p]]].insert(LLi[p]);
          }
        }
      }
    }
  }
}

static int st_update_neigh(int* n_arg, int indloc, int* n_tab, int* r_tab)
{
  int n;

  n = *n_arg;

  // Look for already registered index value

  for (int i = 0; i < n; i++)
  {
    if (r_tab[i] != indloc) continue;
    n_tab[i]++;
    return (0);
  }

  // New index to be registered

  if (n > MAX_NEIGH) return (1);
  r_tab[n] = indloc;
  n_tab[n] = 1;
  n++;
  *n_arg = n;
  return (0);
}

static int st_coarse_type0(const cs* Q,
                           int* indUd,
                           int* indFi,
                           int* indCo,
                           cs** Lret,
                           cs** Ltret)
{
  int *lambda, *Qp, *Qi, n, ip, error;
  cs *L, *Lt, *Ltriplet;
  double *Qx, u, minu;
  double eps = 0.25;

  /* Initializations */

  error = 1;
  Qp = Qi = lambda = nullptr;
  L = Lt = Ltriplet = nullptr;
  Qx                = nullptr;
  n                 = cs_getncol(Q);
  Qp                = Q->p;
  Qx                = Q->x;
  Qi                = Q->i;

  /* Core allocation */

  Ltriplet = cs_spalloc(0, 0, 1, 1, 1);
  if (Ltriplet == nullptr) goto label_end;
  lambda = (int*)mem_alloc(sizeof(int) * n, 0);
  if (lambda == nullptr) goto label_end;

  for (int i = 0; i < n; i++)
    lambda[i] = 0;

  /* Building Ltriplet sparse matrix */

  for (int j = 0; j < n; j++)
  {
    minu = 0.;
    for (int p = Qp[j]; p < Qp[j + 1]; p++)
    {
      u = Qx[p];
      if (u < minu) minu = u;
    }
    for (int p = Qp[j]; p < Qp[j + 1]; p++)
    {
      u = Qx[p];
      if (u < eps * minu)
      {
        ip = Qi[p];
        if (!cs_entry(Ltriplet, j, ip, 1.)) goto label_end;
        lambda[ip]++;
      }
    }
  }

  /* Convert from triplet to sparse matrix */

  L        = cs_triplet(Ltriplet);
  Ltriplet = cs_spfree(Ltriplet);

  // Transpose L

  Lt = cs_transpose(L, 1);
  if (Lt == nullptr) goto label_end;

  // Transform triplet into sparse matrix

  st_get_FiCo(L, Lt, lambda, indUd, indFi, indCo);

  // Set the error return code

  error  = 0;
  *Lret  = L;
  *Ltret = Lt;

label_end:
  mem_free((char*)lambda);
  cs_spfree(Ltriplet);
  if (error)
  {
    cs_spfree(L);
    cs_spfree(Lt);
  }
  return (error);
}

static int st_coarse_typen(cs* /*L*/,
                           cs* Lt,
                           int type,
                           int* indUd,
                           int* indFi,
                           int* indCo)
{
  cs *Lout, *Loutt, *Ltriplet;
  double* Lx;
  int *lambda, *Lp, *Li, n, ip, iq, error, n_neigh, nb;
  int n_neigh_tab[MAX_NEIGH], r_neigh_tab[MAX_NEIGH];

  /* Initializations */

  error = 1;
  Lp = Li = lambda = nullptr;
  Lout = Loutt = Ltriplet = nullptr;
  Lx                      = nullptr;
  n                       = cs_getncol(Lt);
  Lp                      = Lt->p;
  Lx                      = Lt->x;
  Li                      = Lt->i;

  /* Core allocation */

  Ltriplet = cs_spalloc(0, 0, 1, 1, 1);
  if (Ltriplet == nullptr) goto label_end;
  lambda = (int*)mem_alloc(sizeof(int) * n, 0);
  if (lambda == nullptr) goto label_end;

  for (int i = 0; i < n; i++)
    lambda[i] = 0;

  /* Building Ltriplet sparse matrix */

  for (int j = 0; j < n; j++)
  {
    // Skip if 'j' is fine
    if (indCo[j] != 1) continue;
    n_neigh = 0;

    // Loop on the nodes connected to row 'j'
    for (int p = Lp[j]; p < Lp[j + 1]; p++)
    {
      // Skip the non connected node
      if (Lx[p] == 0) continue;
      // Get the column index
      ip = Li[p];

      // Skip if 'ip' is coarse
      if (indCo[ip] == 1) continue;

      // Loop on the neighbors of 'iq'
      for (int q = Lp[ip]; q < Lp[ip + 1]; q++)
      {
        // Skip the non connected node
        if (Lx[q] == 0) continue;
        iq = Li[q];
        if (iq == j) continue;
        if (indCo[iq] == 0) continue;
        if (st_update_neigh(&n_neigh, iq, n_neigh_tab, r_neigh_tab))
          goto label_end;
      }
    }

    for (int k = 0; k < n_neigh; k++)
    {
      iq = r_neigh_tab[k];
      nb = n_neigh_tab[k];
      if (nb >= type)
      {
        if (!cs_entry(Ltriplet, j, iq, 1.)) goto label_end;
        lambda[iq]++;
      }
    }
  }

  /* Convert from triplet to sparse matrix */

  Lout     = cs_triplet(Ltriplet);
  Ltriplet = cs_spfree(Ltriplet);

  // Transpose

  Loutt = cs_transpose(Lout, 1);

  // Transform triplet into sparse matrix

  st_get_FiCo(Lout, Loutt, lambda, indUd, indFi, indCo);

  /* Set the error return code */

  error = 0;

label_end:
  mem_free((char*)lambda);
  cs_spfree(Lout);
  cs_spfree(Loutt);
  cs_spfree(Ltriplet);
  return (error);
}

//
// L:    List of samples strongly connected:
//    for i Qij<0 and ABS(Qij)>eps * max(Qij<0)
// type: 0 for standard coarsening
//    1 for aggressive (type 1)
//    2 for aggressive (type 2)
//
int cs_coarsening(const cs* Q, int type, int** indCo_ret, cs** L_ret)
{
  int *indUd, *indCo, *indFi;
  int n, error;
  cs *L, *Lt;

  // Initializations

  error = 1;
  L = Lt = nullptr;
  indUd = indCo = indFi = nullptr;
  n                     = cs_getncol(Q);

  // Core allocation

  indUd = (int*)mem_alloc(sizeof(int) * n, 0);
  if (indUd == nullptr) goto label_end;
  indCo = (int*)mem_alloc(sizeof(int) * n, 0);
  if (indCo == nullptr) goto label_end;
  indFi = (int*)mem_alloc(sizeof(int) * n, 0);
  if (indFi == nullptr) goto label_end;

  // Construct L (for type = 0)

  for (int i = 0; i < n; i++)
  {
    indUd[i] = 1;
    indFi[i] = 0;
    indCo[i] = 0;
  }
  if (st_coarse_type0(Q, indUd, indFi, indCo, &L, &Lt)) goto label_end;

  if (type == 0)
  {
    error = 0;
    goto label_end;
  }

  /* Construct L for aggressive coarsening */

  for (int i = 0; i < n; i++)
    indUd[i] = indCo[i];

  if (st_coarse_typen(L, Lt, type, indUd, indFi, indCo)) goto label_end;

  /* Set the error return code */

  error = 0;

label_end:
  mem_free((char*)indUd);
  mem_free((char*)indFi);
  cs_spfree(L);
  if (error)
  {
    indCo = (int*)mem_free((char*)indCo);
    Lt    = cs_spfree(Lt);
  }
  *indCo_ret = indCo;
  *L_ret     = Lt;
  return (error);
}

cs* cs_interpolate(const cs* AA, const cs* Lt, const int* Co)
{
  cs *IH, *IHtriplet;
  double *u, *AAx, *Ltx, sunip, sunim, supim, supip, alpha, beta, fact, val;
  int *Nip, *Pip, *AAp, *AAi, *Ltp, *Lti, n, error, rCo, ip, *indCo, s;

  // Initialization
  error = 1;
  IH = IHtriplet = nullptr;
  u              = nullptr;
  Nip = Pip = nullptr;
  AAp       = AA->p;
  AAx       = AA->x;
  AAi       = AA->i;
  n         = cs_getncol(AA);
  indCo     = nullptr;

  // Core allocation */
  u = (double*)mem_alloc(sizeof(double) * n, 0);
  if (u == nullptr) goto label_end;
  indCo = (int*)mem_alloc(sizeof(int) * n, 0);
  if (indCo == nullptr) goto label_end;

  Nip = (int*)mem_alloc(sizeof(int) * n, 0);
  if (Nip == nullptr) goto label_end;
  Pip = (int*)mem_alloc(sizeof(int) * n, 0);
  if (Pip == nullptr) goto label_end;

  // Transpose
  IHtriplet = cs_spalloc(0, 0, 1, 1, 1);
  if (IHtriplet == nullptr) goto label_end;
  Ltp = Lt->p;
  Ltx = Lt->x;
  Lti = Lt->i;

  s = 0;
  for (int i = 0; i < n; i++)
  {
    indCo[i] = 0;
    if (!Co[i]) continue;
    indCo[i] = s;
    s++;
  }

  // Loop on the Fine elements (not Coarse)

  for (int i = 0; i < n; i++)
  {
    if (Co[i] == 1) continue;

    // Initializations
    sunim = supim = sunip = supip = 0.;

    // Connected neighbors of 'i' ('i' excluded)
    for (int p = AAp[i]; p < AAp[i + 1]; p++)
    {
      ip      = AAi[p];
      Nip[ip] = 0;
      val = u[ip] = AAx[p];
      if (val == 0.) continue;
      if (ip != i)
      {
        if (val > 0)
        {
          Nip[ip] = 1;
          sunip += val;
        }
        else
        {
          sunim += val;
        }
      }
    }

    // Strongly connected neighbors on Coarse
    for (int p = AAp[i]; p < AAp[i + 1]; p++)
    {
      ip      = AAi[p];
      Pip[ip] = 0;
    }

    for (int p = Ltp[i]; p < Ltp[i + 1]; p++)
    {
      ip  = Lti[p];
      val = u[ip];

      if (val == 0.) continue;
      if (Co[ip] && Ltx[p] != 0.)
      {
        if (val > 0)
        {
          Pip[ip] = 1;
          supip += val;
        }
        else
        {
          Pip[ip] = -1;
          supim += val;
        }
      }
    }

    alpha = sunim / supim;
    fact  = u[i];
    if (supip > 0.)
    {
      sunip = supip = 0;
      for (int p = AAp[i]; p < AAp[i + 1]; p++)
      {
        ip = AAi[p];
        if (Nip[ip]) sunip += AAx[p];
        if (Pip[ip] > 0) supip += AAx[p];
      }
      beta = sunip / supip;
    }
    else
    {
      fact += sunip;
      beta = 0.;
    }

    // Add the entries relative to 'Pip'

    for (int p = AAp[i]; p < AAp[i + 1]; p++)
    {
      ip = AAi[p];
      if (!Co[ip]) continue;
      if (Pip[ip] != 0.)
      {
        val = (Pip[ip] > 0) ? beta : alpha;
        val *= -u[ip] / fact;
        if (!cs_entry(IHtriplet, indCo[ip], i, val)) goto label_end;
      }
    }
  }

  // Add the entries corresponding to the diagonal

  rCo = -1;
  for (int p = 0; p < n; p++)
  {
    if (!Co[p]) continue;
    rCo++;
    if (!cs_entry(IHtriplet, rCo, p, 1.)) goto label_end;
  }

  // Transform from triplet to sparse
  IH        = cs_triplet(IHtriplet);
  IHtriplet = cs_spfree(IHtriplet);

  // Set the error return code
  error = 0;

label_end:
  mem_free((char*)u);
  mem_free((char*)indCo);
  mem_free((char*)Nip);
  mem_free((char*)Pip);
  cs_spfree(IHtriplet);
  if (error) IH = cs_spfree(IH);
  return (IH);
}

int cs_lsolve_lowtri(const cs* L, const double* x, double* y)
/*
 Purpose:

 CS_LSOLVE solves LowerTri(L)*y=x
 Where LowerTri(L) is the lower triangular part of L (symmetric)
 */
{
  int j, n, p, *Lp, *Li;
  double *Lx, pivot, value;
  if (!L || !x) return (0); /* check inputs */
  n  = cs_getncol(L);
  Lp = L->p;
  Li = L->i;
  Lx = L->x;

  for (int i = 0; i < n; i++)
  {
    value = x[i];
    pivot = 0;
    for (p = Lp[i]; p < Lp[i + 1]; p++)
    {
      j = Li[p];
      if (j < i)
        value -= Lx[p] * y[j];
      else if (i == j)
        pivot = Lx[p];
    }

    y[i] = value / pivot;
  }
  return (1);
}

int cs_lsolve_uptri(const cs* L, const double* x, double* y)
/*
 Purpose:

 CS_LSOLVE solves UpperTri(L)*y=x
 Where UpperTri(L) is the upper triangular part of L (symmetric)
 */
{
  int j, n, p, *Lp, *Li;
  double *Lx, pivot, value;
  if (!L || !x) return (0); /* check inputs */
  n  = cs_getncol(L);
  Lp = L->p;
  Li = L->i;
  Lx = L->x;

  for (int i = n - 1; i >= 0; i--)
  {
    value = x[i];
    pivot = 0;
    for (p = Lp[i]; p < Lp[i + 1]; p++)
    {
      j = Li[p];
      if (j > i)
        value -= Lx[p] * y[j];
      else if (i == j)
        pivot = Lx[p];
    }

    y[i] = value / pivot;
  }
  return (1);
}

/* Operate product of a vector by triangular upper part of sparse matrix */
/* y = A %*% x */
void cs_mulvec_uptri(const cs* A,
                     int nout,
                     const double* x,
                     double* y,
                     int flag_diag)
{
  int *Ap, *Ai, i, n;
  double *Ax, value;

  n  = cs_getncol(A);
  Ap = A->p;
  Ai = A->i;
  Ax = A->x;

  for (int j = 0; j < nout; j++)
    y[j] = 0.;

  for (int j = 0; j < n; j++)
  {
    value = x[j];
    for (int p = Ap[j]; p < Ap[j + 1]; p++)
    {
      i = Ai[p];
      if (!flag_diag && i == j) continue;
      if (j < i) continue;
      y[i] += Ax[p] * value;
    }
  }
}

/* Operate product of a vector by triangular lower part of sparse matrix */
/* y = A %*% x */
void cs_mulvec_lowtri(const cs* A,
                      int nout,
                      const double* x,
                      double* y,
                      int flag_diag)
{
  int *Ap, *Ai, i, n;
  double *Ax, value;

  n  = cs_getncol(A);
  Ap = A->p;
  Ai = A->i;
  Ax = A->x;

  for (int j = 0; j < nout; j++)
    y[j] = 0.;

  for (int j = 0; j < n; j++)
  {
    value = x[j];
    for (int p = Ap[j]; p < Ap[j + 1]; p++)
    {
      i = Ai[p];
      if (!flag_diag && i == j) continue;
      if (j > i) continue;
      y[i] += Ax[p] * value;
    }
  }
}

void cs_print_file(const char* radix, int rank, const cs* A)
{
  FILE* file;
  char filename[100];

  if (A == nullptr) return;

  if (!IFFFF(rank))
    (void)gslSPrintf(filename, "%s-%d", radix, rank);
  else
    (void)gslStrcpy(filename, radix);

  file = gslFopen(filename, "w");
  if (file == nullptr) return;

  NF_Triplet NF_T = csToTriplet(A);

  for (int i = 0; i < NF_T.getNElements(); i++)
  {
    fprintf(file, "%10d %10d %20.10lf\n", NF_T.getCol(i), NF_T.getRow(i), NF_T.getValue(i));
  }

  (void)fclose(file);
}

} // namespace gstlrn