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

#include "LinearOp/InvNuggetOp.hpp"
#include "Basic/AStringable.hpp"
#include "Covariances/CovAniso.hpp"
#include "Db/Db.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "Matrix/NF_Triplet.hpp"
#include "Model/Model.hpp"
#include "geoslib_define.h"

namespace gstlrn
{

  InvNuggetOp::InvNuggetOp(Db* dbin, Model* model, const SPDEParam& params)
    : ASimulable()
    , _invNuggetMatrix(nullptr)
  {
    _buildInvNugget(dbin, model, params);
  }

  int InvNuggetOp::getSize() const
  {
    return _invNuggetMatrix ? _invNuggetMatrix->getNRows() : 0;
  }

  InvNuggetOp::~InvNuggetOp()
  {
    
  }

  int InvNuggetOp::_addSimulateToDest(const constvect whitenoise, vect outv) const
  {
   DECLARE_UNUSED(whitenoise,outv)
   messerr("Not implemented: InvNuggetOp::_addSimulateToDest");
   return 1;
  }

  int InvNuggetOp::_addToDest(const constvect inv, vect outv) const
  {
    return _invNuggetMatrix->addToDest(inv, outv);
  }

  static void _addVerrConstant(MatrixSymmetric& sills, const VectorDouble& verrDef)
  {
    int nverr = (int)verrDef.size();
    if (nverr > 0)
    {
      for (int iverr = 0; iverr < nverr; iverr++)
        sills.updValue(iverr, iverr, EOperator::ADD, verrDef[iverr]);
    }
  }

  static void _checkMinNugget(MatrixSymmetric& sills, const VectorDouble& minNug)
  {
    int nvar = (int)minNug.size();

    // Check that the diagonal of the Sill matrix is large enough
    for (int ivar = 0; ivar < nvar; ivar++)
      sills.setValue(ivar, ivar, MAX(sills.getValue(ivar, ivar), minNug[ivar]));
  }

  static MatrixSymmetric _buildSillPartialMatrix(const MatrixSymmetric& sillsRef,
                                                 int nvar,
                                                 int ndef,
                                                 const VectorInt& identity)
  {
    MatrixSymmetric sills;
    if (ndef == nvar)
      sills = sillsRef;
    else
    {
      sills = MatrixSymmetric(ndef);
      for (int idef = 0; idef < ndef; idef++)
        for (int jdef = 0; jdef <= idef; jdef++)
          sills.setValue(idef, jdef, sillsRef.getValue(identity[idef], identity[jdef]));
    }
    return sills;
  }

  static int _loadPositions(int iech,
                            const VectorVectorInt& index1,
                            const VectorInt& cumul,
                            VectorInt& positions,
                            VectorInt& identity,
                            int* rank_arg)
  {
    int nvar = (int)cumul.size();
    int ndef = 0;
    int rank = 0;
    for (int ivar = 0; ivar < nvar; ivar++)
    {
      rank     = 2 * rank;
      int ipos = VH::whereElement(index1[ivar], iech);
      if (ipos < 0)
        positions[ivar] = -1;
      else
      {
        positions[ivar] = ipos + cumul[ivar];
        identity[ndef]  = ivar;
        ndef++;
        rank += 1;
      }
    }
    *rank_arg = rank;
    return ndef;
  }

  /**
   * Build the inverse of the Nugget Effect matrix
   * It is established for:
   * - the number of variables defined in 'dbin' (and in 'Model')
   * - the active samples of 'dbin'
   * - the samples where Z-variable (and possibly V-variable) is defined
   *
   * @param db Input Db structure
   * @param model Input Model structure
   * @param params A structure for ruling the parameters of SPDE
   */
  void InvNuggetOp::_buildInvNugget(Db* db, Model* model, const SPDEParam& params)
  {
    _invNuggetMatrix = nullptr;
    if (db == nullptr) return;
    int nech = db->getNSample();
    if (model == nullptr) return;
    int nvar = db->getNLoc(ELoc::Z);
    if (nvar != model->getNVar())
    {
      messerr("'db' and 'model' should have the same number of variables");
      return;
    }
    bool hasnugget = false;
    CovAniso* cova = nullptr;

    for (int icov = 0; icov < model->getNCov(); icov++)
    {
      if (model->getCovAniso(icov)->getType() == ECov::NUGGET)
      {
        cova      = model->getCovAniso(icov);
        hasnugget = true;
        break;
      }
    }
    if (!hasnugget)
    {
      MatrixSymmetric sills(model->getNVar());
      cova = CovAniso::createIsotropicMulti(*model->getContext(), ECov::NUGGET, 0, sills);
    }
    VectorInt ivars = VH::sequence(nvar);

    // Get the minimum value for diagonal terms
    double eps = params.getEpsNugget();
    VectorDouble minNug(nvar);
    for (int ivar = 0; ivar < nvar; ivar++)
      minNug[ivar] = eps * model->getTotalSill(ivar, ivar);

    // Play the non-stationarity (if needed)
    bool flag_nostat_sill = cova->isNoStatForVariance();
    if (flag_nostat_sill)
      cova->informDbInForSills(db);

    // Create sets of Vector of valid sample indices per variable (not masked and defined)
    VectorVectorInt index1 = db->getSampleRanks(ivars);
    // 'cumul' counts the number of valid positions for all variables before 'ivar'
    VectorInt cumul = VH::cumulIncrement(index1);

    // Check the various possibilities
    // - flag_verr: True if Variance of Measurement Error variable is defined
    // - flag_isotropic: True in Isotopic case
    // - flag_uniqueVerr: True if the Variance of Measurement Error is constant per variable
    // - flag_nostat: True is some non-stationarity is defined
    int nverr          = db->getNLoc(ELoc::V);
    bool flag_verr     = (nverr > 0);
    bool flag_isotopic = true;
    for (int ivar = 1; ivar < nvar && flag_isotopic; ivar++)
      if (!VH::isEqual(index1[ivar], index1[0])) flag_isotopic = false;
    bool flag_uniqueVerr = true;
    VectorDouble verrDef(nverr, 0.);
    if (flag_verr)
    {
      for (int iverr = 0; iverr < nverr && flag_uniqueVerr; iverr++)
      {
        VectorDouble verr = db->getColumnByLocator(ELoc::V, iverr);
        if ((int)VH::unique(verr).size() > 1) flag_uniqueVerr = false;
        verrDef[iverr] = verr[0];
      }
    }
    bool flag_constant = (!flag_nostat_sill && (!flag_verr || flag_uniqueVerr));

    // Elaborate the Sill matrix for the Nugget Effect component
    MatrixSymmetric sillsRef = cova->getSill();
    int count                = (int)pow(2, nvar);
    _sillsInv.resize(count);

    // Pre-calculate the inverse of the sill matrix (if constant)

    if (flag_constant)
    {
      // In case of (Unique) Variance of measurement error, patch sill matrix
      if (flag_verr) _addVerrConstant(sillsRef, verrDef);

      // Check that the diagonal of the Sill matrix is large enough
      _checkMinNugget(sillsRef, minNug);
    }

    // Constitute the triplet
    NF_Triplet NF_T;

    // Loop on the samples
    int rank;
    int ndef = nvar;
    VectorInt position(nvar);
    VectorInt identity(nvar);
    for (int iech = 0; iech < nech; iech++)
    {
      if (!db->isActive(iech)) continue;

      // Count the number of variables for which current sample is valid
      ndef = _loadPositions(iech, index1, cumul, position, identity, &rank);
      if (ndef <= 0) continue;

      // If all samples are defined, in the stationary case, use the inverted sill matrix
      if (flag_constant)
      {
        if (_sillsInv[rank].empty())
        {
          _sillsInv[rank] = _buildSillPartialMatrix(sillsRef, nvar, ndef, identity);
          if (_sillsInv[rank].invert() != 0) return;
        }
        for (int idef = 0; idef < ndef; idef++)
          for (int jdef = 0; jdef < ndef; jdef++)
            NF_T.add(position[identity[idef]], position[identity[jdef]],
                     _sillsInv[rank].getValue(idef, jdef));
      }
      else
      {
        // Update due to non-stationarity (optional)
        if (flag_nostat_sill)
        {
          cova->updateCovByPoints(1, iech, 1, iech);
          sillsRef = cova->getSill();
        }

        // Establish a local matrix
        MatrixSymmetric local(ndef);
        for (int idef = 0; idef < ndef; idef++)
          for (int jdef = 0; jdef <= idef; jdef++)
          {
            // Load the sill value of the Nugget Effect component
            double value = sillsRef.getValue(identity[idef], identity[jdef]);

            // Patch the diagonal term of the local matrix
            if (idef == jdef)
            {
              // Add the Variance of measurement error (optional)
              if (flag_verr && idef < nverr)
                value += db->getFromLocator(ELoc::V, iech, identity[idef]);

              // Check the minimum values over the diagonal
              value = MAX(value, MAX(local.getValue(idef, idef), minNug[idef]));
            }

            local.setValue(idef, jdef, value);
          }
        if (local.invert() != 0) return;

        for (int idef = 0; idef < ndef; idef++)
          for (int jdef = 0; jdef < ndef; jdef++)
            NF_T.add(position[identity[idef]], position[identity[jdef]],
                     local.getValue(idef, jdef));
      }
    }

    // Convert from triplet to sparse matrix
    _invNuggetMatrix = std::shared_ptr<MatrixSparse>(MatrixSparse::createFromTriplet(NF_T));

    // Free the non-stationary specific allocation
    if (!hasnugget)
      delete cova;
  }

std::shared_ptr<MatrixSparse> buildInvNugget(Db* dbin, Model* model, const SPDEParam& params)
  {
    InvNuggetOp invnugg(dbin, model, params);
    return invnugg.getInvNuggetMatrix();
  }

double InvNuggetOp::computeLogDet(int nMC) const
{
  DECLARE_UNUSED(nMC)
  return TEST;
}

const MatrixSparse* InvNuggetOp::cloneInvNuggetMatrix() const 
{
 return _invNuggetMatrix->clone(); 
}
} // namespace gstlrn