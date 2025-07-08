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

#include "API/SPDEParam.hpp"
#include "LinearOp/ASimulable.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "gstlearn_export.hpp"

#include <memory>
namespace gstlrn
{
  class Db;
  class Model;
  class MatrixSparse;
  class MatrixSymmetric;
    /**
     * invNuggetOp
     *
     * This class is used to represent the inverse covariance matrix operator of a nugget effect.
     * It can be multivariate with anisotropy.
     * It stores the inverse submatrices for each location, allowing efficient computation of 
     * the logdeterminant.
     * It inherits from ASimulable, allowing it to be used in simulations.
     */
  class GSTLEARN_EXPORT InvNuggetOp: public ASimulable
  {
  public:
    InvNuggetOp(Db* dbin = nullptr, Model* model = nullptr, const SPDEParam& params = SPDEParam());
    InvNuggetOp(const InvNuggetOp& m) = delete;
    InvNuggetOp& operator=(const InvNuggetOp& m) = delete;
    virtual ~InvNuggetOp();
    int getSize() const override;
    std::shared_ptr<MatrixSparse> getInvNugget() const { return _invNugget; }
  protected:
   int _addSimulateToDest(const constvect whitenoise, vect outv) const override;
   int _addToDest(constvect inv, vect outv) const override;
   
  private:
   void _buildInvNugget(Db* dbin = nullptr, Model* model = nullptr, const SPDEParam& params = SPDEParam());
  
  private :
    std::shared_ptr<MatrixSparse> _invNugget; // The inverse nugget matrix
    std::vector<MatrixSymmetric> _sillsInv;
  };

GSTLEARN_EXPORT std::shared_ptr<MatrixSparse> buildInvNugget(Db* dbin, Model* model, const SPDEParam& params = SPDEParam());
} // namespace gstlrn

