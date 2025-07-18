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
#include "LinearOp/InvNuggetOp.hpp"
#include "LinearOp/SPDEOp.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "gstlearn_export.hpp"

namespace gstlrn
{
class ProjMultiMatrix;
class MatrixSparse;
class PrecisionOpMultiMatrix;

class GSTLEARN_EXPORT SPDEOpMatrix: public SPDEOp
{
public:
  SPDEOpMatrix(const PrecisionOpMultiMatrix* pop = nullptr,
               const ProjMultiMatrix* A          = nullptr,
               const InvNuggetOp* invNoise       = nullptr,
               const ProjMultiMatrix* projOut    = nullptr);
  virtual ~SPDEOpMatrix();

  double computeLogDetOp(int nbsimu) const override;

  VectorDouble stdev(const VectorDouble& dat, int nMC, int seed) const override;

#ifndef SWIG

private:
  int _addToDest(const constvect inv, vect outv) const override;
  int _solve(const constvect inv, vect outv) const override;
#endif

private:
  mutable std::shared_ptr<MatrixSparse> _QpAinvNoiseAt; // mutable is required to perform the Cholesky decomposition
  mutable CholeskySparse* _chol;                        // when needed, e.g in a const method.
};
} // namespace gstlrn