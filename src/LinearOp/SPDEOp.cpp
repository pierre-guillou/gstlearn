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
#include "LinearOp/SPDEOp.hpp"
#include "LinearOp/ALinearOp.hpp"
#include "LinearOp/ProjMulti.hpp"
#include "LinearOp/PrecisionOpMulti.hpp"

SPDEOp::SPDEOp(const PrecisionOpMulti* pop, const ProjMulti* A, const ALinearOp* invNoise)
: _Q(pop)
, _A(A)
, _invNoise(invNoise)
{
}

SPDEOp::~SPDEOp() {}

int SPDEOp::getSize() const
{ 
  return _Q->getSize(); 
}
/*****************************************************************************/
/*!
**  Evaluate the product (by the SPDEOp) : 'outv' = I * 'inv' = 'inv'
**
** \param[in]  inv     Array of input values
**
** \param[out] outv    Array of output values
**
*****************************************************************************/
int SPDEOp::_addToDest(const Eigen::VectorXd& inv,
                          Eigen::VectorXd& outv) const
{
  for (int i = 0, n = getSize(); i < n; i++)
    outv[i] += inv[i];
  return 0;
}