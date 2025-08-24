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
#include "Covariances/CovExponential.hpp"
#include "Basic/Law.hpp"
#include "Basic/VectorHelper.hpp"
#include "Covariances/CovContext.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Simulation/TurningBandOperate.hpp"

#include <cmath>

namespace gstlrn
{
CovExponential::~CovExponential()
{
}

double CovExponential::getScadef() const
{
  return (2.995732);
}

String CovExponential::getFormula() const
{
  return "C(h)=exp \\left( -\\frac{h}{a_t} \\right)";
}

double CovExponential::simulateTurningBand(double t0, TurningBandOperate& operTB) const
{
  return operTB.spectralOne(t0);
}

MatrixDense CovExponential::simulateSpectralOmega(Id nb) const
{
  auto ndim    = static_cast<Id>(getContext().getNDim());
  double param = 0.5;
  MatrixDense mat(nb, ndim);

  for (Id irow = 0; irow < nb; irow++)
  {
    double scale = sqrt(param / law_gamma(param));
    for (Id icol = 0; icol < ndim; icol++)
      mat.setValue(irow, icol, scale * law_gaussian());
  }
  return mat;
}

double CovExponential::_evaluateCovOnSphere(double alpha,
                                            double scale,
                                            Id degree) const
{
  DECLARE_UNUSED(degree);
  double nu = scale * getScadef();
  return exp(-nu * alpha);
}

VectorDouble CovExponential::_evaluateSpectrumOnSphere(Id n, double scale) const
{
  double nu    = scale * getScadef();
  double nu2   = nu * nu;
  double expnu = exp(-nu * GV_PI);

  VectorDouble sp(1 + n, 0.);
  Id k;

  k     = 0;
  sp[k] = 1. / 2. * (1. + expnu) / (1. + nu2);
  k     = 1;
  sp[k] = 3. / 2. * (1. - expnu) / (4. + nu2);

  while (1)
  {
    k++;
    if (k >= n + 1) break;
    sp[k] = (2. * k + 1.) / (2. * k - 3.) * (nu2 + (k - 2.) * (k - 2.)) / (nu2 + (k + 1.) * (k + 1.)) * sp[k - 2];
  }

  VH::normalize(sp, 1);

  return sp;
}

// double CovExponential::_evaluateCovDerivative(double h) const
// {
//   if (h > MAX_EXP) return (0.);
//   double cov = -exp(-h);
//   return (cov);
// }

} // namespace gstlrn
