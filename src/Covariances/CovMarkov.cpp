/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "Covariances/CovMarkov.hpp"
#include "Covariances/CovContext.hpp"
#include "Basic/MathFunc.hpp"
#include "geoslib_f.h"
#include "math.h"

#define MAXTAB 100

CovMarkov::CovMarkov(const CovContext& ctxt)
: ACovFunc(ECov::MARKOV, ctxt)
{
  setParam(1);
}

CovMarkov::CovMarkov(const CovMarkov &r)
: ACovFunc(r)
{
}

CovMarkov& CovMarkov::operator=(const CovMarkov &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovMarkov::~CovMarkov()
{
}

double CovMarkov::getScadef() const
{
  return getParam();
}

double CovMarkov::_evaluateCov(double /*h*/) const
{
  return TEST;
}

String CovMarkov::getFormula() const
{
  return "C(h)=\\int_{R^d} \\frac{e^{-i\\omega^t.h}}{P(||\\omega||^2)}d\\omega";
}

double CovMarkov::_evaluateCovOnSphere(double scale, int degree) const
{
  double s = 0.;
  int n = (int)_markovCoeffs.size();
  double nnp1 = scale * scale * (double) degree * ((double) degree + 1.);
  for(int i = 0; i< n;i++)
  {
    s += _markovCoeffs[i] * pow(nnp1,i);
  }
  return scale * scale * (2. * degree + 1.) / (4 * GV_PI * s);
}

double CovMarkov::evaluateSpectrum(double freq, double scale, int /*ndim*/) const
{
  double s = 0.;
  int n = (int)_markovCoeffs.size();
  for(int i = 0; i< n;i++)
  {
      s += _markovCoeffs[i] * pow(freq,i);
  }
  return 1. / s;
}
