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
#include "Simulation/SimuSphericalParam.hpp"
#include "Basic/Vector.hpp"

SimuSphericalParam::SimuSphericalParam(int special,
                                       int nbf,
                                       int nfmax,
                                       int degmax,
                                       int ndisc,
                                       double tol)
    : AStringable(),
      _special(special),
      _nbf(nbf),
      _nfmax(nfmax),
      _degmax(degmax),
      _ndisc(ndisc),
      _tol(tol)
{
}

SimuSphericalParam::SimuSphericalParam(const SimuSphericalParam &r)
    : AStringable(r),
      _special(r._special),
      _nbf(r._nbf),
      _nfmax(r._nfmax),
      _degmax(r._degmax),
      _ndisc(r._ndisc),
      _tol(r._tol)
{
}

SimuSphericalParam& SimuSphericalParam::operator=(const SimuSphericalParam &r)
{
  if (this != &r)
  {
    AStringable::operator =(r);
    _special = r._special;
    _nbf = r._nbf;
    _nfmax = r._nfmax;
    _degmax = r._degmax;
    _ndisc = r._ndisc;
    _tol = r._tol;
  }
  return *this;
}

SimuSphericalParam::~SimuSphericalParam()
{
}

String SimuSphericalParam::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  sstr << toTitle(1, "Covariance spectrum in Spherical Coordinates generated");
  if (_special == 1)
    sstr << "Using Chentsov construction";
  else if (_special == 2)
    sstr << "Using Exponential construction";
  else
    sstr << "Number of discretization  = " << _ndisc << std::endl;

  sstr << "Spectrum Tolerance        = " << _tol << std::endl;
  sstr << "Number of basic functions = " << _nbf << std::endl;
  return sstr.str();
}
