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
#include "Fractures/FracFamily.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"
#include "Basic/SerializeHDF5.hpp"

FracFamily::FracFamily(double orient,
                       double dorient,
                       double theta0,
                       double alpha,
                       double ratcst,
                       double prop1,
                       double prop2,
                       double aterm,
                       double bterm,
                       double range)
    : AStringable(),
      ASerializable(),
      _orient(orient),
      _dorient(dorient),
      _theta0(theta0),
      _alpha(alpha),
      _ratcst(ratcst),
      _prop1(prop1),
      _prop2(prop2),
      _aterm(aterm),
      _bterm(bterm),
      _range(range)
{
}

FracFamily::FracFamily(const FracFamily& r)
    : AStringable(r),
      ASerializable(r),
      _orient(r._orient),
      _dorient(r._dorient),
      _theta0(r._theta0),
      _alpha(r._alpha),
      _ratcst(r._ratcst),
      _prop1(r._prop1),
      _prop2(r._prop2),
      _aterm(r._aterm),
      _bterm(r._bterm),
      _range(r._range)
{
}

FracFamily& FracFamily::operator=(const FracFamily& r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    ASerializable::operator=(r);
    _orient = r._orient;
    _dorient =r ._dorient;
    _theta0 = r._theta0;
    _alpha = r._alpha;
    _ratcst = r._ratcst;
    _prop1 = r._prop1;
    _prop2 = r._prop2;
    _aterm = r._aterm;
    _bterm = r._bterm;
    _range = r._range;
  }
  return *this;
}

FracFamily::~FracFamily()
{
}

String FracFamily::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  sstr << "Average Fault Orientation       = " << _orient << " (degree)" << std::endl;
  sstr << "Tolerance for Orientation       = " << _dorient << " (degree)" << std::endl;
  sstr << "Reference Poisson Intensity     = " << _theta0 << std::endl;
  sstr << "Intensity from thick. exponent  = " << _alpha << std::endl;
  sstr << "Intensity Constant/Shaped ratio = " << _ratcst << std::endl;
  sstr << "Survival constant probability   = " << _prop1 << std::endl;
  sstr << "Survival length-dependent proba = " << _prop2 << std::endl;
  sstr << "Survival cumul. length exponent = " << _aterm << std::endl;
  sstr << "Survival thickness exponent     = " << _bterm << std::endl;
  sstr << "Fracture repulsion area Range   = " << _range << std::endl;

  return sstr.str();
}

bool FracFamily::_deserializeAscii(std::istream& is, bool /*verbose*/)
{
  bool ret = true;
  ret = ret && _recordRead<double>(is, "Mean orientation", _orient);
  ret = ret && _recordRead<double>(is, "Tolerance for orientation", _dorient);
  ret = ret && _recordRead<double>(is, "Reference Poisson intensity", _theta0);
  ret = ret && _recordRead<double>(is, "Power dependency between layer and intensity", _alpha);
  ret = ret && _recordRead<double>(is, "Ratio of constant vs. shaped intensity", _ratcst);
  ret = ret && _recordRead<double>(is, "Survival probability (constant term)", _prop1);
  ret = ret && _recordRead<double>(is, "Survival probability (length dependent term)", _prop2);
  ret = ret && _recordRead<double>(is, "Survival probability (cumulative length exponent)", _aterm);
  ret = ret && _recordRead<double>(is, "Survival probability (layer thickness exponent", _bterm);
  ret = ret && _recordRead<double>(is, "Fracture repulsion area Range", _range);
  return ret;
}

bool FracFamily::_serializeAscii(std::ostream& os, bool /*verbose*/) const
{
  bool ret = true;
  ret = ret && _recordWrite<double>(os, "Mean orientation", _orient);
  ret = ret && _recordWrite<double>(os, "Tolerance for orientation", _dorient);
  ret = ret && _recordWrite<double>(os, "Reference Poisson intensity", _theta0);
  ret = ret && _recordWrite<double>(os, "Power dependency between layer and intensity", _alpha);
  ret = ret && _recordWrite<double>(os, "Ratio of constant vs. shaped intensity", _ratcst);
  ret = ret && _recordWrite<double>(os, "Survival probability (constant term)", _prop1);
  ret = ret && _recordWrite<double>(os, "Survival probability (length dependent term)", _prop2);
  ret = ret && _recordWrite<double>(os, "Survival probability (cumulative length exponent)", _aterm);
  ret = ret && _recordWrite<double>(os, "Survival probability (layer thickness exponent)", _bterm);
  ret = ret && _recordWrite<double>(os, "Fracture repulsion area Range", _range);
  return ret;
}

#ifdef HDF5
bool FracFamily::_deserializeH5(H5::Group& grp, [[maybe_unused]] bool verbose)
{
  // Call SerializeHDF5::getGroup to get the subgroup of grp named
  // "FracFamily" with some error handling
  auto famG = SerializeHDF5::getGroup(grp, "FracFamily");
  if (!famG)
  {
    return false;
  }

  /* Read the grid characteristics */
  bool ret = true;

  ret = ret && SerializeHDF5::readValue(*famG, "Orient", _orient);
  ret = ret && SerializeHDF5::readValue(*famG, "TolOrient", _dorient);
  ret = ret && SerializeHDF5::readValue(*famG, "Intensity", _theta0);
  ret = ret && SerializeHDF5::readValue(*famG, "Power", _alpha);
  ret = ret && SerializeHDF5::readValue(*famG, "Ratio", _ratcst);
  ret = ret && SerializeHDF5::readValue(*famG, "Surv_cste", _prop1);
  ret = ret && SerializeHDF5::readValue(*famG, "Surv_length", _prop2);
  ret = ret && SerializeHDF5::readValue(*famG, "Surv_cumExp", _aterm);
  ret = ret && SerializeHDF5::readValue(*famG, "Surv_thickExp", _bterm);
  ret = ret && SerializeHDF5::readValue(*famG, "Repulsion", _range);

  return ret;
}

bool FracFamily::_serializeH5(H5::Group& grp, [[maybe_unused]] bool verbose) const
{
  auto famG = grp.createGroup("FracFamily");

  bool ret = true;
  ret      = ret && SerializeHDF5::writeValue(famG, "Orient", _orient);
  ret      = ret && SerializeHDF5::writeValue(famG, "TolOrient", _dorient);
  ret      = ret && SerializeHDF5::writeValue(famG, "Intensity", _theta0);
  ret      = ret && SerializeHDF5::writeValue(famG, "Power", _alpha);
  ret      = ret && SerializeHDF5::writeValue(famG, "Ratio", _ratcst);
  ret      = ret && SerializeHDF5::writeValue(famG, "Surv_cste", _prop1);
  ret      = ret && SerializeHDF5::writeValue(famG, "Surv_length", _prop2);
  ret      = ret && SerializeHDF5::writeValue(famG, "Surv_cumExp", _aterm);
  ret      = ret && SerializeHDF5::writeValue(famG, "Surv_thickExp", _bterm);
  ret      = ret && SerializeHDF5::writeValue(famG, "Repulsion", _range);

  return ret;
}
#endif
