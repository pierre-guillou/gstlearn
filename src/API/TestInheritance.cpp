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
#include "API/TestInheritance.hpp"

#include <sstream>

TestInheritance::TestInheritance()
 : _iproj(nullptr)
{
}

TestInheritance::~TestInheritance()
{
}

String TestInheritance::toString(const AStringFormat* strfmt) const
{
  DECLARE_UNUSED(strfmt);
  std::stringstream sstr;
  if (nullptr != _iproj)
  {
    sstr << "Number of Apices: " << _iproj->getNApex()  << std::endl;
    sstr << "Number of Points: " << _iproj->getNPoint();
  }
  else
  {
    sstr << "Projection Matrix not yet defined!";
  }
  return sstr.str();
}
