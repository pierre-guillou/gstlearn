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
#pragma once

#include "Basic/AStringable.hpp"
#include "Basic/ICloneable.hpp"
#include "Interfaces/interface_d.hpp"
#include "gstlearn_export.hpp"

#include <memory>

namespace gstlrn
{

class VariableDouble;
class VariableBool;
class VariableInt;
class VariableString;

class GSTLEARN_EXPORT AVariable: public AStringable, public ICloneable
{
public:
  AVariable() = default;
  AVariable(const String& name)
    : AStringable()
    , _name(name)
  {
  }
  ~AVariable() = default;

  virtual size_t getNValues() const = 0;
  virtual VectorDouble getValues() const  = 0;

  static std::unique_ptr<AVariable> createVariable(const String& type, const String& name);

  void setName(const String& name)
  {
    _name = name;
  }

  const String& getName() const
  {
    return _name;
  }

protected:
  String _name {UNDEF_STRING};
};

} // namespace gstlrn
