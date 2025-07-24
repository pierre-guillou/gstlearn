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

#include "Interfaces/AVariable.hpp"
#include "gstlearn_export.hpp"

#include <iostream>
#include <vector>

namespace gstlrn
{

template<class T>
class GSTLEARN_EXPORT AVariableTemplate: public AVariable
{
public:
  AVariableTemplate()
    : AVariable()
    , _values()
  {
  }
  AVariableTemplate(const String& name)
    : AVariable(name)
    , _values()
  {
  }
  AVariableTemplate(const String& name, const std::vector<T>& values)
    : AVariable(name)
    , _values(values)
  {
  }

  size_t getNValues() const override { return _values.size(); }

  virtual void resize(const size_t n, const T& val) { _values.resize(n, val); }

  virtual void setValue(const size_t i, const T& val) { _values.at(i) = val; }

  virtual void setValues(const std::vector<T>& values) { _values = values; }

  T getValueAsType(const size_t i) const { return _values.at(i); }

  virtual const std::vector<T>& getValuesAsType() const { return _values; }

  virtual void display_old() const
  {
    std::cout << "Variable " << _name << " :\n";
    for (size_t i = 0; i < _values.size(); i++)
    {
      std::cout << "  var[" << i << "]=" << _values[i] << " (" << _values[i] << ")\n";
    }
  }

protected:
  std::vector<T> _values;
};

} // namespace gstlrn
