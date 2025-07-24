#pragma once

#include "Interfaces/interface_d.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

#include <iostream>

namespace gstlrn
{

class GSTLEARN_EXPORT Category
{
public:
  Category(int value = UNDEF_CAT_VAL, const String& label = UNDEF_CAT_LABEL);
  virtual ~Category();

  const String& getLabel() const;
  int getValue() const;

  bool operator==(const Category& rhs) const { return this->_label == rhs._label; }

#ifndef SWIG
  // conversion operator
  explicit operator int() const
  {
    return (_value);
  }
  explicit operator double() const
  {
    return (_value);
  }

#endif

private:
  int _value;
  String _label;
};

std::ostream& operator<<(std::ostream& stream, const Category& cat);

} // namespace gstlrn
