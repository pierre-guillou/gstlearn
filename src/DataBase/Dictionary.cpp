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
#include "DataBase/Dictionary.hpp"

namespace gstlrn
{
  bool Dictionary::addCategory(const Id key, const String& val)
  {
    for (const auto& p: this->_data)
    {
      if (p.first == key || p.second == val)
      {
        return false;
      }
    }
    this->_data[key] = val;
    return true;
  }

  bool Dictionary::hasCategory(const Category& cat) const
  {
    const auto it = this->_data.find(cat.first);
    if (it == this->_data.end()) return false;
    if (it->second != cat.second) return false;
    return true;
  }

} // namespace gstlrn
