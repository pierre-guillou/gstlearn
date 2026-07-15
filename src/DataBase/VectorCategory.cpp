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
#include "DataBase/VectorCategory.hpp"

namespace gstlrn
{
  bool VectorCategory::setCategory(const size_t isample, const Category& cat)
  {
    if (isample >= this->_data.size()) return false;
    if (!this->_dict.hasCategory(cat)) return false;
    this->_data[isample] = cat;
    return true;
  }
} // namespace gstlrn
