/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2025) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "geoslib_define.h"

#include <map>
#include <optional>
#include <string_view>
#include <utility>

namespace gstlrn
{
  class Dictionary
  // class GSTLEARN_EXPORT Dictionary
  {
  public:
    using Category = std::pair<Id, std::string_view>;

    Dictionary() = default;

    Dictionary(std::map<Id, String>&& data)
      : _data{std::move(data)}
    {
    }

    bool addCategory(const Id key, const String& val);

    bool hasCategory(const Category& cat) const;

#ifndef SWIG
    std::optional<Category> operator[](const Id key) const
    {
      const auto it = this->_data.find(key);
      if (it == this->_data.end()) return {};
      return *it;
    }
#endif

  private:
    std::map<Id, String> _data;
  };

} // namespace gstlrn
