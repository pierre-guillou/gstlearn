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

#include "DataBase/Dictionary.hpp"
#include <cstddef> // size_t
#include <optional> // std::optional
#include <vector> // std::vector

namespace gstlrn
{
  class VectorCategory
  // class GSTLEARN_EXPORT VectorCategory
  {
  public:
    using Category = Dictionary::Category;
    using value_type = Category;

    VectorCategory(
      const size_t nsample = 0,
      const Dictionary& dict = Dictionary())
      : _data(nsample)
      , _dict{dict}
    {
    }

    void resize(const size_t count) { this->_data.resize(count); }

    size_t size() const { return this->_data.size(); }

    std::optional<Category> getCategory(const size_t isample) const
    {
      if (isample >= this->_data.size()) return {};
      return this->_data[isample];
    }

    bool setCategory(const size_t isample, const Category& cat);

#ifndef SWIG
    std::optional<Category>& operator[](const size_t isample)
    {
      return this->_data[isample];
    }

    const std::optional<Category>& operator[](const size_t isample) const
    {
      return this->_data[isample];
    }
#endif

  private:
    std::vector<std::optional<Category>> _data;
    Dictionary _dict;
  };

} // namespace gstlrn
