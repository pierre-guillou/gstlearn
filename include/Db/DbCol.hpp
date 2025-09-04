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

#include "Basic/VectorNumT.hpp"
#include "Basic/VectorT.hpp"
#include "gstlearn_export.hpp"

#include <functional>
#include <memory>
#include <optional>

namespace gstlrn
{

/**
 * @brief Generic implementation of a Column inside a Db
 */
class GSTLEARN_EXPORT ADbCol
{
public:
  ADbCol(const String& name)
    : _name {name}
  {
  }

  const String& GetName() const { return this->_name; }

private:
  String _name;
};

template<typename VectorType>
class DbColTemplate: public ADbCol
{
  VectorType _values;
};

using DbColDouble = DbColTemplate<VectorDouble>;
using DbColFloat  = DbColTemplate<VectorFloat>;
using DbColId     = DbColTemplate<VectorInt>;
using DbColUchar  = DbColTemplate<VectorUChar>;
using DbColStr    = DbColTemplate<VectorString>;
using DbColBool   = DbColTemplate<VectorBool>;

/**
 * @brief Similar to vtkFieldData
 */
class DbData
{
public:
  template<typename VectorType>
  void AddArray(const VectorType& array)
  {
    if constexpr (std::is_same<VectorType, VectorDouble>::value)
    {
      this->_cols.emplace_back(std::make_unique<DbColDouble>(array));
    }
    else if constexpr (std::is_same<VectorType, VectorFloat>::value)
    {
      this->_cols.emplace_back(std::make_unique<DbColFloat>(array));
    }
    else if constexpr (std::is_same<VectorType, VectorInt>::value)
    {
      this->_cols.emplace_back(std::make_unique<DbColId>(array));
    }
    else if constexpr (std::is_same<VectorType, VectorUChar>::value)
    {
      this->_cols.emplace_back(std::make_unique<DbColUchar>(array));
    }
    else if constexpr (std::is_same<VectorType, VectorString>::value)
    {
      this->_cols.emplace_back(std::make_unique<DbColStr>(array));
    }
    else if constexpr (std::is_same<VectorType, VectorBool>::value)
    {
      this->_cols.emplace_back(std::make_unique<DbColBool>(array));
    }
  }

  void RemoveArray(const String& name)
  {
    const auto col = this->GetArray(name);
    const auto id  = col ? col->second : -1;
    this->RemoveArray(id);
  }

  void RemoveArray(const Id index)
  {
    this->_cols.erase(this->_cols.begin() + index);
  }

  std::optional<std::reference_wrapper<ADbCol>> GetArray(const Id index)
  {
    if (index < 0 || index >= static_cast<Id>(this->_cols.size()))
    {
      return {};
    }
    return {*this->_cols[index]};
  }

  std::optional<std::pair<std::reference_wrapper<ADbCol>, Id>> GetArray(const String& name)
  {
    Id index {};
    for (const auto& col: this->_cols)
    {
      if (col->GetName() == name)
      {
        return {{*col, index}};
      }
      index++;
    }

    return {};
  }

  std::optional<String> GetArrayName(const Id index)
  {
    auto arr = this->GetArray(index);
    return arr ? std::optional<String> {arr->get().GetName()} : std::nullopt;
  }

  bool HasArray(const String& name)
  {
    auto arr = this->GetArray(name);
    return arr.has_value();
  }

private:
  std::vector<std::unique_ptr<ADbCol>> _cols;
};

} // namespace gstlrn
