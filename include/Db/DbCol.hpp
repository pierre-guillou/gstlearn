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
#include <optional>
#include <type_traits>
#include <variant>

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

  template<typename T>
  const T& getValue(const Id i) const
  {
    return std::visit([=](auto&& arg)
                      {
                        using VectorType = std::decay_t<decltype(arg)>;
                        if constexpr (std::is_same<typename VectorType::value_type, T>::value)
                        {
                          return arg[i];
                        }static_assert(false, "non-exhaustive visitor!"); }, this->_data);
  }

  template<typename T>
  void setValue(const Id i, const T& v)
  {
    std::visit([&](auto&& arg)
               {
     using VectorType = std::decay_t<decltype(arg)>;
     if constexpr (std::is_same<typename VectorType::value_type, T>::value)
     {
       arg[i] = v;
     }
     else
     {
       static_assert(false, "non-exhaustive visitor!");
     } }, this->_data);
  }

private:
  String _name;
  std::variant<VectorDouble, VectorFloat, VectorInt, VectorUChar, VectorString, VectorBool> _data;
};

/**
 * @brief Similar to vtkFieldData
 */
class DbData
{
public:
  template<typename VectorType>
  void AddArray(const VectorType& array)
  {
    this->_cols.emplace_back(array);
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
    return {this->_cols[index]};
  }

  std::optional<std::pair<std::reference_wrapper<ADbCol>, Id>> GetArray(const String& name)
  {
    Id index {};
    for (auto& col: this->_cols)
    {
      if (col.GetName() == name)
      {
        return {{col, index}};
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
  std::vector<ADbCol> _cols;
};

} // namespace gstlrn
