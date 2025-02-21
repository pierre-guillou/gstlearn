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

  const String& GetName() const { return _name; }

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
  void AddArray(ADbCol& array);
  void RemoveArray(const String& name);
  void RemoveArray(const Id index);
  ADbCol* GetArray(const Id index)
  {
    if (index < 0 || index >= static_cast<Id>(_cols.size()))
    {
      return {};
    }
    return _cols[index].get();
  }
  ADbCol* GetArray(const String& name, Id& index)
  {
    index = 0;
    for (const auto& col: _cols)
    {
      if (col->GetName() == name)
      {
        return col.get();
      }
      index++;
    }

    index = -1;
    return {};
  }
  ADbCol* GetArray(const String& name)
  {
    Id i;
    return this->GetArray(name, i);
  }
  std::optional<String> GetArrayName(const Id index)
  {
    auto* arr = this->GetArray(index);
    return arr ? std::optional<String> {arr->GetName()} : std::nullopt;
  }
  bool HasArray(const String& name)
  {
    Id i;
    auto* arr = this->GetArray(name, i);
    return arr != nullptr;
  }

private:
  std::vector<std::unique_ptr<ADbCol>> _cols;
};

} // namespace gstlrn
