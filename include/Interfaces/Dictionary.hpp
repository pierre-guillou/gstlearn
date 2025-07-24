#pragma once

#include "Interfaces/Category.hpp"
#include "gstlearn_export.hpp"

#include <vector>

namespace gstlrn
{

class GSTLEARN_EXPORT Dictionary
{
public:
  Dictionary() = default;
  Dictionary(const String& name)
    : _name {name}
  {
  }

  const String& getName() const { return _name; }
  const Category& getUndefCategory() const { return _undefCategory; }

  const Category& getCategory(const int value) const;
  const Category& getCategory(const String& label) const;
  void addCategory(const int value, const String& label);
  bool hasCategory(const int value) const;
  bool hasCategory(const String& label) const;
  bool hasCategory(const Category& cat) const;

private:
  String _name;
  Category _undefCategory;
  std::vector<Category> _categories;
};

} // namespace gstlrn
