#include "Interfaces/Dictionary.hpp"

using namespace gstlrn;

const Category& Dictionary::getCategory(const int value) const
{
  for (const auto& el: _categories)
  {
    if (el.getValue() == value)
    {
      return el;
    }
  }
  return _undefCategory;
}

const Category& Dictionary::getCategory(const String& label) const
{
  for (const auto& el: _categories)
  {
    if (el.getLabel() == label)
    {
      return el;
    }
  }
  return _undefCategory;
}

void Dictionary::addCategory(int value, const String& label)
{
  if (!hasCategory(value) && !hasCategory(label))
  {
    _categories.emplace_back(value, label);
  }
}

bool Dictionary::hasCategory(const int value) const
{
  return getCategory(value) != _undefCategory;
}

bool Dictionary::hasCategory(const String& label) const
{
  return getCategory(label) != _undefCategory;
}

bool Dictionary::hasCategory(const Category& cat) const
{
  if (hasCategory(cat.getValue()))
  {
    return getCategory(cat.getValue()).getLabel() == cat.getLabel();
  }
  return (false);
}
