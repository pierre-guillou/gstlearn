#include "Interfaces/VariableCategorical.hpp"

using namespace gstlrn;

VariableCategorical::VariableCategorical(const String& name,
                                         const Dictionary& dico)
  : AVariableTemplate(name)
  , _dico(dico)
{
}

void VariableCategorical::resize(const size_t n, const Category& val)
{
  AVariableTemplate::resize(
    n,
    _dico.hasCategory(val) ? val : _dico.getUndefCategory());
}

VectorDouble VariableCategorical::getValues() const
{
  VectorDouble vals;
  for (const auto& val: _values)
  {
    vals.push_back(val.getValue());
  }
  return (vals);
}

void VariableCategorical::setValue(const size_t i, const Category& cat)
{
  AVariableTemplate::setValue(
    i,
    _dico.hasCategory(cat) ? cat : _dico.getUndefCategory());
}

void VariableCategorical::setValue(const size_t i, int value)
{
  AVariableTemplate::setValue(
    i,
    _dico.hasCategory(value) ? _dico.getCategory(value) : _dico.getUndefCategory());
}

void VariableCategorical::setValue(const size_t i, const String& label)
{
  AVariableTemplate::setValue(
    i,
    _dico.hasCategory(label) ? _dico.getCategory(label) : _dico.getUndefCategory());
}
