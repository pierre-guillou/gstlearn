#include "Basic/ListParams.hpp"
#include "Basic/AStringable.hpp"
#include "geoslib_define.h"
#include <cstddef>
#include <sstream>

ListParams::ListParams()
  : AStringable()
{
}

void ListParams::addParam(ParamInfo& param)
{
  _params.push_back(param);
}


double ListParams::getValue(int index) const
{
  if (index < 0 || index >= static_cast<int>(_params.size()))
  {
    messerr("Index out of range in ListParams::getValue");
    return TEST;
  }
  return _params[index].get().getValue();
}
void ListParams::setValue(int index, double value)
{
  if (index < 0 || index >= static_cast<int>(_params.size()))
  {
    messerr("Index out of range in ListParams::setValue");
    return;
  }
  _params[index].get().setValue(value);
}


String ListParams::toString(const AStringFormat* strfmt) const
{
  DECLARE_UNUSED(strfmt);
  std::stringstream result;
  result << toTitle(1,"List of Parameters:");
  for (int ipar = 0, npar = (int) _params.size(); ipar < npar; ipar++)
  {
    result << ipar+1 << " - " << _params[ipar].get().toString() << std::endl;
  }
  return result.str();
}

std::vector<double> ListParams::getValues() const
{
  size_t nparam = _params.size();
  std::vector<double> values(nparam);
  for (size_t i = 0; i < nparam; ++i)
  {
    values[i] = _params[i].get().getValue();
  }
  return values;
}

std::vector<double> ListParams::getMinValues() const
{
  size_t nparam = _params.size();
  std::vector<double> values(nparam);
  for (size_t i = 0; i < nparam; ++i)
  {
    values[i] = _params[i].get().getUserMin();
  }
  return values;
}

std::vector<double> ListParams::getMaxValues() const
{
  size_t nparam = _params.size();
  std::vector<double> values(nparam);
  for (size_t i = 0; i < nparam; ++i)
  {
    values[i] = _params[i].get().getUserMax();
  }
  return values;
}

void ListParams::setValues(const std::vector<double>& values)
{
  size_t size = values.size();
  for (size_t i = 0; i < size; i++)
  {
    _params[i].get().setValue(values[i]);
  }
}