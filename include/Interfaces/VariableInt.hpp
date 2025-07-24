/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

#include "Interfaces/AVariableTemplate.hpp"
#include "gstlearn_export.hpp"

namespace gstlrn
{

class GSTLEARN_EXPORT VariableInt: public AVariableTemplate<Id>
{
public:
  VariableInt();
  VariableInt(const String& name);
  VariableInt(const VectorInt& values);
  VariableInt(const VariableInt& ref);
  ~VariableInt() override;
  VariableInt& operator=(const VariableInt& ref);

  /// Cloneable interface
  IMPLEMENT_CLONING(VariableInt)

  VectorDouble getValues() const override;
  // bool                  isUndefined(int i) const override;

#ifdef _USE_NETCDF
  netCDF::NcVar serialize(netCDF::NcFile& file, std::vector<netCDF::NcDim>& dims) const override;
  void deserialize(const netCDF::NcFile& file, const netCDF::NcVar& var) override;
#endif

private:
};

} // namespace gstlrn
