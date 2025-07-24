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
#include "Interfaces/Dictionary.hpp"
#include "gstlearn_export.hpp"

namespace gstlrn
{

/// TODO : to be terminated
class GSTLEARN_EXPORT VariableCategorical: public AVariableTemplate<Category>
{
public:
  VariableCategorical(const String& name, const Dictionary& dico);

  /// Cloneable interface
  IMPLEMENT_CLONING(VariableCategorical)

  void resize(const size_t n, const Category& val) override;
  VectorDouble getValues() const override;
  void setValue(const size_t i, const Category& cat) override;
  void setValue(const size_t i, const int value);
  void setValue(const size_t i, const String& label);

private:
  Dictionary _dico;
};

} // namespace gstlrn
