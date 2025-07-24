#include "Interfaces/AVariable.hpp"
#include "Interfaces/VariableBool.hpp"
#include "Interfaces/VariableDouble.hpp"
#include "Interfaces/VariableInt.hpp"
#include "Interfaces/VariableString.hpp"

using namespace gstlrn;

/**
 * Create a variable with appropriate type
 *
 * @param[in] type   String describing the desired type
 *                   (int, bool, double, string)
 * @param[in] name   name of the created variable.
 *
 * @remark: if type does not match any existing type, an exception is throw
 */
std::unique_ptr<AVariable> AVariable::createVariable(const String& type, const String& name)
{
  if (type == "int")
  {
    return std::make_unique<VariableInt>(name);
  }
  if (type == "bool")
  {
    return std::make_unique<VariableBool>(name);
  }
  if (type == "double")
  {
    return std::make_unique<VariableDouble>(name);
  }
  if (type == "string")
  {
    return std::make_unique<VariableString>(name);
  }

  throw("Wrong variable type");
}
