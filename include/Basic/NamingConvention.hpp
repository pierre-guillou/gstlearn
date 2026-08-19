/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "Basic/AStringable.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

#include "Enum/ELoc.hpp"

namespace gstlrn
{
  class Db;

  /**
   * @brief Naming Convention facility.
   *
   * This class describes the way variables created within the current procedure
   * will be named afterwards and will possibly be assigned a locator.
   *
   * The generic name is generated as follows:
   *      'prefix'.'varname'.'qualifier'|'rank'
   *
   * - prefix: string provided in the constructor of this class
   * - varname: name of the (input) variable on which the procedure is performed
   * - qualifier: type of element stored in the variable
   * - rank: rank of the output variable (if several simulations are generated)
   *
   * The choice of the 'prefix' is done by the user when launching the procedure;
   * the other parameters are usually defined within the calling procedure.
   *
   * For example, when running 'kriging' function with several variables defined
   * in the input Db - say "Pb" and "Zn" (they are assigned a Z-locator),
   * using the following command:
   *    kriging( ... namconv = NamingConvention("MyPrefix") )
   *
   * Then the kriging procedure generates variables such as:
   * - MyPrefix.Pb.estim (estimation of Pb by CoKriging)
   * - MyPrefix.Zn.estim (estimation of Zn by CoKriging)
   * - MyPrefix.Pb.stdev (St. Dev. of estimation error of Pb by CoKriging)
   * - MyPrefix.Zn.stdev (St. Dev. of estimation error of Zn by CoKriging)
   *
   * Then the non-conditional simulation procedure generates variables such as:
   * - MyPrefix.S1 (for first simulation)
   * - MyPrefix.S2 (for second simulation)
   * ...
   *
   * Then the conditional simulation procedure generates variables such as:
   * - MyPrefix.var.S1 (for first simulation)
   * - MyPrefix.var.S2 (for second simulation)
   * ...
   *
   * For multivariate simulations, the setOutputForSimulations method
   * provides consistent naming with explicit Variable and Simulation indicators:
   *
   * Non-conditional multivariate simulations (e.g., 2 variables, 2 simulations):
   * - MyPrefix.V1.S1, MyPrefix.V1.S2, MyPrefix.V2.S1, MyPrefix.V2.S2
   *
   * Conditional multivariate simulations (e.g., variables Fe and Al, 2 simulations):
   * - MyPrefix.Fe.S1, MyPrefix.Fe.S2, MyPrefix.Al.S1, MyPrefix.Al.S2
   *
   * Ultimately, the newly created variables are assigned a locator.
   *
   * Note: the related method getNameEncoded provides a static way to retrieve
   * the variable name based on the same convention (see comments).
   */
  class GSTLEARN_EXPORT NamingConvention: public AStringable
  {
  public:
    /**
     * @brief Constructor.
     *
     * @param prefix Prefix used for naming the output variables.
     * @param flag_varname If true, the variable name is included in the
     * generated name.
     * @param flag_qualifier If true, the qualifier is included in the
     * generated name.
     * @param flag_locator If true, a locator is assigned to the output
     * variables.
     * @param locatorOutType Type of locator assigned to the output variables.
     * @param delim Delimiter used to separate the components of the generated
     * name.
     * @param cleanSameLocator If true, variables with the same locator are
     * cleaned beforehand.
     */
    NamingConvention(
      const String& prefix = "",
      bool flag_varname = true,
      bool flag_qualifier = true,
      bool flag_locator = true,
      const ELoc& locatorOutType = ELoc::fromKey("Z"),
      const String& delim = ".",
      bool cleanSameLocator = true);

    /**
     * @brief Copy constructor.
     *
     * @param m NamingConvention object to copy.
     */
    NamingConvention(const NamingConvention& m);

    /**
     * @brief Assignment operator.
     *
     * @param m NamingConvention object to copy.
     * @return Reference to the current object.
     */
    NamingConvention& operator=(const NamingConvention& m);

    /**
     * @brief Destructor.
     */
    virtual ~NamingConvention();

    /// AStringable Interface
    String toString(const AStringFormat* strfmt = nullptr) const override;

    /**
     * @brief Create a NamingConvention object.
     *
     * This static method is a convenience function for creating a
     * NamingConvention object with the specified naming options.
     *
     * @param prefix Prefix used for naming the output variables.
     * @param flag_varname If true, the variable name is included in the
     * generated name.
     * @param flag_qualifier If true, the qualifier is included in the
     * generated name.
     * @param flag_locator If true, a locator is assigned to the output
     * variables.
     * @param locatorOutType Type of locator assigned to the output variables.
     * @param delim Delimiter used to separate the components of the generated
     * name.
     * @param cleanSameLocator If true, variables with the same locator are
     * cleaned beforehand.
     * @return Pointer to the newly created NamingConvention object.
     */
    static NamingConvention* create(
      const String& prefix = "",
      bool flag_varname = true,
      bool flag_qualifier = true,
      bool flag_locator = true,
      const ELoc& locatorOutType = ELoc::fromKey("Z"),
      const String& delim = ".",
      bool cleanSameLocator = true);

    /**
     * @brief Generate names for output variables.
     *
     * The generated names are assigned to the output variables starting at
     * the specified attribute index. Depending on the options of the
     * NamingConvention object, the input variable name and qualifier can be
     * included in the generated names.
     *
     * The output variables can also be assigned the configured locator.
     *
     * @param names Names of the input variables.
     * @param nvar Number of variables.
     * @param dbout Output Db containing the variables to be named.
     * @param iattout_start Index of the first output variable.
     * @param qualifier Qualifier describing the output variables.
     * @param nitems Number of items generated for each variable.
     * @param flagSetLocator If true, assign the configured locator to the
     * output variables.
     * @param locatorShift Shift applied when assigning the locator.
     */
    void setOutput(
      const VectorString& names,
      Id nvar,
      Db* dbout,
      Id iattout_start,
      const String& qualifier = "",
      Id nitems = 1,
      bool flagSetLocator = true,
      Id locatorShift = 0) const;

    /**
     * @brief Generate names for simulation output variables.
     *
     * This method generates names using both the variable and simulation
     * indices. The order of these two indices is controlled by
     * `flagSimuFirst`.
     *
     * For example, for two variables and two simulations, the generated names
     * can be:
     * - V1.S1, V1.S2, V2.S1, V2.S2 when `flagSimuFirst` is false;
     * - S1.V1, S1.V2, S2.V1, S2.V2 when `flagSimuFirst` is true.
     *
     * @param names Names of the input variables.
     * @param nvar Number of variables.
     * @param dbout Output Db containing the simulation variables.
     * @param iattout_start Index of the first output variable.
     * @param nbsimu Number of simulations.
     * @param flagSimuFirst If true, the simulation index is placed before the
     * variable index.
     * @param flagSetLocator If true, assign the configured locator to the
     * output variables.
     * @param locatorShift Shift applied when assigning the locator.
     */
    void setOutputForSimulations(
      const VectorString& names,
      Id nvar,
      Db* dbout,
      Id iattout_start,
      Id nbsimu,
      bool flagSimuFirst = true,
      bool flagSetLocator = true,
      Id locatorShift = 0) const;

    /**
     * @brief Set the delimiter used in generated names.
     *
     * @param delim Delimiter separating the different components of a name.
     */
    void setDelim(const String& delim) { _delim = delim; }

    /**
     * @brief Set the locator type assigned to output variables.
     *
     * @param l Locator type to assign to the output variables.
     */
    void setLocatorOutType(const ELoc& l) { _locatorOutType = l; }

    /**
     * @brief Set the prefix used in generated names.
     *
     * @param prefix Prefix used for generated names.
     */
    void setPrefix(const String& prefix) { _prefix = prefix; }

    /**
     * @brief Set whether variables with the same locator are cleaned.
     *
     * @param cleanSameLocator If true, variables with the same locator are
     * cleaned beforehand.
     */
    void setFlagClean(bool cleanSameLocator)
    {
      _cleanSameLocator = cleanSameLocator;
    }

    /**
     * @brief Assign the configured locator to output variables.
     *
     * The locator is assigned to a set of output variables starting at
     * `iattout_start`.
     *
     * @param dbout Output Db containing the variables.
     * @param iattout_start Index of the first variable receiving the locator.
     * @param nvar Number of variables.
     * @param nitems Number of items associated with each variable.
     * @param locatorShift Shift applied when assigning the locator.
     */
    void setLocators(
      Db* dbout,
      Id iattout_start,
      Id nvar,
      Id nitems = 1,
      Id locatorShift = 0) const;

    /**
     * @brief Test whether the qualifier is included in generated names.
     *
     * @return True if the qualifier is included.
     */
    bool isFlagQualifier() const { return _flagQualifier; }

    /**
     * @brief Test whether the variable name is included in generated names.
     *
     * @return True if the variable name is included.
     */
    bool isFlagVarname() const { return _flagVarname; }

    /**
     * @brief Return the current prefix.
     *
     * @return Prefix used for generated names.
     */
    String getPrefix() const { return _prefix; }

    /**
     * @brief Return the current delimiter.
     *
     * @return Delimiter used to separate the components of generated names.
     */
    String getDelim() const { return _delim; }

    /**
     * @brief Generate a variable name according to the naming convention.
     *
     * This static method provides a way to generate a variable name without
     * creating a NamingConvention object.
     *
     * Depending on the arguments, the generated name can include the input
     * variable name, variable rank, simulation rank and an extension.
     *
     * @param prefix Prefix used in the generated name.
     * @param db Db containing the input variable names. May be nullptr when
     * the variable name is not required.
     * @param ivar Index of the input variable.
     * @param nvar Number of variables.
     * @param isimu Index of the simulation.
     * @param nbsimu Number of simulations.
     * @param extension Additional extension appended to the generated name.
     * @param delim Delimiter used to separate the different components.
     * @return Generated variable name.
     */
    static String getNameEncoded(
      const String& prefix,
      const Db* db = nullptr,
      Id ivar = 0,
      Id nvar = 0,
      Id isimu = 0,
      Id nbsimu = 0,
      const String& extension = "",
      const String& delim = ".");

    /**
     * @brief Activate or deactivate the old naming convention.
     *
     * This method is provided for compatibility with the historical naming
     * convention.
     *
     * @param status If true, activate the old naming convention.
     */
    static void Naming_Old_Style(bool status);

  private:
    void _setNames(
      Db* dbout,
      Id iattout_start,
      const VectorString& names,
      Id nvar,
      const String& qualifier,
      Id nitems) const;

    VectorString _createNames(
      const VectorString& names,
      Id nvar,
      const String& qualifier = "",
      Id nitems = 1) const;

    VectorString _createSimulationNames(
      const VectorString& names,
      Id nvar,
      Id nbsimu,
      bool flagSimuFirst) const;

    static Id _getNameCount(const VectorString& names, Id nvar);

  private:
    String _prefix; //!< String used as 'prefix'
    String _delim; //!< Character used as 'delimiter'
    bool _flagVarname; //!< When TRUE, add the 'variable name'
    bool _flagQualifier; //!< When TRUE, add the 'qualifier'
    bool _flagLocator; //!< When TRUE, assign a locator to the new variables
    ELoc _locatorOutType; //!< Type of locator assigned ('flagLocator' is TRUE)
    bool _cleanSameLocator; //!< Clean variables with same locator beforehand
  };

  // typedef NamingConvention NC;
  class NC: public NamingConvention
  {
  };

} // namespace gstlrn
