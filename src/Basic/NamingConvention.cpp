/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstLearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "Basic/NamingConvention.hpp"
#include "Basic/String.hpp"
#include "Db/Db.hpp"

#include <string>

namespace gstlrn
{
  // Default value for the Style used for Variable encoding
  bool Old_Style = false;

  void NamingConvention::Naming_Old_Style(bool status)
  {
    Old_Style = status;
  }

  NamingConvention::NamingConvention(
    const String& prefix,
    bool flag_varname,
    bool flag_qualifier,
    bool flag_locator,
    const ELoc& locatorOutType,
    const String& delim,
    bool cleanSameLocator)
    : AStringable()
    , _prefix(prefix)
    , _delim(delim)
    , _flagVarname(flag_varname)
    , _flagQualifier(flag_qualifier)
    , _flagLocator(flag_locator)
    , _locatorOutType(locatorOutType)
    , _cleanSameLocator(cleanSameLocator)
  {
  }

  NamingConvention::NamingConvention(const NamingConvention& m)
    : AStringable(m)
    , _prefix(m._prefix)
    , _delim(m._delim)
    , _flagVarname(m._flagVarname)
    , _flagQualifier(m._flagQualifier)
    , _flagLocator(m._flagLocator)
    , _locatorOutType(m._locatorOutType)
    , _cleanSameLocator(m._cleanSameLocator)
  {
  }

  NamingConvention& NamingConvention::operator=(const NamingConvention& m)
  {
    if (this != &m)
    {
      AStringable::operator=(m);
      _prefix = m._prefix;
      _locatorOutType = m._locatorOutType;
      _flagVarname = m._flagVarname;
      _flagQualifier = m._flagQualifier;
      _flagLocator = m._flagLocator;
      _delim = m._delim;
      _cleanSameLocator = m._cleanSameLocator;
    }
    return *this;
  }

  NamingConvention::~NamingConvention() {}

  NamingConvention* NamingConvention::create(
    const String& prefix,
    bool flag_varname,
    bool flag_qualifier,
    bool flag_locator,
    const ELoc& locatorOutType,
    const String& delim,
    bool cleanSameLocator)
  {
    return new NamingConvention(
      prefix, flag_varname, flag_qualifier, flag_locator, locatorOutType, delim,
      cleanSameLocator);
  }

  void NamingConvention::setOutput(
    const VectorString& names,
    Id nvar,
    Db* dbout,
    Id iattout_start,
    const String& qualifier,
    Id nitems,
    bool flagSetLocator,
    Id locatorShift) const
  {
    if (iattout_start < 0) return;

    auto nameloc = names;
    if (nameloc.empty())
    {
      // 'names' is not provided, 'nvar' prevails (if not defined, it is set to 1)
      if (nvar <= 0) nvar = 1;
    }
    else
    {
      // 'names' is provided
      auto namesize = static_cast<Id>(nameloc.size());
      if (nvar <= 0)
      {
        // If 'nvar' is not defined, argument 'names' prevails
        nvar = namesize;
      }
      else
      {
        // 'names' and 'nvar' are both defined: 'nvar' prevails
        if (namesize == 1 && nvar > 1)
        {
          // Particular case where 'nvar' > 1 but 'names' contains a single name: the name is expanded
          nameloc = generateMultipleNames(names[0], nvar);
        }
        else
        {
          // 'names' and 'nvar' are both defined: 'names' is reset to 'nvar' if needed
          nameloc.resize(nvar);
        }
      }
    }

    _setNames(dbout, iattout_start, nameloc, nvar, qualifier, nitems);

    if (flagSetLocator)
      setLocators(dbout, iattout_start, nvar, nitems, locatorShift);
  }

  void NamingConvention::setOutputForSimulations(
    const VectorString& names,
    Id nvar,
    Db* dbout,
    Id iattout_start,
    Id nbsimu,
    bool flagSimuFirst,
    bool flagSetLocator,
    Id locatorShift) const
  {
    if (iattout_start < 0) return;

    if (names.empty())
    {
      if (nvar <= 0) nvar = 1;
    }
    else
    {
      nvar = static_cast<Id>(names.size());
    }

    // Create simulation names
    VectorString outnames =
      _createSimulationNames(names, nvar, nbsimu, flagSimuFirst);

    // Set the names in the database
    Id ntotal = nvar * nbsimu;
    for (Id i = 0; i < ntotal; i++)
    {
      dbout->setNameByUID(iattout_start + i, outnames[i]);
    }

    if (flagSetLocator)
    {
      if (_flagLocator && _locatorOutType != ELoc::UNDEFINED)
      {
        // Erase already existing locators of the same Type
        if (_cleanSameLocator && locatorShift == 0)
          dbout->clearLocators(_locatorOutType);

        // Set the locator for all variables
        for (Id i = 0; i < ntotal; i++)
          dbout->setLocatorByUID(
            iattout_start + i, _locatorOutType, i + locatorShift);
      }
    }
  }

  void NamingConvention::setLocators(
    Db* dbout,
    Id iattout_start,
    Id nvar,
    Id nitems,
    Id locatorShift) const
  {
    if (!_flagLocator || _locatorOutType == ELoc::UNDEFINED) return;

    // Erase already existing locators of the same Type
    // (this is only done if you are not precisely adding higher order version for given locator)
    if (_cleanSameLocator && locatorShift == 0)
      dbout->clearLocators(_locatorOutType);

    // Set the locator for all variables
    for (Id ecr = 0; ecr < nvar * nitems; ecr++)
      dbout->setLocatorByUID(
        iattout_start + ecr, _locatorOutType, ecr + locatorShift);
  }

  Id NamingConvention::_getNameCount(const VectorString& names, Id nvar)
  {
    if (nvar <= 0)
    {
      // Argument 'nvar' is not defined yet
      if (names.empty()) return 1;
      return static_cast<Id>(names.size());
    }

    // Argument 'nvar' is provided: is it consistent with 'names'
    if (names.empty()) return nvar;

    // Both 'nvar' and 'names' are provided. For safety reasons,
    // the number of variables is the minimum between the two
    return MIN(nvar, static_cast<Id>(names.size()));
  }

  void NamingConvention::_setNames(
    Db* dbout,
    Id iattout_start,
    const VectorString& names,
    Id nvar,
    const String& qualifier,
    Id nitems) const
  {
    auto nloc = _getNameCount(names, nvar);
    VectorString outnames = _createNames(names, nloc, qualifier, nitems);
    correctNamesForDuplicates(outnames, dbout->getAllNames());

    Id ecr = 0;
    for (Id ivar = 0; ivar < nloc; ivar++)
    {
      for (Id item = 0; item < nitems; item++)
      {
        dbout->setNameByUID(iattout_start + ecr, outnames[ecr]);
        ecr++;
      }
    }
  }

  VectorString NamingConvention::_createNames(
    const VectorString& names,
    Id nvar,
    const String& qualifier,
    Id nitems) const
  {
    VectorString outnames;

    for (Id ivar = 0; ivar < nvar; ivar++)
    {
      // Variable 'local' defined for each variable, is:
      // - extracted from the array 'names' (if defined)
      // - generated as the rank of the variable (if several)
      String loc_varname;
      String loc_number;
      if (_flagVarname)
      {
        if (static_cast<Id>(names.size()) == nvar) loc_varname = names[ivar];
        if (loc_varname.empty() && nvar > 1)
        {
          if (Old_Style)
            loc_varname = std::to_string(ivar + 1);
          else
            loc_varname = concatenateString("V", ivar + 1, "");
        }
      }
      else
      {
        // Build the rank from the variable number (possibly overwritten by item number)
        if (nvar > 1)
        {
          if (Old_Style)
            loc_number = std::to_string(ivar + 1);
          else
            loc_number = concatenateString("V", ivar + 1, "");
        }
      }

      for (Id item = 0; item < nitems; item++)
      {
        String loc_qualifier;

        if (_flagQualifier)
        {
          loc_qualifier = qualifier;
          if (nitems > 1)
          {
            if (Old_Style)
              loc_number = std::to_string(item + 1);
            else
              loc_number = concatenateString("S", item + 1, "");
          }
        }

        // Compose the variable name
        String name = concatenateStrings(
          _delim, _prefix, loc_varname, loc_qualifier, loc_number);

        if (name.empty()) name = "Dummy";
        outnames.push_back(name);
      }
    }
    return outnames;
  }

  VectorString NamingConvention::_createSimulationNames(
    const VectorString& names,
    Id nvar,
    Id nbsimu,
    bool flagSimuFirst) const
  {
    if (nvar <= 0 || nbsimu <= 0) return {};

    VectorString outnames;

    // Determine variable names
    VectorString varnames;
    bool isConditional = !names.empty();

    if (isConditional)
    {
      // Conditional simulation: use provided variable names
      varnames = names;
      if (static_cast<Id>(varnames.size()) != nvar && nvar > 0)
        varnames.resize(nvar);
    }
    else
    {
      // Non-conditional simulation: use V1, V2, ... format
      for (Id ivar = 0; ivar < nvar; ivar++)
      {
        String varname = 'V' + std::to_string(ivar + 1);
        varnames.push_back(varname);
      }
    }

    // Create names based on storage order
    if (flagSimuFirst)
    {
      // Simulation varies first: V1.S1, V1.S2, ..., V1.Sn, V2.S1, V2.S2, ...
      for (Id ivar = 0; ivar < nvar; ivar++)
      {
        for (Id isimu = 0; isimu < nbsimu; isimu++)
        {
          String simuname = 'S' + std::to_string(isimu + 1);
          String name =
            concatenateStrings(_delim, _prefix, varnames[ivar], simuname);
          if (name.empty()) name = "Dummy";
          outnames.push_back(name);
        }
      }
    }
    else
    {
      // Variable varies first: V1.S1, V2.S1, ..., Vn.S1, V1.S2, V2.S2, ...
      for (Id isimu = 0; isimu < nbsimu; isimu++)
      {
        for (Id ivar = 0; ivar < nvar; ivar++)
        {
          String simuname = 'S' + std::to_string(isimu + 1);
          String name =
            concatenateStrings(_delim, _prefix, varnames[ivar], simuname);
          if (name.empty()) name = "Dummy";
          outnames.push_back(name);
        }
      }
    }

    return outnames;
  }

  String NamingConvention::toString(const AStringFormat* /*strfmt*/) const
  {
    std::stringstream sstr;

    sstr << toStrTitle(0, "Naming Convention");
    sstr << "- Prefix  = " << _prefix << std::endl;
    sstr << "- Delimitor = '" << _delim << "'" << std::endl;
    sstr << "- Add the Variable Name = " << _flagVarname << std::endl;
    sstr << "- Add the Qualifier     = " << _flagQualifier << std::endl;
    sstr << "- Assign a Locator to output variables = " << _flagLocator
         << std::endl;
    sstr << "- Type of assigned locator = " << _locatorOutType.getDescr()
         << std::endl;
    sstr << "- Clean any other similar locator = " << _cleanSameLocator
         << std::endl;

    return sstr.str();
  }

  String NamingConvention::getNameEncoded(
    const String& prefix,
    const Db* db,
    Id ivar,
    Id nvar,
    Id isimu,
    Id nbsimu,
    const String& extension,
    const String& delim)
  {
    String loc_varname;
    if (db != nullptr)
    {
      if (db->getNLoc(ELoc::Z) > 0)
        loc_varname = db->getNameByLocator(ELoc::Z, ivar);
    }
    else
    {
      if (ivar < 0)
      {
        loc_varname = "*";
      }
      else
      {
        if (nvar > 1 && ivar > 0)
          loc_varname = concatenateString("V", ivar, "");
      }
    }

    String loc_qualifier;
    if (!extension.empty())
    {
      loc_qualifier = extension;
    }
    else
    {
      if (isimu < 0)
      {
        loc_qualifier = "*";
      }
      else
      {
        if (nbsimu > 1 && isimu > 0)
          loc_qualifier = concatenateString("S", isimu, "");
      }
    }

    String name = concatenateStrings(delim, prefix, loc_varname, loc_qualifier);
    if (name.empty()) name = "Dummy";

    return name;
  }

} // namespace gstlrn
