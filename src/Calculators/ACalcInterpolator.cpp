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
#include "geoslib_f.h"

#include "Calculators/ACalcInterpolator.hpp"
#include "Calculators/ACalculator.hpp"
#include "Db/Db.hpp"
#include "Model/Model.hpp"
#include "Neigh/ANeighParam.hpp"

ACalcInterpolator::ACalcInterpolator()
    : ACalculator(),
      _dbin(nullptr),
      _dbout(nullptr),
      _model(nullptr),
      _neighparam(nullptr),
      _namconv(),
      _listVariablePermDbIn(),
      _listVariablePermDbOut(),
      _listVariableTempDbIn(),
      _listVariableTempDbOut()
{
}

ACalcInterpolator::~ACalcInterpolator()
{
}

int ACalcInterpolator::_getNDim() const
{
  int ndim = 0;
  if (_dbin != nullptr)
  {
    if (ndim > 0)
    {
      if (ndim != _dbin->getNDim()) return -1;
    }
    else
    {
      ndim = _dbin->getNDim();
    }
  }

  if (_dbout!= nullptr)
  {
    if (ndim > 0)
    {
      if (ndim != _dbout->getNDim()) return -1;
    }
    else
    {
      ndim = _dbout->getNDim();
    }
  }

  if (_model != nullptr)
  {
    if (ndim > 0)
    {
      if (ndim != _model->getDimensionNumber()) return -1;
    }
    else
    {
      ndim = _model->getDimensionNumber();
    }
  }
  return ndim;
}

int ACalcInterpolator::_getNCova() const
{
  if (_model == nullptr) return -1;
  return getModel()->getCovaNumber();
}

int ACalcInterpolator::_getNVar() const
{
  if (_model == nullptr) return -1;
  return getModel()->getVariableNumber();
}

bool ACalcInterpolator::_check() const
{
  /**************************************************/
  /* Cross-checking the Space Dimension consistency */
  /**************************************************/

  int ndim = 0;
  if (_dbin != nullptr)
  {
    if (ndim > 0)
    {
      if (ndim != _dbin->getNDim())
      {
        messerr("Inconsistent Space dimension:");
        messerr("- Current dimension = %d",ndim);
        messerr("- Space Dimension of 'dbin' = %d",_dbin->getNDim());
        return false;
      }
    }
    else
    {
      ndim = _dbin->getNDim();
    }
  }

  if (_dbout!= nullptr)
  {
    if (ndim > 0)
    {
      if (ndim != _dbout->getNDim())
      {
        messerr("Inconsistent Space dimension:");
        messerr("- Current dimension = %d",ndim);
        messerr("- Space Dimension of 'dbout' = %d",_dbout->getNDim());
        return false;
      }
    }
    else
    {
      ndim = _dbout->getNDim();
    }
  }

  if (_model != nullptr)
  {
    if (ndim > 0)
    {
      if (ndim != _model->getDimensionNumber())
      {
        messerr("Inconsistent Space dimension:");
        messerr("- Current dimension = %d",ndim);
        messerr("- Space Dimension of 'model' = %d",_model->getDimensionNumber());
        return false;
      }
    }
    else
    {
      ndim = _model->getDimensionNumber();
    }
  }

  if (_neighparam != nullptr)
  {
    if (ndim > 0)
    {
      if (ndim != _neighparam->getNDim())
      {
        messerr("Inconsistent Space dimension:");
        messerr("- Current dimension = %d",ndim);
        messerr("- Space Dimension of '_neighparam' = %d",_neighparam->getNDim());
        return false;
      }
    }
    else
    {
      ndim = _neighparam->getNDim();
    }
  }

  /**************************************************/
  /* Cross-Checking the Variable Number consistency */
  /**************************************************/

  int nvar = 0;
  if (_dbin != nullptr)
  {
    if (nvar > 0)
    {
      if (nvar != _dbin->getVariableNumber())
      {
        messerr("Inconsistent the Variable Number:");
        messerr("- Current number = %d",nvar);
        messerr("- Number of variables in 'dbin' = %d",_dbin->getVariableNumber());
        return false;
      }
    }
    else
    {
      nvar = _dbin->getVariableNumber();
    }
  }

  if (_model != nullptr)
  {
    if (nvar > 0)
    {
      if (nvar != _model->getVariableNumber())
      {
        messerr("Inconsistent the Variable Number:");
        messerr("- Current number = %d",nvar);
        messerr("- Number of variables in 'model' = %d",_model->getVariableNumber());
        return false;
      }
    }
    else
    {
      nvar = _model->getVariableNumber();
    }
  }

  /************************************************************/
  /* Cross-Checking the number of External Drifts consistency */
  /************************************************************/

  int nfex = 0;
  if (_model != nullptr)
  {
    nfex = _model->getExternalDriftNumber();
    if (nfex > 0)
    {
      if (_dbin != nullptr)
      {
        if (_dbin->getExternalDriftNumber() != nfex)
        {
          messerr("The _model requires %d external drift(s)", nfex);
          messerr("but the input Db refers to %d external drift variables",
                  _dbin->getExternalDriftNumber());
          return false;
        }
      }

      if (_dbout != nullptr)
      {
        if (_dbout->getExternalDriftNumber() != nfex)
        {
          messerr("The _model requires %d external drift(s)", nfex);
          messerr("but the output Db refers to %d external drift variables",
                  _dbout->getExternalDriftNumber());
          return false;
        }
      }
    }
  }

  /**************************************/
  /* Checking the Validity of the _model */
  /**************************************/

  if (_model != nullptr)
  {
    if (_model->getCovaNumber() <= 0)
    {
      messerr("The number of covariance must be positive");
      return false;
    }
  }

  /*********************************/
  /* Calculate the field extension */
  /*********************************/

  if (_model != nullptr)
  {
    VectorDouble db_mini(ndim, TEST);
    VectorDouble db_maxi(ndim, TEST);
    if (_dbin != nullptr) _dbin->getExtensionInPlace(db_mini, db_maxi);
    if (_dbout != nullptr) _dbout->getExtensionInPlace(db_mini, db_maxi);
    _model->setField(ut_vector_extension_diagonal(db_mini, db_maxi));
  }
  return true;
}

/*****************************************************************************/
/*!
 **  Derive some information from the Output db (if Grid) to the Input Db
 **
 ** \return  Error return code
 **
 ** \param[in]  mode        1 for allocation; -1 for deallocation
 ** \param[in]  locatorType Type of the pointer (ELoc)
 **
 ** \remark This function is only valid when the Output Db is a grid
 ** \remark However, in case of a Point output Db, this function should not
 ** \remark be used: the external drift functions should already be present
 ** \remark in the output Db.
 ** \remark If this is not the case, an error is issued.
 **
 ** \remark When called with mode=1, the new variables are added
 ** \remark When called with mode=-1, the variables are deleted (by type)
 **
 *****************************************************************************/
int ACalcInterpolator::_expandInformation(int mode, const ELoc &locatorType)
{
  if (_dbin == nullptr || _dbout == nullptr) return 0;
  int nechin = _dbin->getSampleNumber();
  int ninfo;
  if (_dbout->isGrid() && locatorType == ELoc::X)
    ninfo = _dbout->getNDim();
  else
    ninfo = _dbout->getFromLocatorNumber(locatorType);
  if (ninfo <= 0) return 0;

  /* Case when the Output Db is not a grid */

  if (! _dbout->isGrid())
  {
    if (_dbin->getFromLocatorNumber(locatorType) == ninfo) return 0;
    messerr("The Output Db is not a Grid file");
    messerr("The Input Db does not contain the correct number of External Drifts");
    return 1;
  }
  DbGrid* dbgrid = dynamic_cast<DbGrid*>(_dbout);

  /* Dispatch */

  if (mode > 0)
  {
    VectorDouble tab(nechin, 0.);
    VectorInt iatts = dbgrid->getUIDsByLocator(locatorType);
    _storeInVariableList(1, 2, iatts);
    if (migrateByAttribute(dbgrid, _dbin, iatts)) return 1;
  }
  else
  {
    _dbin->deleteColumnsByLocator(locatorType);
  }
  return 0;
}

/**
 * Returns a pointer to the relevant Db and issue a message if not defined
 * @param whichDb 1 for 'dbin'  and 2 for 'dbout'
 * @return A pointer to the Db or nullptr
 */
Db* ACalcInterpolator::_whichDb(int whichDb)
{
  Db* db;
  if (whichDb == 1)
    db = _dbin;
  else
    db = _dbout;
  if (db == nullptr)
  {
    messerr("Impossible to add variables in non-defined Db");
  }
  return db;
}

/**
 * Store the IUID of the new variable in the relevant internal list
 * @param whichDb 1 for variable in 'dbin'; 2 for variable in 'dbout'
 * @param status  1 for variables to be stored; 2 for Temporary variable
 * @param iuid    Vector of UIDs of the new variable
 */
void ACalcInterpolator::_storeInVariableList(int whichDb, int status, const VectorInt& iuids)
{
  int number = (int) iuids.size();
  if (number <= 0) return;

  if (whichDb == 1)
  {
    if (status == 1)
    {
      for (int i = 0; i < number; i++)
        _listVariablePermDbIn.push_back(iuids[i]);
    }
    else
    {
      for (int i = 0; i < number; i++)
        _listVariableTempDbIn.push_back(iuids[i]);
    }
  }
  else
  {
    if (status == 1)
    {
      for (int i = 0; i < number; i++)
        _listVariablePermDbOut.push_back(iuids[i]);
    }
    else
    {
      for (int i = 0; i < number; i++)
        _listVariableTempDbOut.push_back(iuids[i]);
    }
  }
}
int ACalcInterpolator::_addVariableDb(int whichDb,
                                      int status,
                                      const ELoc &locatorType,
                                      int number,
                                      double valinit)
{
  Db* db = _whichDb(whichDb);
  if (db == nullptr) return -1;
  int iuid = db->addColumnsByConstant(number, valinit, String(), locatorType);
  if (iuid < 0) return -1;
  VectorInt iuids = ut_ivector_sequence(number, iuid);
  _storeInVariableList(whichDb, status, iuids);
  return iuid;
}

void ACalcInterpolator::_renameVariable(const ELoc &locatorType,
                                        int nvar,
                                        int iptr,
                                        const String &name,
                                        int count)
{
  _namconv.setNamesAndLocators(_dbin, locatorType, nvar, _dbout, iptr, name, count);
}

void ACalcInterpolator::_cleanVariableDb(int status)
{
  // Dispatch

  if (status == 1)
  {
    // In 'dbin'
    if (!_listVariablePermDbIn.empty())
    {
      for (int i = 0; i < (int) _listVariablePermDbIn.size(); i++)
        _dbin->deleteColumnByUID(_listVariablePermDbIn[i]);
    }
    _listVariablePermDbIn.clear();

    // In 'dbout'
    if (!_listVariablePermDbOut.empty())
    {
      for (int i = 0; i < (int) _listVariablePermDbOut.size(); i++)
        _dbout->deleteColumnByUID(_listVariablePermDbOut[i]);
    }
    _listVariablePermDbOut.clear();
  }
  else
  {
    // In 'dbin'
    if (!_listVariableTempDbIn.empty())
    {
      for (int i = 0; i < (int) _listVariableTempDbIn.size(); i++)
        _dbin->deleteColumnByUID(_listVariableTempDbIn[i]);
    }
    _listVariableTempDbIn.clear();

    // In 'dbout'
    if (!_listVariableTempDbOut.empty())
    {
      for (int i = 0; i < (int) _listVariableTempDbOut.size(); i++)
        _dbout->deleteColumnByUID(_listVariableTempDbOut[i]);
    }
    _listVariableTempDbOut.clear();
  }

}