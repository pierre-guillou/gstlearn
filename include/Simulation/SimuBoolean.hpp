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

#include "gstlearn_export.hpp"
#include "geoslib_define.h"

#include "Boolean/ModelBoolean.hpp"
#include "Simulation/ASimulation.hpp"
#include "Basic/AStringable.hpp"

class AShape;
class BooleanObject;
class SimuBooleanParam;
class DbGrid;
class Db;

class GSTLEARN_EXPORT SimuBoolean: public ASimulation, public AStringable
{
public:
  SimuBoolean(int nbsimu = 0, int seed = 4324324);
  SimuBoolean(const SimuBoolean &r);
  SimuBoolean& operator=(const SimuBoolean &r);
  virtual ~SimuBoolean();

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  int simulate(Db *dbin,
               DbGrid *dbout,
               ModelBoolean* tokens,
               const SimuBooleanParam& boolparam,
               int iptr_simu,
               int iptr_rank,
               bool verbose = false);

  VectorDouble extractObjects() const;

private:
  void _clearAllObjects();
  int _getNObjects(int mode = 0) const;
  int _getRankUncovered(const Db* db, int rank);
  int _getObjectRank(int mode, int rank);
  int _deleteObject(int mode, Db* dbin);
  int _getAverageCount(const DbGrid* dbout,
                       const ModelBoolean* tokens,
                       const SimuBooleanParam& boolparam) const;
  int _countConditioningPore(const Db* db);
  int _countConditioningGrain(const Db* db);
  int _generatePrimary(Db* dbin,
                       DbGrid* dbout,
                       const ModelBoolean* tokens,
                       const SimuBooleanParam& boolparam,
                       bool verbose = false);
  int _generateSecondary(Db* dbin,
                         DbGrid* dbout,
                         const ModelBoolean* tokens,
                         const SimuBooleanParam& boolparam,
                         bool verbose = false);
  void _projectToGrid(DbGrid* dbout,
                      const SimuBooleanParam& boolparam,
                      int iptr_simu,
                      int iptr_rank);

private:
  std::vector<BooleanObject*> _objlist;
};