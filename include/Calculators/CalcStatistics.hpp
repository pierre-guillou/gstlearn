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

#include "ACalcDbToDb.hpp"

#include "geoslib_define.h"


class Db;
class DbGrid;
class EStatOption;

class GSTLEARN_EXPORT CalcStatistics: public ACalcDbToDb
{
public:
  CalcStatistics();
  CalcStatistics(const CalcStatistics &r) = delete;
  CalcStatistics& operator=(const CalcStatistics &r) = delete;
  virtual ~CalcStatistics();

  void setFlagStats(bool flagStats) { _flagStats = flagStats; }
  void setRadius(int radius) { _radius = radius; }
  void setOper(const EStatOption &oper) { _oper = oper; }

  void setFlagRegr(bool flagRegr) { _flagRegr = flagRegr; }
  void setFlagCste(bool flagCste) { _flagCste = flagCste; }
  void setNamaux(const VectorString &namaux) { _namaux = namaux; }
  void setRegrMode(int regrMode) { _regrMode = regrMode; }

private:
  virtual bool _check() override;
  virtual bool _preprocess() override;
  virtual bool _run() override;
  virtual bool _postprocess() override;
  virtual void _rollback() override;

private:
  int    _iattOut;

  bool _flagStats;
  EStatOption _oper;
  int _radius;

  bool _flagRegr;
  bool _flagCste;
  int  _regrMode;
  VectorString _namaux;
};

GSTLEARN_EXPORT int dbStatisticsOnGrid(Db *db,
                                       DbGrid *dbgrid,
                                       const EStatOption &oper,
                                       int radius = 0,
                                       const NamingConvention &namconv = NamingConvention(
                                           "Stats"));
GSTLEARN_EXPORT int dbRegression(Db *db1,
                                 Db *db2,
                                 int mode,
                                 const String &name,
                                 const VectorString &namaux,
                                 bool flagCste = true,
                                 const NamingConvention &namconv = NamingConvention(
                                     "Regr"));