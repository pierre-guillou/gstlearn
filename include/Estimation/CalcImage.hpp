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

#include "Calculators/ACalcInterpolator.hpp"
#include "Morpho/EMorpho.hpp"

class DbGrid;
class ANeighParam;

class GSTLEARN_EXPORT CalcImage: public ACalcInterpolator
{
public:
  CalcImage();
  CalcImage(const CalcImage &r) = delete;
  CalcImage& operator=(const CalcImage &r) = delete;
  virtual ~CalcImage();

  void setFlagFilter(bool flagFilter) { _flagFilter = flagFilter; }

  void setFlagMorpho(bool flagMorpho) { _flagMorpho = flagMorpho; }
  void setOper(const EMorpho& oper) { _oper = oper; }
  void setOption(int option) { _option = option; }
  void setRadius(const VectorInt& radius) { _radius = radius; }
  void setVmin(double vmin) { _vmin = vmin; }
  void setVmax(double vmax) { _vmax = vmax; }
  void setVerbose(bool verbose) { _verbose = verbose; }

private:
  virtual bool _check() override;
  virtual bool _preprocess() override;
  virtual bool _run() override;
  virtual bool _postprocess() override;
  virtual void _rollback() override;
  virtual int  _getNVar() const override;

private:
  int    _iattOut;

  bool _flagFilter;

  bool _flagMorpho;
  EMorpho _oper;
  double _vmin;
  double _vmax;
  int _option;
  VectorInt _radius;
  bool _verbose;
};

GSTLEARN_EXPORT int krimage(DbGrid *dbgrid,
                            Model *model,
                            NeighImage *neighparam,
                            const NamingConvention& namconv = NamingConvention("Filtering"));

GSTLEARN_EXPORT int morpho(DbGrid *dbgrid,
                           const EMorpho &oper,
                           double vmin = 0.,
                           double vmax = 1.5,
                           int option = 0,
                           const VectorInt &radius = VectorInt(),
                           bool verbose = false,
                           const NamingConvention &namconv = NamingConvention(
                               "Morpho"));