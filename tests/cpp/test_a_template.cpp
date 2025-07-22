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
#include "Basic/OptCst.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Estimation/CalcKriging.hpp"
#include "Model/Model.hpp"
#include "Neigh/NeighMoving.hpp"
#include "geoslib_define.h"

using namespace gstlrn;

int main(int argc, char* argv[])
{
  DECLARE_UNUSED(argc);
  DECLARE_UNUSED(argv);

  OptCst::define(ECst::NTROW, -1);
  OptCst::define(ECst::NTCOL, -1);

  int ndat = 100;
  Db* db   = Db::createFillRandom(ndat, 2, 1, 0, 5, 0., 0.3);

  double range = 0.3;
  double sill  = 2.;
  Model* model = Model::createFromParam(ECov::SPHERICAL, range, sill);

  NeighMoving* neigh = NeighMoving::create(false, 10, 1., 1, 1, 10);

  neigh->setBallSearch(true, 5);

  (void)xvalid(db, model, neigh, false);

  DbStringFormat* dbfmt = DbStringFormat::create(FLAG_RESUME | FLAG_ARRAY | FLAG_VARS);
  db->display(dbfmt);

  return (0);
}
