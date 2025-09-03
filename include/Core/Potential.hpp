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

#include "gstlearn_export.hpp"
#include "geoslib_define.h"

#include "Basic/VectorNumT.hpp"

namespace gstlrn
{
class Db;
class DbGrid;
class ANeigh;
class Model;

GSTLEARN_EXPORT Id potential_kriging(Db* db,
                                      Db* dbgrd,
                                      Db* dbtgt,
                                      DbGrid* dbout,
                                      Model* model,
                                      ANeigh* neigh,
                                      double nugget_grd   = 0.,
                                      double nugget_tgt   = 0.,
                                      bool flag_pot       = true,
                                      bool flag_grad      = false,
                                      bool flag_trans     = false,
                                      bool flag_save_data = false,
                                      Id opt_part        = 0,
                                      bool verbose        = false);
GSTLEARN_EXPORT Id potential_cov(Model* model,
                                  bool verbose,
                                  Id type1,
                                  const VectorDouble& x10,
                                  const VectorDouble& x1p,
                                  const VectorDouble& tx1,
                                  Id type2,
                                  const VectorDouble& x20,
                                  const VectorDouble& x2p,
                                  const VectorDouble& tx2,
                                  VectorDouble& covtab);
GSTLEARN_EXPORT Id potential_simulate(Db* dbiso,
                                       Db* dbgrd,
                                       Db* dbtgt,
                                       DbGrid* dbout,
                                       Model* model,
                                       ANeigh* neigh,
                                       double nugget_grd   = 0.,
                                       double nugget_tgt   = 0.,
                                       double dist_tempere = TEST,
                                       bool flag_trans     = false,
                                       Id seed            = 135674,
                                       Id nbsimu          = 1,
                                       Id nbtuba          = 100,
                                       bool verbose        = false);
GSTLEARN_EXPORT Id potential_xvalid(Db* dbiso,
                                     Db* dbgrd,
                                     Db* dbtgt,
                                     Model* model,
                                     ANeigh* neigh,
                                     double nugget_grd   = 0.,
                                     double nugget_tgt   = 0.,
                                     bool flag_dist_conv = false,
                                     bool verbose        = false);
}