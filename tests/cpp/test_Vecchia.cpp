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
#include "Anamorphosis/AnamHermite.hpp"
#include "Anamorphosis/CalcAnamTransform.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/File.hpp"
#include "Basic/OptCst.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Estimation/CalcGlobal.hpp"
#include "Estimation/CalcImage.hpp"
#include "Estimation/CalcKriging.hpp"
#include "Estimation/Vecchia.hpp"
#include "Matrix/MatrixT.hpp"
#include "Model/Model.hpp"
#include "Neigh/NeighImage.hpp"
#include "Neigh/NeighMoving.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"
#include "Tree/Ball.hpp"

#include "geoslib_f.h"

using namespace gstlrn;
/****************************************************************************/
/*!
 ** Main Program
 **
 *****************************************************************************/
int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  // Global parameters
  int nb_vecchia = 3;
  int mode       = 0;
  bool verbose   = false;
  OptCst::define(ECst::NTCOL, -1);
  OptCst::define(ECst::NTROW, -1);

  int ndat = 5;
  Db* db   = Db::createFillRandom(ndat, 2, 1);
  // db->setZVariable(3, 0, TEST);
  DbStringFormat dbfmt1(FLAG_ARRAY);
  db->display(&dbfmt1);

  double range = 0.5;
  Model* model = Model::createFromParam(ECov::EXPONENTIAL, range);

  int nx       = 2;
  DbGrid* grid = DbGrid::create({nx, nx}, {1. / nx, 1. / nx});

  if (mode == 0 || mode == 1)
  {
    verbose = true;
    mestitle(0, "Checking Vecchia Class");
    Vecchia V          = Vecchia(model, nb_vecchia, db);
    MatrixT<int> Ranks = findNN(db, nullptr, nb_vecchia + 1, false, verbose);
    (void)V.computeLower(Ranks, verbose);
  }

  if (mode == 0 || mode == 2)
  {
    verbose = true;
    mestitle(0, "Kriging with Vecchia approximation");
    krigingVecchia(db, grid, model, nb_vecchia, verbose);

    // Get some statistics for check printout
    DbStringFormat dbfmt(FLAG_STATS, {"Vecchia*"});
    grid->display(&dbfmt);
  }

  if (mode == 0 || mode == 3)
  {
    delete db;
    db = Db::createFillRandom(20000, 2, 1);

    const double result = logLikelihoodVecchia(db, model, 20, false);
    message("Log-likelihood = %f\n", result);
  }
  // ====================== Free pointers ==================================
  delete db;
  delete grid;
  delete model;

  return (0);
}
