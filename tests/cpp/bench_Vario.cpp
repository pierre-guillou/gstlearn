/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/*                                                                            */
/* This file is meant to demonstrate the process of using PGS                 */
/*                                                                            */
/******************************************************************************/
#include "geoslib_f.h"

#include "Enum/ECalcVario.hpp"
#include "Enum/ECov.hpp"

#include "Variogram/Vario.hpp"
#include "Model/Model.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/File.hpp"
#include "Basic/Timer.hpp"
#include "Db/Db.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"

/****************************************************************************/
/*!
** Main Program for bench marking the variogram calculation
** - on a set of isolated points
** - on a regular grid
**
*****************************************************************************/
int main(int argc, char *argv[])
{
  Timer timer;
  bool verbose = false;

  // Standard output redirection to file
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc <= 1);

  int ndim = 2;
  defineDefaultSpace(ESpaceType::RN, ndim);

  // Creating a Point Data base in the 1x1 square with 'nech' samples
  int nech = 2000;
  Db* db = Db::createFromBox(nech,{0.,0.},{1.,1.}, 3242);

  // Creating a grid covering the same space
  VectorInt nx = { 200, 200 };
  VectorDouble dx = { 0.005, 0.005 };
  DbGrid* grid = DbGrid::create(nx, dx);

  // Creating a Model(s) for simulating a variable
  Model* model = Model::createFromParam(ECov::BESSEL_K,0.2);

  // Perform a non-conditional simulation on the Db and on the Grid
  (void) simtub(nullptr,db,model);
  (void) simtub(nullptr,grid,model);

  // ===============
  // On Data samples
  // ===============

  mestitle(1, "Experimental variogram on Data Samples");
  timer.reset();
  int nlag = 20;
  VarioParam* varioparamP = VarioParam::createMultiple(4, nlag, 0.5 / nlag);
  Vario* varioP = Vario::computeFromDb(varioparamP,db,ECalcVario::VARIOGRAM);
  timer.displayIntervalMilliseconds("Variogram on Isolated Points", 3600);
  if (verbose) varioP->display();

  // ===============
  // On Grid samples
  // ===============

  mestitle(1, "Experimental variogram on Grid");
  timer.reset();
  VarioParam* varioparamG = VarioParam::createMultipleFromGrid(grid, nlag);
  Vario* varioG = Vario::computeFromDb(varioparamG, grid, ECalcVario::VARIOGRAM);
  timer.displayIntervalMilliseconds("Variogram on Regular Grid", 1500);
  if (verbose) varioG->display();

  // ==========================================
  // Calculating Variogram Map on Isolated Data
  // ==========================================

  mestitle(1, "Variogram Map on Isolated Data");
  timer.reset();
  Db* vmapP = db_vmap_compute(db, ECalcVario::VARIOGRAM, {50,50});
  timer.displayIntervalMilliseconds("Variogram Map on Isolated Points", 2400);

  // =================================
  // Calculating Variogram Map on Grid
  // =================================

  mestitle(1, "Variogram Map on Grid");
  timer.reset();
  Db* vmapG = db_vmap_compute(grid, ECalcVario::VARIOGRAM, {100,100});
  timer.displayIntervalMilliseconds("Variogram Map on Regular Grid", 100);

  delete db;
  delete grid;
  delete model;
  delete varioparamP;
  delete varioparamG;
  delete varioP;
  delete varioG;
  delete vmapP;
  delete vmapG;

  return (0);
}