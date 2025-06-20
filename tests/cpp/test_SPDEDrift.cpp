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

#include "Basic/File.hpp"
#include "Model/ConsItem.hpp"
#include "Model/Constraints.hpp"
#include "Model/Model.hpp"
#include "Variogram/Vario.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Db/DbStringFormat.hpp"
#include "API/SPDE.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Estimation/CalcKriging.hpp"
#include "utils.hpp"

int main(int argc, char *argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  bool verbose = true;
  // This crashes under MingGW/Windows due to compatibility issue with getTestData
  //ASerializable::setPrefixName("test_SPDEDrift-");

  String filename;

  filename = gstlearn::getTestData("Scotland","temperatures.ascii");
  Db* temperatures = Db::createFromNF(filename,verbose);
  temperatures->setLocator("January_temp", ELoc::Z, 0);
  temperatures->display();

  filename = gstlearn::getTestData("Scotland","grid.ascii");
  DbGrid* grid = DbGrid::createFromNF(filename,verbose);
  grid->display();

  filename = gstlearn::getTestData("Scotland","model.ascii");
  Model* model = Model::createFromNF(filename,verbose);

  model->display();

  filename = gstlearn::getTestData("Scotland","vario.ascii");
  Vario* vario = Vario::createFromNF(filename,verbose);

  vario->display();

  auto structs = {ECov::NUGGET,ECov::MATERN};
  ConsItem consNug = ConsItem::define(EConsElem::SILL,0,0,0, EConsType::UPPER,0.1);
  ConsItem consParam = ConsItem::define(EConsElem::PARAM,1, 0, 0, EConsType::EQUAL,1.);
  Constraints constraints;
  constraints.addItem(&consNug);
  constraints.addItem(&consParam);

  //OptDbg::define(EDbg::CONVERGE);
  (void)model->fit(vario, structs, constraints);
  model->display();

  NeighUnique* neighU = NeighUnique::create();
  neighU->display();

  //////////////////////
  /// Kriging using SPDE
  int useCholesky = 0;
  SPDE spde(model,grid,temperatures,ESPDECalcMode::KRIGING,nullptr,useCholesky);
  VectorDouble coeffs = spde.getCoeffs();
  spde.compute(grid);

  ///////////////////////
  /// Traditional Kriging
  kriging(temperatures, grid, model, neighU);

  (void) grid->dumpToNF("Grid.ascii",verbose);

  DbStringFormat dbfmt(FLAG_STATS,{"spde*estim"});
  grid->display(&dbfmt);

  delete temperatures;
  delete grid;
  delete model;
  delete neighU;
  delete vario;

  return 0;
}

