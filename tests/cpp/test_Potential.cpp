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

// This test is meant to demonstrate the Potential Model
// through estimation, cross-validation and simulations

#include "Enum/ECov.hpp"
#include "Enum/ESpaceType.hpp"

#include "Space/ASpaceObject.hpp"
#include "Space/SpaceRN.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Model/Model.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Basic/File.hpp"
#include "Basic/OptDbg.hpp"
#include "Core/Potential.hpp"
#include "Covariances/CovAniso.hpp"

/****************************************************************************/
/*!
 ** Main Program
 **
 *****************************************************************************/
int main(int argc, char *argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  ASerializable::setPrefixName("Potential-");

  // Global parameters
  int ndim = 2;
  defineDefaultSpace(ESpaceType::RN, ndim);

  // Iso-Potential file
  VectorDouble tabiso = { 7., 6., 1.,
                          5., 6., 1.,
                          6., 5., 1.,
                          3., 6., 2.,
                          7., 7., 2.,
                          8., 3., 2.,
                          8., 1., 3.,
                          7., 9., 3.,
                          10., 5., 3.,
                          3., 1., 3. };
  Db* dbiso = Db::createFromSamples(10, ELoadBy::SAMPLE, tabiso,
                                    {"x","y","iso"},
                                    {"x1","x2","layer"});

  // Gradient file
  VectorDouble tabgrd = { 1., 6., 1., 0.,
                          9., 2., -1., 1.,
                          7., 8., 0., -1 };
  Db* dbgrd = Db::createFromSamples(3, ELoadBy::SAMPLE, tabgrd,
                                    {"x","y","gx","gy"},
                                    {"x1","x2","g1","g2"});

  // Tangent file
  VectorDouble tabtgt = { 3., 7., 1., 0.,
                          9., 7., 0.5, -0.5 };
  Db* dbtgt = Db::createFromSamples(2, ELoadBy::SAMPLE, tabtgt,
                                    {"x","y","tx","ty"},
                                    {"x1","x2","tangent1","tangent2"});

  // Generate the output grid
  VectorInt nx = {101,101};
  VectorDouble dx = {0.1, 0.1};
  DbGrid* grid = DbGrid::create(nx, dx);

  // Create the model
  Model* model = Model::createFromParam(ECov::CUBIC, 6.);
  model->switchToGradient();

  // Create the Neighborhood (unique)
  NeighUnique* neighU = NeighUnique::create();

  // Launch the Potential estimation
  (void) potential_kriging(dbiso, dbgrd, dbtgt, grid, model, neighU,
                           0., 0., true, false, false, false, 0, true);

  (void) grid->dumpToNF("Grid2D.ascii");

  // ====================== Free pointers ==================================
  delete dbiso;
  delete dbgrd;
  delete dbtgt;
  delete grid;
  delete model;
  delete neighU;

  //============================================================//
  // Exemple in 1-D
  //============================================================//

  mestitle(0,"Working in 1-D");
  ndim = 1;
  defineDefaultSpace(ESpaceType::RN, ndim);

  // Iso-Potential file
  tabiso = { 30., 1.,
             80., 1.,
             40., 2.,
             50., 2.};
  dbiso = Db::createFromSamples(4, ELoadBy::SAMPLE, tabiso,
                                { "x", "iso" },
                                { "x1", "layer" });

  // Gradient file
  tabgrd = { 0., 1.};
  dbgrd = Db::createFromSamples(1, ELoadBy::SAMPLE, tabgrd,
                                { "x", "gx" },
                                { "x1", "g1" });

  // Generate the output grid
  nx = { 101 };
  dx = { 1 };
  grid = DbGrid::create(nx, dx);

  // Create the model
  double range = scale2range(ECov::GAUSSIAN, 20.);
  model = Model::createFromParam(ECov::GAUSSIAN, range);
  model->switchToGradient();

  // Create the Neighborhood (unique)
  auto space = SpaceRN::create(ndim);
  neighU = NeighUnique::create(false, space);

  // Launch the Potential estimation
  // In case we would like to examine the calculation details,
  // set the rank of the target node in the next line
  OptDbg::setReference(-1);
  (void) potential_kriging(dbiso, dbgrd, nullptr, grid, model, neighU,
                           0., 0., true, true, false, true, 0, true);
  OptDbg::setReference(-1);

  // Visualize the results
  (void) grid->dumpToNF("Grid1D.ascii");

  // Cross-validation
  (void) potential_xvalid(dbiso, dbgrd, nullptr, model, neighU,
                          0., 0., true, true);

  delete dbiso;
  delete dbgrd;
  delete grid;
  delete model;
  delete neighU;

  return (0);
}
