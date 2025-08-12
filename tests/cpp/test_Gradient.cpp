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
#include "Estimation/CalcKriging.hpp"
#include "Neigh/NeighUnique.hpp"
#include "geoslib_old_f.h"

#include "Enum/ECst.hpp"
#include "Enum/ESpaceType.hpp"

#include "API/SPDE.hpp"

#include "Basic/File.hpp"
#include "Basic/OptCst.hpp"
#include "Basic/OptDbg.hpp"
#include "Covariances/CovContext.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "LinearOp/PrecisionOpMatrix.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Model/Model.hpp"
#include "Space/ASpaceObject.hpp"

using namespace gstlrn;

/****************************************************************************/
/*!
** This test is meant to check the elaboration of the CovGradient class
** and its use in Kriging (or Depth adn Gradient)
**
*****************************************************************************/
int main(int argc, char* argv[])

{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  /***********************/
  /* 1 - Initializations */
  /***********************/
  int ndim = 2;
  defineDefaultSpace(ESpaceType::RN, ndim);
  ASerializable::setPrefixName("test_Gradient-");

  // Setup constants

  OptDbg::reset();
  OptCst::define(ECst::NTCAR, 10.);
  OptCst::define(ECst::NTDEC, 6.);

  // Define the data set
  Db* data = Db::create();
  data->addColumns({0., 1., 1., 2., 3.}, "x1", ELoc::X, 0);
  data->addColumns({0., 3., 0., 0., 0.}, "x2", ELoc::X, 1);
  data->addColumns({0., 0., 1., 2., 3.}, "z1", ELoc::Z, 0);
  data->addColumns({1., -1., 0., 1., 2.}, "g1", ELoc::G, 0);
  data->addColumns({1., -1., 0., 2., 1.}, "g2", ELoc::G, 1);
  data->display();

  // Define the output Grid
  DbGrid* grid = DbGrid::create({4, 4});
  grid->display();

  // Define the (monovariate) Model
  Model* model = Model::createFromParam(ECov::GAUSSIAN, 2.);
  model->setDriftIRF(1);
  model->display();

  // Define the (unique) Neighborhood
  NeighUnique* neigh = NeighUnique::create();
  neigh->display();

  // Update the Model for Gradients
  double ball_radius = 0.01;
  (void)db_gradient_update(data);
  Model* new_model = model_duplicate_for_gradient(model, ball_radius);
  new_model->display();

  // Perform Kriging
  OptDbg::setReference(3);
  if (kriging(data, grid, new_model, neigh,
              1, 1, 0)) messageAbort("kriging");
  OptDbg::setReference(0);

  // Free memory
  delete data;
  delete grid;
  delete model;
  delete new_model;
  delete neigh;

  return (0);
}
