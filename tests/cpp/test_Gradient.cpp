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
#include "API/SPDE.hpp"
#include "Basic/File.hpp"
#include "Basic/OptCst.hpp"
#include "Basic/OptCustom.hpp"
#include "Basic/OptDbg.hpp"
#include "Covariances/CovContext.hpp"
#include "Covariances/CovGradientGeneric.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Drifts/DriftFactory.hpp"
#include "Enum/ECst.hpp"
#include "Enum/ESpaceType.hpp"
#include "Estimation/CalcKriging.hpp"
#include "LinearOp/PrecisionOpMatrix.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Model/Model.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Space/ASpaceObject.hpp"
#include "geoslib_old_f.h"

using namespace gstlrn;

/****************************************************************************/
/*!
** This test is meant to check the elaboration of the CovGradientGenric class
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
  Id ndim = 2;
  defineDefaultSpace(ESpaceType::RN, ndim);
  ASerializable::setPrefixName("test_Gradient-");
  int mode = 2;

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

  // Options
  OptDbg::setReference(3);
  double ballradius = 0.01;

  // Using the old-style Depth and Gradient fake LMC Model
  if (mode <= 0 || mode == 1)
  {
    // Update the Model for Gradients
    (void)db_gradient_update(data);
    ModelGeneric* new_model = model_duplicate_for_gradient(model, ballradius);

    // Perform Kriging
    (void)kriging(data, grid, new_model, neigh, 1, 1, 0);

    delete new_model;
  }

  // Use the new CovGradientGeneric
  if (mode <= 0 || mode == 2)
  {
    (void)db_gradient_update(data); // Transform ELoc::G into ELoc::Z
    CovContext ctxt(data->getNLoc(ELoc::Z), static_cast<Id>(model->getNDim()));
    auto* new_model = new ModelGeneric(ctxt);
    auto covg       = CovGradientGeneric(*model->getCovAniso(0), ballradius);
    new_model->setCov(&covg);

    DriftList* drifts = DriftFactory::createDriftListForGradients(model->getDriftList(), ctxt);
    new_model->setDriftList(drifts);

    new_model->setContext(ctxt);

    // Perform Kriging
    OptCustom::define("NotOptimSimpleCase", 1.);
    (void)kriging(data, grid, new_model, neigh, 1, 1, 0);

    delete drifts;
    delete new_model;
  }

  // Free memory
  delete data;
  delete grid;
  delete model;
  delete neigh;

  return (0);
}
