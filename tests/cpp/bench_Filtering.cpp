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

// This test is mean to check the Factorial Kriging Analysis

#include "Basic/AStringFormat.hpp"
#include "Basic/NamingConvention.hpp"
#include "Covariances/CovContext.hpp"
#include "Enum/ESpaceType.hpp"

#include "Space/ASpaceObject.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Basic/Law.hpp"
#include "Model/Model.hpp"
#include "Basic/File.hpp"
#include "Basic/OptDbg.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Neigh/NeighMoving.hpp"
#include "Estimation/CalcKriging.hpp"

int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  // Global parameters
  int ndim = 2;
  law_set_random_seed(32131);
  AStringFormat format;
  defineDefaultSpace(ESpaceType::RN, ndim);

  // Parameters
  bool verbose    = true;
  int nech        = 3;
  int nvar        = 2; 
  bool flagSK     = true;

  // Generate the data base
  Db* data = Db::createFillRandom(nech, ndim, nvar, 0);
  data->setLocVariable(ELoc::Z, 1, 0, TEST);
  DbStringFormat* dbfmt = DbStringFormat::create(FLAG_ARRAY);
  data->display(dbfmt);

  // Generate the target file
  Db* target = Db::createFillRandom(1, ndim, 0, 0);
  target->setCoordinate(0, 0, data->getCoordinate(0, 0));
  target->setCoordinate(0, 1, data->getCoordinate(0, 1));

  // Create the Model
  int order = (flagSK) ? -1 : 0;
  Model* model = Model::createFillRandom(ndim, nvar, {ECov::NUGGET, ECov::SPHERICAL},
                                         1., order);
  model->setCovFiltered(0, true);
  model->display();

  // Neighborhood
  ANeigh* neigh = NeighUnique::create();

  // Define the verbose option
  if (verbose) OptDbg::setReference(1);

  // Test on Collocated CoKriging in Unique Neighborhood
  kriging(data, target, model, neigh);
  dbfmt = DbStringFormat::create(FLAG_STATS, {"Kriging.*"});
  target->display(dbfmt);

  // Free pointers

  delete neigh;
  delete data;
  delete target;
  delete model;

  return (0);
}
