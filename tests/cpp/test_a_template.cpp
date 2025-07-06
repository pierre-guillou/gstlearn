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
#include "Basic/ASerializable.hpp"
#include "Basic/OptCst.hpp"
#include "Basic/VectorHelper.hpp"
#include "Db/DbGrid.hpp"
#include "Db/DbStringFormat.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Model/Model.hpp"

/**
 * This file is meant check Takahashi algorithm
 * based on a small set of samples which coincide with some apices
 * of a known Meshing
 */
int main(int argc, char* argv[])
{
  DECLARE_UNUSED(argc);
  DECLARE_UNUSED(argv);
  ASerializable::setPrefixName("test_Template-");
  DbStringFormat* dbfmt = DbStringFormat::createFromFlags(false, true, false,
                                                          false, true);
  OptCst::define(ECst::NTCOL, -1);
  OptCst::define(ECst::NTROW, -1);

  // Create a Grid covering the 1 x 1 square domain
  int ndim         = 2;
  int nx           = 10;
  VectorInt nxs    = {nx, nx};
  VectorDouble dxs = {1. / nx, 1. / nx};
  DbGrid* grid     = DbGrid::create(nxs, dxs);

  // Determine the Mesh adapted to the grid
  MeshETurbo* mesh = MeshETurbo::create(nxs, dxs);
  int npoints      = mesh->getNApices();

  // Sample the Meshing to constitute a Db
  double proportion = 0.1;
  VectorInt ranks   = VH::sampleRanks(npoints, proportion);
  message("Npoints=%d\n", npoints);
  int number = (int)ranks.size();
  VectorDouble coords(2. * number);
  for (int i = 0, ecr = 0; i < number; i++)
    for (int idim = 0; idim < ndim; idim++, ecr++)
      coords[ecr] = mesh->getApexCoor(ranks[i], idim);
  Db* data = Db::createFromSamples(number, ELoadBy::SAMPLE, coords, {"x1", "x2"},
                                   VectorString(), false);
  data->setLocators({"x*"}, ELoc::X);
  data->addColumnsRandom(1, "z", ELoc::Z);
  data->display(dbfmt);

  // Create a Model
  double range = 1. / 3.;
  double sill  = 9.;
  Model* model = Model::createFromParam(ECov::MATERN, range, sill);

  // Perform Kriging using SPDE
  (void)krigingSPDE(data, grid, model, true, true, 1, {mesh});
  grid->display(dbfmt);

  data->dumpToNF("data.NF");
  grid->dumpToNF("grid.NF");
  mesh->dumpToNF("mesh.NF");

  delete grid;
  delete mesh;
  delete data;
  delete model;
  delete dbfmt;
  return (0);
}
