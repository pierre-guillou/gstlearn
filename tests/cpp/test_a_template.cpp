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
#include "Db/DbGrid.hpp"
#include "LinearOp/ProjMatrix.hpp"
#include "LinearOp/ProjMultiMatrix.hpp"
#include "LinearOp/SPDEOpMatrix.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Model/Model.hpp"
/**
 * This file is meant to parametrized the ModelGeneric in terms of ParamInfo
 * and to fit the values of these parameters according to the Maximum LogLikelihood
 * method and using the Vecchia approximation.
 */
using namespace gstlrn;

int main(int argc, char* argv[])
{
  DECLARE_UNUSED(argc);
  DECLARE_UNUSED(argv);

  double rangev = 0.2;
  double sill   = 1.;

  int nx           = 200;
  double dx        = 1. / (nx - 1);
  VectorInt nxs    = {nx, nx};
  VectorDouble dxs = {dx, dx};
  VectorDouble x0s = {0., 0.};

  int nxm           = 100;
  double dxm        = 1.5 / (nxm - 1);
  VectorInt nxms    = {nxm, nxm};
  VectorDouble dxms = {dxm, dxm};
  VectorDouble x0ms = {-0.25, -0.25};

  auto* grid          = DbGrid::create(nxs, dxs, x0s);
  auto* gridExt       = DbGrid::create(nxms, dxms, x0ms);
  auto* mesh          = MeshETurbo::createFromGrid(gridExt);
  VectorMeshes meshes = {mesh};

  auto* model = Model::createFromParam(ECov::MATERN, rangev, sill);

  auto proj      = ProjMatrix(grid, mesh);
  auto projmulti = ProjMultiMatrix({{&proj}});

  bool verbose  = true;
  bool flagSimu = true;
  auto spde     = SPDE(nullptr, grid, model, flagSimu, 1);
  (void)spde.defineMeshes(meshes, VectorMeshes(), verbose);
  (void)spde.defineProjections(&projmulti, nullptr, verbose);
  (void)spde.defineShiftOperator();

  delete model;

  return (0);
}
