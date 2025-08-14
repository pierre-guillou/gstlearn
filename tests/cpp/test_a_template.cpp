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
/**
 * This file is meant to parametrized the ModelGeneric in terms of ParamInfo
 * and to fit the values of these parameters according to the Maximum LogLikelihood
 * method and using the Vecchia approximation.
 */
#include "Basic/VectorHelper.hpp"
#include "Db/Db.hpp"
#include "Enum/ECov.hpp"
#include "Enum/ESpaceType.hpp"
#include "LinearOp/ProjMatrix.hpp"
#include "Mesh/MeshSphericalExt.hpp"
#include "Model/Model.hpp"
#include "Space/ASpaceObject.hpp"
#include "geoslib_define.h"
using namespace gstlrn;

int main(int argc, char* argv[])
{
  DECLARE_UNUSED(argc);
  DECLARE_UNUSED(argv);

  defineDefaultSpace(ESpaceType::SN, EARTH_RADIUS);

  auto mesh0 = MeshSphericalExt();
  (void)mesh0.resetFromDb(nullptr, nullptr, "-r3");

  int nsample = 4000;
  VH::sample(mesh0.getNApices(), nsample);
  VectorInt ind = VH::sampleRanks(mesh0.getNApices(), 0., nsample);

  // Creation of the db
  VectorDouble X = mesh0.getCoordinatesPerApex(0);
  VectorDouble Y = mesh0.getCoordinatesPerApex(1);

  Db dbdat         = Db();
  VectorDouble dbX = VH::sample(X, ind);
  VectorDouble dbY = VH::sample(Y, ind);
  dbdat.addColumns(dbX, "x1", ELoc::X, 0);
  dbdat.addColumns(dbY, "x2", ELoc::X, 1);

  double ratioRange   = 0.2;
  double scale        = EARTH_RADIUS * ratioRange;
  double sill         = 2.;
  VectorDouble coeffs = {1, -1, .5};
  auto* model         = Model::createFromParam(ECov::MARKOV, scale, sill);
  model->setMarkovCoeffs(0, coeffs);

  auto mesh = MeshSphericalExt();
  (void)mesh.resetFromDb(nullptr, nullptr, "-r2");

  auto Aproj = ProjMatrix(&dbdat, &mesh);

  Id ndef = VH::countUndefined(Aproj.getValues());
  message("Number of undef = %ld\n", ndef);

  exit(0);
}
