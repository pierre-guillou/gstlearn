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
#include "Calculators/CalcSimuPostPropByLayer.hpp"
#include "Db/DbGrid.hpp"
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

  auto* g2 = DbGrid::createFillRandom({10, 10}, 1);
  auto* g3 = DbGrid::createFillRandom({10, 10, 10}, 0);
  g2->display();

  VectorInt cells = {1, 2, 3};
  simuPostPropByLayer(g2, g3, {"z"}, false, false,
                      EPostUpscale::MEAN, EPostStat::fromKeys({"MINI", "MAXI", "MED", "MEAN", "STD"}), true, cells, 0);
  g3->display();

  return (0);
}
