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

  const auto* g3 = DbGrid::create({2, 3, 4}, VectorDouble(), VectorDouble(), {10., 20., 30.});
  g3->display();
  const auto* g2 = DbGrid::create({2, 3}, VectorDouble(), VectorDouble(), {10., 20.});
  g2->display();
  return (0);
}
