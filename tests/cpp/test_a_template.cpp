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
#include "Basic/ASerializable.hpp"
#include "Calculators/CalcGridToGrid.hpp"
#include "Db/DbGrid.hpp"
#include "Model/Model.hpp"
#include "Simulation/Simulations.hpp"
#include "geoslib_define.h"

using namespace gstlrn;

int main(int argc, char* argv[])
{
  // Unless you test a specific feature, bring it to its minimal expansion
  // e.g. the few lines below
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);
  ASerializable::setPrefixName("test_a_template-"); // Here set the test name

  auto* grid = DbGrid::create({150, 100});

  auto* model = Model::createFromParam(ECov::CUBIC, 30, 10.0);
  (void)simtub(nullptr, grid, model, nullptr, 2);
  grid->setName("Simu.S1", "Top");
  grid->setName("Simu.S2", "Bot");
  for (Id i = 0; i < grid->getNSample(); i++)
  {
    grid->setValue("Top", i, grid->getValue("Top", i) + 110.);
    grid->setValue("Bot", i, grid->getValue("Bot", i) + 100.);
  }

  model = Model::createFromParam(ECov::CUBIC, 30, 10.);
  (void)simtub(nullptr, grid, model, nullptr, 2);
  grid->setName("Simu", "VBot");

  model = Model::createFromParam(ECov::SPHERICAL, 10, 3.);
  (void)simtub(nullptr, grid, model, nullptr, 2);
  grid->setName("Simu", "VTop");

  auto* g3D = DbGrid::create({150, 100, 3}, {1, 1, 1}, {0, 0, 91});
  grid->setLocators({"VBot", "VTop"}, ELoc::Z);
  (void)dbg2gInterpolate(grid, g3D, {"Top"}, {"Bot"});

  return (0);
}
