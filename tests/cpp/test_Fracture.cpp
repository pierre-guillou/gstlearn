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

// This test is meant to demonstrate the fracture Simulation

#include "Enum/ESpaceType.hpp"

#include "Space/ASpaceObject.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Db/DbStringFormat.hpp"
#include "Basic/Law.hpp"
#include "Basic/File.hpp"
#include "Fractures/FracEnviron.hpp"
#include "Fractures/FracFamily.hpp"
#include "Fractures/FracFault.hpp"
#include "Fractures/FracList.hpp"

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

  ASerializable::setPrefixName("Fractures-");

  // Global parameters
  int ndim = 2;
  bool verbose = true;
  int ndisc = 1000.;
  law_set_random_seed(32131);

  defineDefaultSpace(ESpaceType::RN, ndim);
  DbStringFormat dbfmt(FLAG_STATS);

  // Generate the output grid
  VectorInt nx = {301,101};
  DbGrid* grid = DbGrid::create(nx);
  grid->display();

  // Creating the Fracture Environment
  double xmax = grid->getExtend(0);
  double ymax = grid->getExtend(1);
  double deltax = 0.;
  double deltay = 0.;
  double mean   = 20.;
  double stdev  = 10.;
  FracEnviron env = FracEnviron(xmax, ymax, deltax, deltay, mean, stdev);

  // Creating the Fault Families

  // Family1: orient=0, dorient=20, theta0=0.2, alpha=1, ratcst=1,
  //          prop1 = 0.5, prop2=0.2, aterm=1.2, bterm=2.4, range=5
  FracFamily family1 = FracFamily(0., 20., 0.2, 1., 1., 0.5, 0.2, 1.2, 2.4, 5.);

  // Family2: orient=30, dorient=5., theta0=0.2, alpha=1, ratcst=1,
  //          prop1 = 0.5, prop2=0.2, aterm=1.2, bterm=2.4, range=5
  FracFamily family2 = FracFamily(30., 5., 0.2, 1., 1., 0.5, 0.2, 1.2, 2.4, 5.);

  // Creating the Major Fault
  double coord   = 30.;
  double forient = 0.;
  FracFault fault = FracFault(coord, forient);
  double thetal  = 1.;
  double thetar  = 2.;
  double rangel  = 10.;
  double ranger  = 20.;
  fault.addFaultPerFamily(thetal, thetar, rangel, ranger);
  fault.addFaultPerFamily(thetal, thetar, rangel, ranger);

  // Gluing all elements within the Environment
  env.addFault(fault);
  env.addFamily(family1);
  env.addFamily(family2);
  env.display();

  // Simulating fractures
  FracList flist = FracList();
  int seed = 432431;
  flist.simulate(env, true, true, seed, verbose, VectorDouble());
  flist.display();

  // Plunge the set of fractures on the Grid
  VectorDouble permtab = { 20., 10., 15. };
  double perm_mat   = 0.;
  double perm_bench = 5.;
  (void) flist.fractureToBlock(grid, xmax, permtab, perm_mat, perm_bench,
                               ndisc, verbose);

  grid->display(&dbfmt);
  (void) grid->dumpToNF("Grid.ascii");

  // ====================== Free pointers ==================================
  delete grid;

  return (0);
}
