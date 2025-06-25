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
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Neigh/NeighMoving.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Variogram/VarioParam.hpp"
#include "Variogram/Vario.hpp"
#include "Matrix/Table.hpp"
#include "Model/Model.hpp"
#include "Basic/File.hpp"
#include "Basic/Law.hpp"
#include "Basic/PolyLine2D.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/AStringFormat.hpp"
#include "Polygon/Polygons.hpp"
#include "LithoRule/Rule.hpp"

/****************************************************************************/
/*!
** Main Program for testing serialize/deserialize
** in Neutral or HDF5 format
**
*****************************************************************************/
int main(int argc, char *argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  ASerializable::setPrefixName("TS-");

  // Next flag indicates if the format is NF (true) or H5 (false)
  bool flagNeutral = false;
  bool verbose     = false;
  int mode         = 10;

  // =================
  // Preliminary tasks
  // =================

  int nech          = 20;
  Db* db            = Db::createFromBox(nech, {0., 0.}, {1., 1.}, 32432);
  VectorDouble vec1 = VH::simulateGaussian(nech);
  db->addColumns(vec1, "myvar1", ELoc::Z, 0);
  VectorDouble vec2 = VH::simulateGaussian(nech);
  db->addColumns(vec2, "myvar2");

  DbGrid* dbg = DbGrid::create({12, 10}, {0.1, 0.3}, {0.2, 0.4});
  vec1        = VH::simulateGaussian(dbg->getNSample());
  dbg->addColumns(vec1, "myvar1", ELoc::Z, 0);
  vec2    = VH::simulateGaussian(dbg->getNSample());
  vec2[2] = TEST;
  vec2[5] = TEST;
  dbg->addColumns(vec2, "myvar2", ELoc::Z, 1);

  // ===========
  // Checking Db
  // ===========

  if (mode == 0 || mode == 1)
  {
    Db* db1      = db->clone();
    db1->display();

    // Serialize
    if (flagNeutral)
      (void)db1->dumpToNF("Neutral.Db.ascii");
    else
      (void)db1->dumpToH5("Neutral.Db.h5");

    // Deserialize
    Db* db2 = nullptr;
    if (flagNeutral)
      db2 = Db::createFromNF("Neutral.Db.ascii", verbose);
    else
      db2 = Db::createFromH5("Neutral.Db.h5", verbose);
    db2->display();

    delete db1;
    delete db2;
  }

  // ===============
  // Checking DbGrid
  // ===============

  if (mode == 0 || mode == 2)
  {
    DbGrid* dbg1 = dbg->clone();
    dbg1->display();

    // Serialize
    if (flagNeutral)
      (void)dbg1->dumpToNF("Neutral.Dbg.ascii");
    else
      (void)dbg1->dumpToH5("Neutral.Dbg.h5");

    // Deserialize
    Db* dbg2 = nullptr;
    if (flagNeutral)
      dbg2 = DbGrid::createFromNF("Neutral.Dbg.ascii", verbose);
    else
      dbg2 = DbGrid::createFromH5("Neutral.Dbg.h5", verbose);
    dbg2->display();

    delete dbg1;
    delete dbg2;
  }

  // =================
  // Checking Polygons
  // =================

  if (mode == 0 || mode == 3)
  {
    Polygons* poly1 = new Polygons();
    poly1->resetFromDb(db);
    Polygons* polyb = new Polygons();
    polyb->resetFromDb(dbg);
    poly1->addPolyElem(polyb->getPolyElem(0));
    poly1->display();

    // Serialize
    if (flagNeutral)
      (void)poly1->dumpToNF("Neutral.Polygon.ascii");
    else
      (void)poly1->dumpToH5("Neutral.Polygon.h5");

    // Deserialize
    Polygons* poly2 = nullptr;
    if (flagNeutral)
      poly2 = Polygons::createFromNF("Neutral.Polygon.ascii", verbose);
    else
      poly2 = Polygons::createFromH5("Neutral.Polygon.h5", verbose);
    poly2->display();

    delete polyb;
    delete poly1;
    delete poly2;
  }

  // ==============
  // Checking Vario
  // ==============

  if (mode == 0 || mode == 4)
  {
    VarioParam varioparam;
    DirParam dirparam(10, 0.02);
    varioparam.addDir(dirparam);
    Vario vario1 = Vario(varioparam);
    vario1.compute(db, ECalcVario::VARIOGRAM);
    vario1.display();

    // Serialize
    if (flagNeutral)
      (void)vario1.dumpToNF("Neutral.Vario.ascii");
    else
      (void)vario1.dumpToH5("Neutral.Vario.h5");

    // Deserialize
    Vario* vario2 = nullptr;
    if (flagNeutral)
      vario2 = Vario::createFromNF("Neutral.Vario.ascii", verbose);
    else
      vario2 = Vario::createFromH5("Neutral.Vario.h5", verbose);
    vario2->display();

    delete vario2;
  }

  // ==============
  // Checking Model
  // ==============

  if (mode == 0 || mode == 5)
  {
    Model* model1 = Model::createFromParam(ECov::EXPONENTIAL, 0.3, 0.2, 1.);
    model1->display();

    // Serialize model1
    if (flagNeutral)
      (void)model1->dumpToNF("Neutral.Model.ascii");
    else
      (void)model1->dumpToH5("Neutral.Model.h5");

    // Deserialize model2
    Model* model2 = nullptr;
    if (flagNeutral)
      model2 = Model::createFromNF("Neutral.Model.ascii", verbose);
    else
      model2 = Model::createFromH5("Neutral.Model.h5", verbose);
    model2->display();

    delete model1;
    delete model2;
  }

  // =======================
  // Checking Table
  // =======================

  if (mode == 0 || mode == 6)
  {
    VectorVectorDouble table;
    int ncols     = 3;
    int nrows     = 10;
    Table* table1 = Table::create(nrows, ncols);
    for (int irow = 0; irow < nrows; irow++)
      for (int icol = 0; icol < ncols; icol++)
        table1->setValue(irow, icol, law_uniform());
    table1->display();

    // Serialize table
    if (flagNeutral)
      (void)table1->dumpToNF("Neutral.Table.ascii");
    else
      (void)table1->dumpToH5("Neutral.Table.h5");

    // Deserialize table1
    Table* table2 = nullptr;
    if (flagNeutral)
      table2 = Table::createFromNF("Neutral.Table.ascii", verbose);
    else
      table2 = Table::createFromH5("Neutral.Table.h5", verbose);
    table2->display();

    delete table1;
    delete table2;
  }

  // =============
  // Checking Rule
  // =============

  if (mode == 0 || mode == 7)
  {
    Rule* rule1 = Rule::createFromNames({"S", "F1", "T", "F2", "S", "F3", "F4"});
    rule1->display();

    // Serialize
    if (flagNeutral)
      (void)rule1->dumpToNF("Neutral.Rule.ascii");
    else
      (void)rule1->dumpToH5("Neutral.Rule.h5");

    // Deserialize
    Rule* rule2 = nullptr;
    if (flagNeutral)
      rule2 = Rule::createFromNF("Neutral.Rule.ascii", verbose);
    else
      rule2 = Rule::createFromH5("Neutral.Rule.h5", verbose);
    rule2->display();

    delete rule1;
    delete rule2;
  }

  // ===================
  // Checking PolyLine2D
  // ===================

  if (mode == 0 || mode == 8)
  {
    int npolyline          = 100;
    VectorDouble xpolyline = VH::simulateGaussian(npolyline);
    VectorDouble ypolyline = VH::simulateGaussian(npolyline);
    PolyLine2D* polyline1  = new PolyLine2D(xpolyline, ypolyline);
    AStringFormat afmt(3);
    polyline1->display(&afmt);

    // Serialize
    if (flagNeutral)
      (void)polyline1->dumpToNF("Neutral.Polyline.ascii");
    else
      (void)polyline1->dumpToH5("Neutral.Polyline.h5");

    // Deserialize
    PolyLine2D* polyline2 = nullptr;
    if (flagNeutral)
      polyline2 = PolyLine2D::createFromNF("Neutral.Polyline.ascii", verbose);
    else
      polyline2 = PolyLine2D::createFromH5("Neutral.Polyline.h5", verbose);
    polyline2->display(&afmt);

    delete polyline1;
    delete polyline2;
  }

  // ============================
  // Checking Moving Neighborhood
  // ============================

  if (mode == 0 || mode == 9)
  {
    int nmaxi = 20;
    double radius = 4.;
    int nmini = 2;
    int nsect = 5;
    int nsmax           = 3;
    VectorDouble coeffs = {2., 3.};
    VectorDouble angles = {25., 0.};
    bool useBallTree    = true;

    NeighMoving* neigh1 = NeighMoving::create(false, nmaxi, radius,
                                              nmini, nsect, nsmax, 
                                              coeffs, angles, useBallTree);
    neigh1->display();

    // Serialize
    if (flagNeutral)
      (void)neigh1->dumpToNF("Neutral.NeighMoving.ascii");
    else
      (void)neigh1->dumpToH5("Neutral.NeighMoving.h5");

    // Deserialize
    NeighMoving* neigh2 = nullptr;
    if (flagNeutral)
      neigh2 = NeighMoving::createFromNF("Neutral.NeighMoving.ascii", verbose);
    else
      neigh2 = NeighMoving::createFromH5("Neutral.NeighMoving.h5", verbose);
    neigh2->display();

    delete neigh1;
    delete neigh2;
  }

  // ============================
  // Checking Unique Neighborhood
  // ============================

  if (mode == 0 || mode == 10)
  {
    NeighUnique* neigh1 = NeighUnique::create();
    neigh1->display();
    flagNeutral = true;

    // Serialize
    if (flagNeutral)
      (void)neigh1->dumpToNF("Neutral.NeighUnique.ascii");
    else
      (void)neigh1->dumpToH5("Neutral.NeighUnique.h5");

    // Deserialize
    NeighUnique* neigh2 = nullptr;
    if (flagNeutral)
      neigh2 = NeighUnique::createFromNF("Neutral.NeighUnique.ascii", verbose);
    else
      neigh2 = NeighUnique::createFromH5("Neutral.NeighUnique.h5", verbose);
    neigh2->display();

    delete neigh1;
    delete neigh2;
  }

  // Cleaning procedure
  delete db;
  delete dbg;

  return (0);
}
