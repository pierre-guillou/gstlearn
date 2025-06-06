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
#include "Basic/OptCustom.hpp"
#include "Enum/ECov.hpp"
#include "Enum/ESpaceType.hpp"

#include "Space/ASpaceObject.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Model/Model.hpp"
#include "Model/Constraints.hpp"
#include "Basic/Law.hpp"
#include "Basic/File.hpp"
#include "Basic/OptDbg.hpp"
#include "Neigh/NeighMoving.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Anamorphosis/CalcAnamTransform.hpp"
#include "Variogram/VarioParam.hpp"
#include "Variogram/Vario.hpp"
#include "Estimation/CalcKriging.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"
#include "Calculators/CalcStatistics.hpp"
#include "Estimation/CalcKrigingFactors.hpp"

#include <math.h>

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

  // Global parameters
  int ndim = 2;
  law_set_random_seed(32131);
  defineDefaultSpace(ESpaceType::RN, ndim);
  DbStringFormat dbfmt(FLAG_STATS);

  // General parameters
  bool verbose    = true;

  // Generate initial grid
  DbGrid* grid = DbGrid::create({100,100}, {0.01,0.01});

  // Create grid of (single) Panel
  double dx_P = 0.250;
  double x0_P = 0.375;
  DbGrid* panel = DbGrid::create({1,1}, {dx_P, dx_P}, {x0_P, x0_P});
  panel->display();

  // Discretization with blocs
  int nx_B = 5;
  double x0_B = x0_P + dx_P / nx_B / 2.;
  double dx_B = dx_P / nx_B;
  DbGrid* blocs = DbGrid::create({nx_B,nx_B}, {dx_B,dx_B}, {x0_B,x0_B});
  blocs->display();

  // Simulation of the Data
  Model* model_init = Model::createFromParam(ECov::EXPONENTIAL, 0.1, 1.);
  (void) simtub(nullptr, grid, model_init);
  grid->setName("Simu", "Y");
  grid->display();

  // Non-linear transform
  double m_Z = 1.5;
  double s_Z = 0.5;
  VectorDouble Zval = grid->getColumn("Y");
  for (int i = 0; i < (int) Zval.size(); i++)
    Zval[i] = m_Z * exp(s_Z * Zval[i] - 0.5 * s_Z * s_Z);
  grid->addColumns(Zval, "Z");
  grid->display(&dbfmt);

  // Data extraction
  int np = 500;
  Db* data = Db::createSamplingDb(grid, 0., np, {"x1","x2","Y","Z"});
  data->setLocator("Z", ELoc::Z, 0);
  data->display(&dbfmt);

  // Gaussian Anamorphosis with 10 coefficients
  AnamHermite* anam = AnamHermite::create(20);
  anam->fitFromLocator(data);
  anam->display();

  // Selectivity
  Selectivity* selectivity =
      Selectivity::createByCodes( { ESelectivity::Q, ESelectivity::T },
                                  { 0., 0.5 }, true, true);

  // Global experimental selectivity
  data->display();
  selectivity->eval(data).display();

  // Selectivity in the model
  selectivity->evalFromAnamorphosis(anam).display();

  // Define the variogram calculation parameters
  VarioParam* varioparam = VarioParam::createOmniDirection(10, 0.025);

  // Calculate the variogram of the raw variable
  Vario* vario_raw = Vario::computeFromDb(*varioparam, data);
  vario_raw->display();

  // Fitting the Model on the Raw variable
  Model* model_raw = new Model(1, ndim);
  (void) model_raw->fit(vario_raw);
  model_raw->display();

  // Transform Data into Gaussian variable
  (void) anam->rawToGaussianByLocator(data);
  data->setName("Y.Z","Gauss.Z");
  data->display();

  // Calculate the variogram of the Gaussian variable
  Vario* vario = Vario::computeFromDb(*varioparam, data);
  vario->display();

  // Fitting the Model on the Gaussian transformed variable
  Model* model = new Model(1, ndim);
  Constraints constraints = Constraints(1.);
  (void) model->fit(vario, {ECov::EXPONENTIAL, ECov::EXPONENTIAL}, constraints);
  model->display();

  // Creating a Moving Neighborhood
  int nmini = 5;
  int nmaxi = 5;
  double radius = 1.;
  NeighMoving* neighM = NeighMoving::create(false, nmaxi, radius, nmini);
  neighM->display();

  // ====================== Conditional Expectation =====================

  // Estimating the Gaussian Variable on the nodes of the Blocks
  data->display();
  (void) kriging(data, blocs, model, neighM, 
                 true, true, false, KrigOpt(),
                 NamingConvention("G_PTS"));

  // Calculating the Conditional Expectation
  (void) ConditionalExpectation(blocs, anam, selectivity, "G_PTS*estim",
                                "G_PTS*stdev", false, TEST, 0,
                                NamingConvention("PTS_Recovery",false));
  blocs->display();

  // ====================== Point Disjunctive Kriging =====================

  // Setting the trace
  if (verbose) OptDbg::setReference(1);

  // Attach the Anamorphosis
  model->setAnam(anam);

  // Computing the Point factors
  int nfactor = 3;
  (void) anam->rawToFactor(data, nfactor);
  data->display();

  // Simple Point Kriging over the blocks
  (void) krigingFactors(data, blocs, model, neighM,
                        true, true, KrigOpt(), NamingConvention("DK_Pts"));
  blocs->display();

  // Simple Block Kriging over the blocks
  VectorInt ndisc_B = { 5, 5 };
  KrigOpt krigopt;
  krigopt.setOptionCalcul(EKrigOpt::BLOCK, ndisc_B);
  (void)krigingFactors(data, blocs, model, neighM,
                       true, true, krigopt, NamingConvention("DK_Blk"));
  blocs->display();

  // Simple Block Kriging over the panel(s)
  VectorInt ndisc_P = { 10, 10 };
  krigopt.setOptionCalcul(EKrigOpt::BLOCK, ndisc_P);
  (void) krigingFactors(data, panel, model, neighM, 
                        true, true, krigopt, NamingConvention("DK_Blk"));
  panel->display();

  // ====================== Uniform Conditioning ==================================

  OptCustom::define("ompthreads",1);
  // Perform the Point Kriging of the Raw Variable
  data->clearLocators(ELoc::Z);
  data->setLocator("Z",ELoc::Z, 0);
  data->display();
  (void) kriging(data, blocs, model, neighM, 
                 true, true, false, KrigOpt(), 
                 NamingConvention("Z_PTS"));
  blocs->display();

  // Perform the Uniform Conditioning over Blocks
  (void) UniformConditioning(blocs, anam, selectivity,
                             "Z_PTS*estim", "Z_PTS*stdev",
                             NamingConvention("UC",false));
  blocs->display();
  data->setLocator("Gauss.Z",ELoc::Z, 0);

  // ====================== Block Disjunctive Kriging (DGM-1) =====================

  // Calculate the change of support coefficient
  double r1 = anam->evalSupportCoefficient(1, model, blocs->getDXs(), ndisc_B);

  // Update the Model with Block anamorphosis
  AnamHermite* anam_b1 = anam->clone();
  anam_b1->setRCoef(r1);

  // Regularization of the point model by the block support
  Vario* vario_b1_Z = Vario::createRegularizeFromModel(*model, *varioparam,
                                                       blocs->getDXs(),
                                                       ndisc_B, blocs->getAngles());
  Vario* vario_b1_Y = Vario::createTransformZToY(*vario_b1_Z, anam);

  // Fitting the regularized model on the point Gaussian variable
  Model* model_b1_Y = new Model(1, ndim);
  constraints.setConstantSillValue(1);
  (void) model_b1_Y->fit(vario_b1_Y, { ECov::CUBIC, ECov::EXPONENTIAL }, constraints);

  // Update the Model with Block Anamorphosis
  model_b1_Y->setAnam(anam_b1);
  model_b1_Y->display();

  // Simple Point Kriging over the blocs(s) with Model with Change of Support
  (void) krigingFactors(data, blocs, model_b1_Y, neighM, 
                        true, true, KrigOpt(), NamingConvention("DK_DGM1"));
  blocs->display();

  // Simple Point Kriging over the panel(s) with Model with Change of Support
  VectorInt ndisc_Bc = { nx_B, nx_B };
  krigopt.setOptionCalcul(EKrigOpt::BLOCK, ndisc_Bc);
  (void) krigingFactors(data, panel, model_b1_Y, neighM, 
                        true, true, krigopt, NamingConvention("DK_DGM1"));
  panel->display();

  // ====================== Block Disjunctive Kriging (DGM-2) =====================

  // Calculate the change of support coefficient
  double r2 = anam->evalSupportCoefficient(2, model, blocs->getDXs(), ndisc_B);

  // Update the Model with Block anamorphosis
  AnamHermite* anam_b2 = anam->clone();
  anam_b2->setRCoef(r2);

  // Regularization of the point model by the block support
  Vario* vario_b2_Y = Vario::createRegularizeFromModel(*model, *varioparam,
                                                       blocs->getDXs(),
                                                       ndisc_B, blocs->getAngles());

  // Fitting the regularized model on the point Gaussian variable
  Model* model_b2_Y = new Model(1, ndim);
  constraints.setConstantSillValue(r2 * r2);
  (void) model_b2_Y->fit(vario_b2_Y, { ECov::CUBIC, ECov::EXPONENTIAL }, constraints);

  // Normalization of the block model to a total sill equal to 1.0
  model_b2_Y->normalize(1.0);
  model_b2_Y->display();

  // Update the Model with Block anamorphosis
  model_b2_Y->setAnam(anam_b2);
  model_b2_Y->display();

  // Simple Point Kriging over the blocs(s) with Model with Change of Support
  (void) krigingFactors(data, blocs, model_b2_Y, neighM, 
                        true, true, KrigOpt(), NamingConvention("DK_DGM2"));
  blocs->display();

  // Simple Point Kriging over the panel(s) with Model with Change of Support
  krigopt.setOptionCalcul(EKrigOpt::BLOCK, ndisc_Bc);
  (void) krigingFactors(data, panel, model_b2_Y, neighM, 
                        true, true, krigopt, NamingConvention("DK_DGM2"));
  panel->display();

  // ====================== Selectivity Function ==================================

  DisjunctiveKriging(blocs, anam, selectivity,
                      blocs->getName("DK_Pts*estim"),
                      blocs->getName("DK_Pts*stdev"),
                      NamingConvention("QT",false));
  blocs->display();

  // ====================== Free pointers ==================================

  delete data;
  delete grid;
  delete model_raw;
  delete model_init;
  delete model;
  delete model_b1_Y;
  delete model_b2_Y;
  delete panel;
  delete blocs;
  delete anam;
  delete anam_b1;
  delete anam_b2;
  delete varioparam;
  delete vario;
  delete vario_raw;
  delete vario_b1_Z;
  delete vario_b1_Y;
  delete vario_b2_Y;
  delete neighM;
  delete selectivity;

  return (0);
}
