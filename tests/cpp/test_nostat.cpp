#include "Basic/AException.hpp"
#include "Basic/Vector.hpp"
#include "Basic/Law.hpp"
#include "Covariances/CovAniso.hpp"
#include "Db/Db.hpp"
#include "geoslib_e.h"

#include "geoslib_f.h"

#include "Model/Model.hpp"
#include "Model/NoStatArray.hpp"
#include "Model/NoStatFunctional.hpp"
#include "Basic/FunctionalSpirale.hpp"
#include "LinearOp/PrecisionOp.hpp"
#include "LinearOp/ShiftOpCs.hpp"
#include "Mesh/MeshETurbo.hpp"

#include "Matrix/MatrixRectangular.hpp"

#include <algorithm>
#include <math.h>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#define __USE_MATH_DEFINES
#include <cmath>

/****************************************************************************/
/*!
 ** Main Program
 **
 *****************************************************************************/
int main(int argc, char *argv[])

{
  int seed = 10355;
  law_set_random_seed(seed);

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("Nostat-");

  // Creating the 2-D Db
  auto nx = { 101, 101 };
  Db workingDbc(nx);

  // Creating the Non-stationary Model
  Model model = Model(&workingDbc);
  CovAniso cova = CovAniso(ECov::BESSEL_K,model.getContext());
  cova.setRanges({45,10});
  model.addCova(&cova);

  FunctionalSpirale spirale(0., -1.4, 1., 1., 50., 50.);
  VectorDouble angle = spirale.getFunctionValues(&workingDbc);
  MatrixSquareGeneral hh(2);

  int nech = workingDbc.getSampleNumber();
  VectorDouble h11(nech);
  VectorDouble h12(nech);
  VectorDouble h22(nech);
  for (int i = 0; i < nech; i++)
  {
    VectorDouble coor = workingDbc.getSampleCoordinates(i);
    VectorVectorDouble dirs = spirale.getFunctionVectors(coor);

    MatrixSquareGeneral rotmat(2);
    rotmat.setValue(0,0,dirs[0][0]);
    rotmat.setValue(1,0,dirs[1][0]);
    rotmat.setValue(0,1,dirs[0][1]);
    rotmat.setValue(1,1,dirs[1][1]);

    VectorDouble diag = ut_vector_power(cova.getScales(), 2.);
    MatrixSquareSymmetric temp(2);
    temp.setDiagonal(diag);
    hh.normMatrix(temp, rotmat);

    h11[i] = hh.getValue(0,0);
    h12[i] = hh.getValue(0,1);
    h22[i] = hh.getValue(1,1);
  }
  workingDbc.addFields(h11,"H1-1",ELoc::NOSTAT,0);
  workingDbc.addFields(h12,"H1-2",ELoc::NOSTAT,1);
  workingDbc.addFields(h22,"H2-2",ELoc::NOSTAT,2);
  workingDbc.display(1);

  // Inquiry the value of the Non-stationary parameters at a given sample
  int target = 1000;
  VectorDouble vect = workingDbc.getSampleAttributes(ELoc::NOSTAT,target);
  ut_vector_display("Non-stationary parameters at sample", vect);

  NoStatArray NoStat({"H1-1","H1-2","H2-2"},&workingDbc);
  model.addNoStat(&NoStat);
  model.display(1);

  message("Test performed successfully\n");

  MeshETurbo mesh(workingDbc);
  ShiftOpCs S(&mesh, &model, &workingDbc);
  PrecisionOp Qsimu(&S, &cova, EPowerPT::MINUSHALF, false);

  int nvertex = Qsimu.getSize();
  VectorDouble vectnew = ut_vector_simulate_gaussian(nvertex);

  VectorDouble result(Qsimu.getSize());
  Qsimu.eval(vectnew,result);
  workingDbc.addFields(result,"Simu",ELoc::Z);

  workingDbc.serialize("spirale.ascii");

  return 0;
}