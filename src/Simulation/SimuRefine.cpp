/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "geoslib_f.h"
#include "geoslib_old_f.h"
#include "Db/DbGrid.hpp"
#include "Db/Db.hpp"
#include "Model/Model.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Simulation/ASimulation.hpp"
#include "Simulation/SimuRefine.hpp"
#include "Simulation/SimuRefineParam.hpp"
#include "Basic/Law.hpp"

#include <math.h>

#define LHS(i,j) (lhs[(i) * neq + (j)])
#define RHS(i)   (rhs[(i)])

SimuRefine::SimuRefine(int nbsimu, int seed)
    : ASimulation(nbsimu, seed),
      _param(),
      _model(nullptr),
      _ndim(0),
      _nx1(3),
      _dx1(3),
      _x01(3),
      _nx2(3),
      _dx2(3),
      _x02(3)
{
}

SimuRefine::~SimuRefine()
{
}

/****************************************************************************/
/*!
 **  Refine the simulation
 **
 ** \return  Refined Grid
 **
 ** \param[in]  dbin       Input grid Db structure
 ** \param[in]  model      Model structure
 ** \param[in]  param      SimuRefineParam structure
 **
 *****************************************************************************/
DbGrid* SimuRefine::simulate(DbGrid *dbin, Model* model, const SimuRefineParam& param)
{
  DbGrid *db1, *db2;

  /* Initializations */

  db1 = db2 = nullptr;
  db1 = dbin;
  law_set_random_seed(getSeed());
  _ndim = dbin->getNDim();
  _param = param;
  _model = model;

  /* Store information from the input grid */

  int iatt1 = db_attribute_identify(dbin, ELoc::Z, 0);
  if (iatt1 <= 0) return nullptr;

  /* Loop on the refinement factors */

  for (int imult = 0; imult < param.getNmult(); imult++)
  {

    /* Create the output grid */

    _dim_1_to_2(db1);
    db2 = db_create_grid(0, _ndim, 1, ELoadBy::SAMPLE, 1, _nx2, _x02, _dx2,
                         dbin->getGrid().getRotAngles());
    int iatt2 = db2->addColumnsByConstant(1, TEST);

    /* Establish the Kriging system */

    if (_kriging_define())
    {
      if (db1 != dbin) delete db1;
      delete db2;
      return nullptr;
    }

    /* Copy the initial data */

    _merge_data(db1, iatt1, db2, iatt2);

    /* Perform the simulation */

    _simulate_nodes(db2, iatt2);

    /* Create the new input file (for next step) */

    if (db1 != dbin) db1 = db_delete(db1);
    _dim_2_to_1(db2);
    db1 = db_create_grid(0, _ndim, 1, ELoadBy::SAMPLE, 1, _nx1, _x01, _dx1,
                         dbin->getGrid().getRotAngles());
    iatt1 = db1->addColumnsByConstant(1, TEST);

    /* Truncate the output grid for next step */

    _truncate_result(db2, iatt2, db1, iatt1);

    /* Delete the output file */

    db2 = db_delete(db2);
  }

  return db1;
}
//
/****************************************************************************/
/*!
 **  Define the characteristics from Db1 to Db2
 **
 ** \param[in]  db  Staring grid Db structure
 **
 *****************************************************************************/
void SimuRefine::_dim_1_to_2(DbGrid *db)

{

  /* Input file */

  _nx1[0] = (_ndim >= 1) ? db->getNX(0) : 1;
  _nx1[1] = (_ndim >= 2) ? db->getNX(1) : 1;
  _nx1[2] = (_ndim >= 3) ? db->getNX(2) : 1;
  _dx1[0] = (_ndim >= 1) ? db->getDX(0) : 1.;
  _dx1[1] = (_ndim >= 2) ? db->getDX(1) : 1.;
  _dx1[2] = (_ndim >= 3) ? db->getDX(2) : 1.;
  _x01[0] = (_ndim >= 1) ? db->getX0(0) : 0.;
  _x01[1] = (_ndim >= 2) ? db->getX0(1) : 0.;
  _x01[2] = (_ndim >= 3) ? db->getX0(2) : 0.;

  /* Output file */

  _nx2[0] = _nx1[0] * 2 + 1;
  _nx2[1] = _nx1[1] * 2 + 1;
  _nx2[2] = _nx1[2];
  _dx2[0] = _dx1[0] / 2.;
  _dx2[1] = _dx1[1] / 2.;
  _dx2[2] = _dx1[2];
  _x02[0] = _x01[0] - _dx2[0];
  _x02[1] = _x01[1] - _dx2[1];
  _x02[2] = _x01[2];
}

/****************************************************************************/
/*!
 **  Define the characteristics from Db2 to Db1
 **
 ** \param[in]  db  Starting grid Db structure
 **
 *****************************************************************************/
void SimuRefine::_dim_2_to_1(DbGrid *db)

{

  /* Input file */

  _nx2[0] = (_ndim >= 1) ? db->getNX(0) : 1;
  _nx2[1] = (_ndim >= 2) ? db->getNX(1) : 1;
  _nx2[2] = (_ndim >= 3) ? db->getNX(2) : 1;
  _dx2[0] = (_ndim >= 1) ? db->getDX(0) : 1.;
  _dx2[1] = (_ndim >= 2) ? db->getDX(1) : 1.;
  _dx2[2] = (_ndim >= 3) ? db->getDX(2) : 1.;
  _x02[0] = (_ndim >= 1) ? db->getX0(0) : 0.;
  _x02[1] = (_ndim >= 2) ? db->getX0(1) : 0.;
  _x02[2] = (_ndim >= 3) ? db->getX0(2) : 0.;

  /* Output file */

  _nx1[0] = _nx2[0] - 2;
  _nx1[1] = _nx2[1] - 2;
  _nx1[2] = _nx2[2];
  _dx1[0] = _dx2[0];
  _dx1[1] = _dx2[1];
  _dx1[2] = _dx2[2];
  _x01[0] = _x02[0] + _dx2[0];
  _x01[1] = _x02[1] + _dx2[1];
  _x01[2] = _x02[2];
}

/****************************************************************************/
/*!
 **  Establish and solve the different Kriging systems
 **
 ** \return  Error return code
 **
 *****************************************************************************/
int SimuRefine::_kriging_define()

{
  /* Define the kriging system for the cell centers */

  _neigh_simfine(0, 0, -1, -1,  0);
  _neigh_simfine(0, 1,  1, -1,  0);
  _neigh_simfine(0, 2,  1,  1,  0);
  _neigh_simfine(0, 3, -1,  1,  0);
  _neigh_simfine(0, 4,  0,  0, -1);

  if (_kriging_solve(0, 0, 4)) return (1);
  if (_kriging_solve(0, 1, 5)) return (1);

  /* Define the Kriging system for the mid-vertices */

  _neigh_simfine(1, 0, -1,  0,  0);
  _neigh_simfine(1, 1,  0, -1,  0);
  _neigh_simfine(1, 2,  1,  0,  0);
  _neigh_simfine(1, 3,  0,  1,  0);
  _neigh_simfine(1, 4,  0,  0, -1);

  if (_kriging_solve(1, 0, 4)) return (1);
  if (_kriging_solve(1, 1, 5)) return (1);

  return (0);
}

/****************************************************************************/
/*!
 **  Define the location of a neighborhood point
 **
 ** \param[in]  type   Type of kriging
 ** \param[in]  rank   Rank of the neighboring data
 ** \param[in]  idx    Shift along X
 ** \param[in]  idy    Shift along Y
 ** \param[in]  idz    Shift along Z
 **
 *****************************************************************************/
void SimuRefine::_neigh_simfine(int type, int rank, int idx, int idy, int idz)
{
  _IX[type][rank] = idx;
  _IY[type][rank] = idy;
  _IZ[type][rank] = idz;
  _XN[type][rank] = idx * _dx2[0];
  _YN[type][rank] = idy * _dx2[1];
  _ZN[type][rank] = idz * _dx2[2];
}

/****************************************************************************/
/*!
 **  Copy the initial information from db1 into db2
 **
 ** \param[in]  db1     Input grid Db structure
 ** \param[in]  iatt1   Rank of the attribute to be read from db1
 ** \param[in]  db2     Output grid Db structure
 ** \param[in]  iatt2   Rank of the attribute to be written into db2
 **
 *****************************************************************************/
void SimuRefine::_merge_data(DbGrid *db1, int iatt1, DbGrid *db2, int iatt2)
{
  for (int ix1 = 0; ix1 < _nx1[0]; ix1++)
    for (int iy1 = 0; iy1 < _nx1[1]; iy1++)
      for (int iz1 = 0; iz1 < _nx1[2]; iz1++)
      {
        int ix2 = 1 + 2 * ix1;
        int iy2 = 1 + 2 * iy1;
        int iz2 = iz1;
        double value = _read(db1, iatt1, ix1, iy1, iz1, 0, 0, 0);
        _write(db2, iatt2, ix2, iy2, iz2, value);
      }
}

/****************************************************************************/
/*!
 **  Read a value in a Db
 **
 ** \return  The value read
 **
 ** \param[in]  db     Grid Db structure
 ** \param[in]  iatt   Rank of the attribute to be read
 ** \param[in]  ix0    Index of the target along X
 ** \param[in]  iy0    Index of the target along Y
 ** \param[in]  iz0    Index of the target along Z
 ** \param[in]  idx    Shift along X
 ** \param[in]  idy    Shift along Y
 ** \param[in]  idz    Shift along Z
 **
 *****************************************************************************/
double SimuRefine::_read(DbGrid *db,
                         int iatt,
                         int ix0,
                         int iy0,
                         int iz0,
                         int idx,
                         int idy,
                         int idz)
{
  VectorInt ind(3);
  int ix = ix0 + idx;
  if (ix < 0 || ix >= db->getNX(0)) ix = ix0 - idx;
  int iy = iy0 + idy;
  if (iy < 0 || iy >= db->getNX(1)) iy = iy0 - idy;
  int iz = iz0 + idz;
  if (iz < 0 || iz >= db->getNX(2)) iz = iz0 - idz;
  ind[0] = ix;
  ind[1] = iy;
  ind[2] = iz;
  int iad = db->indiceToRank(ind);
  return db->getArray(iad, iatt);
}

/****************************************************************************/
/*!
 **  Write a value in a Db
 **
 ** \param[in]  db     Grid Db structure
 ** \param[in]  iatt   Rank of the attribute to be written
 ** \param[in]  ix0    Index of the target along X
 ** \param[in]  iy0    Index of the target along Y
 ** \param[in]  iz0    Index of the target along Z
 ** \param[in]  value  Value to be written
 **
 *****************************************************************************/
void SimuRefine::_write(DbGrid *db, int iatt, int ix0, int iy0, int iz0, double value)
{
  VectorInt ind(3);
  ind[0] = ix0;
  ind[1] = iy0;
  ind[2] = iz0;
  int iad = db->indiceToRank(ind);
  db->setArray(iad, iatt, value);
}

/****************************************************************************/
/*!
 **  Truncate the resulting information from db2 into db1
 **
 ** \param[in]  db2     Input grid Db structure
 ** \param[in]  iatt2   Rank of the attribute to be read from db2
 ** \param[in]  db1     Output grid Db structure
 ** \param[in]  iatt1   Rank of the attribute to be written into db1
 **
 *****************************************************************************/
void SimuRefine::_truncate_result(DbGrid *db2, int iatt2, DbGrid *db1, int iatt1)
{
  for (int ix = 0; ix < _nx1[0]; ix++)
    for (int iy = 0; iy < _nx1[1]; iy++)
      for (int iz = 0; iz < _nx1[2]; iz++)
      {
        double value = _read(db2, iatt2, ix, iy, iz, 1, 1, 0);
        _write(db1, iatt1, ix, iy, iz, value);
      }
}

/****************************************************************************/
/*!
 **  Solve the kriging system
 **
 ** \return  Error return code
 **
 ** \param[in]  type    Type of kriging
 ** \param[in]  rank    Rank of the neighboring data
 ** \param[in]  nb      Number of equations
 ** \param[in]  verbose Verbose flag
 **
 *****************************************************************************/
int SimuRefine::_kriging_solve(int type,
                               int rank,
                               int nb,
                               bool verbose)
{
  int neq = (_param.isFlagSK()) ? nb : nb + 1;
  VectorDouble d1(3);
  VectorDouble lhs(36);
  VectorDouble rhs(6);

  CovCalcMode mode;
  mode.setMember(ECalcMember::RHS);

  /* Establish the kriging L.H.S. */

  for (int i = 0; i < nb; i++)
    for (int j = 0; j < nb; j++)
    {
      d1[0] = _XN[type][i] - _XN[type][j];
      d1[1] = _YN[type][i] - _YN[type][j];
      d1[2] = _ZN[type][i] - _ZN[type][j];
      model_calcul_cov(NULL, _model, mode, 1, 1., d1, lhs.data());
    }

  /* Establish the kriging R.H.S. */

  for (int i = 0; i < nb; i++)
  {
    d1[0] = _XN[type][i];
    d1[1] = _YN[type][i];
    d1[2] = _ZN[type][i];
    model_calcul_cov(NULL, _model, mode, 1, 1., d1, rhs.data());
  }

  /* Add the Universality condition (optional) */

  if (! _param.isFlagSK())
  {
    for (int i = 0; i < nb; i++)
    {
      LHS(i,nb) = 1.;
      LHS(nb,i) = 1.;
    }
    LHS(nb,nb) = 0;
    RHS(nb) = 1.;
  }

  /* Derive the Kriging weights */

  if (matrix_invert(lhs.data(), neq, -1))
  {
    messerr("Kriging matrix inversion failed");
    messerr("Check the consistency between the model and the SK/OK option");
    return (1);
  }
  matrix_product(neq, neq, 1, lhs.data(), rhs.data(), _WGT[type][rank]);

  /* Calculate the variance */

  mode.setMember(ECalcMember::VAR);
  for (int i = 0; i < 3; i++)
    d1[i] = 0.;
  double var[2];
  model_calcul_cov(NULL,_model, mode, 1, 1., d1, &var[0]);
  matrix_product(1, neq, 1, rhs.data(),_WGT[type][rank], &var[1]);
  double variance = var[0] - var[1];
  _STDV[type][rank] = (variance > 0) ? sqrt(variance) : 0.;

  /* Printout of the weights */

  if (verbose)
  {
    message("\nDisplay of the Kriging weights\n");
    for (int i = 0; i < nb; i++)
      message("X=%10.3lf Y=%10.3lf Z=%10.3lf W=%10.6lf\n",
              _XN[type][i], _YN[type][i], _ZN[type][i], _WGT[type][rank][i]);
    message("Variance of error           = %10.6lf\n", variance);
    message("Standard deviation of error = %10.6lf\n", _STDV[type][rank]);
  }

  return (0);
}

/****************************************************************************/
/*!
 **  Simulate the missing nodes
 **
 ** \param[in]  db     Grid Db structure
 ** \param[in]  iatt   Rank of the column
 **
 *****************************************************************************/
void SimuRefine::_simulate_nodes(DbGrid *db, int iatt)
{
  for (int iz = 0; iz < _nx2[2]; iz++)
    for (int ix = 0; ix < _nx2[0]; ix++)
      for (int iy = 0; iy < _nx2[1]; iy++)
        if ((ix % 2 == 0) && (iy % 2 == 0))
          _simulate_target(db, 0, iatt, ix, iy, iz);

  /* Perform the cell mid-vertices */

  for (int iz = 0; iz < _nx2[2]; iz++)
    for (int ix = 0; ix < _nx2[0]; ix++)
      for (int iy = 0; iy < _nx2[1]; iy++)
        if (((ix % 2 == 0) && (iy % 2 == 1)) || ((ix % 2 == 1) && (iy % 2 == 0)))
          _simulate_target(db, 1, iatt, ix, iy, iz);
}

/****************************************************************************/
/*!
 **  Simulate the target cell
 **
 ** \param[in]  db     Grid Db structure
 ** \param[in]  type   Type of kriging
 ** \param[in]  iatt   Rank of the attribute to be read
 ** \param[in]  ix0    Index of the target along X
 ** \param[in]  iy0    Index of the target along Y
 ** \param[in]  iz0    Index of the target along Z
 **
 *****************************************************************************/
void SimuRefine::_simulate_target(DbGrid *db,
                                  int type,
                                  int iatt,
                                  int ix0,
                                  int iy0,
                                  int iz0)
{
  double value = 0.;
  if (iz0 == 0)
  {

    /* Case of the first layer */

    for (int i = 0; i < 4; i++)
      value += (_WGT[type][0][i]
          * _read(db, iatt, ix0, iy0, iz0, _IX[type][i], _IY[type][i], _IZ[type][i]));
    value += _STDV[type][0] * law_gaussian();
  }
  else
  {

    /* Case of a subsequent layer */

    for (int i = 0; i < 5; i++)
      value += (_WGT[type][1][i]
          * _read(db, iatt, ix0, iy0, iz0, _IX[type][i], _IY[type][i], _IZ[type][i]));
    value += _STDV[type][1] * law_gaussian();
  }

  _write(db, iatt, ix0, iy0, iz0, value);
}
