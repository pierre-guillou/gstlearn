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
#include "LinearOp/ShiftOpStencil.hpp"

#include "Basic/AStringable.hpp"
#include "Basic/Indirection.hpp"
#include "Basic/VectorNumT.hpp"
#include "LinearOp/AShiftOp.hpp"
#include "LinearOp/ShiftOpMatrix.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Covariances/CovAniso.hpp"
#include "Basic/Grid.hpp"

#include "geoslib_define.h"

ShiftOpStencil::ShiftOpStencil(const MeshETurbo* mesh,
                               const CovAniso* cova,
                               bool verbose)
  : AShiftOp()
  , _relativeShifts()
  , _absoluteShifts()
  , _weights()
  , _isInside()
  , _lambdaVal(0.)
  , _useLambdaSingleVal(true)
  , _useModifiedShift(false)
  , _mesh()
{
  if (_buildInternal(mesh, cova, verbose)) return;
}

ShiftOpStencil::ShiftOpStencil(const ShiftOpStencil& shift)
  : AShiftOp(shift)
  , _relativeShifts(shift._relativeShifts)
  , _absoluteShifts(shift._absoluteShifts)
  , _weights(shift._weights)
  , _weightsSimu(shift._weightsSimu)
  , _isInside(shift._isInside)
  , _lambdaVal(shift._lambdaVal)
  , _useLambdaSingleVal(shift._useLambdaSingleVal)
  , _useModifiedShift(shift._useModifiedShift)
  , _mesh(shift._mesh)
{
}

ShiftOpStencil& ShiftOpStencil::operator=(const ShiftOpStencil& shift)
{
  AShiftOp::operator=(shift);
  if (this != &shift)
  {
    _relativeShifts = shift._relativeShifts;
    _absoluteShifts = shift._absoluteShifts;
    _weights        = shift._weights;
    _weightsSimu    = shift._weightsSimu;
    _isInside       = shift._isInside;
    _lambdaVal      = shift._lambdaVal;
    _useLambdaSingleVal = shift._useLambdaSingleVal;
    _useModifiedShift   = shift._useModifiedShift;
    _mesh           = shift._mesh;
  }
  return *this;
}

ShiftOpStencil::~ShiftOpStencil() {}

int ShiftOpStencil::_addToDest(const constvect inv, vect outv) const
{
  const VectorDouble* currentWeights;
  if (_useModifiedShift)
  {
    currentWeights = &_weightsSimu;
  }
  else 
  {
    currentWeights = &_weights;
  }

  int nw = _getNWeights();
  int size = (int) inv.size();
  const Indirection& indirect = _mesh->getGridIndirect();

  double total;
  if (!indirect.isDefined())
  {
    // Use the fast option when no selection is defined on the Grid
    for (int ic = 0; ic < size; ic++)
    {
      total = 0.;
      if (_isInside[ic])
      {
        for (int iw = 0; iw < nw; iw++)
        {
          int iabs = ic + _absoluteShifts[iw];
          double value = _isInside[iabs] ? inv[iabs] : 0.;
          total += (*currentWeights)[iw] * value;
        }
      }
      outv[ic] = total;
    }
  }
  else
  {
    const Grid& grid = _mesh->getGrid();
    int ndim         = _mesh->getNDim();

    VectorInt center(ndim);
    VectorInt local(ndim);
    for (int ic = 0; ic < size; ic++)
    {
      total = 0.;

      // Check if the target point is not on the edge and not masked
      if (_isInside[ic] && indirect.getAToR(ic) >= 0)
      {
        grid.rankToIndice(ic, center);
        for (int iw = 0; iw < nw; iw++)
        {
          local = center;
          VH::addInPlace(local, _relativeShifts[iw]);
          int ie = grid.indiceToRank(local);
          if (indirect.getAToR(ie) >= 0) total += (*currentWeights)[iw] * inv[ie];
        }
      }
      outv[ic] = total;
    }
  }
  return 0;
}

void ShiftOpStencil::resetModif()
{
  _useModifiedShift = false;
}

void ShiftOpStencil::normalizeLambdaBySills(const AMesh* mesh)
{
  if (_cova->isNoStatForVariance())
  {
    _Lambda.resize(_napices);
    for (auto &e : _Lambda)
    {
      e = _lambdaVal;
    }
    AShiftOp::normalizeLambdaBySills(mesh);
    _useLambdaSingleVal = false;
  }
  else 
  {
    _lambdaVal /= sqrt(_cova->getSill(0,0));
  }
}

double ShiftOpStencil::getMaxEigenValue() const
{
  double s = 0.;
  for (const auto &e : _weights)
  {
    s += ABS(e);
  }
  return s;
}

double ShiftOpStencil::getLambda(int iapex) const
{
  if (_useLambdaSingleVal) return _lambdaVal;
  return AShiftOp::getLambda(iapex);
}

void ShiftOpStencil::multiplyByValueAndAddDiagonal(double v1, double v2) 
{
  _weightsSimu = VectorDouble(_weights.size());
  for (int i = 0; i < (int)_weights.size(); i++)
    _weightsSimu[i] = v1 * _weights[i];
  int center = (int) _relativeShifts.size() / 2;
  _weightsSimu[center] += v2;
  _useModifiedShift = true;
}

int ShiftOpStencil::_buildInternal(const MeshETurbo* mesh,
                                   const CovAniso* cova,
                                   bool verbose)
{
  // Preliminary checks

  if (cova == nullptr)
  {
    messerr("The argument 'cova' must be provided");
    return 1;
  }
  if (mesh == nullptr)
  {
    messerr("The argument 'mesh' must be provided");
    return 1;
  }

  _mesh = mesh;
  _napices = mesh->getNApices();
  int ndim        = _mesh->getNDim();

  _setCovAniso(cova);

  if (_cova->isNoStatForAnisotropy())
  {
    messerr("The Shiftop as a Stencil is incompatible with non-stationarity");
    return 1;
  }

  // Create a local Turbo Meshing starting from a DbGrid (Dimension 5 sgould be sufficient)
  const Grid& grid = mesh->getGrid();
  VectorInt NXs    = grid.getNXs();
  VectorInt nxlocal(ndim, 5);

  MeshETurbo localMesh(nxlocal, grid.getDXs(), grid.getX0s(), grid.getRotAngles());
  ShiftOpMatrix shiftMat(&localMesh, cova, nullptr, verbose);

  // Display the vector of the 'S' matrix for the center Apex
  MatrixSparse* S           = shiftMat.getS();
  int centerApex            = localMesh.getNApices() / 2;
  VectorDouble centerColumn = S->getColumn(centerApex);

  // Fill lambda
  _lambdaVal = shiftMat.getLambda(centerApex);

  // Get the indices of the centerApex
  VectorInt center(ndim);
  VectorInt other(ndim);
  localMesh.getApexIndicesInPlace(centerApex, center);

  // Get the non-zero elements of the center column
  _relativeShifts.clear();
  _weights.clear();
  for (int i = 0, n = (int)centerColumn.size(); i < n; i++)
  {
    double weight = centerColumn[i];
    if (ABS(weight) < EPSILON6) continue;
    localMesh.getApexIndicesInPlace(i, other);
    VH::subtractInPlace(other, center);
    _relativeShifts.push_back(other);
    _weights.push_back(weight);
  }
  int nw = _getNWeights();

  // Calculate the shifts (from the center cell) for each weight
  // This is calculated for a reference pixel (center of the grid)
  _absoluteShifts.fill(0, nw);
  for (int idim = 0; idim < ndim; idim++) center[idim] = NXs[idim] / 2;
  int iorigin = grid.indiceToRank(center);

  for (int iw = 0; iw < nw; iw++)
  {
    VectorInt local = center;
    VH::addInPlace(local, _relativeShifts[iw]);
    int iabs = grid.indiceToRank(local);
    _absoluteShifts[iw] = iabs - iorigin;
  }

  // Delineate the border of the grid (not to be treated)
  int size = _mesh->getNApices();
  _isInside.fill(true, size);

  for (int i = 0; i < size; i++)
  {
    grid.rankToIndice(i, center);
    bool flagInside = true;
    for (int idim = 0; idim < ndim && flagInside; idim++)
    {
      int ival = center[idim];
      if (ival <= 0 || ival >= NXs[idim] - 1) flagInside = false;
    }
    _isInside[i] = flagInside;
  }

  // Print the contents of non-zero elements
  if (verbose) _printStencil();

  return 0;
}

void ShiftOpStencil::_printStencil() const
{
  int nweight = _getNWeights();
  int ndim    = _mesh->getNDim();
  int size    = _mesh->getNApices();
  mestitle(0, "Stencil contents");
  for (int i = 0; i < nweight; i++)
  {
    message("Weight %d/%d - Relative (", i + 1, nweight);
    for (int idim = 0; idim < ndim; idim++)
      message("%2d ", _relativeShifts[i][idim]);
    message(") - Absolute (%4d)", _absoluteShifts[i]);
    message(" : %lf\n", _weights[i]);
  }

  int ntreated = 0;
  for (int i = 0; i < size; i++)
    if (_isInside[i]) ntreated++;
  message("Number of pixels inside the grid (no edge effect) = %d/%d\n", ntreated, size);
}