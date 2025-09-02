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
#include "Mesh/MeshETurbo.hpp"
#include "Basic/Grid.hpp"
#include "Basic/SerializeHDF5.hpp"
#include "Covariances/CovAniso.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Geometry/Rotation.hpp"
#include "LinearOp/ProjMatrix.hpp"
#include "Matrix/MatrixSquare.hpp"
#include "Matrix/NF_Triplet.hpp"
#include "Mesh/AMesh.hpp"
#include "Mesh/Delaunay.hpp"

#include <cmath>

namespace gstlrn
{
MeshETurbo::MeshETurbo(Id mode)
  : AMesh()
  , _grid()
  , _nPerCell(0)
  , _isPolarized(false)
  , _meshIndirect(mode)
  , _gridIndirect(mode)
{
}

MeshETurbo::MeshETurbo(const VectorInt& nx,
                       const VectorDouble& dx,
                       const VectorDouble& x0,
                       const VectorDouble& angles,
                       bool flag_polarized,
                       bool verbose,
                       Id mode)
  : AMesh()
  , _grid()
  , _nPerCell(0)
  , _isPolarized(flag_polarized)
  , _meshIndirect(mode)
  , _gridIndirect(mode)
{
  (void)initFromGridByAngles(nx, dx, x0, angles, VectorDouble(), flag_polarized, verbose);
}

MeshETurbo::MeshETurbo(const DbGrid* dbgrid,
                       bool flag_polarized,
                       bool verbose,
                       Id mode)
  : AMesh()
  , _grid()
  , _nPerCell(0)
  , _isPolarized(flag_polarized)
  , _meshIndirect(mode)
  , _gridIndirect(mode)
{
  if (!dbgrid->isGrid()) return;
  VectorDouble sel = dbgrid->getSelections();
  (void)initFromGridByMatrix(dbgrid->getNXs(), dbgrid->getDXs(),
                             dbgrid->getX0s(), dbgrid->getRotMat(), sel,
                             flag_polarized, verbose);
}

MeshETurbo::MeshETurbo(const MeshETurbo& r)
  : AMesh(r)
  , _grid()
  , _nPerCell(0)
  , _isPolarized(r._isPolarized)
  , _meshIndirect(r._meshIndirect)
  , _gridIndirect(r._gridIndirect)
{
  _grid = r._grid;
}

MeshETurbo& MeshETurbo::operator=(const MeshETurbo& r)
{
  _grid         = r._grid;
  _nPerCell     = r._nPerCell;
  _isPolarized  = r._isPolarized;
  _meshIndirect = r._meshIndirect;
  _gridIndirect = r._gridIndirect;

  return *this;
}

MeshETurbo::~MeshETurbo()
{
}

/**
 * Returns the total number of apices of the whole grid
 * (not accounting for possible mask on meshes)
 * @return
 */
Id MeshETurbo::getNApices() const
{
  if (_gridIndirect.isDefined()) return _gridIndirect.getRelSize();
  return _grid.getNTotal();
}

/**
 * Returns the total number of possible meshes built using the whole grid
 * (not accounting for possible mask on triangles)
 * @return
 */
Id MeshETurbo::_nmeshInCompleteGrid() const
{
  Id nmesh = 1;
  for (Id idim = 0, ndim = getNDim(); idim < ndim; idim++)
    nmesh *= (_grid.getNX(idim) - 1);
  nmesh *= _nPerCell;
  // _meshIndirect.getRelSize();
  return nmesh;
}

/**
 * Actual number of (active) meshes
 * @return
 */
Id MeshETurbo::getNMeshes() const
{
  if (_meshIndirect.isDefined()) return _meshIndirect.getRelSize();
  return _nmeshInCompleteGrid();
}

double MeshETurbo::getMeshSize(Id /*imesh*/) const
{
  double size = 1.;
  for (Id idim = 0, ndim = getNDim(); idim < ndim; idim++)
    size *= _grid.getDX(idim);
  size /= _nPerCell;
  return size;
}

/****************************************************************************/
/*!
** Returns the Apex 'rank' of the Mesh 'imesh'
**
** \returns The rank of the target apex
**
** \param[in]  imesh    Rank of active Mesh (starting from 0)
** \param[in]  rank     Rank of Apex within a Mesh (from 0 to _nApexPerMesh-1)
**
*****************************************************************************/
Id MeshETurbo::getApex(Id imesh, Id rank) const
{
  Id node, icas;
  auto ndim = getNDim();
  _indg.resize(ndim);

  Id jmesh = _meshIndirect.getRToA(imesh);
  _getGridFromMesh(jmesh, &node, &icas);
  _grid.rankToIndice(node, _indg);
  auto ipol = _getPolarized(_indg);

  for (Id idim = 0; idim < ndim; idim++)
    _indg[idim] += MSS(ndim, ipol, icas, rank, idim);

  Id igrid = _grid.indiceToRank(_indg);

  Id irel = _gridIndirect.getAToR(igrid);
  if (irel < 0)
  {
    messerr("Problem for mesh=%d rank=%d grid=%d -> Mesh relative rank is negative",
            imesh + 1, rank + 1, igrid + 1);
  }
  return irel;
}

double MeshETurbo::getCoor(Id imesh, Id rank, Id idim) const
{
  _indg.resize(getNDim());

  auto irel = getApex(imesh, rank);
  Id iabs = _gridIndirect.getRToA(irel);
  _grid.rankToIndice(iabs, _indg);
  return _grid.indiceToCoordinate(idim, _indg);
}

void MeshETurbo::getCoordinatesPerMeshInPlace(Id imesh, Id rank, VectorDouble& coords) const
{
  _indg.resize(getNDim());

  auto irel = getApex(imesh, rank);
  Id iabs = _gridIndirect.getRToA(irel);
  _grid.rankToCoordinatesInPlace(iabs, coords);
}

double MeshETurbo::getApexCoor(Id i, Id idim) const
{
  // _meshIndirect.getRelSize();
  _indg.resize(getNDim());

  Id iabs = _gridIndirect.getRToA(i);
  _grid.rankToIndice(iabs, _indg);
  return _grid.indiceToCoordinate(idim, _indg);
}

void MeshETurbo::getApexIndicesInPlace(Id i, VectorInt& indg) const
{
  // _meshIndirect.getRelSize();
  indg.resize(getNDim());

  Id iabs = _gridIndirect.getRToA(i);
  _grid.rankToIndice(iabs, indg);
}

void MeshETurbo::getApexCoordinatesInPlace(Id i, VectorDouble& coords) const
{
  _indg.resize(getNDim());

  Id iabs = _gridIndirect.getRToA(i);
  _grid.rankToIndice(iabs, _indg);

  for (Id idim = 0; idim < getNDim(); idim++)
    coords[idim] = _grid.indiceToCoordinate(idim, _indg);
}

Id MeshETurbo::initFromGridByAngles(const VectorInt& nx,
                                     const VectorDouble& dx,
                                     const VectorDouble& x0,
                                     const VectorDouble& angles,
                                     const VectorDouble& sel,
                                     bool flag_polarized,
                                     bool verbose)
{
  Id ndim = static_cast<Id>(nx.size());
  _setNDim(ndim);

  /* Create the internal (rotated) grid by Angles */

  if (_grid.resetFromVector(nx, dx, x0)) return 1;
  _grid.setRotationByAngles(angles);

  return _initFromGridInternal(sel, flag_polarized, verbose);
}

Id MeshETurbo::initFromGridByMatrix(const VectorInt& nx,
                                     const VectorDouble& dx,
                                     const VectorDouble& x0,
                                     const VectorDouble& rotmat,
                                     const VectorDouble& sel,
                                     bool flag_polarized,
                                     bool verbose)
{
  Id ndim = static_cast<Id>(nx.size());
  _setNDim(ndim);

  /* Create the internal (rotated) grid by Rotation Matrix */

  if (_grid.resetFromVector(nx, dx, x0)) return 1;
  _grid.setRotationByVector(rotmat);

  return _initFromGridInternal(sel, flag_polarized, verbose);
}

Id MeshETurbo::_initFromGridInternal(const VectorDouble& sel,
                                      bool flag_polarized,
                                      bool verbose)
{
  auto ndim = getNDim();

  // Get grid extension
  // TODO: the grid extension should be calculated in Grid and take
  // case of a possible rotation
  VectorDouble extendmin(ndim);
  VectorDouble extendmax(ndim);
  for (Id idim = 0; idim < ndim; idim++)
  {
    extendmin[idim] = _grid.getX0(idim);
    extendmax[idim] =
      _grid.getX0(idim) + (_grid.getNX(idim) - 1) * _grid.getDX(idim);
  }
  if (_setExtend(extendmin, extendmax)) return (1);

  // Define the number of Elements per Cell

  _setNElementPerCell();

  // Set polarization

  _isPolarized = flag_polarized;

  // Convert optional selection array into mask on meshes

  _buildMaskInMeshing(sel);

  // Optional printout

  if (verbose) messageFlush(toString());

  return 0;
}

void MeshETurbo::_buildMaskInMeshing(const VectorDouble& sel)
{
  Id node, icas;
  std::map<Id, Id> map;

  // If no selection is defined on the grid, the vector of Meshing Mask is cancelled
  if (sel.empty()) return;

  // Creating the Masking information for Meshing 'meshActiveToAbsolute'
  // which gives the Absolute meshing index from its Active index

  auto ndim    = getNDim();
  Id nmesh   = _nmeshInCompleteGrid();
  auto ncorner = getNApexPerMesh();
  VectorInt indg0(ndim);
  _indg.resize(ndim);

  // Loop on all possible meshes
  Id meshNactive = 0;
  for (Id imesh = 0; imesh < nmesh; imesh++)
  {
    _getGridFromMesh(imesh, &node, &icas);
    _grid.rankToIndice(node, indg0);
    auto ipol = _getPolarized(indg0);

    // Loop on the corners of the mesh (polarization is taken into account)
    bool flagMasked = false;
    for (Id icorner = 0; icorner < ncorner && !flagMasked; icorner++)
    {

      // Generate the indices of the mesh apex
      for (Id idim = 0; idim < ndim; idim++)
        _indg[idim] = indg0[idim] + MSS(ndim, ipol, icas, icorner, idim);
      Id iad = _grid.indiceToRank(_indg);
      if (sel[iad] == 0.) flagMasked = true;
    }

    // The triangle is not masked, store its index
    if (flagMasked) continue;
    map[imesh] = meshNactive;
    meshNactive++;
  }

  // Creating the Indirection information for Meshing

  _meshIndirect.buildFromMap(map, nmesh);

  // Creating the selection of the active grid nodes
  // It is at least equal to 'sel'. In addition, non active grid nodes are also discarded

  VectorDouble selbis(sel.size(), 0.);
  for (Id imesh = 0; imesh < meshNactive; imesh++)
  {
    Id jmesh = _meshIndirect.getRToA(imesh);
    _getGridFromMesh(jmesh, &node, &icas);
    _grid.rankToIndice(node, indg0);
    auto ipol = _getPolarized(indg0);
    for (Id icorner = 0; icorner < ncorner; icorner++)
    {
      for (Id idim = 0; idim < ndim; idim++)
        _indg[idim] = indg0[idim] + MSS(ndim, ipol, icas, icorner, idim);
      Id iad     = _grid.indiceToRank(_indg);
      selbis[iad] = 1;
    }
  }

  _gridIndirect.buildFromSel(selbis);
}

/**
 * Create a MeshETurbo by loading the contents of a Neutral File
 *
 * @param NFFilename Name of the Neutral File (MeshEStandard format)
 * @param verbose    Verbose
 */
MeshETurbo* MeshETurbo::createFromNF(const String& NFFilename, bool verbose)
{
  MeshETurbo* mesh = new MeshETurbo;
  if (mesh->_fileOpenAndDeserialize(NFFilename, verbose)) return mesh;
  delete mesh;
  return nullptr;
}

MeshETurbo* MeshETurbo::create(const VectorInt& nx,
                               const VectorDouble& dx,
                               const VectorDouble& x0,
                               const VectorDouble& angles,
                               bool flag_polarized,
                               bool verbose)
{
  auto* mesh = new MeshETurbo(nx, dx, x0, angles, flag_polarized, verbose);
  return mesh;
}

MeshETurbo* MeshETurbo::createFromGrid(const DbGrid* dbgrid,
                                       bool flag_polarized,
                                       bool verbose,
                                       Id mode)
{
  auto* mesh = new MeshETurbo(dbgrid, flag_polarized, verbose, mode);
  return mesh;
}

MeshETurbo* MeshETurbo::createFromGridInfo(const Grid* grid,
                                           bool flag_polarized,
                                           bool verbose,
                                           Id mode)
{
  auto* mesh = new MeshETurbo(mode);
  if (mesh->initFromGridByMatrix(grid->getNXs(), grid->getDXs(), grid->getX0s(),
                                 grid->getRotMat(), VectorDouble(), flag_polarized,
                                 verbose))
    return nullptr;
  return mesh;
}

MeshETurbo* MeshETurbo::createFromCova(const CovAniso& cova,
                                       const Db* field,
                                       double ratio,
                                       Id nbExt,
                                       bool isPolarized,
                                       bool useSel,
                                       Id nxmax,
                                       bool verbose)
{
  auto* mesh = new MeshETurbo();
  if (mesh->initFromCova(cova, field, ratio, nbExt, isPolarized, useSel,
                         nxmax, verbose))
    return nullptr;
  return mesh;
}

/****************************************************************************/
/*!
** Create the meshing
**
** \param[in]  extendmin       Minimum of the dilated rotated bounding box
** \param[in]  extendmax       Minimum of the dilated rotated bounding box
** \param[in]  cellsize        Array giving the cell size
** \param[in]  rotmat          Rotation matrix (optional)
** \param[in]  flag_polarized  Switching ON/OFF the polarization
** \param[in]  verbose         Verbose flag
**
*****************************************************************************/
Id MeshETurbo::initFromExtend(const VectorDouble& extendmin,
                               const VectorDouble& extendmax,
                               const VectorDouble& cellsize,
                               const VectorDouble& rotmat,
                               bool flag_polarized,
                               bool verbose)
{
  Id ndim = static_cast<Id>(extendmin.size());
  _setNDim(ndim);
  if (_setExtend(extendmin, extendmax)) return (1);

  /* Create the internal (rotated) grid */

  if (_defineGrid(cellsize)) return (1);
  _grid.setRotationByVector(rotmat);

  // Define the number of Elements per Cell

  _setNElementPerCell();

  // Set polarization

  _isPolarized = flag_polarized;

  // Optional printout

  if (verbose) messageFlush(toString());

  return 0;
}

bool MeshETurbo::_addElementToTriplet(NF_Triplet& NF_T,
                                      Id iech,
                                      const VectorDouble& coor,
                                      const VectorInt& indg0,
                                      bool verbose) const
{
  auto ncorner = getNApexPerMesh();
  _indices.resize(ncorner);
  _lambdas.resize(ncorner);

  for (Id icas = 0; icas < _nPerCell; icas++)
  {
    if (_addWeights(icas, indg0, coor, _indices, _lambdas, verbose) == 0)
    {
      for (Id icorner = 0; icorner < ncorner; icorner++)
        NF_T.add(iech, _indices[icorner], _lambdas[icorner]);
      return true;
    }
  }
  return false;
}

/**
 * @brief Given the coordinates of a point, return the corresponding mesh index
 * and updates the apex indices
 *
 * @param coor Input coordinates
 * @param indices Returned vector of apex indices
 * @param lambdas Returned vector of weights (barycenter coordinates)
 * @return Rank of the mesh (-1 if point does not belong to the meshing)
 */
Id MeshETurbo::getMeshFromCoordinates(const VectorDouble& coor,
                                       VectorInt& indices,
                                       VectorDouble& lambdas) const
{
  auto ndim = getNDim();
  VectorInt indg0(ndim);
  if (_grid.coordinateToIndicesInPlace(coor, indg0) != 0)
  {
    messerr("The target coordinate does not belong to the Meshing");
    return -1;
  }
  auto ncorner = getNApexPerMesh();
  indices.resize(ncorner);
  lambdas.resize(ncorner);

  VectorInt nx = _grid.getNXs();
  VH::addConstant(nx, -1);

  Id iref = MAX(0, indg0[ndim - 1]);
  for (Id idim = ndim - 2; idim >= 0; idim--)
    iref = iref * nx[idim] + MAX(0, indg0[idim]);
  iref *= _nPerCell;

  for (Id icas = 0; icas < _nPerCell; icas++)
  {
    if (_addWeights(icas, indg0, coor, indices, lambdas) == 0)
      return iref + icas;
  }
  return -1;
}

/****************************************************************************/
/*!
** Returns the Sparse Matrix used to project a Db onto the Meshing
**
** \param[out] m         Projection matrix to be initialized
** \param[in]  db        Db structure
** \param[in]  rankZ     Rank of the Z-locator to be tested (see remarks)
** \param[in]  verbose   Verbose flag
**
** \remarks If rankZ>=0, a sample is only considered if:
**          - the value of the corresponding variable is defined
**          - the sample is covered by the grid of the Turbo Meshing
*****************************************************************************/
void MeshETurbo::resetProjFromDb(ProjMatrix* m,
                                 const Db* db,
                                 Id rankZ,
                                 bool verbose) const
{
  auto ndim = getNDim();
  VectorInt indg0(ndim);
  VectorDouble coor(ndim);
  _grid.initThread();
  // Preliminary checks

  if (isCompatibleDb(db)) return;

  // Core allocation

  NF_Triplet NF_T;

  /* Optional title */

  if (verbose) mestitle(0, "Mesh Barycenter");

  /* Loop on the samples */

  Id iech   = 0;
  Id nout   = 0;
  Id nvalid = 0;
  for (Id jech = 0; jech < db->getNSample(); jech++)
  {
    if (!db->isActive(jech)) continue;
    if (rankZ >= 0)
    {
      if (FFFF(db->getFromLocator(ELoc::Z, jech, rankZ))) continue;
    }
    nvalid++;

    /* Identification */

    db->getCoordinatesInPlace(coor, jech);

    /* Calculate the grid indices */

    if (_grid.coordinateToIndicesInPlace(coor, indg0) != 0)
    {
      messerr("Sample #%d does not belong to the Turbo Meshing internal Grid",
              jech + 1);
      iech++;
      continue;
    }

    /* Optional printout */

    if (verbose)
      message("Sample %4d assigned to Grid Node %4d :",
              jech + 1, _grid.indiceToRank(indg0) + 1);

    // Finding the active mesh to which the sample belongs
    bool found = _addElementToTriplet(NF_T, iech, coor, indg0, verbose);

    // In the case the target coordinate is on the edge of the grid
    // try to shift the point down by one node
    if (!found)
    {
      bool flag_correct = false;
      for (Id idim = 0; idim < ndim; idim++)
      {
        if (indg0[idim] != _grid.getNX(idim) - 1) continue;
        indg0[idim] -= 1;
        flag_correct = true;
      }
      if (flag_correct)
        found = _addElementToTriplet(NF_T, iech, coor, indg0, verbose);
    }

    // The point does not belong to any active mesh, issue a message (optional)
    if (!found)
    {
      nout++;
      if (verbose)
        messerr("Sample #%d does not belong to the meshing", jech + 1);
    }
    iech++;
  }

  /* Add the extreme value to force dimension */

  NF_T.force(nvalid, getNApices());

  /* Convert the triplet into a sparse matrix */

  if (verbose && nout > 0)
    messerr("%d / %d samples which do not belong to the Meshing",
            nout, db->getNSample(true));

  return m->resetFromTriplet(NF_T);
}

/****************************************************************************/
/*!
** Print the contents of the meshing
**
** \param[in] strfmt    Format for printout
**
*****************************************************************************/
String MeshETurbo::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;
  sstr << toTitle(0, "Turbo Meshing");
  if (_isPolarized) sstr << "Diamond construction is activated" << std::endl;
  sstr << _grid.toString(strfmt);
  sstr << AMesh::toString(strfmt);

  if (_meshIndirect.isDefined())
  {
    sstr << toTitle(2, "Mask Information");
    sstr << "Mesh Masking Indexing" << std::endl;
    sstr << _meshIndirect.toString(strfmt) << std::endl;
    sstr << "Grid Masking Indexing" << std::endl;
    sstr << _gridIndirect.toString(strfmt) << std::endl;
  }

  return sstr.str();
}

/****************************************************************************/
/*!
** Define the internal grid
**
** \param[in]  cellsize  Array giving the cell size (see details)
**
*****************************************************************************/
Id MeshETurbo::_defineGrid(const VectorDouble& cellsize)

{
  Id ndim;

  // Initializations

  if (cellsize.empty())
  {
    messerr("The argument 'cellsize' must be provided");
    return (1);
  }
  ndim = getNDim();

  // Create the grid internal structure

  _grid.resetFromSpaceDimension(ndim);

  // Copy the grid main characteristics

  for (Id idim = 0; idim < ndim; idim++)
  {
    _grid.setX0(idim, getExtendMin(idim));
    _grid.setDX(idim, cellsize[idim]);
    _grid.setNX(idim, static_cast<Id>(ceil((getExtendMax(idim) - getExtendMin(idim)) /
                                            cellsize[idim]) +
                                       1));
  }

  return 0;
}

void MeshETurbo::_setNElementPerCell()
{
  auto ndim = getNDim();

  if (ndim == 1)
    _nPerCell = 1;
  else if (ndim == 2)
    _nPerCell = 2;
  else if (ndim == 3)
    _nPerCell = 6;
}

/**
 * Return the weights assigned to the corners
 * @param icas   Corner indication
 * @param indg0  Indices of the starting grid node
 * @param coor   Coordinates of the targte point
 * @param indices Grid indices of the target (in active ranks)
 * @param lambda  Weights
 * @param verbose Verbose flag
 * @return
 *
 * @remark The function returns 1 if:
 * @remark - the grid node corresponding to a mesh apex is outside the grid
 * @remark - the grid node corresponding to a mesh apex is not active
 */
Id MeshETurbo::_addWeights(Id icas,
                            const constvectint indg0,
                            const constvect coor,
                            const vectint indices,
                            const vect lambda,
                            bool verbose) const
{
  auto ndim    = getNDim();
  auto ncorner = getNApexPerMesh();
  auto ipol    = _getPolarized(indg0);
  MatrixSquare lhs;
  _rhs.resize(ncorner);
  _indgg.resize(ndim);

  // Build the LHS matrix

  lhs.reset(ncorner, ncorner);
  for (Id icorner = 0; icorner < ncorner; icorner++)
  {
    // Generate the indices of the mesh apex
    for (Id idim = 0; idim < ndim; idim++)
      _indgg[idim] = indg0[idim] + MSS(ndim, ipol, icas, icorner, idim);
    Id igrid = _grid.indiceToRank(_indgg);
    if (igrid < 0) return 1; // grid node outside grid

    indices[icorner] = _gridIndirect.getAToR(igrid);
    if (indices[icorner] < 0) return 1; // grid node not active

    // Update the LHS matrix
    for (Id idim = 0; idim < ndim; idim++)
      lhs.setValue(idim, icorner, _grid.indiceToCoordinate(idim, _indgg));
    lhs.setValue(ndim, icorner, 1.);
  }

  // Generate the right-hand side vector
  for (Id idim = 0; idim < ndim; idim++)
    _rhs[idim] = coor[idim];
  _rhs[ndim] = 1;

  // Invert the matrix
  if (lhs.invert()) return 1;

  // Calculate the weights
  lhs.prodMatVecInPlaceC(_rhs, lambda);

  // Check that all weights are positive
  for (Id icorner = 0; icorner < ncorner; icorner++)
  {
    if (lambda[icorner] < -EPSILON6) return 1;
    if (lambda[icorner] < 0) lambda[icorner] = 0.;
    if (lambda[icorner] > 1 + EPSILON6) return 1;
    if (lambda[icorner] > 1) lambda[icorner] = 1.;
  }

  // Optional printout
  if (verbose)
  {
    for (Id icorner = 0; icorner < ncorner; icorner++)
      message(" %4d (%4.2lf)", indices[icorner], lambda[icorner]);
    message("\n");
  }

  return 0;
}

Id MeshETurbo::_getPolarized(const constvectint indg) const
{
  auto ndim = getNDim();
  if (!_isPolarized) return (0);

  // Polarization has only been coded for the 2-D case
  if (ndim != 2) return (0);
  if ((indg[0] + indg[1]) % 2 == 1)
    return (0);
  return (1);
}

/**
 * Returns the (starting) grid node, given the absolute rank of the mesh
 * @param imesh Absolute Rank of the  mesh
 * @param node  Starting grid node
 * @param icas  Sorting used for reviewing grid meshes (takes polarization into account)
 */
void MeshETurbo::_getGridFromMesh(Id imesh, Id* node, Id* icas) const
{
  _indg.resize(getNDim());
  Id ncas = _nPerCell;
  Id rank = imesh / ncas;
  *icas    = imesh - rank * ncas;
  _grid.rankToIndice(rank, _indg, true);
  *node = _grid.indiceToRank(_indg);
}

Id MeshETurbo::initFromCova(const CovAniso& cova,
                             const Db* field,
                             double ratio,
                             Id nbExt,
                             bool isPolarized,
                             bool useSel,
                             Id nxmax,
                             bool verbose)
{
  // Initializations
  auto ndim = cova.getNDim();
  Id nval = static_cast<Id>(pow(2., ndim));

  // Get the rotation linked to the covariance
  const Rotation& rot = cova.getAnisoRotation();

  // Project the corners of the grid
  VectorDouble extendMinRot(ndim, TEST);
  VectorDouble extendMaxRot(ndim, TEST);
  VectorDouble cornerRot(ndim);
  VectorInt ic(ndim, 0);
  VectorVectorDouble extremesData;
  VectorDouble cornerRef(ndim, 0.);
  const DbGrid* fieldGrid = nullptr;
  const auto isGrid = field->isGrid();
  if (isGrid)
  {
    fieldGrid = dynamic_cast<const DbGrid*>(field);
    cornerRef = fieldGrid->getGrid().getCoordinatesByCorner(ic);
  }
  else
  {
    cornerRef    = field->getCoorMinimum(useSel);
    extremesData = field->getExtremas(useSel);
  }

  for (Id icorner = 0; icorner < nval; icorner++)
  {

    // Get one corner
    Id jcorner = icorner;
    for (Id idim = 0; idim < ndim; idim++)
    {
      ic[idim] = jcorner % 2;
      jcorner /= 2;
    }
    VectorDouble corner1(ndim, 0.);
    if (isGrid)
    {
      corner1 = fieldGrid->getGrid().getCoordinatesByCorner(ic);
    }
    else
    {
      for (Id idim = 0; idim < ndim; idim++)
        corner1[idim] = (ic[idim] == 0) ? extremesData[idim][0] : extremesData[idim][1];
    }
    VH::subtractInPlace(corner1, cornerRef);

    // Rotate this corner in the Covariance Rotation system
    rot.rotateInverse(corner1, cornerRot);

    // Calculate the minimum and maximum in the Covariance rotated system
    for (Id idim = 0; idim < ndim; idim++)
    {
      if (FFFF(extendMinRot[idim]) || cornerRot[idim] < extendMinRot[idim])
        extendMinRot[idim] = cornerRot[idim];
      if (FFFF(extendMaxRot[idim]) || cornerRot[idim] > extendMaxRot[idim])
        extendMaxRot[idim] = cornerRot[idim];
    }
  }

  // Calculating the Mesh of the Grid
  VectorInt nx(ndim);
  VectorDouble dx(ndim);
  VectorDouble x0(ndim);

  double dxmin = MAXIMUM_BIG;
  for (Id idim = 0; idim < ndim; idim++)
  {
    double delta = extendMaxRot[idim] - extendMinRot[idim];
    dx[idim]     = cova.getRange(idim) / ratio;
    nx[idim]     = static_cast<Id>(ceil(delta / dx[idim])) + 2 * nbExt + 1;
    // Adapt the number of nodes if too large (compared to 'nxmax' if defined)
    if (nxmax > 0 && nx[idim] > nxmax)
    {
      nx[idim] = nxmax;
      dx[idim] = delta / (nxmax - 2 * nbExt);
    }
    if (dx[idim] < dxmin) dxmin = dx[idim];
  }

  if (cova.isNoStatForAnisotropy())
  {
    // In case of non-stationarity on the anisotropy rotation angle
    // use the minimum mesh (for internal grid)
    for (Id idim = 0; idim < ndim; idim++)
      dx[idim] = dxmin;

    for (Id idim = 0; idim < ndim; idim++)
    {
      double delta = extendMaxRot[idim] - extendMinRot[idim];
      nx[idim]     = static_cast<Id>(ceil(delta / dx[idim])) + 2 * nbExt + 1;
    }
  }

  for (Id idim = 0; idim < ndim; idim++)
  {
    extendMinRot[idim] -= nbExt * dx[idim];
    extendMaxRot[idim] += nbExt * dx[idim];
  }

  // Get the rotated Bounding Box in the initial system
  rot.rotateDirect(extendMinRot, x0);
  VH::addInPlace(x0, cornerRef);

  initFromGridByMatrix(nx, dx, x0, rot.getMatrixDirectVec(), VectorDouble(), isPolarized, verbose);
  return 0;
}

bool MeshETurbo::_deserializeAscii(std::istream& is, bool /*verbose*/)
{
  Id ndim = 0;
  VectorInt nx;
  VectorDouble dx;
  VectorDouble x0;
  VectorDouble rotmat;
  VectorInt relranks;
  Id flag_polarized = 0;
  Id nmesh_active   = 0;
  Id ngrid_active   = 0;
  Id nmesh_mask     = 0;
  Id ngrid_mask     = 0;
  Id mode           = 0;

  bool ret = true;
  ret      = ret && _recordRead<Id>(is, "Space Dimension", ndim);
  ret      = ret && _recordReadVec<Id>(is, "NX", nx, ndim);
  ret      = ret && _recordReadVec<double>(is, "DX", dx, ndim);
  ret      = ret && _recordReadVec<double>(is, "X0", x0, ndim);
  ret      = ret && _recordReadVec<double>(is, "Rotation", rotmat, ndim * ndim);
  ret      = ret && _recordRead<Id>(is, "Polarization", flag_polarized);
  ret      = ret && _recordRead<Id>(is, "Storing Mode", mode);

  if (ret)
  {
    _meshIndirect.setMode(mode);
    _gridIndirect.setMode(mode);
  }

  if (ret)
    (void)initFromGridByMatrix(nx, dx, x0, rotmat, VectorDouble(), static_cast<bool>(flag_polarized), 0);

  ret = ret && _recordRead<Id>(is, "Mesh Active Count", nmesh_active);
  ret = ret && _recordRead<Id>(is, "Mesh Masking Count", nmesh_mask);
  if (ret && nmesh_mask > 0)
  {
    ret = ret && _recordReadVec<Id>(is, "Mesh Masking", relranks, nmesh_active);
    if (ret) _meshIndirect.buildFromRankRInA(relranks, _nmeshInCompleteGrid());
  }

  ret = ret && _recordRead<Id>(is, "Grid Active Count", ngrid_active);
  ret = ret && _recordRead<Id>(is, "Mesh Masking Count", ngrid_mask);
  if (ret && ngrid_mask > 0)
  {
    ret = ret && _recordReadVec<Id>(is, "Grid Masking", relranks, ngrid_active);
    if (ret) _gridIndirect.buildFromRankRInA(relranks, _grid.getNTotal());
  }
  return ret;
}

bool MeshETurbo::_serializeAscii(std::ostream& os, bool /*verbose*/) const
{
  bool ret = true;
  ret      = ret && _recordWrite<Id>(os, "Space Dimension", getNDim());
  ret      = ret && _recordWriteVec<Id>(os, "NX", _grid.getNXs());
  ret      = ret && _recordWriteVec<double>(os, "DX", _grid.getDXs());
  ret      = ret && _recordWriteVec<double>(os, "X0", _grid.getX0s());
  ret      = ret && _recordWriteVec<double>(os, "Rotation", _grid.getRotMat());
  ret      = ret && _recordWrite<Id>(os, "Polarization", _isPolarized);
  ret      = ret && _recordWrite<Id>(os, "Storing Mode", _meshIndirect.getMode());

  // Dumping the Mesh Masking map
  Id nmesh_mask = static_cast<Id>(_meshIndirect.getRelRanks().size());
  ret            = ret && _recordWrite<Id>(os, "Mesh Active Count", getNMeshes());
  ret            = ret && _recordWrite<Id>(os, "Mesh Masking Count", nmesh_mask);
  if (nmesh_mask > 0)
    ret = ret && _recordWriteVec<Id>(os, "Mesh Masking", _meshIndirect.getRelRanks());

  // Dumping the Grid Masking map
  Id ngrid_mask = static_cast<Id>(_gridIndirect.getRelRanks().size());
  ret            = ret && _recordWrite<Id>(os, "Grid Active Count", getNApices());
  ret            = ret && _recordWrite<Id>(os, "Grid Masking Count", ngrid_mask);
  if (ngrid_mask > 0)
    ret = ret && _recordWriteVec<Id>(os, "Grid Masking", _gridIndirect.getRelRanks());

  return ret;
}

/**
 * @brief Check if a series of Meshes (included in 'meshes') are Turbo
 *
 * @param meshes
 * @return True if ALL meshes are TURBO
 */
bool isTurbo(const VectorMeshes& meshes)
{
  if (meshes.empty()) return false;
  for (Id imesh = 0, nmesh = static_cast<Id>(meshes.size()); imesh < nmesh; imesh++)
  {
    const auto* mTurbo = dynamic_cast<const MeshETurbo*>(meshes[imesh]);
    if (mTurbo == nullptr) return false;
  }
  return true;
}
#ifdef HDF5
bool MeshETurbo::_deserializeH5(H5::Group& grp, [[maybe_unused]] bool verbose)
{
  auto meshG = SerializeHDF5::getGroup(grp, "MeshETurbo");
  if (!meshG)
  {
    return false;
  }

  /* Read the grid characteristics */
  bool ret           = true;
  Id ndim           = 0;
  Id flag_polarized = 0;
  Id mode           = 0;
  VectorInt nx;
  VectorDouble dx;
  VectorDouble x0;
  VectorDouble rotmat;

  ret = ret && SerializeHDF5::readValue(*meshG, "NDim", ndim);
  ret = ret && SerializeHDF5::readVec(*meshG, "NX", nx);
  ret = ret && SerializeHDF5::readVec(*meshG, "DX", dx);
  ret = ret && SerializeHDF5::readVec(*meshG, "X0", x0);
  ret = ret && SerializeHDF5::readVec(*meshG, "Rotation", rotmat);
  ret = ret && SerializeHDF5::readValue(*meshG, "Polarization", _isPolarized);
  ret = ret && SerializeHDF5::readValue(*meshG, "Mode", mode);

  if (ret)
  {
    _meshIndirect.setMode(mode);
    _gridIndirect.setMode(mode);
    (void)initFromGridByMatrix(nx, dx, x0, rotmat, VectorDouble(), static_cast<bool>(flag_polarized), 0);
  }

  Id nmesh_active = 0;
  Id nmesh_mask   = 0;
  VectorInt nmesh_ranks;
  ret = ret && SerializeHDF5::readValue(*meshG, "MeshCountActive", nmesh_active);
  ret = ret && SerializeHDF5::readValue(*meshG, "MeshCountMask", nmesh_mask);
  if (nmesh_mask > 0)
  {
    ret = ret && SerializeHDF5::readVec(*meshG, "MeshMask", nmesh_ranks);
    if (ret) _meshIndirect.buildFromRankRInA(nmesh_ranks, _nmeshInCompleteGrid());
  }

  Id ngrid_active = 0;
  Id ngrid_mask   = 0;
  VectorInt ngrid_ranks;

  ret = ret && SerializeHDF5::readValue(*meshG, "GridCountActive", ngrid_active);
  ret = ret && SerializeHDF5::readValue(*meshG, "GridCountMask", ngrid_mask);
  if (ngrid_mask > 0)
  {
    ret = ret && SerializeHDF5::readVec(*meshG, "GridMask", ngrid_ranks);
    if (ret) _gridIndirect.buildFromRankRInA(ngrid_ranks, _grid.getNTotal());
  }

  return ret;
}

bool MeshETurbo::_serializeH5(H5::Group& grp, [[maybe_unused]] bool verbose) const
{
  auto meshG = grp.createGroup("MeshETurbo");

  bool ret = true;

  ret = ret && SerializeHDF5::writeValue(meshG, "NDim", getNDim());
  ret = ret && SerializeHDF5::writeVec(meshG, "NX", _grid.getNXs());
  ret = ret && SerializeHDF5::writeVec(meshG, "DX", _grid.getDXs());
  ret = ret && SerializeHDF5::writeVec(meshG, "X0", _grid.getX0s());
  ret = ret && SerializeHDF5::writeVec(meshG, "Rotation", _grid.getRotMat());
  ret = ret && SerializeHDF5::writeValue(meshG, "Polarization", _isPolarized);
  ret = ret && SerializeHDF5::writeValue(meshG, "Mode", _meshIndirect.getMode());

  // Dumping the Mesh Masking map
  Id nmesh_mask = static_cast<Id>(_meshIndirect.getRelRanks().size());

  ret = ret && SerializeHDF5::writeValue(meshG, "MeshCountActive", getNMeshes());
  ret = ret && SerializeHDF5::writeValue(meshG, "MeshCountMask", nmesh_mask);
  if (nmesh_mask > 0)
    ret = ret && SerializeHDF5::writeVec(meshG, "MeshMask", _meshIndirect.getRelRanks());

  // Dumping the Grid Masking map
  Id ngrid_mask = static_cast<Id>(_gridIndirect.getRelRanks().size());

  ret = ret && SerializeHDF5::writeValue(meshG, "GridCountActive", getNApices());
  ret = ret && SerializeHDF5::writeValue(meshG, "GridCountMask", ngrid_mask);
  if (ngrid_mask > 0)
    ret = ret && SerializeHDF5::writeVec(meshG, "GridMask", _gridIndirect.getRelRanks());

  return ret;
}
#endif
} // namespace gstlrn
