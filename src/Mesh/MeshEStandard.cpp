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
#include "Matrix/MatrixDense.hpp"
#include "LinearOp/ProjMatrix.hpp"
#include "Mesh/MeshEStandard.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Basic/AException.hpp"

MeshEStandard::MeshEStandard()
  : AMesh()
  , _apices()
  , _meshes()
{
}

MeshEStandard::MeshEStandard(const MeshEStandard &m) 
  : AMesh(m)
{
  _recopy(m);
}

MeshEStandard& MeshEStandard::operator= (const MeshEStandard &m)
{
  _recopy(m);
  return *this;
}

MeshEStandard::~MeshEStandard()
{
  _deallocate();
}

/****************************************************************************/
/*!
** Returns the number of Apices
**
** \returns Number of apices
**
*****************************************************************************/
int MeshEStandard::getNApices() const 
{
  return _apices.getNRows();
}

/****************************************************************************/
/*!
** Returns the number of Meshes
**
** \returns Number of meshes
**
*****************************************************************************/
int MeshEStandard::getNMeshes() const
{
  return (static_cast<int> (_meshes.size()) / getNApexPerMesh());
}

/****************************************************************************/
/*!
** Returns the size of the Mesh 'imesh'
**
** \returns mesh dimension
**
** \param[in]  imesh    Rank of the Mesh (from 0 to _nMeshes-1))
**
*****************************************************************************/
double MeshEStandard::getMeshSize(int imesh) const
{
  VectorVectorDouble corners = getCoordinatesPerMesh(imesh);
  return _getMeshUnit(corners);
}

int MeshEStandard::resetFromTurbo(const MeshETurbo& turbo, bool verbose)
{
  int ndim     = turbo.getNDim();
  int napices  = turbo.getNApices();
  int nmeshes  = turbo.getNMeshes();
  int npermesh = turbo.getNApexPerMesh();

  // Dimension the members

  _apices = MatrixDense(napices, ndim);
  _meshes = MatrixInt(nmeshes, npermesh);

  // Load the apices;
  VectorDouble local(ndim);
  for (int ip = 0; ip < napices; ip++)
  {
    turbo.getApexCoordinatesInPlace(ip, local);
    for (int idim = 0; idim < ndim; idim++)
      _apices.setValue(ip, idim, local[idim]);
  }

  // Load the meshes
  for (int imesh = 0; imesh < nmeshes; imesh++)
    for (int rank = 0; rank < npermesh; rank++)
      _meshes.setValue(imesh,  rank, turbo.getApex(imesh, rank));

  // Check consistency

  _checkConsistency();

  // Define and store the Bounding Box extension

  _defineBoundingBox();

  // Optional printout

  if (verbose) messageFlush(toString());

  return 0;
}

/****************************************************************************/
/*!
** Create the meshing (from mesh information)
**
** \param[in]  apices          Matrix of Apices
** \param[in]  meshes          Array of mesh indices
** \param[in]  verbose         Verbose flag
**
*****************************************************************************/
int MeshEStandard::reset(const MatrixDense& apices,
                         const MatrixInt& meshes,
                         bool verbose)
{
  int ndim = apices.getNCols();
  _setNDim(ndim);

  // Core allocation

  _apices = apices;
  _meshes = meshes;

  // Check consistency

  _checkConsistency();

  // Define and store the Bounding Box extension

  _defineBoundingBox();

  // Optional printout

  if (verbose) messageFlush(toString());

  return(0);
}

/****************************************************************************/
/*!
** Create the meshing (from mesh information)
**
** \param[in]  ndim            Space Dimension
** \param[in]  napexpermesh    Number of apices per mesh
** \param[in]  apices          Vector of Apex information
** \param[in]  meshes          Vector of mesh indices
** \param[in]  byCol           true for Column major; false for Row Major
** \param[in]  verbose         Verbose flag
**
** \remark The argument 'byCol' concerns 'apices' and 'meshes'
**
*****************************************************************************/
int MeshEStandard::reset(int ndim,
                         int napexpermesh,
                         const VectorDouble &apices,
                         const VectorInt &meshes,
                         bool byCol,
                         bool verbose)
{
  _setNDim(ndim);
  int npoints = static_cast<int> (apices.size()) / ndim;
  int nmeshes = static_cast<int> (meshes.size()) / napexpermesh;

  // Core allocation

  _apices.reset(npoints,ndim);
  _apices.setValues(apices, byCol);
  _meshes.reset(nmeshes,napexpermesh);
  _meshes.setValues(meshes, byCol);

  // Check consistency

  _checkConsistency();

  // Define and store the Bounding Box extension

  _defineBoundingBox();

  // Optional printout

  if (verbose) messageFlush(toString());

  return(0);
}

/****************************************************************************/
/*!
** Returns the rank of the Apex 'rank' of the Mesh 'imesh'
**
** \returns The rank of the target apex
**
** \param[in]  imesh    Rank of Mesh (from 0 to _nMeshes-1))
** \param[in]  rank     Rank of Apex within a Mesh (from 0 to _nApexPerMesh)
**
*****************************************************************************/
int MeshEStandard::getApex(int imesh, int rank) const
{
  return _meshes.getValue(imesh,rank);
}

void MeshEStandard::_setApex(int imesh, int rank, int value)
{
  _meshes.setValue(imesh,rank,value);
}

/****************************************************************************/
/*!
** Returns the coordinate 'idim' of the Apex 'rank' of the Mesh 'imesh'
**
** \returns The coordinate of the target apex
**
** \param[in]  imesh    Rank of the Mesh (from 0 to _nMeshes-1))
** \param[in]  rank     Rank of Apex within a Mesh (from 0 to _nApexPerMesh-1)
** \param[in]  idim     Rank of the coordinate (from 0 to _nDim-1)
**
*****************************************************************************/
double MeshEStandard::getCoor(int imesh,
                              int rank,
                              int idim) const
{
  return _apices(getApex(imesh,rank),idim);
}

/****************************************************************************/
/*!
** Print the contents of the meshing
**
** \param[in] strfmt    Format for printout
**
*****************************************************************************/
String MeshEStandard::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;
  sstr << toTitle(0,"Standard Meshing");
  sstr << AMesh::toString(strfmt);
  return sstr.str();
}

/**
 * Create a MeshEStandard by loading the contents of a Neutral File
 *
 * @param neutralFilename Name of the Neutral File (MeshEStandard format)
 * @param verbose         Verbose
 */
MeshEStandard* MeshEStandard::createFromNF(const String& neutralFilename, bool verbose)
{
  MeshEStandard* mesh = nullptr;
  std::ifstream is;
  mesh = new MeshEStandard;
  bool success = false;
  if (mesh->_fileOpenRead(neutralFilename, is, verbose))
  {
    success =  mesh->deserialize(is, verbose);
  }
  if (! success)
  {
    delete mesh;
    mesh = nullptr;
  }
  return mesh;
}

MeshEStandard* MeshEStandard::createFromExternal(const MatrixDense &apices,
                                                 const MatrixInt &meshes,
                                                 bool verbose)
{
  MeshEStandard* mesh = new MeshEStandard;
  mesh->reset(apices, meshes, verbose);
  return mesh;
}

/**
 * Returns the list of mesh vertex information
 * This List is organized as a single Vector of Double
 * It is dimensioned to Nrows=getNApices() and Ncols=getNDim()
 * @param byCol true if the values must be sorted by column
 * @return
 */
VectorDouble MeshEStandard::getPointList(bool byCol) const
{
  VectorDouble list;

  if (byCol)
  {
    for (int idim = 0; idim < getNDim(); idim++)
      for (int irow = 0; irow < getNApices(); irow++)
      {
        list.push_back(getApexCoor(irow, idim));
      }
  }
  else
  {
    for (int irow = 0; irow < getNApices(); irow++)
      for (int idim = 0; idim < getNDim(); idim++)
      {
        list.push_back(getApexCoor(irow, idim));
      }
  }
  return list;
}

double MeshEStandard::getApexCoor(int i, int idim) const
{
  return _apices(i, idim);
}

/****************************************************************************/
/*!
**  Check if a point, defined by its coordinates, belong to a Mesh
**
** \return true if the point belongs to the Mesh; false otherwise
**
** \param[in]  coor      Array of target coordinates
** \param[in]  imesh     Mesh Index
** \param[in]  meshsize  Dimension of the mesh
**
** \param[out] weights   Array of barycentric weights (Dim: NApexPerMesh)
** 
** \remarks The argument 'meshsize' is used to speed the calculations
**
*****************************************************************************/
bool MeshEStandard::_coorInMesh(const VectorDouble& coor,
                                int imesh,
                                double meshsize,
                                VectorDouble& weights) const
{
  VectorVectorDouble corners = getCoordinatesPerMesh(imesh);
  return _weightsInMesh(coor, corners, meshsize, weights);
}

void MeshEStandard::_deallocate()
{

}

int MeshEStandard::_recopy(const MeshEStandard &m)
{
  _deallocate();

  _apices = m._apices;
  _meshes = m._meshes;
  AMesh::_recopy(m);
  return(0);
}

/**
 * This function checks the consistency between the number of points
 * and the vertices indication
 */
void MeshEStandard::_checkConsistency() const
{
  for (int imesh = 0; imesh < getNMeshes(); imesh++)
    for (int ic = 0; ic < getNApexPerMesh(); ic++)
    {
      int apex = getApex(imesh, ic);
      if (apex < 0 || apex >= getNApices())
      {
        my_throw("Mesh indices are not compatible with the Points");
      }
    }
}

bool MeshEStandard::_deserialize(std::istream& is, bool verbose)
{
  DECLARE_UNUSED(verbose);
  int ndim = 0;
  int napices = 0;
  int napexpermesh = 0;
  int nmeshes = 0;

  bool ret = true;
  ret = ret && _recordRead<int>(is, "Space Dimension", ndim);
  ret = ret && _recordRead<int>(is, "Napices", napices);
  ret = ret && _recordRead<int>(is, "Number of Apices per Mesh", napexpermesh);
  ret = ret && _recordRead<int>(is, "Number of Meshes", nmeshes);

  if (ret)
  {
    VectorDouble apices_local;
    ret = ret && _recordReadVec<double>(is, "Apices", apices_local, napices * ndim);
    _apices = MatrixDense(napices, ndim);
    _apices.setValues(apices_local);
  }

  if (ret)
  {
    VectorInt meshes_local;
    ret = ret && _recordReadVec<int>(is, "Meshes", meshes_local, nmeshes * napexpermesh);
    _meshes = MatrixInt(nmeshes, napexpermesh);
    _meshes.setValues(meshes_local);
  }
  return ret;
}

bool MeshEStandard::_serialize(std::ostream& os, bool /*verbose*/) const
{
  bool ret = true;

  ret = ret && _recordWrite<int>(os, "Space Dimension", getNDim());
  ret = ret && _recordWrite<int>(os, "Napices", getNApices());
  ret = ret && _recordWrite<int>(os, "Number of Apices per Mesh", getNApexPerMesh());
  ret = ret && _recordWrite<int>(os, "Number of Meshes", getNMeshes());
  ret = ret && _recordWriteVec<double>(os, "Apices", _apices.getValues());
  ret = ret && _recordWriteVec<int>(os, "Meshes", _meshes.getValues());
  return ret;
}

/**
 * Validate Meshing, in particular when transiting from Old Meshing to New one
 *
 * @remark For safety and considering the rank of Apices stored in 'meshes',
 * @remark the minimum value is considered
 * @remark If not equal to 0 or 1, a fatal error is issued
 * @remark If equal to 1 (old style), values are decreased by 1.
 */
void MeshEStandard::_validate()
{
  // For safety, the minimum value of the array meshes is considered
  // If equal to , the mesh numbering must be modified in order to start from 0
  // This code is temporary and there to cope with different numberings

  int nmeshes = getNMeshes();
  int ncorner = getNApexPerMesh();
  int shift_min = 1000;
  for (int imesh = 0; imesh < nmeshes; imesh++)
    for (int ic = 0; ic < ncorner; ic++)
    {
      int ipos = getApex(imesh, ic);
      if (ipos < shift_min) shift_min = ipos;
    }

  if (shift_min != 0 && shift_min != 1)
    my_throw("Wrong minimum shift: it should be 0 or 1");

  if (shift_min == 1)
  {
    for (int imesh = 0; imesh < nmeshes; imesh++)
      for (int ic = 0; ic < ncorner; ic++)
        _setApex(imesh, ic, getApex(imesh, ic) - shift_min);
  }
}

void MeshEStandard::_defineBoundingBox(void)
{
  VectorDouble extendmin, extendmax;
  double coor,mini,maxi;
  int ndim = getNDim();

  // Initializations
  extendmin.resize(ndim);
  extendmax.resize(ndim);

  // Loop on the Space dimensions
  for (int idim=0; idim<ndim; idim++)
  {
    mini =  1.e30;
    maxi = -1.e30;

    // Loop on the apices
    for (int i=0; i<getNApices(); i++)
    {
      coor = getApexCoor(i,idim);
      if (coor < mini) mini = coor;
      if (coor > maxi) maxi = coor;
    }
    extendmin[idim] = mini;
    extendmax[idim] = maxi;
  }

  // Store the Bounding Box extension
  (void) _setExtend(extendmin,extendmax);
}

