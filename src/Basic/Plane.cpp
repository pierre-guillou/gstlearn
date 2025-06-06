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
#include "Basic/Plane.hpp"
#include "Basic/Law.hpp"
#include "Db/DbGrid.hpp"

#include <math.h>

Plane::Plane()
    : AStringable(),
      _coor(3),
      _intercept(0),
      _value(0.),
      _rndval(0.)
{
}

Plane::Plane(const Plane &m)
    : AStringable(m),
      _coor(m._coor),
      _intercept(m._intercept),
      _value(m._value),
      _rndval(m._rndval)
{
}

Plane& Plane::operator=(const Plane &m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
    _coor = m._coor;
    _intercept = m._intercept;
    _value = m._value;
    _rndval = m._rndval;
  }
  return *this;
}

Plane::~Plane()
{

}

String Plane::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  return sstr.str();
}

/****************************************************************************/
/*!
 **  Generate the Poisson planes that cover the grid
 **
 ** \param[in]  dbgrid   Db corresponding to the target grid
 ** \param[in]  np       Number of planes
 **
 ** \remarks  The array 'planes' contains successively a,b,c,d such that
 ** \remarks  ax + by + cz + d = 0
 ** \remarks  The valuation of each line is assigned a uniform value [0,1]
 **
 *****************************************************************************/
std::vector<Plane> Plane::poissonPlanesGenerate(DbGrid *dbgrid, int np)
{
  double ap[3];

  VectorDouble center = dbgrid->getCenters();
  center.resize(3,0.);
  double diagonal = dbgrid->getExtensionDiagonal();
  std::vector<Plane> planes;
  planes.resize(np);

  /* Loop on the planes to be generated */

  for (int ip = 0; ip < np; ip++)
  {
    double d0 = diagonal * law_uniform(-1., 1.) / 2.;
    double u = 0.;
    for (int idim = 0; idim < 3; idim++)
    {
      ap[idim] = law_gaussian();
      u += ap[idim] * ap[idim];
    }
    u = sqrt(u);
    for (int idim = 0; idim < 3; idim++)
    {
      ap[idim] /= u;
    }
    // Check position of the Center (in its OWN space dimension)
    for (int idim = 0; idim < (int) center.size(); idim++)
    {
      d0 -= ap[idim] * center[idim];
    }
    if (d0 < 0)
    {
      for (int idim = 0; idim < 3; idim++)
        ap[idim] = -ap[idim];
      d0 = -d0;
    }

    /* Storing the plane */

    for (int idim = 0; idim < 3; idim++)
      planes[ip].setCoor(idim, ap[idim]);
    planes[ip].setIntercept(d0);
    planes[ip].setRndval(law_uniform(0., 1.));
  }

  return planes;
}

double Plane::getCoor(int idim) const
{
  if (idim < (int) _coor.size())
    return _coor[idim];
  return 0.;
}

void Plane::setCoor(int idim, double value)
{
  if (idim < (int) _coor.size())
    _coor[idim] = value;
}
