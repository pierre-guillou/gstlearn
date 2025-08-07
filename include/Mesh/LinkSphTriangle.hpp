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
#pragma once

#include "Basic/VectorNumT.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

namespace gstlrn
{
class Db;

class SphTriangle
{
public:
  int n_nodes; /* Number of nodes */
  int sph_size; /* Size of arrays sph_list and sph_lptr */
  VectorDouble sph_x; /* Array of X-coordinates for nodes */
  VectorDouble sph_y; /* Array of Y-coordinates for nodes */
  VectorDouble sph_z; /* Array of Z-coordinates for nodes */
  std::vector<int> sph_list; /* Set of nodal indexes */
  std::vector<int> sph_lptr; /* Set of pointers (sph_list indexes) */
  std::vector<int> sph_lend; /* Set of pointers to adjacency lists */
};



GSTLEARN_EXPORT void meshes_2D_sph_init(SphTriangle* t);
GSTLEARN_EXPORT void meshes_2D_sph_free(SphTriangle* t, Id mode);
GSTLEARN_EXPORT Id meshes_2D_sph_from_db(Db* db, SphTriangle* t);
GSTLEARN_EXPORT Id meshes_2D_sph_from_points(Id nech,
                                              double* x,
                                              double* y,
                                              SphTriangle* t);
GSTLEARN_EXPORT Id meshes_2D_sph_from_auxiliary(const String& triswitch,
                                                 SphTriangle* t);
GSTLEARN_EXPORT void meshes_2D_sph_print(SphTriangle* t, Id brief);
GSTLEARN_EXPORT Id meshes_2D_sph_create(Id verbose, SphTriangle* t);

} // namespace gstlrn
