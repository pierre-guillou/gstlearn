/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Basic/VectorNumT.hpp"
#include "Tree/ball_algorithm.h"

namespace gstlrn
{

class GSTLEARN_EXPORT KNN
{
public:
  KNN();

  void setNNeighbors(int n_neighbors) { _n_neighbors = n_neighbors; }
  void setNSamples(int n_samples) { _n_samples = n_samples; }

  int btree_query(const t_btree& tree,
                  const MatrixT<double>& x,
                  int n_samples,
                  int n_features,
                  int n_neigh);
  int btree_query_inPlace(const t_btree& tree,
                          const MatrixT<double>& x,
                          int n_samples,
                          int n_features,
                          int n_neigh,
                          int rank,
                          VectorInt& indices,
                          VectorDouble& distances);
#ifndef SWIG
  constvectint getIndices(int rank = 0) const;
#endif // SWIG
  int getIndex(int rank = 0, int ineigh = 0) const;
  constvect getDistances(int rank = 0) const;
  double getDistance(int rank = 0, int ineigh = 0) const;

private:
  bool _query(const t_btree& tree,
              const MatrixT<double>& x,
              int n_samples,
              int n_features,
              int n_neigh);

private:
  t_nheap _heap;
  int _n_samples;
  int _n_neighbors;
};
} // namespace gstlrn
