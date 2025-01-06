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
#include "Tree/Ball.hpp"
#include "Basic/VectorNumT.hpp"
#include "Tree/ball_algorithm.h"
#include "Db/Db.hpp"
#include "Space/SpacePoint.hpp"

Ball::Ball(const VectorVectorDouble& data,
           int leaf_size)
{
  int n_samples     = data.empty() ? 0 : (int)data[0].size();
  int n_features    = (int)data.size();
  _tree = t_btree(data, n_samples, n_features, leaf_size);
}

Ball::Ball(const Db* db,
           int leaf_size,
           bool useSel)
{
  VectorVectorDouble data = db->getAllCoordinates(useSel);
  int n_samples           = (int)data[0].size();
  int n_features          = (int)data.size();
  _tree = t_btree(data, n_samples, n_features, leaf_size);
}

void Ball::init(const Db* db,
                int leaf_size,
                bool useSel)
{
  VectorVectorDouble data = db->getAllCoordinates(useSel);
  int n_samples           = (int)data[0].size();
  int n_features          = (int)data.size();
    _tree = t_btree(data, n_samples, n_features, leaf_size);
}

Ball::~Ball()
{
}

KNN Ball::query(const std::vector<SpacePoint> &test, int n_samples, int n_features, int n_neighbors)
{
  KNN knn;
  (void) knn.btree_query(_tree, test, n_samples, n_features, n_neighbors);
  return knn;
}

KNN Ball::queryAsVVD(const std::vector<SpacePoint>& test, int n_neighbors)
{
  KNN knn;
  if (test.empty()) return knn;
  int n_samples = (int) test[0].getNDim();
  int n_features = (int) test.size();
  (void) knn.btree_query(_tree, test, n_samples, n_features, n_neighbors);
  return knn;
}

KNN Ball::queryOneAsVDFromSP(const SpacePoint& Pt, int n_neighbors)
{
  KNN knn;
  int n_features = Pt.getNDim();
  const VectorDouble internal = {Pt.getCoords().begin(), Pt.getCoords().end()};
  (void)knn.btree_query(_tree, {internal}, 1, n_features,
                        n_neighbors);
  return knn;
}

VectorInt Ball::getIndices(const SpacePoint& Pt, int n_neighbors)
{
  KNN knn = queryOneAsVDFromSP(Pt, n_neighbors);
  return knn.getIndices(0);
}

int Ball::queryClosest(const VectorDouble& test)
{
  KNN knn;
  int n_features = (int) test.size();
  if (knn.btree_query(_tree, {test}, 1, n_features, 1)) return ITEST;
  return knn.getIndex(0, 0);
}

int Ball::queryOneInPlace(const VectorDouble& test,
                          int n_neighbors,
                          VectorInt& indices,
                          VectorDouble& distances,
                          int rank)
{
  KNN knn;
  int n_features         = (int)test.size();
  return knn.btree_query_inPlace(_tree, {test}, 1,
                                 n_features, n_neighbors, rank, indices,
                                 distances);
}

void Ball::display(int level) const
{
  _tree.display(level);
}
