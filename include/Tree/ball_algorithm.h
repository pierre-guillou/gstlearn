/*
                                      ball

Original Author: Eung Bum Lee
Website: https://42.fr
License: MIT
*/

/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   ball.h                                             :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: elee <elee@student.42.us.org>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2017/06/28 10:55:27 by elee              #+#    #+#             */
/*   Updated: 2017/06/28 20:56:18 by elee             ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

/*
Modified by MINES Paris / ARMINES (2023)
Authors: gstlearn Team
Website: https://gstlearn.org
License: BSD 3-clause
*/

#pragma once

#include "Basic/VectorNumT.hpp"
#include "Matrix/MatrixT.hpp"

namespace gstlrn
{

struct t_btree;

struct t_nheap
{
  MatrixT<double> distances;
  MatrixT<int> indices;
  int n_pts;
  int n_nbrs;

  t_nheap() = default;
  void resize(int n_pts, int n_nbrs);

  /*
  ** neighbors_heap.c
  */
  double largest(int row) const;
  int push(int row, double val, int i_val);
  void sort();
  void load(const t_btree& b, const MatrixT<double>& x);
};

struct t_nodedata
{
  int idx_start;
  int idx_end;
  double radius;
  bool is_leaf;
};

struct t_btree
{
  MatrixT<double> data;
  VectorBool accept;
  VectorInt idx_array;
  std::vector<t_nodedata> node_data;
  MatrixT<double> node_bounds;

  int n_samples;
  int n_features;

  int leaf_size;
  int n_levels;
  int n_nodes;
  int default_distance_function {1};

  t_btree() = default;
  t_btree(MatrixT<double>&& data,
          int n_samples,
          int n_features,
          bool has_constraints,
          int leaf_size,
          int default_distance_function);

  void display(int level = -1) const;
  int init_node(int i_node, int idx_start, int idx_end);
  void recursive_build(int i_node, int idx_start, int idx_end);
  double min_dist(int i_node, const constvect pt) const;
  int query_depth_first(int i_node, const constvect pt, int i_pt, t_nheap& heap, double dist) const;
};

/*
** metrics.c
*/
double manhattan_distance(const double* x1, const double* x2, int n_features);
double euclidean_distance(const double* x1, const double* x2, int n_features);

} // namespace gstlrn
