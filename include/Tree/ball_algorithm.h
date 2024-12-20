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

#include <Basic/VectorNumT.hpp>
#include <Space/SpacePoint.hpp>

struct t_btree;

struct t_nheap
{
  VectorVectorDouble distances;
  VectorVectorInt indices;
  int n_pts{};
  int n_nbrs{};

  t_nheap() = default;
  t_nheap(int n_pts, int n_nbrs);

  double largest(const int row) const;
  int    push(int row, double val, int i_val);
  void   sort();
  void   load(t_btree &b, const std::vector<SpacePoint> &x);
};

struct t_nodedata
{
  int idx_start{};
  int idx_end{};
  double radius{};
  bool is_leaf{};
};

struct t_btree
{
  std::vector<SpacePoint> data;
  VectorInt idx_array;
  std::vector<t_nodedata> node_data;
  std::vector<SpacePoint> node_bounds;

  int n_samples{};
  int n_features{};

  int leaf_size{40};
  int n_levels{};
  int n_nodes{};

  t_btree() = default;
  t_btree(const VectorVectorDouble &data,
          int n_samples,
          int n_features,
          int leaf_size);

  void display(int level = -1) const;
  int init_node(int i_node, int idx_start, int idx_end);
  void recursive_build(int i_node, int idx_start, int idx_end);
  double min_dist(int i_node, const SpacePoint &pt);
  int query_depth_first(int i_node, const SpacePoint &pt, int i_pt, t_nheap &heap, double dist);
};
