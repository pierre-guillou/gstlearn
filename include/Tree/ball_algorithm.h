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

#include "gstlearn_export.hpp"
#include <Basic/VectorNumT.hpp>

struct t_nheap
{
  VectorVectorDouble distances;
  VectorVectorInt indices;
  int n_pts{};
  int n_nbrs{};

  t_nheap() = default;
  t_nheap(int n_pts, int n_nbrs);
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
  VectorVectorDouble data;
  VectorInt idx_array;
  std::vector<t_nodedata> node_data;
  VectorVectorDouble node_bounds;

  int n_samples{};
  int n_features{};

  int leaf_size{40};
  int n_levels{};
  int n_nodes{};

  t_btree() = default;
  t_btree(const VectorVectorDouble &data,
          int n_samples,
          int n_features,
          double (*dist_function)(const double* x1,
                                  const double* x2,
                                  int size),
          int leaf_size,
          int default_distance_function);

  void display(int level = -1) const;
};

/*
** metrics.c
*/

double manhattan_distance(const double* x1, const double* x2, int size);
double euclidean_distance(const double* x1, const double* x2, int size);

/*
** neighbors_heap.c
*/

double   nheap_largest(const t_nheap &h, const int row);
int		 nheap_push(t_nheap &h, int row, double val, int i_val);
void     nheap_sort(t_nheap &h);
void     nheap_load(t_nheap &heap, t_btree &b, const VectorVectorDouble &x);
double   min_dist(t_btree &tree, int i_node, const double *pt);
int      query_depth_first(t_btree &b, int i_node, const double *pt, int i_pt, t_nheap &heap, double dist);
