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
#include "Tree/KNN.hpp"
#include "Basic/VectorNumT.hpp"

# include <stdlib.h>
# include <limits.h>
# include <math.h>
# include <float.h>
# include <sys/stat.h>

# define TRUE 1
# define FALSE 0

struct t_nheap
{
  double **distances;
  int **indices;
  int n_pts;
  int n_nbrs;
};

typedef struct
{
  int idx_start;
  int idx_end;
  int is_leaf;
  double radius;
} t_nodedata;

struct t_btree
{
  double **data;
  bool *accept;
  int *idx_array;
  t_nodedata *node_data;
  double ***node_bounds;

  int n_samples;
  int n_features;

  int leaf_size;
  int n_levels;
  int n_nodes;
};

/*
** ball.c
*/
GSTLEARN_EXPORT double** copy_double_arrAsVVD(const VectorVectorDouble& arr);
GSTLEARN_EXPORT double** copy_double_arr(const double** arr, int row, int col);
GSTLEARN_EXPORT bool* init_bool_arr(int col, bool status);
GSTLEARN_EXPORT VectorVectorDouble copy_double_toVVD(const double** arr,
                                                     int row,
                                                     int col);
GSTLEARN_EXPORT VectorVectorInt copy_int_toVVI(const int** arr,
                                               int row,
                                               int col);
GSTLEARN_EXPORT int** copy_int_arr(const int** arr, int row, int col);

GSTLEARN_EXPORT t_btree* btree_init(const double** data,
                                    int n_samples,
                                    int n_features,
                                    bool has_constraints,
                                    double (*dist_function)(const double* x1,
                                                            const double* x2,
                                                            int n_features),
                                    int leaf_size,
                                    int default_distance_function);
GSTLEARN_EXPORT void free_2d_double(double **arr, int row);
GSTLEARN_EXPORT void free_2d_int(int **arr, int row);
GSTLEARN_EXPORT void free_tree(t_btree *tree);
GSTLEARN_EXPORT void finalize_tree(t_btree *tree); // to be suppressed when bug is corrected
GSTLEARN_EXPORT void btree_display(const t_btree* tree, int level = -1);

/*
** metrics.c
*/
double manhattan_distance(const double* x1, const double* x2, int n_features);
double euclidean_distance(const double* x1, const double* x2, int n_features);

/*
** neighbors_heap.c
*/
t_nheap* nheap_init(int n_pts, int n_nbrs);
t_nheap* nheap_free(t_nheap* heap);
double   nheap_largest(t_nheap* h, int row);
int		   nheap_push(t_nheap *h, int row, double val, int i_val);
void     nheap_sort(t_nheap *h);
void     nheap_load(t_nheap *heap, t_btree *b, const double **x);
double   min_dist(const t_btree *tree, int i_node, const double *pt);
int      query_depth_first(const t_btree *b, int i_node, const double *pt, int i_pt, t_nheap *heap, double dist);
