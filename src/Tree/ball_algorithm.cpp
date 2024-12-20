/*
                                      ball

Original Author: Eung Bum Lee
Website: https://42.fr
License: MIT
*/

/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   ball.c                                             :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: elee <elee@student.42.us.org>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2017/06/28 10:45:02 by elee              #+#    #+#             */
/*   Updated: 2017/06/28 20:56:01 by elee             ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

/*
Modified by MINES Paris / ARMINES (2023)
Authors: gstlearn Team
Website: https://gstlearn.org
License: BSD 3-clause
*/

#include "Tree/ball_algorithm.h"
#include "Tree/KNN.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/VectorHelper.hpp"
#include "Space/SpacePoint.hpp"

static double (*st_distance_function)(const double*, const double*, int) = euclidean_distance;

int t_btree::init_node(int i_node, int idx_start, int idx_end)
{
  int n_features = this->n_features;
  int n_points = idx_end - idx_start;
  double* centroid = this->node_bounds[i_node].data();

  for (int j = 0; j < n_features; j++)
    centroid[j] = 0.0;

  for (int i = idx_start; i < idx_end; i++)
    for (int j = 0; j < n_features; j++)
      centroid[j] += this->data[this->idx_array[i]][j];

  for (int j = 0; j < n_features; j++)
    centroid[j] /= n_points;

  double radius = 0.0;
  for (int i = idx_start; i < idx_end; i++)
    radius = fmax(radius, st_distance_function(centroid, this->data[this->idx_array[i]].data(), n_features));

  this->node_data[i_node].radius = radius;
  this->node_data[i_node].idx_start = idx_start;
  this->node_data[i_node].idx_end = idx_end;
  return (0);
}

int find_node_split_dim(VectorVectorDouble &data, const std::vector<int> &node_indices, int n_features, int n_points)
{
	double	min_val, max_val, val, spread;

	int j_max = 0;
	double max_spread = 0;
	for (int j = 0; j < n_features; j++)
	{
		max_val = data[node_indices[0]][j];
		min_val = max_val;
		for (int i = 1; i < n_points; i++)
		{
			val = data[node_indices[i]][j];
			max_val = fmax(max_val, val);
			min_val = fmin(min_val, val);
		}
		spread = max_val - min_val;
		if (spread > max_spread)
		{
			max_spread = spread;
			j_max = j;
		}
	}
	return (j_max);
}

int partition_node_indices(VectorVectorDouble &data, int *node_indices, int split_dim, int n_points, int split_index)
{
  int   midindex;
  double  d1, d2;

  int left = 0;
  int right = n_points - 1;

  while (true)
  {
    midindex = left;
    for (int i = left; i < right; i++)
    {
      d1 = data[node_indices[i]][split_dim];
      d2 = data[node_indices[right]][split_dim];
      if (d1 < d2)
      {
        std::swap(node_indices[i], node_indices[midindex]);
        midindex++;
      }
    }
    std::swap(node_indices[midindex], node_indices[right]);
    if (midindex == split_index)
      break ;
    if (midindex < split_index)
      left = midindex + 1;
    else
      right = midindex - 1;
  }

  return (0);
}

void t_btree::recursive_build(int i_node, int idx_start, int idx_end)
{
	int	imax;
	int n_features = this->n_features;
	int n_points = idx_end - idx_start;
	int n_mid = n_points / 2;

	//initialize the node data
	this->init_node(i_node, idx_start, idx_end);

	if (2 * i_node + 1 >= this->n_nodes)
	{
		this->node_data[i_node].is_leaf = true;
		if (idx_end - idx_start > 2 * this->leaf_size)
			messerr("Memory layout is flawed: not enough nodes allocated");
	}
	else if (idx_end - idx_start < 2)
	{
		messerr("Memory layout is flawed: too many nodes allocated");
		this->node_data[i_node].is_leaf = true;
	}
	else
	{
		this->node_data[i_node].is_leaf = false;
		imax = find_node_split_dim(this->data, this->idx_array, n_features, n_points);
		partition_node_indices(this->data, &this->idx_array[idx_start], imax, n_points, n_mid);
		this->recursive_build(2 * i_node + 1, idx_start, idx_start + n_mid);
		this->recursive_build(2 * i_node + 2, idx_start + n_mid, idx_end);
	}
}

void define_dist_function(double (*dist_function)(const double* x1,
                                                  const double* x2,
                                                  int size),
                          int default_distance_function)
{
  if (dist_function != nullptr)
  {
    st_distance_function = dist_function;
  }
  else
  {
    if (default_distance_function == 1)
      st_distance_function = euclidean_distance;
    if (default_distance_function == 2)
      st_distance_function = manhattan_distance;
  }
}

t_btree::t_btree(const VectorVectorDouble &data,
                 int n_samples,
                 int n_features,
                 double (*dist_function)(const double* x1,
                                         const double* x2,
                                         int size),
                 int leaf_size,
                 int default_distance_function)
{
    this->data = data;
	this->leaf_size = leaf_size;
	
	if (leaf_size < 1)
	{
		messerr("leaf_size must be greater than or equal to 1\n");
		return;
	}

  // Define the relevant distance function
  define_dist_function(dist_function, default_distance_function);

	this->n_samples = n_samples;
	this->n_features = n_features;

	this->n_levels = log2(fmax(1, (this->n_samples - 1) / this->leaf_size)) + 1;
	this->n_nodes = pow(2.0, this->n_levels) - 1;

	this->idx_array.resize(this->n_samples);
	for (int i = 0; i < this->n_samples; i++)
		this->idx_array[i] = i;
	this->node_data.resize(this->n_nodes);
	this->node_bounds.resize(this->n_nodes);
	for (int i = 0; i < this->n_nodes; i++)
	{
      this->node_bounds[i].resize(this->n_features);
	}
	this->recursive_build(0, 0, this->n_samples);
}

double t_btree::min_dist(int i_node, const double *pt)
{
  double dist_pt = st_distance_function(pt, this->node_bounds[i_node].data(), this->n_features);
  return (fmax(0.0, dist_pt - this->node_data[i_node].radius));
}

int t_btree::query_depth_first(int i_node, const double *pt, int i_pt, t_nheap &heap, double dist)
{
  t_nodedata node_info = this->node_data[i_node];
  double dist_pt, dist1, dist2;
  int i1, i2;

  // case 1: query point is outside node radius: trim it from the query
  if (dist > heap.largest(i_pt))
  {
    ;
  }
  // case 2: this is a leaf node. Update set of nearby points
  else if (node_info.is_leaf)
  {
    for (int i = node_info.idx_start; i < node_info.idx_end; i++)
    {
      dist_pt = st_distance_function(pt, this->data[this->idx_array[i]].data(), this->n_features);
      if (dist_pt < heap.largest(i_pt))
        heap.push(i_pt, dist_pt, this->idx_array[i]);
    }
  }
  // case 3: Node is not a leaf, Recursively query sub-nodes starting with the
  // closest
  else
  {
    i1    = 2 * i_node + 1;
    i2    = i1 + 1;
    dist1 = this->min_dist(i1, pt); // implement min_rdist
    dist2 = this->min_dist(i2, pt);
    if (dist1 <= dist2)
    {
      this->query_depth_first(i1, pt, i_pt, heap, dist1);
      this->query_depth_first(i2, pt, i_pt, heap, dist2);
    }
    else
    {
      this->query_depth_first(i2, pt, i_pt, heap, dist2);
      this->query_depth_first(i1, pt, i_pt, heap, dist1);
    }
  }
  return (0);
}

void t_btree::display(int level) const
{
  message("- Number of samples = %d\n", this->n_samples);
  message("- Number of Features = %d\n", this->n_features);
  message("- Number of levels = %d\n", this->n_levels);
  message("- Number of nodes = %d\n", this->n_nodes);
  message("- Size of leaf = %d\n", this->leaf_size);
  if (level < 0) return;

  // Loop on the nodes

  for (int i_node = 0; i_node < this->n_nodes; i_node++)
  {
    const t_nodedata* info = &this->node_data[i_node];
    const auto &centroid = this->node_bounds[i_node];

    message("Node #%3d/%3d - Indices [%5d; %5d[ - Radius = %lf",
            i_node, this->n_nodes, info->idx_start, info->idx_end, info->radius);
    if (info->is_leaf)
      message(" - Terminal Leaf\n");
    else
      message("\n");

    if (level > 0)
    {
      VH::display("Centroid = ", centroid, 0);

      if (level > 1 && info->is_leaf)
      {
        message("  Sample indices = ");
        for (int is = info->idx_start; is < info->idx_end; is++)
          message(" %d", this->idx_array[is]);
        message("\n");
      }
    }
  }
}

/**
 * Returns the Manhattan distance between two points
 * @param x1 Vector of coordinates for the first point
 * @param x2 Vector of coordinates for the second point
 * @param size Number of coordinates
 * @return
 */
double manhattan_distance(const double* x1, const double* x2, int size)
{
  double delta;
  double d1 = 0.;
  for (int i = 0; i < size; i++)
  {
    delta = fabs(x1[i] - x2[i]);
    d1 += delta;
  }
  return (d1);
}

/**
 * Returns the Standard distance between two points
 * @param x1 Vector of coordinates for the first point
 * @param x2 Vector of coordinates for the second point
 * @param size Number of coordinates
 * @return
 */
double euclidean_distance(const double* x1, const double* x2, int size)
{
  SpacePoint p1;
  SpacePoint p2;
  p1.setCoords(x1, size);
  p2.setCoords(x2, size);
  return p1.getDistance(p2);
}
