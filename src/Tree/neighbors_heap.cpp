/*
                                      ball

Original Author: Eung Bum Lee
Website: https://42.fr
License: MIT
*/

/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   neighbors_heap.c                                   :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: elee <elee@student.42.us.org>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2017/06/28 10:45:06 by elee              #+#    #+#             */
/*   Updated: 2017/06/28 20:17:58 by elee             ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

/*
Modified by MINES Paris / ARMINES (2023)
Authors: gstlearn Team
Website: https://gstlearn.org
License: BSD 3-clause
*/

#include "Tree/ball_algorithm.h"

void dual_swap(double* darr, int* iarr, int i1, int i2)
{
  double dtmp = darr[i1];
  darr[i1]    = darr[i2];
  darr[i2]    = dtmp;
  int itmp    = iarr[i1];
  iarr[i1]    = iarr[i2];
  iarr[i2]    = itmp;
}

t_nheap::t_nheap(int n_pts, int n_nbrs) : n_pts{n_pts}, n_nbrs{n_nbrs}
{
  this->distances.resize(n_pts);
  for (int i = 0; i < n_pts; i++)
  {
    this->distances[i].resize(n_nbrs, INFINITY);
  }
  this->indices.resize(n_pts);
  for (int i = 0; i < n_pts; i++)
    this->indices[i].resize(n_nbrs);
}

void t_nheap::load(t_btree &b, const std::vector<SpacePoint> &x)
{
  double dist;
  for (int i = 0; i < this->n_pts; i++)
  {
    dist = b.min_dist(0, x[i]);
    b.query_depth_first(0, x[i], i, *this, dist);
  }
}

double t_nheap::largest(const int row) const
{
  return this->distances[row][0];
}

int t_nheap::push(int row, double val, int i_val)
{
  int ic1, ic2, i_swap;

  int size         = this->n_nbrs;
  double* dist_arr = this->distances[row].data();
  int* ind_arr     = this->indices[row].data();

  // if distance is already greater than the furthest element, don't push
  if (val > dist_arr[0]) return (0);

  // insert the values at position 0
  dist_arr[0] = val;
  ind_arr[0]  = i_val;

  // descend the heap, swapping values until the max heap criterion is met
  int i = 0;
  while (true)
  {
    ic1 = 2 * i + 1;
    ic2 = ic1 + 1;

    if (ic1 >= size)
      break;
    if (ic2 >= size)
    {
      if (dist_arr[ic1] > val)
        i_swap = ic1;
      else
        break;
    }
    else if (dist_arr[ic1] >= dist_arr[ic2])
    {
      if (val < dist_arr[ic1])
        i_swap = ic1;
      else
        break;
    }
    else
    {
      if (val < dist_arr[ic2])
        i_swap = ic2;
      else
        break;
    }
    dist_arr[i] = dist_arr[i_swap];
    ind_arr[i]  = ind_arr[i_swap];
    i           = i_swap;
  }

  dist_arr[i] = val;
  ind_arr[i]  = i_val;

  return (0);
}

void simultaneous_sort(double* dist, int* idx, int size)
{
  int pivot_idx, store_idx;
  double pivot_val;

  if (size <= 1)
    ;
  else if (size == 2)
  {
    if (dist[0] > dist[1]) dual_swap(dist, idx, 0, 1);
  }
  else if (size == 3)
  {
    if (dist[0] > dist[1]) dual_swap(dist, idx, 0, 1);
    if (dist[1] > dist[2])
    {
      dual_swap(dist, idx, 1, 2);
      if (dist[0] > dist[1]) dual_swap(dist, idx, 0, 1);
    }
  }
  else
  {
    pivot_idx = size / 2;
    if (dist[0] > dist[size - 1]) dual_swap(dist, idx, 0, size - 1);
    if (dist[size - 1] > dist[pivot_idx])
    {
      dual_swap(dist, idx, size - 1, pivot_idx);
      if (dist[0] > dist[size - 1]) dual_swap(dist, idx, 0, size - 1);
    }
    pivot_val = dist[size - 1];

    store_idx = 0;
    for (int i = 0; i < size - 1; i++)
    {
      if (dist[i] < pivot_val)
      {
        dual_swap(dist, idx, i, store_idx);
        store_idx++;
      }
    }
    dual_swap(dist, idx, store_idx, size - 1);
    pivot_idx = store_idx;
    if (pivot_idx > 1) simultaneous_sort(dist, idx, pivot_idx);
    if (pivot_idx * 2 < size)
      simultaneous_sort(dist + pivot_idx + 1, idx + pivot_idx + 1,
                        size - pivot_idx - 1);
  }
}

void t_nheap::sort()
{
  for (int row = 0; row < this->n_pts; row++)
    simultaneous_sort(this->distances[row].data(), this->indices[row].data(), this->n_nbrs);
}
