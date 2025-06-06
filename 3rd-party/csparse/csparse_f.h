/*
                                      csparse

Original Author: Timothy Davis
Website: https://people.math.sc.edu/Burkardt/c_src/csparse/csparse.html
License: LGPL v2.1
*/

/*
Author: Timothy Davis

License:

CSPARSE: a Concise Sparse matrix package.
Copyright (c) 2006, Timothy A. Davis.
http://www.cise.ufl.edu/research/sparse/CSparse

CSPARSE is free software; you can redistribute it and/or modify it under the terms of
the GNU Lesser General Public License as published by the Free Software Foundation;
either version 2.1 of the License, or (at your option) any later version.

CSPARSE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with
this Module; if not, write to the Free Software Foundation,
Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
*/

/*
Modified by MINES Paris / ARMINES (2023)
Authors: gstlearn Team
Website: https://gstlearn.org
*/
#ifndef _CS_F_H
#define _CS_F_H

#include <cstdio>
#include <stdlib.h>
#include <stddef.h>
#include "csparse_d.h"

cs     *cs_add (const cs *A, const cs *B, double alpha, double beta) ;
int     cs_cholsol (const cs *A, double *b, int order) ;
int     cs_dupl (cs *A) ;
int     cs_entry (cs *T, int i, int j, double x) ;
int     cs_lusol (const cs *A, double *b, int order, double tol) ;
int     cs_gaxpy (const cs *A, const double *x, double *y) ;
cs     *cs_multiply (const cs *A, const cs *B) ;
int     cs_qrsol (const cs *A, double *b, int order) ;
cs     *cs_transpose (const cs *A, int values) ;
cs     *cs_triplet (const cs *T) ;
double  cs_norm (const cs *A) ;
int     cs_print (const cs *A, int brief) ;
cs     *cs_load (FILE *f) ;

/* utilities */
void   *cs_calloc  (int n, size_t size) ;
void   *cs_free    (void *p) ;
void   *cs_realloc (void *p, int n, size_t size, int *ok) ;
cs     *cs_spalloc (int m, int n, int nzmax, int values, int triplet) ;
cs     *cs_spfree  (cs *A) ;
int     cs_sprealloc (cs *A, int nzmax) ;
void   *cs_malloc  (int n, size_t size) ;

/* --- secondary CSparse routines and data structures ----------------------- */

int    *cs_amd (const cs *A, int order) ;
csn    *cs_chol (const cs *A, const css *S) ;
csd    *cs_dmperm (const cs *A) ;
int     cs_droptol (cs *A, double tol) ;
int     cs_dropzeros (cs *A) ;
int     cs_happly (const cs *V, int i, double beta, double *x) ;
int     cs_ipvec (int n, const int *P, const double *b, double *x) ;
int     cs_lsolve (const cs *L, double *x) ;
int     cs_ltsolve (const cs *L, double *x) ;
csn    *cs_lu (const cs *A, const css *S, double tol) ;
cs     *cs_permute (const cs *A, const int *P, const int *Q, int values) ;
int    *cs_pinv (const int *P, int n) ;
int     add_cs_pvec(int n, const int *P, const double *b, double *x);
int     cs_pvec (int n, const int *P, const double *b, double *x) ;
csn    *cs_qr (const cs *A, const css *S) ;
css    *cs_schol (const cs *A, int order) ;
css    *cs_sqr (const cs *A, int order, int qr) ;
cs     *cs_symperm (const cs *A, const int *Pinv, int values) ;
int     cs_usolve (const cs *U, double *x) ;
int     cs_utsolve (const cs *U, double *x) ;
int     cs_updown (cs *L, int sigma, const cs *C, const int *parent) ;
css    *cs_sfree (css *S) ;
csn    *cs_nfree (csn *N) ;
csd    *cs_dfree (csd *D) ;

/* --- tertiary CSparse routines -------------------------------------------- */
int    *cs_counts (const cs *A, const int *parent, const int *post, int ata) ;
int     cs_cumsum (int *p, int *c, int n) ;
int     cs_dfs (int j, cs *L, int top, int *xi, int *pstack, const int *Pinv) ;
int    *cs_etree (const cs *A, int ata) ;
int     cs_fkeep (cs *A, int (*fkeep) (int, int, double, void *), void *other) ;
double  cs_house (double *x, double *beta, int n) ;
int    *cs_maxtrans (const cs *A) ;
int    *cs_post (int n, const int *parent) ;
int     cs_reach (cs *L, const cs *B, int k, int *xi, const int *Pinv) ;
csd    *cs_scc (cs *A) ;
int     cs_scatter (const cs *A, int j, double beta, int *w, double *x,
                                    int mark, cs *C, int nz) ;
int     cs_splsolve (cs *L, const cs *B, int k, int *xi, double *x,
                                     const int *Pinv) ;
int     cs_tdfs (int j, int k, int *head, const int *next, int *post,
                                 int *stack) ;

/* utilities */
csd    *cs_dalloc (int m, int n) ;
cs     *cs_done (cs *C, void *w, void *x, int ok) ;
int    *cs_idone (int *p, cs *C, void *w, int ok) ;
csn    *cs_ndone (csn *N, cs *C, void *w, void *x, int ok) ;
csd    *cs_ddone (csd *D, cs *C, void *w, int ok) ;

/* accessors */
int     cs_getncol(const cs* mat);
int     cs_getnrow(const cs* mat);

#define CS_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define CS_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define CS_FLIP(i) (-(i)-2)
#define CS_UNFLIP(i) (((i) < 0) ? CS_FLIP(i) : (i))
#define CS_MARKED(Ap,j) (Ap [j] < 0)
#define CS_MARK(Ap,j) { Ap [j] = CS_FLIP (Ap [j]) ; }
#define CS_OVERFLOW(n,size) (n > INT_MAX / (int) size)

#endif
