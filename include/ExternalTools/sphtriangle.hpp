///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// SphTriangle                                                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
#ifndef sphtriangleH
#define sphtriangleH

typedef struct
{
  int n_nodes; /* Number of nodes */
  int sph_size; /* Size of arrays sph_list and sph_lptr */
  double *sph_x; /* Array of X-coordinates for nodes */
  double *sph_y; /* Array of Y-coordinates for nodes */
  double *sph_z; /* Array of Z-coordinates for nodes */
  int *sph_list; /* Set of nodal indexes */
  int *sph_lptr; /* Set of pointers (sph_list indexes) */
  int *sph_lend; /* Set of pointers to adjacency lists */
} SphTriangle;

int trmesh_(int *n,
            double *x,
            double *y,
            double *z__,
            int *list,
            int *lptr,
            int *lend,
            int *lnew,
            int *near__,
            int *next,
            double *dist,
            int *ier);
int trlist_(int *n,
            int *list,
            int *lptr,
            int *lend,
            int *nrow,
            int *nt,
            int *ltri,
            int *ier);

#endif // #ifndef sphtriangleH