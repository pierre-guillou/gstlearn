/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

#include "Enum/ERotation.hpp"

#include "gstlearn_export.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"

class GSTLEARN_EXPORT GeometryHelper
{
public:
  static void rotationGetSinCos(double angle, double *cosa, double *sina);

  static void rotationIdentity(int ndim, double* rot);
  static void rotationInit(double angle, double* rot);
  static void rotationInit(double alpha,
                           double beta,
                           double gamma,
                           double* rot);
  static void rotationInit(int ndim, const double* angles, double* rot);
  static VectorDouble rotationInit(int ndim, const VectorDouble &angles);
  static MatrixSquareGeneral EulerToRotation(const VectorDouble &angles,
                                             const ERotation &convrot = ERotation::fromKey("SXYZ"));

  static void rotationGetDirection(double ct, double st, double* a, double *codir);
  static void rotationGetDirection(int ndim,
                                   int ndir,
                                   const VectorDouble &angles,
                                   VectorDouble &codir);
  static void rotationCopy(int ndim, const double* rotin, double* rotout);

  static int rotationGetAngles(int  ndim, const double* rot, double* angles);
  static void rotationGetAngles(const VectorDouble &codir,
                                VectorDouble &angles);
  static VectorDouble rotationGetAngles(const VectorDouble& codir);
  static VectorDouble rotationToEuler(const MatrixSquareGeneral &mat,
                                      const ERotation &convrot = ERotation::fromKey("SXYZ"),
                                      double eps = EPSILON10);

  static bool rotationIsIdentity(int ndim, double* rot, double eps = EPSILON10);

  static void mergeBoxes(VectorDouble &mini1,
                         VectorDouble &maxi1,
                         VectorDouble &mini2,
                         VectorDouble &maxi2);

  static double distancePointToSegment(double x0,
                                       double y0,
                                       double x1,
                                       double y1,
                                       double x2,
                                       double y2,
                                       double* xd,
                                       double* yd,
                                       int *nint);
  static int segmentIntersect(double xd1,
                              double yd1,
                              double xe1,
                              double ye1,
                              double xd2,
                              double yd2,
                              double xe2,
                              double ye2,
                              double* xint,
                              double* yint);
  static bool isSegmentIntersect(double xd1,
                                 double yd1,
                                 double xe1,
                                 double ye1,
                                 double xd2,
                                 double yd2,
                                 double xe2,
                                 double ye2);

  static double geodeticAngularDistance(double long1,
                                        double lat1,
                                        double long2,
                                        double lat2,
                                        double radius = 1.);
  static void geodeticAngles(double long1,
                             double lat1,
                             double long2,
                             double lat2,
                             double long3,
                             double lat3,
                             double* a,
                             double* b,
                             double* c,
                             double* A,
                             double* B,
                             double* C);
  static double geodeticTrianglePerimeter(double long1,
                                          double lat1,
                                          double long2,
                                          double lat2,
                                          double long3,
                                          double lat3);
  static double geodeticTriangleSurface(double long1,
                                        double lat1,
                                        double long2,
                                        double lat2,
                                        double long3,
                                        double lat3);
  static bool isInSphericalTriangle(double* coor,
                                    double surface,
                                    double *pts1,
                                    double *pts2,
                                    double *pts3,
                                    double *wgts,
                                    double eps = EPSILON6);
  static bool isInSphericalTriangleOptimized(const double *coor,
                                             double* ptsa,
                                             double* ptsb,
                                             double* ptsc,
                                             double* wgts,
                                             double eps = EPSILON6);
  static VectorVectorDouble convertLongLat(const VectorDouble &longitude,
                                           const VectorDouble &latitude,
                                           double dilate = 1.,
                                           double radius_arg = 1.);
  static void convertCart2Sph(double x,
                              double y,
                              double z,
                              double* rlong,
                              double* rlat,
                              double radius_arg = 1.);
  static void convertSph2Cart(double rlong,
                              double rlat,
                              double* x,
                              double* y,
                              double* z,
                              double radius_arg = 1.);
  static MatrixSquareGeneral gradXYToRotmat(double dzoverdx, double dzoverdy);

private:
  static void _decodeConvRot(const ERotation &convrot,
                             int* firstaxis,
                             int* parity,
                             int* repetition,
                             int* frame);
};

//typedef GeometryHelper GH;
class GH: public GeometryHelper {};