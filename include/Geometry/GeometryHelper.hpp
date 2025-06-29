/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "Basic/VectorNumT.hpp"
#include "Enum/ERotation.hpp"

#include "gstlearn_export.hpp"
#include "Matrix/MatrixSquare.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Matrix/MatrixInt.hpp"
#include <vector>

class GSTLEARN_EXPORT GeometryHelper
{
public:
  static void rotationGetSinCos(double angle, double* cosa, double* sina);

  static void rotationMatrixIdentityInPlace(int ndim, VectorDouble& rot);
  static void rotation2DMatrixInPlace(double angle, VectorDouble& rot);
  static void rotation3DMatrixInPlace(double alpha,
                                      double beta,
                                      double gamma,
                                      VectorDouble& rot);
  static void rotationMatrixInPlace(int ndim,
                                    const VectorDouble& angles,
                                    VectorDouble& rot);
  static VectorDouble rotationMatrix(int ndim, const VectorDouble& angles);
  static void rotationGetAnglesInPlace(int ndim,
                                       const double* rot,
                                       double* angles);
  static void rotationGetAnglesInPlace(const VectorDouble& rot,
                                       VectorDouble& angles);
  static void rotationCopy(int ndim, const double* rotin, double* rotout);
  static bool rotationIsIdentity(int ndim, const double* rot, double eps = EPSILON10);
  static void rotationMatrixDerivativesInPlace(int ndim,
                                               VectorDouble& angles,
                                               std::vector<MatrixSquare>& dR);
  static void rotation2DMatrixDerivativesInPlace(double angle, MatrixSquare& dR);
  static void rotation3DMatrixDerivativesInPlace(VectorDouble& angles,
                                                 std::vector<MatrixSquare>& dR);

  static MatrixSquare EulerToRotation(const VectorDouble& angles,
                                      const ERotation& convrot = ERotation::fromKey("SXYZ"));
  static std::vector<MatrixSquare> EulerToRotationDerivatives(const VectorDouble& angles,
                                                              const ERotation& convrot = ERotation::fromKey("SXYZ"));
  static void EulerToRotationDerivativesInPlace(std::vector<MatrixSquare>& dR,
                                                const VectorDouble& angles,
                                                const ERotation& convrot = ERotation::fromKey("SXYZ"));

  static void rotationGetRandomDirection(double ct,
                                         double st,
                                         const double* a,
                                         double* codir);
  static void rotationGetDirection2D(const VectorDouble& angles,
                                     VectorDouble& codir);
  static void rotationGetDirectionDefault(int ndim, VectorDouble& codir);
  static void rotationGetAnglesFromCodirInPlace(const VectorDouble& codir,
                                                VectorDouble& angles);
  static VectorDouble rotationGetAngles(const VectorDouble& codir,
                                        bool flagResize = false);
  static VectorDouble rotationFromIncrements(const VectorDouble& incr, bool flagDegree = false);

  static VectorDouble rotationToEuler(const MatrixSquare& mat,
                                      const ERotation& convrot = ERotation::fromKey(
                                        "SXYZ"),
                                      double eps = EPSILON10);

  static double distancePointToSegment(double x0,
                                       double y0,
                                       double x1,
                                       double y1,
                                       double x2,
                                       double y2,
                                       double* xd,
                                       double* yd,
                                       int* nint);
  static bool segmentIntersect(double xd1,
                               double yd1,
                               double xe1,
                               double ye1,
                               double xd2,
                               double yd2,
                               double xe2,
                               double ye2,
                               double* xint,
                               double* yint);

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
                                    double* pts1,
                                    double* pts2,
                                    double* pts3,
                                    double* wgts,
                                    double eps = EPSILON6);
  static bool isInSphericalTriangleOptimized(const double* coor,
                                             const double* ptsa,
                                             const double* ptsb,
                                             const double* ptsc,
                                             double* wgts,
                                             double eps = EPSILON6);
  static VectorVectorDouble convertLongLatTo3D(const VectorDouble& longitude,
                                               const VectorDouble& latitude,
                                               double dilate     = 1.,
                                               double radius_arg = 1.);
  static VectorVectorDouble convert3DToLongLat(const VectorDouble& x,
                                               const VectorDouble& y,
                                               const VectorDouble& z,
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
  static MatrixSquare gradXYToRotmat(double dzoverdx, double dzoverdy);
  static MatrixDense* getDirectionsInR3(const MatrixDense* U);
  static MatrixDense* getDirectionsInRn(const MatrixDense* U);

  static double formatAngle(double anglein, double basis = 360.);
  static VectorDouble formatAngles(const VectorDouble& anglesin, double basis = 360.);

  static VectorDouble rayTriangleIntersect(const VectorDouble& dir,
                                           const VectorDouble& v0,
                                           const VectorDouble& v1,
                                           const VectorDouble& v2);
  static VectorVectorDouble sphBarCoord(const VectorVectorDouble& sphPts,
                                        const MatrixDense& apices,
                                        const MatrixInt& meshes);

  static double getCosineAngularTolerance(double tolang);

  static VectorVectorDouble getEllipse(const VectorDouble& center,
                                       double rx,
                                       double ry,
                                       double theta,
                                       int count = 360);

private:
  static void _decodeConvRot(const ERotation& convrot,
                             int* firstaxis,
                             int* parity,
                             int* repetition,
                             int* frame);
};

// typedef GeometryHelper GH;
class GH: public GeometryHelper
{
};
