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

#include "gstlearn_export.hpp"
#include "geoslib_define.h"
#include "Basic/VectorNumT.hpp"
#include <vector>

namespace gstlrn
{
class GSTLEARN_EXPORT VectorHelper
{
public:
  static VectorInt          initVInt(Id nval, Id value = 0.);
  static VectorDouble       initVDouble(Id nval, double value = 0.);
  static VectorVectorDouble initVVDouble(Id nval1, Id nval2, double value = 0.);
  static VectorVectorInt    initVVInt(Id nval1, Id nval2, Id value = 0);

  static VectorInt          initVInt(const Id* values, Id number);
  static VectorDouble       initVDouble(const double* values, Id number);
  static VectorVectorDouble initVVDouble(const double* value, Id n1, Id n2);

  static VectorString initVString(Id ntab, char** names);
  
  static void dump(const String& title, const VectorVectorInt& vect, bool skipLine = true);
  static void dump(const String& title, const VectorVectorDouble& vect, bool skipLine = true);
  static void dump(const String& title, const VectorDouble& vect, bool skipLine = true);
  static void dump(const String& title, const VectorString &vect, bool skipLine = true);
  static void dump(const String& title, const VectorInt &vect, bool skipLine = true);
  static String toStringAsSpan(constvect vec);
  static String toStringAsVD(const VectorDouble& vec); // TODO rename
  static String toStringAsVVD(const VectorVectorDouble& vec);
  static String toStringAsVVI(const VectorVectorInt& vec);
  static String toStringAsVS(const VectorString& vec);
  static String toStringAsVI(const VectorInt& vec);

  #ifndef SWIG
  static void dumpStats(const String& title, constvect vect, Id nmax = -1);
  static void dumpRange(const String& title, constvect vect, Id nmax = -1);
#endif
  static void dumpStats(const String& title, const VectorDouble& vectin, Id nmax = -1);
  static void dumpRange(const String &title, const VectorDouble& vectin, Id nmax = -1);
  static void dumpRange(const String &title, const VectorInt &vect);
  static void dumpNNZ(const String &title, const VectorDouble &vect, Id nclass = 10);

  static Id maximum(const VectorInt &vec, bool flagAbs = false);
  static Id minimum(const VectorInt &vec, bool flagAbs = false);
  static double maximum(const VectorDouble &vec, bool flagAbs = false, const VectorDouble& aux = VectorDouble(), Id mode=0);
  static double minimum(const VectorDouble &vec, bool flagAbs = false, const VectorDouble& aux = VectorDouble(), Id mode=0);
  static double maximum(const VectorVectorDouble &vec, bool flagAbs = false);
  static double maximum(const std::vector<std::vector<double>> &vec, bool flagAbs = false);

  static double minimum(const VectorVectorDouble &vec, bool flagAbs = false);
  static Id product(const VectorInt& vec);
  static double product(const VectorDouble& vec);
  static Id countUndefined(const VectorDouble& vec);
  static Id countDefined(const VectorDouble& vec);
  static bool hasUndefined(const VectorDouble& vec);
  static double extensionDiagonal(const VectorDouble& mini, const VectorDouble& maxi);

  static Id    count(const VectorVectorInt& vec);
  static Id    cumul(const VectorInt& vec);
  static Id    cumul(const VectorVectorInt& vec);
  static double cumul(const VectorDouble &vec);
  static double cumulLog(const VectorDouble &vec);
  static VectorInt cumulIncrement(const VectorVectorInt& vec);
  static double mean(const VectorDouble& vec);
  static double variance(const VectorDouble &vec, bool scaleByN = false);
  static double stdv(const VectorDouble &vec, bool scaleByN = false);
  static double norm(const VectorDouble &vec);
  static double normL1(const VectorDouble &vec);
  static double norminf(const VectorDouble &vec);
  static double median(const VectorDouble& vec);
  static double normDistance(const VectorDouble& veca, const VectorDouble& vecb);
  static double correlation(const VectorDouble &veca, const VectorDouble &vecb);
  static VectorDouble quantiles(const VectorDouble& vec, const VectorDouble& probas);

  static bool isConstant(const VectorDouble& vect, double refval = TEST);
  static bool isConstant(const VectorInt& vect, Id refval = ITEST);
  static bool isEqual(const VectorDouble &v1,
                     const VectorDouble &v2,
                     double eps = EPSILON10);
  static bool isEqual(const VectorInt& v1, const VectorInt& v2);
  static bool isEqualExtended(const VectorDouble& v1,
                              const VectorDouble& v2,
                              double eps           = EPSILON10,
                              bool flagRelative    = true,
                              bool flagAbsolute    = false,
                              const String& string = "");

  static void sequenceInPlace(Id n, VectorInt& vec);
  static VectorInt sequence(Id number, Id ideb = 0, Id step = 1);
  static VectorDouble sequence(double valFrom,
                               double valTo,
                               double valStep = 1.,
                               double ratio   = 1.);
  static void fill(VectorDouble& vec, double v, Id size = 0);
  static void fill(VectorInt& vec, Id v, Id size = 0);
  static void fill(VectorVectorDouble &vec, double value);
  static void fillUndef(VectorDouble& vec, double repl);
  
#ifndef SWIG
  static void addMultiplyConstantInPlace(double val1,
                                         const constvect in,
                                         vect out,
                                         Id iad);
  static double innerProduct(const constvect veca, const constvect vecb);

  static void addMultiplyVectVectInPlace(const constvect in1,
                                         const constvect in2,
                                         vect out,
                                         Id iad);

  static void addInPlace(constvect in, vect dest);

#endif
  static VectorDouble add(const VectorDouble &veca, const VectorDouble &vecb);
  static void addInPlace(VectorDouble& dest, const VectorDouble& src);
  static void addInPlace(VectorInt& dest, const VectorInt& src);
  static void addInPlace(std::vector<double>& dest, const std::vector<double> &src);
  static void addInPlace(const VectorDouble& veca,
                         const VectorDouble& vecb,
                         VectorDouble& res,
                         Id size = 0);
  static void addInPlace(const VectorInt& veca,
                         const VectorInt& vecb,
                         VectorInt& res,
                         Id size = 0);
  static void addInPlace(const double *veca,
                         const double *vecb,
                         double *res,
                         Id size);
  static void addInPlace(const VectorVectorDouble& in1,
                         const VectorVectorDouble& in2,
                         VectorVectorDouble& outv);
  static void addInPlace(const std::vector<std::vector<double>>& in1,
                         const std::vector<std::vector<double>>& in2,
                         std::vector<std::vector<double>>& outv);
  
  static void addSquareInPlace(VectorDouble &dest, const VectorDouble &src);
  static VectorDouble subtract(const VectorDouble& veca, const VectorDouble& vecb);
  static VectorDouble subtract(constvect veca,constvect vecb);
  static VectorInt    subtract(const VectorInt& veca, const VectorInt& vecb);
  static void   subtractInPlace(const constvect in1,
                                const constvect in2,
                                vect  outv);
  static void subtractInPlace(VectorDouble &dest, const VectorDouble &src);
  static void subtractInPlace(VectorInt &dest, const VectorInt &src);
  static void subtractInPlace(const VectorVectorDouble& in1,
                              const VectorVectorDouble& in2,
                              VectorVectorDouble& outv);
  static void subtractInPlace(const std::vector<std::vector<double>>& in1,
                              const std::vector<std::vector<double>>& in2,
                              std::vector<std::vector<double>>& outv);

  static VectorDouble multiply(const VectorDouble& veca, const VectorDouble& vecb);
  static void multiplyInPlace(const VectorDouble& veca, const VectorDouble& vecb, VectorDouble& res);
  static void multiplyInPlace(VectorDouble& vec, const VectorDouble& v);
  static VectorDouble divide(const VectorDouble& veca, const VectorDouble& vecb);
  static void divideInPlace(const VectorDouble& veca, const VectorDouble& vecb, VectorDouble& res);
  static void divideInPlace(VectorDouble& vec, const VectorDouble& v);
  static void divideInPlace(std::vector<double>& vec, const std::vector<double>& v);
  static void multiplyComplexInPlace(const VectorDouble& vecaRe,
                                     const VectorDouble& vecaIm,
                                     const VectorDouble& vecbRe,
                                     const VectorDouble& vecbIm,
                                     VectorDouble& resRe,
                                     VectorDouble& resIm);

  static void multiplyConstant(VectorDouble& vec, double v);
  static void multiplyConstantInPlace(const VectorDouble& vec, double v, VectorDouble& vecout);
  static void multiplyConstantSelfInPlace(VectorDouble &vec, double v);
  static void addMultiplyConstantInPlace(double val1,
                                         const VectorDouble &in1,
                                         VectorDouble &outv,
                                         Id iad);
  static void addMultiplyConstantInPlace(double val1,
                                         const VectorVectorDouble &in1,
                                         VectorVectorDouble &outv);
  static void divideConstant(VectorDouble& vec, double v);
  static void copy(const VectorDouble& vecin, VectorDouble& vecout, Id size = -1);
  static void copy(const VectorInt &vecin, VectorInt &vecout, Id size = -1);
  static void copy(const VectorVectorDouble &inv, VectorVectorDouble &outv);
  static void copy(const std::vector<std::vector<double>> &inv, std::vector<std::vector<double>> &outv);

  static void addConstant(VectorDouble& vec, double v);
  static void addConstant(VectorInt& vec, Id v);
  static void mean1AndMean2ToStdev(const VectorDouble &mean1,
                                   const VectorDouble &mean2,
                                   VectorDouble &std,
                                   Id number);

  static void normalize(VectorDouble& vec, Id norm=2);
  static void normalize(double *tab, Id ntab);
  static void normalizeFromGaussianDistribution(VectorDouble &vec,
                                                double mini = 0.,
                                                double maxi = 1.);
  static VectorDouble normalScore(const VectorDouble& data,
                                  const VectorDouble& wt = VectorDouble());
  static VectorDouble qnormVec(const VectorDouble& vec);
  static VectorDouble pnormVec(const VectorDouble& vec);
  static VectorDouble concatenate(const VectorDouble& veca, const VectorDouble& vecb);
  static void concatenateInPlace(VectorDouble& veca, const VectorDouble& vecb);
  static VectorDouble power(const VectorDouble& vec, double power);
  static VectorDouble inverse(const VectorDouble& vec);
#ifndef SWIG
  static void power(VectorDouble& res, const constvect vec, double power);
  static void inverse(VectorDouble& res, const constvect vec);
#endif // !SWIG

  static double innerProduct(const VectorDouble &veca, const VectorDouble &vecb, Id size = -1);
  static double innerProduct(const double* veca, const double* vecb, Id size);
  static double innerProduct(const VectorVectorDouble &x,
                             const VectorVectorDouble &y);
  static double innerProduct(const std::vector<double> &veca, const std::vector<double> &vecb, Id size = -1);

  static VectorDouble crossProduct3D(const VectorDouble &veca, const VectorDouble &vecb);
  static void crossProduct3DInPlace(const double *a, const double *b, double *v);

  static VectorDouble cumsum(const VectorDouble &vecin, bool flagAddZero, bool revert=false);
  static void cumulateInPlace(VectorDouble &vec);
  static void cumulate(VectorDouble &veca, const VectorDouble &vecb, double coeff = 1., double addval = 0.);
  static void getMostSignificant(const VectorDouble& vec, double tol = EPSILON6, Id nmax = -1);

  static VectorDouble simulateUniform(Id n = 1,
                                      double mini = 0.,
                                      double maxi = 1.);
  static VectorDouble simulateBernoulli(Id n = 1,
                                        double proba = 0.5,
                                        double vone = 1.,
                                        double velse = 0.);
  static VectorDouble simulateGaussian(Id n = 1,
                                       double mean = 0.,
                                       double sigma = 1.);
  static void simulateGaussianInPlace(VectorDouble &vec,
                                      double mean = 0.,
                                      double sigma = 1.);
  static VectorInt sampleRanks(Id ntotal,
                               double proportion = 0.,
                               Id number = 0,
                               Id seed = 242141,
                               Id optSort = 0);
  static void normalizeCodir(Id ndim, VectorDouble &codir);

  static bool         isInList(const VectorInt& vec, Id item);
  static VectorInt    sort(const VectorInt& vecin, bool ascending = true, Id size = -1);
  static VectorDouble sort(const VectorDouble& vecin, bool ascending = true, Id size = -1);
  static void         sortInPlace(VectorInt& vecin, bool ascending = true, Id size = -1);
  static void         sortInPlace(VectorDouble& vecin, bool ascending = true, Id size = -1);
  static bool         isSorted(const VectorDouble& vec, bool ascending = true);
  static VectorDouble unique(const VectorDouble& vecin, Id size = -1);
  static VectorInt    unique(const VectorInt& vecin, Id size = -1);
  static VectorInt    orderRanks(const VectorInt& vecin, bool ascending = true, Id size = -1);
  static VectorInt    orderRanks(const VectorDouble& vecin, bool ascending = true, Id size = -1);
  static VectorInt    sortRanks(const VectorDouble& vecin, bool ascending = true, Id size = -1);
  static VectorInt    reorder(const VectorInt& vecin, const VectorInt& order, Id size = -1);
  static VectorDouble reorder(const VectorDouble& vecin, const VectorInt& order, Id size = -1);
  static VectorDouble revert(const VectorDouble& vecin);
  static VectorInt    revert(const VectorInt& vecin);
  static VectorDouble sample(const VectorDouble& vecin, const VectorInt& indKeep);
  
  static void arrangeInPlace(Id safe,
                             VectorInt& ranks,
                             VectorDouble& values,
                             bool ascending = true,
                             Id size       = -1);
  static void arrangeInPlace(Id safe,
                             VectorInt &ranks,
                             VectorInt &values,
                             bool ascending = true,
                             Id size = -1);
  static VectorInt filter(const VectorInt &vecin,
                          Id vmin = ITEST,
                          Id vmax = ITEST,
                          bool ascending = true);
  static VectorInt complement(const VectorInt& vec, const VectorInt& sel);

  static std::pair<double,double> rangeVals(const VectorDouble& vec);
  static void unflattenInPlace(const std::vector<double>& vd, std::vector<std::vector<double>>& vvd);
  static void flattenInPlace(const std::vector<std::vector<double>>& vvd, std::vector<double>& vd);
  static VectorDouble flatten(const VectorVectorDouble& vvd);
  static VectorVectorDouble unflatten(const VectorDouble& vd, const VectorInt& sizes);
  static std::vector<double> flatten(const std::vector<std::vector<double>>& vvd);
  static std::vector<std::vector<double>> unflatten(const std::vector<double>& vd, const VectorInt& sizes);
  static void flattenInPlace(const VectorVectorDouble& vvd, VectorDouble& vd);
  static void linearCombinationInPlace(double val1,
                                       const VectorDouble &vd1,
                                       double val2,
                                       const VectorDouble &vd2,
                                       VectorDouble &outv);
  static void linearCombinationVVDInPlace(double val1,
                                          const VectorVectorDouble &vvd1,
                                          double val2,
                                          const VectorVectorDouble &vvd2,
                                          VectorVectorDouble &outv);
  static double innerProduct(const std::vector<std::vector<double>> &x,
                             const std::vector<std::vector<double>> &y);
  static void linearCombinationVVDInPlace(double val1,
                                          const std::vector<std::vector<double>> &vvd1,
                                          double val2,
                                          const std::vector<std::vector<double>> &vvd2,
                                          std::vector<std::vector<double>> &outv);

  static VectorDouble suppressTest(const VectorDouble& vecin);
  static void extractInPlace(const VectorDouble& vecin, VectorDouble& vecout, Id start);
  static void mergeInPlace(const VectorDouble& vecin, VectorDouble& vecout, Id start);

  static void transformVD(VectorDouble& tab, Id oper_choice = 1);

  static void squeezeAndStretchInPlaceForward(const VectorDouble &vecin,
                                              VectorDouble &vecout,
                                              double origin,
                                              double mesh,
                                              double top,
                                              double bot);
  static void squeezeAndStretchInPlaceBackward(const VectorDouble &vecin,
                                               VectorDouble &vecout,
                                               double origin,
                                               double mesh,
                                               double top,
                                               double bot);

  static Id whereMinimum(const VectorDouble& tab);
  static Id whereMaximum(const VectorDouble& tab);
  static Id whereElement(const VectorInt& tab, Id target);
  static double norm(const std::vector<double>& vec);
  static bool isIsotropic(const VectorVectorInt& sampleRanks);

  static VectorDouble reduceOne(const VectorDouble &vecin, Id index);
  static VectorDouble reduce(const VectorDouble &vecin, const VectorInt& vindex);
  static VectorDouble compress(const VectorDouble &vecin, const VectorInt& vindex);
  static void truncateDecimalsInPlace(VectorDouble& vec, Id ndec);
  static void truncateDigitsInPlace(VectorDouble& vec, Id ndec);
  static void simulateGaussianInPlace(std::vector<double> &vec,
                                           double mean = 0.,
                                           double sigma = 1.);
};

//typedef VectorHelper VH;
class VH: public VectorHelper {};

}