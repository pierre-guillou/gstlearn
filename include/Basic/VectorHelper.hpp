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

#include "gstlearn_export.hpp"
#include "geoslib_define.h"

#include "Basic/VectorNumT.hpp"

class GSTLEARN_EXPORT VectorHelper
{
public:
  static VectorInt          initVInt(int nval, int value = 0.);
  static VectorDouble       initVDouble(int nval, double value = 0.);
  static VectorVectorDouble initVVDouble(int nval1, int nval2, double value = 0.);
  static VectorVectorInt    initVVInt(int nval1, int nval2, int value = 0);

  static VectorInt          initVInt(int* values, int number);
  static VectorDouble       initVDouble(double* values, int number);
  static VectorVectorDouble initVVDouble(double* value, int n1, int n2);

  static void display(const String &title, const VectorDouble &vect); // TODO rename
  static void display(const String &title, const VectorVectorDouble &vect);
  static void display(const String &title, const VectorString &vect);
  static void display(const String &title, const VectorInt &vect);

  static String toString(const VectorDouble& vec); // TODO rename
  static String toString(const VectorVectorDouble& vec);
  static String toString(const VectorString& vec);
  static String toString(const VectorInt& vec);

  static void displayStats(const String &title, const VectorDouble &vect);
  static void displayRange(const String &title, const VectorDouble &vect);
  static void displayRange(const String &title, const VectorInt &vect);

  static int maximum(const VectorInt &vec);
  static int minimum(const VectorInt &vec);
  static double maximum(const VectorDouble &vec, bool flagAbs = false);
  static double minimum(const VectorDouble &vec, bool flagAbs = false);
  static double maximum(const VectorVectorDouble &vec, bool flagAbs = false);
  static double minimum(const VectorVectorDouble &vec, bool flagAbs = false);
  static int product(const VectorInt& nx);
  static double product(const VectorDouble& nx);
  static int countUndefined(const VectorDouble& vec);
  static int countDefined(const VectorDouble& vec);
  static double extensionDiagonal(const VectorDouble& mini, const VectorDouble& maxi);

  static double cumul(const VectorDouble &vec);
  static double mean(const VectorDouble &vec);
  static double variance(const VectorDouble &vec);
  static double stdv(const VectorDouble &vec);
  static double norm(const VectorDouble &vec);
  static double correlation(const VectorDouble &veca, const VectorDouble &vecb);

  static bool isConstant(const VectorDouble& vect, double refval = TEST);
  static bool isConstant(const VectorInt& vect, int refval = ITEST);
  static bool isSame(const VectorDouble &v1,
                     const VectorDouble &v2,
                     double eps = EPSILON10);
  static bool isSame(const VectorInt &v1, const VectorInt &v2);

  static VectorInt sequence(int number, int ideb = 0);
  static VectorDouble sequence(double valFrom, double valTo, double valStep);
  static void fill(VectorDouble& vec, double v, int size = 0);
  static void fill(VectorInt& vec, int v, int size = 0);
  static void fill(VectorVectorDouble &vec, double value);

  static VectorDouble add(const VectorDouble &veca, const VectorDouble &vecb);
  static void addInPlace(VectorDouble &dest, const VectorDouble &src);
  static void addInPlace(const VectorDouble &veca,
                         const VectorDouble &vecb,
                         VectorDouble &res);
  static VectorDouble subtract(const VectorDouble& veca, const VectorDouble& vecb);
  static void subtractInPlace(VectorDouble &dest, const VectorDouble &src);
  static void multiplyInPlace(VectorDouble& vec, const VectorDouble& v);
  static void divideInPlace(VectorDouble& vec, const VectorDouble& v);

  static void multiplyConstant(VectorDouble& vec, double v);
  static void divideConstant(VectorDouble& vec, double v);
  static void copy(VectorDouble& veca, const VectorDouble& vecb);
  static void addConstant(VectorDouble& vec, double v);
  static void addConstant(VectorInt& vec, int v);

  static void normalize(VectorDouble& vec);
  static VectorDouble concatenate(const VectorDouble &veca, const VectorDouble &vecb);
  static VectorDouble power(const VectorDouble& vec, double power);
  static VectorDouble inverse(const VectorDouble& vec);

  static double innerProduct(const VectorDouble &veca, const VectorDouble &vecb);
  static VectorDouble crossProduct(const VectorDouble &veca, const VectorDouble &vecb);

  static void cumulate(VectorDouble &veca,
                       const VectorDouble &vecb,
                       double coeff = 1.,
                       double addval = 0.);

  static VectorDouble simulateUniform(int n,
                                      double mini = 0.,
                                      double maxi = 1.);
  static VectorDouble simulateBernoulli(int n,
                                        double proba,
                                        double vone = 1.,
                                        double velse = 0.);
  static VectorDouble simulateGaussian(int n,
                                       double mean = 0.,
                                       double sigma = 1.);
  static void simulateGaussianInPlace(VectorDouble &vect,
                                      double mean = 0.,
                                      double sigma = 1.);
  static VectorInt sampleRanks(int ntotal,
                               double proportion = 0.,
                               int number = 0,
                               int seed = 242141);

  static VectorInt    sort(const VectorInt& vecin, bool ascending = true);
  static VectorDouble sort(const VectorDouble& vecin, bool ascending = true);
  static VectorInt    sortRanks(const VectorDouble& vecin);

  static std::pair<double,double> rangeVals(const VectorDouble& vec);

private:
  static int _myrandom (int i);
};

//typedef VectorHelper VH;
class VH: public VectorHelper {};
