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
#include "Matrix/MatrixSquare.hpp"
#include "Enum/EOperator.hpp"

#include <map>
#include <cmath>
#include <math.h>

typedef struct
{
  int number;
  int nvalid;
  double mini;
  double maxi;
  double delta;
  double mean;
  double stdv;
} StatResults;

GSTLEARN_EXPORT bool   isInteger(double value, double eps = EPSILON10);
GSTLEARN_EXPORT int    getClosestInteger(double value);
GSTLEARN_EXPORT bool   isMultiple(int nbig, int nsmall);
GSTLEARN_EXPORT bool   isOdd(int number);
GSTLEARN_EXPORT bool   isEven(int number);
GSTLEARN_EXPORT bool   isZero(double value, double eps = EPSILON10);
GSTLEARN_EXPORT bool   isOne(double value, double eps = EPSILON10);
GSTLEARN_EXPORT bool   isEqual(double v1, double v2, double eps = EPSILON10);
GSTLEARN_EXPORT double getMin(double val1, double val2);
GSTLEARN_EXPORT double getMax(double val1, double val2);
GSTLEARN_EXPORT double ut_deg2rad(double angle);
GSTLEARN_EXPORT double ut_rad2deg(double angle);

GSTLEARN_EXPORT bool isEqualExtended(double v1,
                                     double v2,
                                     double eps           = EPSILON10,
                                     bool flagRelative    = true,
                                     bool flagAbsolute    = false,
                                     const String& string = "");

// No need this stuff through SWIG (because we use target language NAs)
#ifndef SWIG

GSTLEARN_EXPORT bool   FFFF(double value); // TODO isNA<double>
GSTLEARN_EXPORT bool   IFFFF(int value);   // TODO isNA<int>
GSTLEARN_EXPORT double getTEST();  // TODO getNA<double>
GSTLEARN_EXPORT int    getITEST(); // TODO getNA<int>

#define DOUBLE_NA  TEST
#define    INT_NA  ITEST
#define STRING_NA  "NA"
#define  FLOAT_NA  static_cast<float>(TEST)   // 1.234e30 is ok for 4 bytes but needs a cast for Windows

template <typename T> inline T getNA();
template <> inline double getNA() { return DOUBLE_NA; }
template <> inline int    getNA() { return INT_NA; }
template <> inline String getNA() { return STRING_NA; }
template <> inline float  getNA() { return FLOAT_NA; }

template <typename T> inline bool isNA(const T& v);
template <> inline bool isNA(const double& v) { return (v == getNA<double>() || std::isnan(v) || std::isinf(v)); }
template <> inline bool isNA(const int& v)    { return (v == getNA<int>()); }
template <> inline bool isNA(const String& v) { return (v == getNA<String>()); }
template <> inline bool isNA(const float& v)  { return (v == getNA<float>()  || std::isnan(v) || std::isinf(v)); }

#endif // SWIG

// Other Utility functions

GSTLEARN_EXPORT void ut_sort_double(int safe, int nech, int *ind, double *value);
GSTLEARN_EXPORT StatResults ut_statistics(int nech,
                                          const double *tab,
                                          const double *sel = nullptr,
                                          const double *wgt = nullptr);
GSTLEARN_EXPORT void ut_stats_mima_print(const char *title, int nech, double *tab, double *sel);
GSTLEARN_EXPORT void ut_facies_statistics(int nech,
                                          double *tab,
                                          double *sel,
                                          int *nval,
                                          int *mini,
                                          int *maxi);
GSTLEARN_EXPORT void ut_classify(int nech,
                                 const double *tab,
                                 double *sel,
                                 int nclass,
                                 double start,
                                 double pas,
                                 int *nmask,
                                 int *ntest,
                                 int *nout,
                                 int *classe);
GSTLEARN_EXPORT double ut_median(double *tab, int ntab);
GSTLEARN_EXPORT double ut_cnp(int n, int k);
GSTLEARN_EXPORT MatrixSquare ut_pascal(int ndim);
GSTLEARN_EXPORT int* ut_combinations(int n, int maxk, int *ncomb);
GSTLEARN_EXPORT void ut_shuffle_array(int nrow, int ncol, double *tab);

GSTLEARN_EXPORT VectorInt getListActiveToAbsolute(const VectorDouble &sel);
GSTLEARN_EXPORT std::map<int, int> getMapAbsoluteToRelative(const VectorDouble &sel,
                                                            bool verbose = false);
GSTLEARN_EXPORT int getRankMapAbsoluteToRelative(const std::map<int, int>& map, int iabs);
GSTLEARN_EXPORT int getRankMapRelativeToAbsolute(const std::map<int, int>& map, int irel);

typedef double (*operate_function)(double);
GSTLEARN_EXPORT operate_function operate_Identify( int oper );
GSTLEARN_EXPORT double operate_Identity(double x);
GSTLEARN_EXPORT double operate_Inverse(double x);
GSTLEARN_EXPORT double operate_Square(double x);
GSTLEARN_EXPORT double operate_InverseSquare(double x);
GSTLEARN_EXPORT double operate_Sqrt(double x);
GSTLEARN_EXPORT double operate_InverseSqrt(double x);
GSTLEARN_EXPORT double modifyOperator(const EOperator& oper, double oldval, double value);

GSTLEARN_EXPORT double roundZero(double value, double eps = EPSILON6);

GSTLEARN_EXPORT double truncateDecimals(double value, int ndec = 0);
GSTLEARN_EXPORT double truncateDigits(double value, int ndigits);

GSTLEARN_EXPORT void setInternalDebug(bool status);
GSTLEARN_EXPORT bool isInternalDebug();

GSTLEARN_EXPORT void print_range(const char* title, int ntab, const double* tab, const double* sel);

GSTLEARN_EXPORT void convertIndptrToIndices(int ncumul, const int* cumul, int* tab);
