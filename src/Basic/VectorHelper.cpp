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
#include "geoslib_old_f.h"

#include "Basic/VectorHelper.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/AException.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/Law.hpp"
#include "Basic/OptCustom.hpp"
#include <string.h>
#include <algorithm>
#include <iomanip>
#include <ctime>
#include <cstdlib>
#include <math.h>
#include "../../include/Geometry/GeometryHelper.hpp"

VectorInt VectorHelper::initVInt(int nval, int value)
{
  VectorInt tab(nval, value);
  return tab;
}

VectorDouble VectorHelper::initVDouble(int nval, double value)
{
  VectorDouble tab(nval, value);
  return tab;
}

VectorVectorDouble VectorHelper::initVVDouble(int nval1, int nval2, double value)
{
  VectorVectorDouble tab(nval1, VectorDouble(nval2, value));
  return tab;
}

VectorVectorInt VectorHelper::initVVInt(int nval1, int nval2, int value)
{
  VectorVectorInt tab(nval1, VectorInt(nval2, value));
  return tab;
}

VectorInt VectorHelper::initVInt(int* values, int number)
{
  if (values == nullptr) return VectorInt();
  VectorInt vec(number);
  for (int i = 0; i < number; i++) vec[i] = values[i];
  return vec;
}

VectorDouble VectorHelper::initVDouble(double* values, int number)
{
  if (values == nullptr) return VectorDouble();
  VectorDouble vec;
  for (int i = 0; i < number; i++) vec[i] = values[i];
  return vec;
}

VectorVectorDouble VectorHelper::initVVDouble(double* value, int n1, int n2)
{
  if (value == nullptr) return VectorVectorDouble();
  VectorVectorDouble vec;
  vec.resize(n1);
  for (int i1 = 0; i1 < n1; i1++) vec[i1].resize(n2);

  int lec = 0;
  for (int i1 = 0; i1 < n1; i1++)
    for (int i2 = 0; i2 < n2; i2++)
      vec[i1][i2] = value[lec++];
  return vec;
}

void VectorHelper::display(const String &title, const VectorDouble &vect)
{
  if (!title.empty()) message("%s\n", title.c_str());
  messageFlush(VH::toString(vect));
}

void VectorHelper::display(const String &title, const VectorString &vect)
{
  if (!title.empty()) message("%s\n", title.c_str());
  messageFlush(VH::toString(vect));
}

void VectorHelper::display(const String &title, const VectorVectorDouble &vect)
{
  if (!title.empty()) message("%s\n", title.c_str());
  messageFlush(VH::toString(vect));
}

void VectorHelper::display(const String &title, const VectorInt &vect)
{
  if (!title.empty()) message("%s\n", title.c_str());
  messageFlush(VH::toString(vect));
}

String VectorHelper::toString(const VectorDouble &vec)
{
  return toVector(String(), vec);
}

String VectorHelper::toString(const VectorVectorDouble &vec)
{
  return toVector(String(), vec);
}

String VectorHelper::toString(const VectorString& vec)
{
  return toVector(String(), vec);
}

String VectorHelper::toString(const VectorInt &vec)
{
  return toVector(String(), vec);
}

void VectorHelper::displayStats(const String &title, const VectorDouble &vect)
{
  int ntotal = (int) vect.size();
  int number = 0;
  double mean = 0.;
  double stdv = 0.;
  double mini = 1.e30;
  double maxi = -1.e30;

  for (int i = 0; i < ntotal; i++)
  {
    double value = vect[i];
    if (FFFF(value)) continue;
    number++;
    mean += value;
    stdv += value * value;
    if (value < mini) mini = value;
    if (value > maxi) maxi = value;
  }

  if (!title.empty()) message("%s\n", title.c_str());
  if (number > 0)
  {
    mean /= (double) number;
    stdv = stdv / (double) number - mean * mean;
    stdv = (stdv > 0.) ? sqrt(stdv) : 0.;

    message("- Number of samples = %d / %d\n", number, ntotal);
    message("- Minimum  = %s\n", toDouble(mini).c_str());
    message("- Maximum  = %s\n", toDouble(maxi).c_str());
    message("- Mean     = %s\n", toDouble(mean).c_str());
    message("- St. Dev. = %s\n", toDouble(stdv).c_str());
  }
  else
  {
    message("No value defined\n");
  }
}

void VectorHelper::displayRange(const String &title, const VectorDouble &vect)
{
  int ntotal = (int) vect.size();
  int number = 0;
  double mini = 1.e30;
  double maxi = -1.e30;

  for (int i = 0; i < ntotal; i++)
  {
    double value = vect[i];
    if (FFFF(value)) continue;
    number++;
    if (value < mini) mini = value;
    if (value > maxi) maxi = value;
  }

  if (!title.empty()) message("%s\n", title.c_str());
  if (number > 0)
  {
    message("- Number of samples = %d / %d\n", number, ntotal);
    message("- Minimum  = %lf\n", mini);
    message("- Maximum  = %lf\n", maxi);
  }
  else
  {
    message("No value defined\n");
  }
}

void VectorHelper::displayRange(const String &title, const VectorInt &vect)
{
  int ntotal = (int) vect.size();
  int number = 0;
  int mini =  100000000;
  int maxi = -100000000;

  for (int i = 0; i < ntotal; i++)
  {
    int value = vect[i];
    if (FFFF(value)) continue;
    number++;
    if (value < mini) mini = value;
    if (value > maxi) maxi = value;
  }

  if (!title.empty()) message("%s\n", title.c_str());
  if (number > 0)
  {
    message("- Number of samples = %d / %d\n", number, ntotal);
    message("- Minimum  = %d\n", mini);
    message("- Maximum  = %d\n", maxi);
  }
  else
  {
    message("No value defined\n");
  }
}

double VectorHelper::maximum(const VectorDouble &vec, bool flagAbs)
{
  if (vec.size() <= 0) return 0.;
  double max = -1.e30;
  for (auto v : vec)
  {
    if (flagAbs) v = ABS(v);
    if (v > max) max = v;
  }
  return (max);
}

int VectorHelper::maximum(const VectorInt &vec)
{
  if (vec.size() <= 0) return 0;
  int max = -10000000;
  for (auto v : vec)
  {
    if (v > max) max = v;
  }
  return (max);
}

double VectorHelper::maximum(const VectorVectorDouble& vect, bool flagAbs)
{
  double val = VH::maximum(vect[0]);
  for (int i = 1; i < (int)vect.size(); i++)
  {
    val = MAX(val,  VH::maximum(vect[i], flagAbs));
  }
  return val;
}

int VectorHelper::minimum(const VectorInt &vec)
{
  if (vec.size() <= 0) return 0;
  int min = 10000000;
  for (auto v : vec)
  {
    if (v < min) min = v;
  }
  return (min);
}

double VectorHelper::minimum(const VectorDouble &vec, bool flagAbs)
{
  if (vec.size() <= 0) return 0.;
  double min = 1.e30;
  for (auto v : vec)
  {
    if (flagAbs) v = ABS(v);
    if (v < min) min = v;
  }
  return (min);
}

double VectorHelper::minimum(const VectorVectorDouble& vect, bool flagAbs)
{
  double val = VH::minimum(vect[0]);
  for (int i = 1; i < (int)vect.size(); i++)
  {
    val = MAX(val,  VH::minimum(vect[i], flagAbs));
  }
  return val;
}

double VectorHelper::mean(const VectorDouble &vec)
{
  if (vec.size() <= 0) return 0.;
  double mean = 0.;
  int number = 0;
  for (auto v : vec)
  {
    if (FFFF(v)) continue;
    mean += v;
    number++;
  }
  if (number > 0)
    mean /= (double) number;
  else
    mean = TEST;
  return (mean);
}

double VectorHelper::cumul(const VectorDouble& vec)
{
  double total = 0.;
  for (auto v : vec)
  {
    if (FFFF(v)) continue;
    total += v;
  }
  return total;
}

double VectorHelper::variance(const VectorDouble &vec)
{
  if (vec.size() <= 0) return 0.;
  double mean = 0.;
  double var = 0.;
  int number = 0;
  for (auto v : vec)
  {
    if (FFFF(v)) continue;
    var += v * v;
    mean += v;
    number++;
  }
  if (number <= 0) return TEST;
  mean /= (double) number;
  var = var / (double) number - mean * mean;
  return var;
}

double VectorHelper::correlation(const VectorDouble &veca, const VectorDouble& vecb)
{
  if (veca.size() <= 0 || vecb.size() <= 0 || veca.size() != vecb.size()) return 0.;

  double m1  = 0.;
  double m2  = 0.;
  double v11 = 0.;
  double v22 = 0.;
  double v12 = 0.;
  int number = 0;
  for (int i = 0; i < (int) veca.size(); i++)
  {
    double z1 = veca[i];
    double z2 = vecb[i];
    if (FFFF(z1) || FFFF(z2)) continue;
    v11 += z1 * z1;
    v22 += z2 * z2;
    v12 += z1 * z2;
    m1 += z1;
    m2 += z2;
    number++;
  }
  if (number <= 0) return TEST;
  m1 /= (double) number;
  m2 /= (double) number;
  v11 = v11 / (double) number - m1 * m1;
  v22 = v22 / (double) number - m2 * m2;
  v12 = v12 / (double) number - m1 * m2;
  if (v11 <= 0.) return TEST;
  if (v22 <= 0.) return TEST;
  double corr = v12 / sqrt(v11 * v22);
  return corr;
}

double VectorHelper::stdv(const VectorDouble &vec)
{
  double var = variance(vec);
  if (!FFFF(var))
    return (sqrt(var));
  else
    return TEST;
}

double VectorHelper::norm(const VectorDouble &vec)
{
  double ip = static_cast<double>(innerProduct(vec, vec));
  return sqrt(ip);
}

int VectorHelper::product(const VectorInt& vec)
{
  if (vec.empty()) return 0;
  int nprod = 1;
  for (int i = 0; i < (int) vec.size(); i++)
    nprod *= vec[i];
  return nprod;
}

double VectorHelper::product(const VectorDouble& vec)
{
  if (vec.empty()) return 0;
  double nprod = 1.;
  for (int i = 0; i < (int) vec.size(); i++)
    nprod *= vec[i];
  return nprod;
}

void VectorHelper::normalize(VectorDouble &vec)
{
  double ratio = VH::norm(vec);
  if (ratio <= 0.) return;
  for (auto v : vec)
     v /= ratio;
}

/**
 * Check if the contents of a vector is constant (equal to 'refval' is defined)
 * @param vect   Input vector
 * @param refval Reference value (TEST if not defined)
 * @return
 */
bool VectorHelper::isConstant(const VectorDouble& vect, double refval)
{
  if (vect.empty()) return false;
  if (FFFF(refval)) refval = vect[0];
  for (int i = 1; i < (int) vect.size(); i++)
    if (vect[i] != refval) return false;
  return true;
}

/**
 * Check if the contents of a vector is constant (equal to 'refval' is defined)
 * @param vect   Input vector
 * @param refval Reference value (ITEST if not defined)
 * @return
 */
bool VectorHelper::isConstant(const VectorInt& vect, int refval)
{
  if (vect.empty()) return false;
  if (IFFFF(refval)) refval = vect[0];
  for (int i = 1; i < (int) vect.size(); i++)
    if (vect[i] != refval) return false;
  return true;
}

bool VectorHelper::isSame(const VectorDouble &v1, const VectorDouble &v2, double eps)
{
  if (v1.size() != v2.size()) return false;
  for (int i = 0, n = static_cast<int>(v1.size()); i < n; i++)
    if (ABS(v1.at(i) - v2.at(i)) > eps) return false;
  return true;
}

bool VectorHelper::isSame(const VectorInt &v1, const VectorInt &v2)
{
  if (v1.size() != v2.size()) return false;
  for (int i = 0, n = static_cast<int>(v1.size()); i < n; i++)
    if (ABS(v1.at(i) - v2.at(i)) > 0) return false;
  return true;
}

void VectorHelper::fill(VectorDouble &vec, double value, int size)
{
  if (size > 0) vec.resize(size);
  std::fill(vec.begin(), vec.end(), value);
}

void VectorHelper::fill(VectorInt &vec, int value, int size)
{
  if (size > 0) vec.resize(size);
  std::fill(vec.begin(), vec.end(), value);
}

void VectorHelper::fill(VectorVectorDouble &vec, double value)
{
  for (auto &e : vec)
  {
    VH::fill(e, value);
  }
}

/**
 * Create an output vector containing the 'number' consecutive numbers starting from 'ideb'
 *
 * @param number  Length of the output vector
 * @param ideb    Index of the first element of the output vector
 */
VectorInt VectorHelper::sequence(int number, int ideb)
{
  VectorInt vec(number);

  for (int i = 0; i < number; i++)
    vec[i] = ideb + i;
  return vec;
}

/**
 * Create an output vector going from 'valFrom' to 'ValTo' by step of 'valStep'
 */
VectorDouble VectorHelper::sequence(double valFrom, double valTo, double valStep)
{
  VectorDouble vec;

  double value = valFrom;
  while (value <= valTo)
  {
    vec.push_back(value);
    value = value + valStep;
  }
  return vec;
}

VectorDouble VectorHelper::simulateUniform(int n, double mini, double maxi)
{
  VectorDouble vec(n);
  for (int i = 0; i < n; i++)
    vec[i] = law_uniform(mini, maxi);
  return vec;
}

VectorDouble VectorHelper::simulateBernoulli(int n, double proba, double vone, double velse)
{
  VectorDouble vec(n);
  for (int i = 0; i < n; i++)
  {
    double rand = law_uniform(0., 1.);
    if (rand < proba)
      vec[i] = vone;
    else
      vec[i] = velse;
  }
  return vec;
}

VectorDouble VectorHelper::simulateGaussian(int n, double mean, double sigma)
{
  VectorDouble vec(n);
  simulateGaussianInPlace(vec,mean,sigma);
  return vec;
}

void VectorHelper::simulateGaussianInPlace(VectorDouble &vect,
                                           double mean,
                                           double sigma)
{
  int n = (int) vect.size();
  for (int i = 0; i < n; i++)
    vect[i] = mean + sigma * law_gaussian();

}

VectorDouble VectorHelper::concatenate(const VectorDouble &veca,
                                       const VectorDouble &vecb)
{
  VectorDouble res = veca;
  for (auto &e: vecb)
    res.push_back(e);
  return res;
}

void VectorHelper::cumulate(VectorDouble &veca,
                            const VectorDouble &vecb,
                            double coeff,
                            double addval)
{
  if (veca.size() != vecb.size())
  my_throw("Wrong size");
  for (int i = 0, n = static_cast<int>(veca.size()); i < n; i++)
    veca[i] += coeff * vecb[i] + addval;
}

// random generator function:
int VectorHelper::_myrandom(int i)
{
  return std::rand()%i;
}

/**
 * Sample a set of 'ntotal' ranks
 * @param ntotal      Dimension to be sampled
 * @param proportion  Proportion of elected samples (in [0,1])
 * @param number      Number of elected samples
 * @param seed        Seed used for the random number generator
 * @return A vector of indices lying between 0 and ntotal-1. No duplicate.
 *
 * @remark If 'proportion' and 'number' are not specified,
 * @remark the output vector has dimension equal to 'ntotal'
 */
VectorInt VectorHelper::sampleRanks(int ntotal, double proportion, int number, int seed)
{
  if (proportion <= 0. && number <= 0) return VectorInt();

  // Find the number of expected values
  int count;
  if (proportion <= 0. && number <= 0)
    count = ntotal;
  else if (proportion > 0.)
    count = (int) (ntotal * proportion);
  else
    count = number;
  count = MIN(ntotal, MAX(1, count));

  // Define the Seed
  std::srand(seed);

  VectorInt ranks;
  for (int i = 0; i < ntotal; i++) ranks.push_back(i);

  std::random_shuffle ( ranks.begin(), ranks.end(), _myrandom);

  ranks.resize(count);
  return ranks;
}

VectorDouble VectorHelper::add(const VectorDouble &veca, const VectorDouble &vecb)
{
  VectorDouble res;
  if (veca.size() != vecb.size())
  my_throw("Wrong size");
  for (int i = 0, n = static_cast<int>(veca.size()); i < n; i++)
    res.push_back(veca.at(i) + vecb.at(i));
  return res;
}

/**
 * Performs: veca += vecb
 * @param dest Input/Output vector
 * @param src Auxiliary vector
 */
void VectorHelper::addInPlace(VectorDouble &dest, const VectorDouble &src)
{
  VectorDouble res;
  if (dest.size() != src.size())
  my_throw("Wrong size");
  for (int i = 0, n = static_cast<int>(dest.size()); i < n; i++)
    dest[i] += src[i];
}

void VectorHelper::addInPlace(const VectorDouble &veca,
                              const VectorDouble &vecb,
                              VectorDouble &res)
{
  if (veca.size() != vecb.size())
  {
    my_throw("Wrong size");
  }
  int n = (int) veca.size();
  if ((int) res.size() != n)
    res.resize(n);
  for (int i = 0; i < n; i++)
    res[i] = veca[i] + vecb[i];
}

/**
 * Return a vector containing vecb - veca
 * @param veca Input Vector
 * @param vecb Input Vector
 * @return
 */
VectorDouble VectorHelper::subtract(const VectorDouble &veca,
                                    const VectorDouble &vecb)
{
  VectorDouble res;
  if (veca.size() != vecb.size())
  my_throw("Wrong size");
  for (int i = 0, n = static_cast<int>(veca.size()); i < n; i++)
    res.push_back(vecb.at(i) - veca.at(i));
  return res;
}

/**
 * Performs: veca -= vecb
 * @param dest Input/Output vector
 * @param src Auxiliary vector
 */
void VectorHelper::subtractInPlace(VectorDouble &dest, const VectorDouble &src)
{
  VectorDouble res;
  if (dest.size() != src.size())
    my_throw("Wrong size");
  for (int i = 0, n = static_cast<int>(dest.size()); i < n; i++)
    dest[i] -= src[i];
}

void VectorHelper::multiplyInPlace(VectorDouble &vec, const VectorDouble &v)
{
  if (vec.size() != v.size())
    my_throw("Arguments 'vec' and 'v' should have same dimension");
  for (unsigned int i = 0; i < vec.size(); i++)
  {
    vec[i] *= v[i];
  }
}

void VectorHelper::divideInPlace(VectorDouble &vec, const VectorDouble &v)
{
  if (vec.size() != v.size())
    my_throw("Arguments 'vec' and 'v' should have same dimension");
  for (unsigned int i = 0; i < vec.size(); i++)
  {
    if (ABS(v[i]) < EPSILON20)
      my_throw("division by 0");
    vec[i] /= v[i];
  }
}

void VectorHelper::multiplyConstant(VectorDouble &vec, double v)
{
  std::for_each(vec.begin(), vec.end(), [v](double &d)
  { d *= v;});
}

void VectorHelper::divideConstant(VectorDouble &vec, double v)
{
  if (ABS(v) < EPSILON10)
  my_throw("division by 0");
  std::for_each(vec.begin(), vec.end(), [v](double &d)
  { d /= v;});
}

void VectorHelper::copy(VectorDouble &veca, const VectorDouble &vecb)
{
  if (veca.size() != vecb.size())
  my_throw("Wrong size");
  for (int i = 0, n = static_cast<int>(veca.size()); i < n; i++)
    veca[i] = vecb[i];
}

void VectorHelper::addConstant(VectorDouble &vec, double v)
{
  std::for_each(vec.begin(), vec.end(), [v](double &d)
  { d += v;});
}

void VectorHelper::addConstant(VectorInt &vec, int v)
{
  std::for_each(vec.begin(), vec.end(), [v](int &d)
  { d += v;});
}

VectorDouble VectorHelper::power(const VectorDouble &vec, double power)
{
  int size = static_cast<int>(vec.size());
  VectorDouble res(size);
  for (int i = 0; i < size; i++)
    res[i] = pow(vec[i], power);
  return res;
}

VectorDouble VectorHelper::inverse(const VectorDouble& vec)
{
  VectorDouble inv(vec.size());
  for (int i = 0; i < (int)vec.size();i++)
  {
    inv[i] = 1. / vec[i];
  }
  return inv;
}

int VectorHelper::countUndefined(const VectorDouble &vec)
{
  int count = 0;
  for (int i = 0; i < (int) vec.size(); i++)
  {
    if (FFFF(vec[i])) count++;
  }
  return count;
}

int VectorHelper::countDefined(const VectorDouble &vec)
{
  int count = 0;
  for (int i = 0; i < (int) vec.size(); i++)
  {
    if (! FFFF(vec[i])) count++;
  }
  return count;
}

/**
 * Calculate the diagonal of the box extension
 * @param mini Array of lower coordinates of the box
 * @param maxi Array of upper coordinates of the box
 * @return
 * @remark If one coordinate is undefined, TEST is returned.
 */
double VectorHelper::extensionDiagonal(const VectorDouble &mini,
                                       const VectorDouble &maxi)
{
  double diag = 0.;
  VectorDouble delta = VH::subtract(mini, maxi);
  int ndim = (int) delta.size();
  for (int idim = 0; idim < ndim; idim++)
  {
    double dval = delta[idim];
    if (FFFF(dval)) return TEST;
    diag += dval * dval;
  }
  diag = sqrt(diag);
  return diag;
}

VectorInt VectorHelper::sort(const VectorInt& vecin, bool ascending)
{
  if (vecin.empty()) return VectorInt();

  VectorInt vecout = vecin;
  std::sort(vecout.begin(), vecout.end());
  if (! ascending)
    std::reverse(vecout.begin(), vecout.end());
  return vecout;
}

VectorDouble VectorHelper::sort(const VectorDouble& vecin, bool ascending)
{
  if (vecin.empty()) return VectorDouble();

  VectorDouble vecout = vecin;
  std::sort(vecout.begin(), vecout.end());
  if (! ascending)
    std::reverse(vecout.begin(), vecout.end());
  return vecout;
}

VectorInt VectorHelper::sortRanks(const VectorDouble& vecin)
{
  if (vecin.empty()) return VectorInt();

  VectorInt idx(vecin.size());
  for (int i = 0; i < (int) vecin.size(); i++) idx[i] = i;

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values
  stable_sort(idx.begin(), idx.end(),
       [&vecin](size_t i1, size_t i2) {return vecin[i1] < vecin[i2];});

  return idx;
}

std::pair<double,double> VectorHelper::rangeVals(const VectorDouble& vec)
{
  std::pair<double,double> res(vec[0],vec[0]);
  for (int i = 1; i < (int)vec.size(); i++)
  {
    res.first  = MIN(res.first ,vec[i]);
    res.second = MAX(res.second,vec[i]);
  }
  return res;
}

double VectorHelper::innerProduct(const VectorDouble &veca,
                                  const VectorDouble &vecb)
{
  if (veca.size() != vecb.size())
  my_throw("Wrong size");
  double prod = 0.;
  for (int i = 0, n = static_cast<int>(veca.size()); i < n; i++)
    prod += veca.at(i) * vecb.at(i);
  return prod;
}

/**
 * Cross product (limited to 3D)
 * @param veca First vector
 * @param vecb Second Vector
 * @return
 */
VectorDouble VectorHelper::crossProduct(const VectorDouble &veca,
                                     const VectorDouble &vecb)
{
  if (veca.size() != vecb.size())
  my_throw("Wrong size");
  VectorDouble res;
  res.push_back(veca[1] * vecb[2] - veca[2] * vecb[1]);
  res.push_back(veca[2] * vecb[0] - veca[0] * vecb[2]);
  res.push_back(veca[0] * vecb[1] - veca[1] * vecb[0]);
  return res;
}