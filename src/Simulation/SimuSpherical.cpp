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

#include "Db/DbGrid.hpp"
#include "Db/Db.hpp"
#include "Mesh/MeshSpherical.hpp"
#include "Model/Model.hpp"
#include "Simulation/ASimulation.hpp"
#include "Simulation/SimuSpherical.hpp"
#include "Simulation/SimuSphericalParam.hpp"
#include "Basic/Vector.hpp"
#include "Basic/Law.hpp"

#include <math.h>

#define IPTR(ix,iy)        ((iy) * nx + (ix))
#define DISCRET(idisc)     (GV_PI * (0.5 + (idisc)) / ((double) ndisc))

SimuSpherical::SimuSpherical(int nbsimu, int seed)
    : ASimulation(nbsimu, seed),
      AStringable()
{
}

SimuSpherical::SimuSpherical(const SimuSpherical &r)
    : ASimulation(r),
      AStringable(r)
{
}

SimuSpherical& SimuSpherical::operator=(const SimuSpherical &r)
{
  if (this != &r)
  {
    ASimulation::operator=(r);
    AStringable::operator =(r);
  }
  return *this;
}

SimuSpherical::~SimuSpherical()
{
}

String SimuSpherical::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  return sstr.str();
}

int SimuSpherical::simulate(DbGrid *db,
                            Model *model,
                            const SimuSphericalParam& sphepar,
                            int iptr,
                            bool verbose)
{
  int nbf = sphepar.getNbf();
  int special = sphepar.getSpecial();
  int degmax = sphepar.getDegmax();
  int nx = db->getNX(0);
  int ny = db->getNX(1);
  int nech = db->getSampleNumber();
  law_set_random_seed(getSeed());

  /* Core allocation */

  VectorDouble phase(nbf);
  VectorInt degree(nbf);
  VectorInt order(nbf);

  /* Define the spectrum */

  VectorDouble freqs;
  if (special == 1)
    freqs = _spectrum_chentsov(sphepar);
  else if (special == 2)
    freqs = _spectrum_exponential(model, sphepar);
  else
    freqs = _spectrum_any(model, sphepar);
  if (freqs.empty()) return 1;

  /* Optional printout */

  if (verbose)
  {
    sphepar.display();
    message("Random generation seed    = %d\n", law_get_random_seed());
    message("Number of frequencies     = %d\n", (int) freqs.size());
  }
  _spectrum_normalize(verbose, freqs);

  /* Get the model ingredients (generated by separate flows in order to  */
  /* avoid intermeshing dependencies) */

  for (int ibf = 0; ibf < nbf; ibf++)
  {
    degree[ibf] = _gdiscrete(freqs);
    if (degmax > 0) degree[ibf] = MIN(degmax, degree[ibf]);
  }
  for (int ibf = 0; ibf < nbf; ibf++)
    order[ibf] = law_int_uniform(-degree[ibf], degree[ibf]);
  for (int ibf = 0; ibf < nbf; ibf++)
    phase[ibf] = law_uniform(0., 2. * GV_PI);
  if (_check_degree_order(sphepar, freqs, degree, order, verbose)) return 1;

  /* Loop on the samples of the Data Base */
  /* We benefit in writing this as a double loop on coordinates */
  /* as some calculations can be factorized as they only concern latitude */

  int ecr = 0;
  int ntot = nbf * nech;
  for (int iy = 0; iy < ny; iy++)
  {
    double theta = ut_deg2rad(db->getCoordinate(IPTR(0, iy), 1) + 90.); // Latitude[-90,90]
    for (int ibf = 0; ibf < nbf; ibf++)
    {
      double degree_loc = degree[ibf];
      double order_loc = order[ibf];
      double phase_loc = phase[ibf];
      double t1 = ut_flegendre(1, (int) degree_loc, (int) order_loc, theta);
      for (int ix = 0; ix < nx; ix++, ecr++)
      {
        int jech = IPTR(ix, iy);
        mes_process("Simulation on Sphere", ntot, ecr);
        if (!db->isActive(jech)) continue;
        double phi = ut_deg2rad(db->getCoordinate(jech, 0));       // Longitude [0,360]
        double t2 = cos(phi * order_loc + phase_loc);
        db->updArray(jech, iptr, 0, t1 * t2);
      }
    }
  }

  /* Final normalization */

  double val = 2. / sqrt((double) nbf);
  for (int iech = 0; iech < nech; iech++)
  {
    if (db->isActive(iech)) db->updArray(iech, iptr, 3, val);
  }
  return 0;
}

VectorDouble SimuSpherical::simulate_mesh(MeshSpherical *mesh,
                                          Model *model,
                                          const SimuSphericalParam &sphepar,
                                          bool verbose)
{
  int nbf = sphepar.getNbf();
  int special = sphepar.getSpecial();
  int degmax = sphepar.getDegmax();
  int nech = mesh->getNApices();
  law_set_random_seed(getSeed());

  /* Core allocation */

  VectorDouble phase(nbf);
  VectorInt degree(nbf);
  VectorInt order(nbf);
  VectorDouble simu;

  /* Define the spectrum */

  VectorDouble freqs;
  if (special == 1)
    freqs = _spectrum_chentsov(sphepar);
  else if (special == 2)
    freqs = _spectrum_exponential(model, sphepar);
  else
    freqs = _spectrum_any(model, sphepar);
  if (freqs.empty()) return simu;

  /* Optional printout */

  if (verbose)
  {
    sphepar.display();
    message("Random generation seed    = %d\n", law_get_random_seed());
    message("Number of frequencies     = %d\n", (int) freqs.size());
  }
  _spectrum_normalize(verbose, freqs);

  /* Get the model ingredients (generated by separate flows in order to  */
  /* avoid intermeshing dependencies) */

  for (int ibf = 0; ibf < nbf; ibf++)
  {
    degree[ibf] = _gdiscrete(freqs);
    if (degmax > 0) degree[ibf] = MIN(degmax, degree[ibf]);
  }
  for (int ibf = 0; ibf < nbf; ibf++)
    order[ibf] = law_int_uniform(-degree[ibf], degree[ibf]);
  for (int ibf = 0; ibf < nbf; ibf++)
    phase[ibf] = law_uniform(0., 2. * GV_PI);
  if (_check_degree_order(sphepar, freqs, degree, order, verbose)) return simu;

  /* Loop on the samples of the Data Base */
  /* We benefit in writing this as a double loop on coordinates */
  /* as some calculations can be factorized as they only concern latitude */

  simu.resize(nech, 0.);
  for (int iech = 0; iech < nech; iech++)
  {
    double theta = ut_deg2rad(mesh->getApexCoor(iech, 1) + 90.); // Latitude [-90,90]
    for (int ibf = 0; ibf < nbf; ibf++)
    {
      double degree_loc = degree[ibf];
      double order_loc = order[ibf];
      double phase_loc = phase[ibf];
      double t1 = ut_flegendre(1, (int) degree_loc, (int) order_loc, theta);

      double phi = ut_deg2rad(mesh->getApexCoor(iech, 0));  // Longitude [0,360]
      double t2 = cos(phi * order_loc + phase_loc);
      simu[iech] += t1 * t2;
    }
  }

  /* Final normalization */

  double val = 2. / sqrt((double) nbf);
  for (int iech = 0; iech < nech; iech++)
    simu[iech] /= val;

  return simu;
}

/*****************************************************************************/
/*!
 **  Generate the spectrum for Chentsov
 **
 ** \returns The array 'freqs' or NULL
 **
 ** \param[in]  sphepar SimuSphericalParam structure
 **
 *****************************************************************************/
VectorDouble SimuSpherical::_spectrum_chentsov(const SimuSphericalParam& sphepar)
{
  VectorDouble freqs;
  int ifreq = 0;
  double total = 0.;

  /* Loop on the spectrum items */

  freqs.push_back(0.);            // ifreq = 0
  total += freqs[ifreq];
  ifreq++;

  freqs.push_back(0.75);          // ifreq = 1
  total += freqs[ifreq];
  ifreq++;

  while (1)
  {
    freqs.push_back(0.);
    ifreq++;

    double ratio = ((double) (ifreq - 2.)) / ((double) (ifreq + 1.));
    double value = freqs[ifreq - 2] * (2. * ifreq + 1.) / (2. * ifreq - 3.);
    value *= ratio * ratio;
    freqs.push_back(value);
    total += freqs[ifreq];
    ifreq++;

    if (ABS(1. - total) < sphepar.getTol()) break;
    if (sphepar.getNfmax() > 0 && ifreq >= sphepar.getNfmax()) break;
  }
  return (freqs);
}

/*****************************************************************************/
/*!
 **  Generate the spectrum for Exponential
 **
 ** \returns The array 'freqs' or NULL
 **
 ** \param[in]  model   Model (used for its range)
 ** \param[in]  sphepar SimuSphericalParam structure
 **
 *****************************************************************************/
VectorDouble SimuSpherical::_spectrum_exponential(Model *model,
                                                  const SimuSphericalParam& sphepar)
{
  VectorDouble freqs;
  int ifreq = 0;
  double total = 0.;
  double fcs = 1. / model->getCova(0)->getTheoretical();
  double fcs2 = fcs * fcs;
  double expfc = exp(-fcs * GV_PI);

  /* Core allocation */

  double value = (1. + expfc) / (fcs2 + 1.) / 2.;
  if (value < 0.) value = 0.;
  freqs.push_back(value);
  total += freqs[ifreq];
  ifreq++;

  value = 3. * (1. - expfc) / (fcs2 + 4.) / 2.;
  if (value < 0.) value = 0.;
  freqs.push_back(value);
  total += freqs[ifreq];
  ifreq++;

  /* Loop on the spectrum items */

  while (1)
  {
    double r1 = ifreq + 1.;
    double r2 = ifreq - 2.;
    value  = freqs[ifreq - 2] * (2. * ifreq + 1.) / (2. * ifreq - 3.);
    value *= (fcs2 + r2 * r2) / (fcs2 + r1 * r1);
    freqs.push_back(value);
    total += freqs[ifreq];
    ifreq++;

    if (ABS(1. - total) < sphepar.getTol()) break;
    if (sphepar.getNfmax() > 0 && ifreq >= sphepar.getNfmax()) break;
  }

  return (freqs);
}

/*****************************************************************************/
/*!
 **  Generate the spectrum for any covariance
 **
 ** \returns The array 'freqs' or NULL
 **
 ** \param[in]  model   Model (defined in Euclidean space) to be used
 ** \param[in]  sphepar SimuSphericalParam structure
 **
 *****************************************************************************/
VectorDouble SimuSpherical::_spectrum_any(Model *model,
                                          const SimuSphericalParam& sphepar)
{
  int ndisc = sphepar.getNdisc();
  VectorDouble freqs;
  VectorDouble dd(2);
  dd[0] = dd[1] = 0.;
  int ifreq = 0;
  double dincr = GV_PI / ((double) sphepar.getNdisc());
  VectorDouble covs(ndisc);

  /* Calculate the discretized covariance values */

  for (int idisc = 0; idisc < ndisc; idisc++)
  {
    double alpha = DISCRET(idisc);
    dd[0] = 2. * sin(alpha / 2.);
    double ca = 0.;
    for (int icova = 0; icova < model->getCovaNumber(); icova++)
      ca += model_calcul_basic(model, icova, ECalcMember::LHS, dd);
    covs[idisc] = ca;
  }

  /* Loop on the frequencies */

  double total = 0.;
  while (1)
  {
    /* Discretization of the frequency item */

    double an = 0.;
    for (int idisc = 0; idisc < ndisc; idisc++)
    {
      double alpha = DISCRET(idisc);
      double cosa = cos(alpha);
      double sina = sin(alpha);
      an += covs[idisc] * sina * ut_legendre(1, ifreq, cosa);
    }
    double value = an * dincr * sqrt((2. * ifreq + 1) / 2.);
    freqs.push_back(value);
    total += freqs[ifreq];
    ifreq++;

    if (total > 1.)
    {
      ifreq--;
      break;
    }
    if (ABS(1. - total) < sphepar.getTol()) break;
    if (sphepar.getNfmax() > 0 && ifreq >= sphepar.getNfmax()) break;
  }
  return freqs;
}

/*****************************************************************************/
/*!
 **  Normalize the spectrum
 **
 ** \param[in]  verbose Verbose flag
 **
 ** \param[out] freqs   Array of frequencies
 **
 *****************************************************************************/
void SimuSpherical::_spectrum_normalize(int verbose, VectorDouble& freqs)
{
  int nfreq = freqs.size();
  double totpos = 0.;
  double totneg = 0.;
  for (int ifreq = 0; ifreq < nfreq; ifreq++)
  {
    if (freqs[ifreq] < 0)
    {
      totneg -= freqs[ifreq];
      freqs[ifreq] = 0.;
    }
    else
    {
      totpos += freqs[ifreq];
    }
  }

  for (int ifreq = 0; ifreq < nfreq; ifreq++)
    freqs[ifreq] /= totpos;

  /* Printout (optional) */

  if (verbose)
  {
    message("Cumulated Spectrum        = %lf\n", totpos);
    message("Sum of negative weights   = %lf\n", totneg);
  }
}

/*****************************************************************************/
/*!
 **  Simulates the discrete distribution on 0,1,...,n
 **
 ** \param[in]  freqs  Array of frequencies (which add up to 1)
 **
 *****************************************************************************/
int SimuSpherical::_gdiscrete(VectorDouble& freqs)
{
  int nfreq = freqs.size();
  double u = law_uniform(0., 1.);

  double partvec = 0.;
  for (int ifreq = 0; ifreq < nfreq; ifreq++)
  {
    partvec += freqs[ifreq];
    if (u < partvec) return (ifreq);
  }
  return (nfreq - 1);
}

/*****************************************************************************/
/*!
 **  Check the degrees and orders
 **
 ** \return Error code
 **
 ** \param[in]  sphepar SimuSphericalParm structure
 ** \param[in]  freqs   Vector of frequencies
 ** \param[in]  degree  Array of degrees
 ** \param[in]  order   Array of orders
 ** \param[in]  verbose Verbose flag
 **
 *****************************************************************************/
int SimuSpherical::_check_degree_order(const SimuSphericalParam& sphepar,
                                       const VectorDouble& freqs,
                                       VectorInt& degree,
                                       VectorInt& order,
                                       int verbose)
{
  int nbf = sphepar.getNbf();
  int degmax = 0;
  int nfreq  = (int) freqs.size();
  int ordmin =  nfreq;
  int ordmax = -nfreq;

  for (int ibf = 0; ibf < nbf; ibf++)
  {
    if (degree[ibf] > degmax) degmax = degree[ibf];
    if (order[ibf] < ordmin) ordmin = order[ibf];
    if (order[ibf] > ordmax) ordmax = order[ibf];
    if (order[ibf] < -degree[ibf] || order[ibf] > +degree[ibf])
    {
      messerr("Order(%d) must lie in [-degree;+degree] where degree=%d",
              order[ibf], degree[ibf]);
      return (1);
    }
  }

  if (verbose)
  {
    message("Maximum degree            = %d\n", degmax);
    message("Minimum order             = %d\n", ordmin);
    message("Maximum order             = %d\n", ordmax);
  }
  return (0);
}
