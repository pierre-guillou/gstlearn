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

#include "Simulation/ASimulation.hpp"
#include "Basic/AStringable.hpp"

class SimuSphericalParam;
class MeshSpherical;

class GSTLEARN_EXPORT SimuSpherical: public ASimulation, public AStringable
{
public:
  SimuSpherical(int nbsimu = 1, int seed = 4324324);
  SimuSpherical(const SimuSpherical &r);
  SimuSpherical& operator=(const SimuSpherical &r);
  virtual ~SimuSpherical();

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  int simulate(DbGrid *db,
               Model *model,
               const SimuSphericalParam& sphepar,
               int iptr,
               bool verbose = false);

  VectorDouble simulate_mesh(MeshSpherical *mesh,
                             Model *model,
                             const SimuSphericalParam &sphepar,
                             bool verbose = false);

private:
  VectorDouble _spectrum_chentsov(const SimuSphericalParam& sphepar);
  VectorDouble _spectrum_exponential(Model *model, const SimuSphericalParam& sphepar);
  VectorDouble _spectrum_any(Model *model, const SimuSphericalParam& sphepar);
  void _spectrum_normalize(int verbose, VectorDouble& freqs);
  int _gdiscrete(VectorDouble& freqs);
  int _check_degree_order(const SimuSphericalParam& sphepar,
                          const VectorDouble& freqs,
                          VectorInt& degree,
                          VectorInt& order,
                          int verbose);
};