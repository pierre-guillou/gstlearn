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

#include "LinearOp/CGParam.hpp"

class ALinearOp;

/**
 * @brief Definition of the parameters used within SPDE
 *
 * refineK Discretization factor used for Kriging
 * refineS Discretization factor used for Simulation
 * border  Border size
 * nxmax   Maximum number of vertices in the internal mesh (0 : no limit)
 * epsNugget Nugget effect
 * useStencil Default option for no Cholesky (can only be used in stationary case for Turbo Meshing)
 * nMC Number of Monte-Carlo simulations (used for Variance and logdet)
 * seedMC Seed for the random number generator (used for Variance and logdet)
 *
 * cgparams Parameters for the Conjugate Gradient method
 */
class GSTLEARN_EXPORT SPDEParam {

public:
  SPDEParam(int refineK             = 11,
            int refineS             = 18,
            int border              = 8,
            bool flag_polarized     = true,
            int nxmax               = 300,
            double epsNugget        = EPSILON2,
            bool useStencil         = true,
            int nMC                 = 10,
            int seedMC              = 134341,
            const CGParam& cgparams = CGParam());
  SPDEParam(const SPDEParam& m);
  SPDEParam& operator=(const SPDEParam& m);
  virtual ~SPDEParam();

  static SPDEParam* create(int refineK             = 11,
                           int refineS             = 18,
                           int border              = 8,
                           bool flag_polarized     = true,
                           int nxmax               = 300,
                           double epsNugget        = EPSILON2,
                           bool useStencil         = true,
                           int nMC                 = 10,
                           int seedMC              = 134341,
                           const CGParam& cgparams = CGParam());

  int     getBorder() const { return _border; }
  CGParam getCGparams() const { return _CGparams; }
  double  getEpsNugget() const { return _epsNugget; }
  int     getRefineK() const { return _refineK; }
  int     getRefineS() const { return _refineS; }
  bool    isPolarized() const { return _flagPolarized; }
  void    setPolarized(bool flagPolarized) { _flagPolarized = flagPolarized; }
  int     getNxMax() const { return _nxmax; }
  bool    getUseStencil() const { return _useStencil; }
  int     getNMC() const { return _nMC; }
  int     getSeedMC() const { return _seedMC; }

  void setBorder(int border) { _border = border; }
  void setCGparams(const CGParam& CGparams) { _CGparams = CGparams; }
  void setEpsNugget(double epsNugget) { _epsNugget = epsNugget; }
  void setRefineK(int refineK) { _refineK = refineK; }
  void setRefineS(int refineS) { _refineS = refineS; }
  void setNxMax(int nxmax) { _nxmax = nxmax; }
  void setUseStencil(bool useStencil) { _useStencil = useStencil; }
  void setNMC(int nMC) { _nMC = nMC; }
  void setSeedMC(int seedMC) { _seedMC = seedMC; }

private:
  int     _refineK;
  int     _refineS;
  int     _border;
  bool    _flagPolarized;
  int     _nxmax;
  double  _epsNugget;
  bool    _useStencil;
  int     _nMC; // Number of Monte-Carlo simulations
  int     _seedMC; // Seed for the random number generator
  CGParam _CGparams;
};
