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

#include "../Simulation/BooleanObject.hpp"
#include "Basic/Vector.hpp"
#include "Basic/AStringable.hpp"

class AShape;

class GSTLEARN_EXPORT ModelBoolean: public AStringable
{
public:
  ModelBoolean(double thetaCst = 1., bool flagStat = true);
  ModelBoolean(const ModelBoolean &r);
  ModelBoolean& operator=(const ModelBoolean &r);
  virtual ~ModelBoolean();

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  int getNbTokens() const { return (int) _shapes.size(); }
  void addToken(const AShape& token);
  void normalizeProportions();
  BooleanObject* generateObject(int ndim) const;
  const AShape* getToken(int itok) const { return _shapes[itok]; }

  bool   isFlagStat() const { return _flagStat; }
  double getThetaCst() const { return _thetaCst; }

  void setFlagStat(bool flagStat) { _flagStat = flagStat; }
  void setThetaCst(double thetaCst) { _thetaCst = thetaCst; }

private:
  bool   _flagStat;
  double _thetaCst;
  std::vector<AShape*> _shapes; // List of the Token
};