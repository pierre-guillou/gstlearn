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

#include "Basic/Vector.hpp"
#include "Boolean/ETShape.hpp"
#include "Boolean/AToken.hpp"
#include "Boolean/TokenParameter.hpp"

class Object;

class GSTLEARN_EXPORT TokenHalfSinusoid: public AToken
{
public:
  TokenHalfSinusoid(double proportion = 1.,
                    double period = 10.,
                    double amplitude = 1.,
                    double thickness = 1.,
                    double xext = 1.,
                    double zext = 1.,
                    double theta = 0.);
  TokenHalfSinusoid(const TokenHalfSinusoid &r);
  TokenHalfSinusoid& operator=(const TokenHalfSinusoid &r);
  virtual ~TokenHalfSinusoid();

  /// Interface for Iclonable
  virtual IClonable* clone() const override { return new TokenHalfSinusoid(*this); };

  ETShape getType() const override { return ETShape::HALFSINUSOID; }
  int  getNParams() const override { return 6; }
  bool getFlagCutZ() const override { return true; }
  Object* generateObject(int ndim = 3) override;
  bool belongObject(const VectorDouble& coor, const Object* object) const override;
};