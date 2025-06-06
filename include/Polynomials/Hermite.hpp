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

#include "Basic/VectorNumT.hpp"
#include "Matrix/MatrixSquare.hpp"

GSTLEARN_EXPORT VectorDouble hermitePolynomials(double y,
                                                double r,
                                                int nbpoly);
GSTLEARN_EXPORT VectorDouble hermitePolynomials(double y,
                                                double r,
                                                const VectorInt& ifacs);
GSTLEARN_EXPORT VectorDouble hermiteCoefIndicator(double yc, int nbpoly);
GSTLEARN_EXPORT VectorDouble hermiteCoefMetal(double yc,
                                              const VectorDouble &phi);
GSTLEARN_EXPORT VectorDouble hermiteCoefLower(double y, int nbpoly);
GSTLEARN_EXPORT VectorDouble hermiteIndicatorLower(double y, int nbpoly);
GSTLEARN_EXPORT MatrixSquare hermiteIncompleteIntegral(double yc,
                                                              int nbpoly);
GSTLEARN_EXPORT VectorDouble hermiteLognormal(double mean,
                                              double sigma,
                                              int nbpoly);
GSTLEARN_EXPORT double hermiteSeries(const VectorDouble &an,
                                     const VectorDouble &hn);

GSTLEARN_EXPORT VectorDouble hermiteIndicator(double yc,
                                              VectorDouble krigest,
                                              VectorDouble krigstd);
GSTLEARN_EXPORT double hermiteIndicatorElement(double yc,
                                               double krigest,
                                               double krigstd);
GSTLEARN_EXPORT VectorDouble hermiteIndicatorStd(double yc,
                                                 VectorDouble krigest,
                                                 VectorDouble krigstd);
GSTLEARN_EXPORT double hermiteIndicatorStdElement(double yc,
                                                  double krigest,
                                                  double krigstd);
GSTLEARN_EXPORT VectorDouble hermiteMetal(double yc,
                                          VectorDouble krigest,
                                          VectorDouble krigstd,
                                          const VectorDouble &phi);
GSTLEARN_EXPORT double hermiteMetalElement(double yc,
                                           double krigest,
                                           double krigstd,
                                           const VectorDouble &phi);
GSTLEARN_EXPORT VectorDouble hermiteMetalStd(double yc,
                                             VectorDouble krigest,
                                             VectorDouble krigstd,
                                             const VectorDouble &phi);
GSTLEARN_EXPORT double hermiteMetalStdElement(double yc,
                                              double krigest,
                                              double krigstd,
                                              const VectorDouble &phi);
GSTLEARN_EXPORT VectorDouble hermiteCondExp(VectorDouble krigest,
                                            VectorDouble krigstd,
                                            const VectorDouble &phi);
GSTLEARN_EXPORT double hermiteCondExpElement(double krigest,
                                             double krigstd,
                                             const VectorDouble &phi);
GSTLEARN_EXPORT VectorDouble hermiteCondStd(VectorDouble krigest,
                                            VectorDouble krigstd,
                                            const VectorDouble &phi);
GSTLEARN_EXPORT double hermiteCondStdElement(double krigest,
                                             double krigstd,
                                             const VectorDouble &phi);
