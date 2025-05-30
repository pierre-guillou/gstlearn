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

#include "Enum/ELoc.hpp"


/**
 * Gives the Rank of the Attribute for:
 * - a given pointer type
 * - a given locatorIndex rank
 *
 * The dimension of the internal vector is equal to the number of items for a
 * given pointer type
 */
class GSTLEARN_EXPORT PtrGeos {
public:
  bool isLocatorIndexValid(int locatorIndex) const;
  int  getLocatorByIndex(int locatorIndex) const { return _r[locatorIndex]; }
  void setLocatorByIndex(int locatorIndex, int value) { _r[locatorIndex] = value; }
  int  getNLoc() const { return static_cast<int>(_r.size()); }
  bool hasLocator() const { return ! _r.empty(); }
  int  findUIDInLocator(int iuid) const;
  void erase(int locatorIndex);
  void clear();
  void resize(int count) { _r.resize(count,0); }
  String dumpLocator(int rank, const ELoc& locatorType) const;

private:
  VectorInt _r;    /* Rank of the attribute */
};

GSTLEARN_EXPORT int getLocatorTypeFromName(const String& name_type);
GSTLEARN_EXPORT int locatorIdentify(String string,
                                    ELoc* ret_locatorType,
                                    int* ret_locatorIndex,
                                    int* ret_mult);
GSTLEARN_EXPORT bool isLocatorTypeValid(const ELoc& locatorType,
                                        bool unknownValid = false);
GSTLEARN_EXPORT String getLocatorName(const ELoc& locatorType, int locatorIndex=1);
GSTLEARN_EXPORT void   printLocatorList();
GSTLEARN_EXPORT VectorString getLocatorNames();
GSTLEARN_EXPORT VectorInt    getLocatorMultiples();
