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

#include "geoslib_define.h"
#include "gstlearn_export.hpp"

namespace gstlrn
{
GSTLEARN_EXPORT void ascii_study_define(const char* study);
GSTLEARN_EXPORT void ascii_environ_read(char* file_name, bool verbose);
GSTLEARN_EXPORT void ascii_filename(const char* type, Id rank, Id mode, char* filename);
GSTLEARN_EXPORT void ascii_simu_read(char* file_name, bool verbose, Id* nbsimu, Id* nbtuba, Id* seed);
GSTLEARN_EXPORT Id ascii_option_defined(const char* file_name,
                                         bool verbose,
                                         const char* option_name,
                                         Id type,
                                         void* answer);
} // namespace gstlrn
