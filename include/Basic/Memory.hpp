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

#include <cstdlib>

#define mem_free(tab)          mem_free_(__FILE__, __LINE__, tab)
#define mem_alloc(a, b)        mem_alloc_(__FILE__, __LINE__, a, b)
#define mem_realloc(tab, a, b) mem_realloc_(__FILE__, __LINE__, tab, a, b)

namespace gstlrn
{
GSTLEARN_EXPORT char* mem_alloc_(const char* call_file,
                                 size_t call_line,
                                 Id size,
                                 Id flag_fatal);
GSTLEARN_EXPORT char* mem_realloc_(const char* call_file,
                                   size_t call_line,
                                   char* tab,
                                   Id size,
                                   Id flag_fatal);
GSTLEARN_EXPORT char* mem_free_(const char* call_file,
                                size_t call_line,
                                char* tab);
GSTLEARN_EXPORT unsigned long long getTotalSystemMemory();
} // namespace gstlrn
