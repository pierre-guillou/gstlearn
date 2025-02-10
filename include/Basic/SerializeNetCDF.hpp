/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2025) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "geoslib_define.h"

#include "Basic/AStringable.hpp"
#include "Basic/VectorT.hpp"

#include <ncDim.h>
#include <ncFile.h>
#include <ncGroup.h>
#include <ncType.h>
#include <ncVar.h>

#include <optional>

namespace SerializeNetCDF
{
  /**
   * @brief Map values to corresponding netCDF types
   */
  inline netCDF::NcType getNetCDFType([[maybe_unused]] const int a)
  {
    return netCDF::NcType::nc_INT;
  }
  inline netCDF::NcType getNetCDFType([[maybe_unused]] const double a)
  {
    return netCDF::NcType::nc_DOUBLE;
  }
  inline netCDF::NcType getNetCDFType([[maybe_unused]] const long a)
  {
    return netCDF::NcType::nc_INT64;
  }
  inline netCDF::NcType getNetCDFType([[maybe_unused]] const std::string& a)
  {
    return netCDF::NcType::nc_STRING;
  }

  /**
   * @brief Open netCDF file in read mode, check metadata
   */
  inline netCDF::NcFile _fileOpenRead(const String& fname)
  {
    netCDF::NcFile file {fname, netCDF::NcFile::FileMode::read};
    auto metadata = file.getGroup("gstlearn metadata");
    if (metadata.isNull())
    {
      messerr("File %s doesn't contain Gstlearn metadata…", fname.c_str());
    }
    else
    {
      const auto att = metadata.getAtt("Format version");
      std::string version;
      att.getValues(version);
      if (version != "1.0.0")
      {
        messerr("File %s has format version %s, expected 1.0.0", fname.c_str(),
                version.c_str());
      }
    }
    file.close();
    // cannot return file directly as NcFile doesn't have a move constructor
    return netCDF::NcFile {fname, netCDF::NcFile::FileMode::read};
  }

  /**
   * @brief Open netCDF file in write mode, write metadata
   */
  inline netCDF::NcFile _fileOpenWrite(const String& fname)
  {
    netCDF::NcFile file {fname, netCDF::NcFile::FileMode::replace};
    auto metadata = file.addGroup("gstlearn metadata");
    metadata.putAtt("Description",
                    "This file is used to serialize gstlearn's internal data structures");
    metadata.putAtt("Format version", "1.0.0");
    // cannot return file directly as NcFile doesn't have a move constructor
    file.close();
    return netCDF::NcFile {fname, netCDF::NcFile::FileMode::write};
  }

  /**
   * @brief Read a named netCDF variable into a generic VectorT
   *
   * @param[in] grp netCDF group containing the variable to read
   * @param[in] title Name of netCDF variable to read
   * @param[out] vec Vector to be filled with netCDF data (will be resized)
   * @return true if success
   */
  template<typename T>
  bool _readVec(const netCDF::NcGroup& grp, const String& title, VectorT<T>& vec);
  template<>

  /**
   * @brief Read a named netCDF string variable into a VectorString
   *
   * @param[in] grp netCDF group containing the string variable to read
   * @param[in] title Name of netCDF string variable to read
   * @param[out] vec VectorString to be filled with netCDF data (will be resized)
   * @return true if success
   */
  inline bool
  _readVec(const netCDF::NcGroup& grp, const String& title, VectorString& vec);

  /**
   * @brief Write a generic VectorT into a named netCDF variable
   *
   * @param[in,out] grp netCDF group containing the variable to write
   * @param[in] title Name of netCDF variable to write
   * @param[in] vec Vector to write into netCDF
   * @param[in] dim netCDF dimension of the vector to write
   * @return true if success
   */
  template<typename T>
  bool _writeVec(netCDF::NcGroup& grp,
                 const String& title,
                 const VectorT<T>& vec,
                 const netCDF::NcDim& dim);

  /**
   * @brief Write a VectorString into a named netCDF string variable
   *
   * @param[in,out] grp netCDF group containing the string variable to write
   * @param[in] title Name of netCDF string variable to write
   * @param[in] vec VectorString to write to netCDF
   * @param[in] dim netCDF dimension of the vector to write
   * @return true if success
   */
  template<>
  inline bool _writeVec(netCDF::NcGroup& grp,
                        const String& title,
                        const VectorString& vec,
                        const netCDF::NcDim& dim);

  /**
   * @brief Extract a group inside a parent group
   *
   * @param[in] parent Parent group
   * @param[in] name Name of the group to find in parent
   * @return Group if found else nullopt
   */
  inline std::optional<netCDF::NcGroup>
  getGroup(const netCDF::NcGroup& parent, const String& name)
  {
    auto grp = parent.getGroup(name);
    if (grp.isNull())
    {
      messerr("Cannot find group %s in parent group %s", name.c_str(),
              parent.getName().c_str());
      return std::nullopt;
    }
    return grp;
  }

}; // namespace SerializeNetCDF

template<typename T>
bool SerializeNetCDF::_readVec(const netCDF::NcGroup& grp,
                               const String& title,
                               VectorT<T>& vec)
{

  if (grp.isNull())
  {
    messerr("Cannot read vector %s: container group is null", title.c_str());
    return false;
  }

  const auto var = grp.getVar(title);
  if (var.isNull())
  {
    messerr("Cannot read NetCDF Variable of name %s in group %s", title.c_str(),
            grp.getName().c_str());
    return false;
  }

  // Assume variable is of dim 1
  if (var.getDimCount() != 1)
  {
    messerr("NetCDF Variable of name %s in group %s has %ul dims, but we expect only 1",
            title.c_str(), grp.getName().c_str(), var.getDimCount());
    return false;
  }

  const auto dim = var.getDim(0);
  vec.resize(dim.getSize());

  var.getVar(vec.data());

  return true;
}

template<>
bool SerializeNetCDF::_readVec(const netCDF::NcGroup& grp,
                               const String& title,
                               VectorString& vec)
{

  if (grp.isNull())
  {
    messerr("Cannot read string vector %s: container group is null", title.c_str());
    return false;
  }

  const auto var = grp.getVar(title);
  if (var.isNull())
  {
    messerr("Cannot read NetCDF String Variable of name %s in group %s", title.c_str(),
            grp.getName().c_str());
    return false;
  }

  // Assume variable is of dim 1
  if (var.getDimCount() != 1)
  {
    messerr(
      "NetCDF String Variable of name %s in group %s has %ul dims, but we expect only 1",
      title.c_str(), grp.getName().c_str(), var.getDimCount());
    return false;
  }

  const auto dim = var.getDim(0);
  vec.resize(dim.getSize());

  // Use a vector of char* managed by netCDF to read string data
  std::vector<char*> data_ptr(dim.getSize());
  var.getVar(data_ptr.data());

  // copy char pointers into gstlearn managed string vector
  for (size_t i = 0; i < data_ptr.size(); ++i)
  {
    vec[i] = data_ptr[i];
  }

  return true;
}

template<typename T>
bool SerializeNetCDF::_writeVec(netCDF::NcGroup& grp,
                                const String& title,
                                const VectorT<T>& vec,
                                const netCDF::NcDim& dim)
{
  if (vec.size() < dim.getSize())
  {
    messerr("Vector %s has size %ul, which is lesser than requested netCDF dim %ul",
            title.c_str(), vec.size(), dim.getSize());
    return false;
  }

  if (grp.isNull())
  {
    messerr("Cannot write vector %s: container group is null", title.c_str());
    return false;
  }

  // assume vectors of dim 1 (only one NcDim)

  const auto var = grp.addVar(title, getNetCDFType(vec[0]), dim);
  var.putVar(vec.constData());
  return true;
}

template<>
bool SerializeNetCDF::_writeVec(netCDF::NcGroup& grp,
                                const String& title,
                                const VectorString& vec,
                                const netCDF::NcDim& dim)
{
  if (vec.size() < dim.getSize())
  {
    messerr(
      "String Vector %s has size %ul, which is lesser than requested netCDF dim %ul",
      title.c_str(), vec.size(), dim.getSize());
    return false;
  }

  if (grp.isNull())
  {
    messerr("Cannot write string vector %s: container group is null", title.c_str());
    return false;
  }

  // generate a vector of char * to feed netCDF
  std::vector<const char*> data_ptr(vec.size());
  for (size_t i = 0; i < vec.size(); ++i)
  {
    data_ptr[i] = vec[i].c_str();
  }

  // assume vectors of dim 1 (only one NcDim)

  const auto var = grp.addVar(title, getNetCDFType(vec[0]), dim);
  var.putVar(data_ptr.data());
  return true;
}
