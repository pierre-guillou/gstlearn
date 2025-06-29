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
#include "Basic/ASerializable.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/SerializeHDF5.hpp"
#include "Basic/SerializeNeutralFile.hpp"
#include "Basic/File.hpp"
#include "Basic/String.hpp"
#include "Enum/EFormatNF.hpp"

#include <iostream>
#include <filesystem>
#include <fstream>

#if defined(_WIN32) || defined(_WIN64)
#include <windows.h> // for CreateDirectory
#else
#include <unistd.h> // for readlink
#endif

#ifdef __APPLE__
#include <mach-o/dyld.h> // for _NSGetExecutablePath
#endif

#include <sys/stat.h>
#include <sys/types.h>

String ASerializable::_myPrefixName = String();
EFormatNF ASerializable::_defaultFormatNF = EFormatNF::ASCII;

ASerializable::ASerializable()                                    = default;
ASerializable::ASerializable(const ASerializable&)                = default;
ASerializable& ASerializable::operator=(const ASerializable&)     = default;
ASerializable::ASerializable(ASerializable&&) noexcept            = default;
ASerializable& ASerializable::operator=(ASerializable&&) noexcept = default;
ASerializable::~ASerializable()                                   = default;

void ASerializable::defineDefaultFormatNF(const EFormatNF& format)
{
  _defaultFormatNF = format;
}

/**
 * @brief Dump the contents of an object into an Output File
 * using a given Output NF Format
 *
 * @param NFFilename Name of the Output File
 * @param format Choice of the format (see remarks)
 * @param verbose Verbose flag
 * @return true or false
 *
 * @remarks In the argument 'format', the user can select the format for encoding
 * the contents of the output file.
 * If the value DEFAULT is used, the package uses the Format currently defined
 * as the defaulted one. This default value can be updated using the method
 * ASerializable::defineDefaultFormatNF()
 */
bool ASerializable::dumpToNF(const String& NFFilename,
                             const EFormatNF& format, 
                             bool verbose) const
{
  bool ret = true;

  EFormatNF formatLocal = format;
  if (format == EFormatNF::DEFAULT)
    formatLocal = _defaultFormatNF;

  // Check if H5 format is available: otherwise force ASCII
  bool canUseH5 = false;
  #ifdef HDF5
    canUseH5 = true;
  #endif
  if (! canUseH5) formatLocal = EFormatNF::ASCII;
  // TODO: message temporarily printed to help debugging on Windows
  message("Dump to NF is using the Format %s\n", formatLocal.getKey());

  if (formatLocal == EFormatNF::ASCII)
  {
    std::ofstream os;
    if (SerializeNeutralFile::fileOpenWrite(*this, NFFilename, os, true))
    {
      ret = _serializeAscii(os, verbose);
      if (! ret)
        messerr("Problem writing in the Neutral File.");
      os.close();
    }
    return ret;
  }

#ifdef HDF5
  if (formatLocal == EFormatNF::H5)
  {
    auto file = SerializeHDF5::fileOpenWrite(*this, NFFilename);
    bool ret  = _serializeH5(file, verbose);
    if (!ret)
      messerr("Problem writing in the HDF5 file.");
    return ret;
  }
#endif

  messerr("Aserializable::dumpToNF : No Format is defined");
  return false;
}

bool ASerializable::_fileOpenWrite(const String& filename,
                                   std::ofstream& os,
                                   bool verbose) const
{
  return SerializeNeutralFile::fileOpenWrite(*this, filename, os, verbose);
}

bool ASerializable::_fileOpenRead(const String& filename,
                                  std::ifstream& is,
                                  bool verbose) const
{
  return SerializeNeutralFile::fileOpenRead(*this, filename, is, verbose);
}

bool ASerializable::_commentWrite(std::ostream& os, const String& comment)
{
  return SerializeNeutralFile::commentWrite(os, comment);
}

bool ASerializable::_tableWrite(std::ostream& os,
                                const String& string,
                                int ntab,
                                const VectorDouble& tab)
{
  return SerializeNeutralFile::tableWrite(os, string, ntab, tab);
}

bool ASerializable::_tableRead(std::istream &is,
                               const String &string,
                               int ntab,
                               double *tab)
{
  return SerializeNeutralFile::tableRead(is, string, ntab, tab);
}

/**
 * Build a standard filename for Read or Write operation
 * @param status 1 for Read and 2 for Write
 * @param filename Name of the filename (see remark)
 * @param ensureDirExist When TRUE, the Directory is created if not already existing
 * @return
 */
String ASerializable::buildFileName(int status, const String& filename, bool ensureDirExist)
{
  // In the case of Output File (2), 'filename' is appended after the 'containerName' and 'prefixName'
  // In the case of Input file (1), the process depends on the contents of 'filename':
  // - if 'filename' is absolute (starts with '/' or second character is ':'): do nothing
  // - otherwise, add the 'containerName' and 'prefixName' (if defined)

  std::filesystem::path fileLocal {filename};

  if (status == 1 && fileLocal.is_absolute())
  {
    return fileLocal.string();
  }

  fileLocal.clear();

  // container name: first search for the GSTLEARN_OUTPUT_DIR
  // environment variable, then if empty create a `gstlearn_dir'
  // folder in the current directory
  const auto output_dir = gslGetEnv("GSTLEARN_OUTPUT_DIR");

  if (!output_dir.empty())
  {
    fileLocal = output_dir;
  }
  else
  {
    fileLocal = std::filesystem::current_path() / "gstlearn_dir";
  }

  if (ensureDirExist)
  {
    std::filesystem::create_directory(fileLocal);
  }

  const auto fname = _myPrefixName + filename;

  return (fileLocal / fname).string();
}

/**
 * Returns the Identity of a Neutral File which allows knowing its type
 * @param filename Name of the Neutral File
 * @param verbose Verbose flag
 * @return
 */
String ASerializable::getFileIdentity(const String& filename, bool verbose)
{
  // Preliminary check (no message if string is empty ... even in verbose)
  if (filename.empty()) return String();
  if (verbose)
    message("Input File Name = %s\n", filename.c_str());

  // Build the multi-platform filename
  const auto filepath = ASerializable::buildFileName(1, filename, true);
  if (verbose)
    message("Input File Path = %s\n", filepath.c_str());

  // Open the file according to various formats
  int ret_type = -1;
  String classType;

#ifdef HDF5
  if (ret_type < 0)
  {
    // Check if the file is written at format HDF5
    if (H5::H5File::isHdf5(filepath))
    {

      // Attempt to open according to a H5 format
      H5::H5File file {filepath, H5F_ACC_RDONLY};

      if (!file.nameExists("gstlearn metadata"))
      {
        messerr("File %s doesn't contain Gstlearn metadata…", filepath.c_str());
        return String();
      }
      ret_type = 1;

      // Read the class type
      classType = SerializeHDF5::getFileClass(filepath, verbose);
    }
  }
#endif

  if (ret_type < 0)
  {

    // Attempt to open according to the ASCII format
    std::ifstream file(filepath);
    if (!file.is_open())
    {
      if (verbose) messerr("Could not open the Neutral File %s", filepath.c_str());
      return String();
    }
    ret_type = 0;

    // Read the Class Type
    gslSafeGetline(file, classType);
    classType = trimRight(classType);

    // Close the file
    file.clear();
  }

  if (verbose) message("Decoded Type = %s\n", classType.c_str());

  return classType;
}

void ASerializable::setPrefixName(const String& prefixName)
{
  _myPrefixName = prefixName;
}

void ASerializable::unsetPrefixName(void)
{
  _myPrefixName.clear();
}

const String& ASerializable::getPrefixName()
{
  return _myPrefixName;
}
