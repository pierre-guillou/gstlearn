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

#include "Basic/String.hpp"
#include <stdarg.h>

class ASerializable
{
public:
  ASerializable();
  virtual ~ASerializable() {};

  virtual int deSerialize(const String& filename, bool verbose = false) = 0;
  virtual int serialize(const String& filename, bool verbose = false) = 0;

  const String& getCurrentRecord() const { return _currentRecord; }
  const FILE*   getFile() const { return _file; }
  const String& getFileName() const { return _fileName; }
  const String& getFileType() const { return _fileType; }

protected:
  int  _fileOpen(const String& filename, const String& filetype, const String& mode);
  int  _fileClose();
  int  _recordRead(const String& title, const String& format, ...);
  void _recordWrite(const String& format, ...);
  int  _fileRead(const String& format, va_list ap);
  void _fileWrite(const String& format, va_list ap);
  bool _onlyBlanks(char *string);

private:
  String _fileName;
  String _fileType;
  FILE*  _file;
  String _currentRecord;
};