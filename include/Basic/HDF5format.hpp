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
#include "Basic/String.hpp"
#include "Basic/AStringable.hpp"
#include "H5Cpp.h"
#include "typeinfo"

#if H5_VERSION_GE(1,8,20)
#define EXCEPTION_PRINT_ERROR(e) e.printErrorStack();
#else
#define EXCEPTION_PRINT_ERROR(e) e.printError();
#endif

class GSTLEARN_EXPORT HDF5format
{
public:
  HDF5format(const String& filename = String(), const String& varname = String());
  HDF5format(const HDF5format &r);
  HDF5format& operator=(const HDF5format &r);
  virtual ~HDF5format();

public:
  void* readRegular(int flag_compress,
                    hsize_t *start,
                    hsize_t *stride,
                    hsize_t *count,
                    hsize_t *block,
                    hsize_t *dimout);
  int writeRegular(hsize_t *start,
                   hsize_t *stride,
                   hsize_t *count,
                   hsize_t *block,
                   void *wdata);
  int deleteFile();
  void* allocArray(H5::DataType type, int ndim, hsize_t *dims);

  void setFileName(const String& filename) { _filename = filename; }
  void setVarName(const String& varname) { _varname = varname; }

  int displayNames() const;

  void openFile(const String& filename = String()) const;
  void openNewFile(const String& filename) const;
  void openDataSet(const String& varname = String()) const;
  void openNewDataSet(const String& varname,
                      int ndim,
                      hsize_t *dims,
                      H5::DataType type) const;
  void closeFile() const;
  void closeDataSet() const;

  // Functions to be overloaded

  template<typename T>
  void writeData(const T&);
  template<typename T>
  void writeData(const std::vector<T>&);
  template<typename T>
  void writeData(const std::vector<std::vector<T> >&);

  int getDataInt() const;
  float getDataFloat() const;
  double getDataDouble() const;
  VectorInt getDataVInt() const;
  VectorFloat getDataVFloat() const;
  VectorDouble getDataVDouble() const;
  VectorVectorInt getDataVVInt() const;
  VectorVectorFloat getDataVVFloat() const;
  VectorVectorDouble getDataVVDouble() const;

  VectorDouble getDataDoublePartial(int myrank) const;
  int writeDataDoublePartial(int myrank, const VectorDouble& data);

  // Return the size of the data
  // Note that for multi-dim arrays that it gets the total size and not the size of a single row.
  int getSize() const;

  // We now make a proxy class so that we can overload the return type and use a single
  // function to get data whether int or float or double.
  class Proxy
  {
  private:
    HDF5format const* myOwner;
  public:
    Proxy(const HDF5format* owner)
        : myOwner(owner)
    {
    }
    operator int() const
    {
      return myOwner->getDataInt();
    }
    operator float() const
    {
      return myOwner->getDataFloat();
    }
    operator double() const
    {
      return myOwner->getDataDouble();
    }
    operator std::vector<int>() const
    {
      return myOwner->getDataVInt();
    }
    operator std::vector<float>() const
    {
      return myOwner->getDataVFloat();
    }
    operator std::vector<double>() const
    {
      return myOwner->getDataVDouble();
    }
    operator std::vector<std::vector<int> >() const
    {
      return myOwner->getDataVVInt();
    }
    operator std::vector<std::vector<float> >() const
    {
      return myOwner->getDataVVFloat();
    }
    operator std::vector<std::vector<double> >() const
    {
      return myOwner->getDataVVDouble();
    }
  };
  // Here we use the Proxy class to have a single getData function
  Proxy getData() const
  {
    return Proxy(this);
  }

private:
  int _getNDim() const;
  hsize_t* _getDims() const;
  void _getOrderSize(H5T_order_t* order, size_t* size, bool* big_endian) const;
  void _readInt(int *data,
                const H5::DataSpace& memspace = H5::DataSpace::ALL,
                const H5::DataSpace& dataspace = H5::DataSpace::ALL) const;
  void _readFloat(float *data,
                  const H5::DataSpace&  = H5::DataSpace::ALL,
                  const H5::DataSpace& dataspace = H5::DataSpace::ALL) const;
  void _readDouble(double *data,
                   const H5::DataSpace& memspace = H5::DataSpace::ALL,
                   const H5::DataSpace& dataspace = H5::DataSpace::ALL) const;
  void _writeDouble(double *data,
                    const H5::DataSpace& memspace = H5::DataSpace::ALL,
                    const H5::DataSpace& dataspace = H5::DataSpace::ALL) const;
  int _checkClass(int value) const;
  void _writeAll(const char* type, void* a);

public:
  mutable String        _filename;
  mutable String        _varname;
  mutable H5::H5File    _datafile;
  mutable H5::DataSet   _dataset;
  mutable H5::DataType  _datatype;
  mutable H5::DataSpace _dataspace;
};


/**
 * Numeric implementation of our write data function
 * Only accepts numerical values. Integers, floats, or doubles
 * @param data
 */
template<typename T>
void HDF5format::writeData(const T &data)
{
  H5::Exception::dontPrint();
  char* type = (char*) (typeid(T).name());
  auto *a = new T { data };
  _writeAll(type, (void*) a);
  delete a;
}

template<typename T>
void HDF5format::writeData(const std::vector<T> &data)
{
  H5::Exception::dontPrint();
  size_t npts = data.size();
  auto *a = new T[npts];
  char* type = (char*) (typeid(a[0]).name());
  for (size_t i = 0; i < npts; ++i)
    a[i] = data[i];
  _writeAll(type, (void*) a);
  delete[] a;
 }

template<typename T>
void HDF5format::writeData(const std::vector<std::vector<T> > &data)
{
  H5::Exception::dontPrint();
  size_t dim1 = data.size();
  size_t dim2 = data[0].size();
  auto a = new T[dim1 * dim2];
  auto md = new T*[dim1];
  for (size_t i = 0; i < dim1; ++i)
    md[i] = a + i * dim2;
  for (size_t i = 0; i < dim1; ++i)
    for (size_t j = 0; j < dim2; ++j)
      md[i][j] = data[i][j];
  char* type = (char*) (typeid(a[0]).name());
  _writeAll(type, (void* ) a);
  delete[] md;
  delete a;
}