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
#include "Matrix/Table.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/String.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/AException.hpp"
#include "Basic/ASerializable.hpp"
#include "Basic/AStringable.hpp"

Table::Table(int nrow, int ncol)
  : MatrixRectangular(nrow, ncol),
    ASerializable(),
    _title(),
    _rowNames(),
    _colNames()
{
  init(nrow, ncol);
}

Table::Table(const Table &m)
    : MatrixRectangular(m),
      ASerializable(m),
      _title(m._title),
      _rowNames(m._rowNames),
      _colNames(m._colNames)
{

}

Table& Table::operator=(const Table &m)
{
  if (this != &m)
  {
    MatrixRectangular::operator=(m);
    ASerializable::operator=(m);
    _title = m._title;
    _rowNames = m._rowNames;
    _colNames = m._colNames;
  }
  return *this;
}

Table::~Table()
{
}

void Table::_clearContents()
{
  _rowNames.clear();
  _colNames.clear();
  return;
}

Table* Table::create(int nrow, int ncol)
{
  return new Table(nrow, ncol);
}

Table* Table::createFromNames(const VectorString &rownames,
                              const VectorString &colnames)
{
  int nrow = (int) rownames.size();
  int ncol = (int) colnames.size();
  Table* table = new Table(nrow, ncol);
  table->setRowNames(rownames);
  table->setColumnNames(colnames);
  return table;
}

Table* Table::createFromNF(const String& neutralFilename, bool verbose)
{
  Table* table = nullptr;
  std::ifstream is;
  table = new Table();
  bool success = false;
  if (table->_fileOpenRead(neutralFilename, is, verbose))
  {
    success =  table->deserialize(is, verbose);
  }
  if (! success)
  {
    delete table;
    table = nullptr;
  }
  return table;
}

VectorDouble Table::getRange(int icol) const
{
  VectorDouble vec = getColumn(icol);
  if (vec.empty()) return VectorDouble();
  VectorDouble limits(2);
  limits[0] = VH::minimum(vec);
  limits[1] = VH::maximum(vec);
  return limits;
}

VectorDouble Table::getAllRange() const
{
  int ncols = getNCols();
  VectorDouble limits(2);
  limits[0] =  1.e30;
  limits[1] = -1.e30;
  for (int icol = 0; icol < ncols; icol++)
  {
    VectorDouble local = getRange(icol);
    if (local[0] < limits[0]) limits[0] = local[0];
    if (local[1] > limits[1]) limits[1] = local[1];
  }
  return limits;
}

bool Table::_serialize(std::ostream& os, bool /*verbose*/) const
{
  bool ret = true;
  ret = ret && _recordWrite<int>(os, "Number of Columns", getNCols());
  ret = ret && _recordWrite<int>(os, "Number of Rows", getNRows());

  /* Writing the tail of the file */

  for (int irow = 0; ret && irow < getNRows(); irow++)
  {
    for (int icol = 0; ret && icol < getNCols(); icol++)
    {
      ret = ret && _recordWrite<double>(os, "", getValue(irow, icol));
    }
    ret = ret && _commentWrite(os, "");
  }
  return ret;
}

bool Table::_deserialize(std::istream& is, bool /*verbose*/)
{
  int nrows = 0;
  int ncols = 0;
  double value = 0.;

  bool ret = true;
  ret = ret && _recordRead<int>(is, "Number of Columns", ncols);
  ret = ret && _recordRead<int>(is, "Number of Rows", nrows);
  if (! ret) return false;

  reset(nrows, ncols);

  /* Loop on the lines */

  for (int irow = 0; ret && irow < nrows; irow++)
  {
    for (int icol = 0; ret && icol < ncols; icol++)
    {
      ret = ret && _recordRead<double>(is, "Numerical value", value);
      if (ret) setValue(irow, icol, value);
    }
  }
  return ret;
}

/**
 * Print the contents of the statistics
 */
String Table::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  if (isEmpty()) return sstr.str();

  if (_title.empty())
    sstr << toTitle(1, "Table contents");
  else
    sstr << toTitle(1, _title.c_str());

  int ncols = getNCols();
  int nrows = getNRows();
  sstr << "- Number of Rows    = " << nrows << std::endl;
  sstr << "- Number of Columns = " << ncols << std::endl;
  sstr << std::endl;

  // Print optional header (using Column names if defined)

  if (! _colNames.empty())
  {
    if (! _rowNames.empty()) sstr << toStr(" ");
    for (int icol = 0; icol < ncols; icol++)
      sstr << " " << toStr(_colNames[icol]);
    sstr << std::endl;
  }

  // Print the contents of the table
  for (int irow = 0; irow < nrows; irow++)
  {
    if (! _rowNames.empty()) sstr << toStr(_rowNames[irow]);
    for (int icol = 0; icol < ncols; icol++)
    {
      sstr << " " << toDouble(getValue(irow, icol));
    }
    sstr << std::endl;
  }
  return sstr.str();
}

/**
 * Plot the contents of the statistics
 */
void Table::plot(int isimu) const
{
  if (isEmpty()) return;
  String filename = incrementStringVersion("TableStats",isimu+1);
  (void) dumpToNF(filename);
}

void Table::setColumnNames(const VectorString &colNames)
{
  if (getNCols() != (int) colNames.size())
  {
    messerr("The size of 'colNames' (%d) does not match the number of columns (%d)",
            (int) colNames.size(), getNCols());
    return;
  }
  _colNames = colNames;
}

void Table::setColumnName(int icol, const String& name)
{
  if (! _isColumnValid(icol)) return;
  int ncols = getNCols();
  if (_colNames.empty())
    _colNames.resize(ncols, "  ");
  _colNames[icol] = name;
}

void Table::setRowName(int irow, const String& name)
{
  if (! _isRowValid(irow)) return;
  int nrows = getNRows();
  if (_rowNames.empty())
    _rowNames.resize(nrows, "  ");
  _rowNames[irow] = name;
}

void Table::setRowNames(const VectorString &rowNames)
{
  if (getNRows() != (int) rowNames.size())
  {
    messerr("The size of 'rowNames' (%d) does not match the number of rows (%d)",
            (int) rowNames.size(), getNRows());
    return;
  }
  _rowNames = rowNames;
}

String Table::getColumnName(int icol) const
{
  if (! _isColumnValid(icol)) return String();
  return _colNames[icol];
}

String Table::getRowName(int irow) const
{
  if (! _isRowValid(irow)) return String();
  return _rowNames[irow];
}