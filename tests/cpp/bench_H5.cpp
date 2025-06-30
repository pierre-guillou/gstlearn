#include "Enum/EFormatNF.hpp"
#include <Basic/Timer.hpp>
#include <Db/DbGrid.hpp>

int test_NF(const DbGrid& db)
{
  Timer tm;
  db.dumpToNF("dbgrid.NF.ascii", EFormatNF::ASCII);
  auto db2 = std::unique_ptr<DbGrid>(DbGrid::createFromNF("dbgrid.NF.ascii"));
  if (db2 != nullptr)
  {
    db2->dumpToNF("dbgrid2.NF.ascii", EFormatNF::ASCII);
  }
  else
  {
    messerr("Cannot Deserialize `dbgrid.NF.ascii'");
    return 1;
  }
  tm.displayIntervalMilliseconds("Serialize + Deserialize + Serialize Neutral File",
                                 2500);
  return 0;
}

int test_HDF5(const DbGrid& db)
{
#ifdef HDF5
  Timer tm;

  db.dumpToNF("dbgrid.NF.h5", EFormatNF::H5);
  auto db2 = std::unique_ptr<DbGrid>(DbGrid::createFromNF("dbgrid.NF.h5"));
  if (db2 != nullptr)
  {
    db2->dumpToNF("dbgrid2.NF.h5", EFormatNF::H5);
  }
  else
  {
    messerr("Cannot Deserialize `dbgrid.NF.h5'");
    return 1;
  }
  tm.displayIntervalMilliseconds("Serialize + Deserialize + Serialize HDF5", 80);
#endif
  return 0;
}

int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  const VectorInt dims {100, 100, 100};
  auto db = std::unique_ptr<DbGrid>(DbGrid::create(dims));
  int ret {};
  ret += test_NF(*db);
  ret += test_HDF5(*db);
  return ret;
}
