#include <Basic/Timer.hpp>
#include <Db/DbGrid.hpp>

#include <array>

int test(const std::array<int, 3>& dims)
{
  auto db = std::unique_ptr<DbGrid>(DbGrid::create(VectorInt {dims.begin(), dims.end()}));
  {
    Timer tm;
    db->dumpToNF("dbgrid.nf");
    auto db2 = std::unique_ptr<DbGrid>(DbGrid::createFromNF("dbgrid.nf"));
    if (db2 != nullptr)
    {
      db2->dumpToNF("dbgrid2.nf");
    }
    else
    {
      messerr("Cannot deserialize `dbgrid.nf'");
    }
    tm.displayIntervalMilliseconds();
  }
  {
    Timer tm;
    db->dumpToNC("dbgrid.nc");
    auto db2 = std::unique_ptr<DbGrid>(DbGrid::createFromNC("dbgrid.nc"));
    if (db2 != nullptr)
    {
      db2->dumpToNC("dbgrid2.nc");
    }
    else
    {
      messerr("Cannot deserialize `dbgrid.nc'");
    }
    tm.displayIntervalMilliseconds();
  }
  return 0;
}

int main()
{
  std::array<int, 3> dims {100, 100, 100};
  return test(dims);
}
