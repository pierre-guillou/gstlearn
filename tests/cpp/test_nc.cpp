#include <ncFile.h>
#include <ncDim.h>
#include <ncVar.h>

#include <array>
#include <numeric>

int run(const std::array<int, 3>& dims)
{
  netCDF::NcFile ds {"titi.nc", netCDF::NcFile::FileMode::replace};
  auto db = ds.addGroup("DbGrid");

  auto grch = db.addGroup("Grid characteristics");
  auto sp   = grch.addDim("Space Dimension", dims.size());
  std::vector<double> zeros(dims.size());
  std::vector<double> ones(dims.size(), 1.);
  auto nx = grch.addVar("NX", netCDF::NcType::nc_INT, sp);
  nx.putVar(dims.data());
  auto x0 = grch.addVar("X0", netCDF::NcType::nc_FLOAT, sp);
  x0.putVar(zeros.data());
  auto dx = grch.addVar("DX", netCDF::NcType::nc_FLOAT, sp);
  dx.putVar(ones.data());
  auto ag = grch.addVar("ANGLE", netCDF::NcType::nc_FLOAT, sp);
  ag.putVar(zeros.data());

  auto data = db.addGroup("Data");
  std::vector<netCDF::NcDim> data_dims {
    data.addDim("x", dims[0]),
    data.addDim("y", dims[1]),
    data.addDim("z", dims[2]),
  };

  std::vector<double> indices0(dims[0] * dims[1] * dims[2]);
  std::vector<double> indices1(dims[0] * dims[1] * dims[2]);
  std::vector<double> indices2(dims[0] * dims[1] * dims[2]);
  std::vector<double> vals(dims[0] * dims[1] * dims[2]);
  std::iota(vals.begin(), vals.end(), 1.);
  for (int i = 0; i < dims[0]; ++i)
  {
    for (int j = 0; j < dims[1]; ++j)
    {
      for (int k = 0; k < dims[2]; ++k)
      {
        indices0[(i * dims[1] * dims[2]) + (j * dims[2]) + k] = k;
        indices1[(i * dims[1] * dims[2]) + (j * dims[2]) + k] = j;
        indices2[(i * dims[1] * dims[2]) + (j * dims[2]) + k] = i;
      }
    }
  }

  auto x1 = data.addVar("x1", netCDF::NcType::nc_INT, data_dims);
  x1.putAtt("Locators", "x1");
  x1.putVar(indices0.data());
  auto x2 = data.addVar("x2", netCDF::NcType::nc_INT, data_dims);
  x2.putAtt("Locators", "x2");
  x2.putVar(indices1.data());
  auto x3 = data.addVar("x3", netCDF::NcType::nc_INT, data_dims);
  x3.putAtt("Locators", "x3");
  x3.putVar(indices2.data());

  auto rank = data.addVar("rank", netCDF::NcType::nc_FLOAT, data_dims);
  rank.putAtt("Locators", "NA");
  rank.putVar(vals.data());

  return 0;
}

int main()
{
  std::array<int, 3> dims {10, 10, 10};
  run(dims);
  return 0;
}
