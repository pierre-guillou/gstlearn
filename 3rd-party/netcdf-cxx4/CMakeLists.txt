find_package(HDF5 REQUIRED COMPONENTS C HL)
find_package(netCDF REQUIRED COMPONENTS C)

if(netCDF_FOUND)
  message(STATUS "Found netCDF version ${netCDF_VERSION} from ${netCDF_CONFIG}")
else()
  message(FATAL_ERROR "netCDF not found")
endif()

add_library(
  netcdf-cxx4 STATIC
  src/ncAtt.cpp
  src/ncByte.cpp
  src/ncChar.cpp
  src/ncCheck.cpp
  src/ncCompoundType.cpp
  src/ncDim.cpp
  src/ncDouble.cpp
  src/ncEnumType.cpp
  src/ncException.cpp
  src/ncFile.cpp
  src/ncFill.cpp
  src/ncFilter.cpp
  src/ncFloat.cpp
  src/ncGroup.cpp
  src/ncGroupAtt.cpp
  src/ncInt.cpp
  src/ncInt64.cpp
  src/ncOpaqueType.cpp
  src/ncShort.cpp
  src/ncString.cpp
  src/ncType.cpp
  src/ncUbyte.cpp
  src/ncUint.cpp
  src/ncUint64.cpp
  src/ncUshort.cpp
  src/ncVar.cpp
  src/ncVarAtt.cpp
  src/ncVlenType.cpp)
target_include_directories(netcdf-cxx4 PUBLIC include)
target_link_libraries(netcdf-cxx4 PUBLIC hdf5::hdf5_hl netCDF::netcdf)

# Any client who links the library needs -fPIC
set_target_properties(netcdf-cxx4 PROPERTIES POSITION_INDEPENDENT_CODE 1)

export(
  TARGETS netcdf-cxx4
  NAMESPACE ${PROJECT_NAME}::
  FILE ${GSTLEARN_CMAKE_FILE}
  APPEND)
