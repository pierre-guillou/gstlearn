# use ncxx4-config to find location of netcdf-cxx4 headers and libraries

find_program(
  NCXX4_CONFIG
  NAMES ncxx4-config
  HINTS $ENV{NetCDF-CXX4_DIR}/bin)
mark_as_advanced(NCXX4_CONFIG)

execute_process(
  COMMAND ${NCXX4_CONFIG} --cflags
  OUTPUT_VARIABLE netcdf-cxx4_CFLAGS
  OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(
  COMMAND ${NCXX4_CONFIG} --libs
  OUTPUT_VARIABLE netcdf-cxx4_LIBRARIES
  OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(
  COMMAND ${NCXX4_CONFIG} --version
  OUTPUT_VARIABLE netcdf-cxx4_VERSION
  OUTPUT_STRIP_TRAILING_WHITESPACE)
string(SUBSTRING ${netcdf-cxx4_VERSION} 12 -1 netcdf-cxx4_VERSION)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  NetCDF-CXX4 REQUIRED_VARS netcdf-cxx4_CFLAGS netcdf-cxx4_LIBRARIES)

if(NetCDF-CXX4_FOUND)
  if(NOT TARGET NetCDF-CXX4)
    add_library(NetCDF-CXX4 INTERFACE IMPORTED)
    set_target_properties(
      NetCDF-CXX4
      PROPERTIES INTERFACE_COMPILE_OPTIONS ${netcdf-cxx4_CFLAGS}
                 INTERFACE_LINK_LIBRARIES ${netcdf-cxx4_LIBRARIES}
                 VERSION ${netcdf-cxx4_VERSION})
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
  endif()
endif()
