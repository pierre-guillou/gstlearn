# Make Release version the default (only for single configuration generators)
# TODO : Differentiate build directories for Debug and Release
if(NOT IS_MULTI_CONFIG)
  if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting build type to 'Release' as none was specified")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
  endif()
  # Show current configuration
  message(STATUS "BUILD_TYPE=" ${CMAKE_BUILD_TYPE})
endif()

# Add c++20 support whatever the compiler
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Warning fiesta!
# https://cmake.org/cmake/help/latest/command/add_compile_options.html
if (MSVC)
  # Warning level 4 (4 = maximum, 0 = none)
  add_compile_options(/bigobj /W4 /wd4251 /wd4244 /wd4127 /wd4250) # Except those few ones
  # Silence MSVC warnings about unsafe C standard library functions
  add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
else()
  # Lots of warnings (-Wall = add some warnings, -Wextra = add a ton of warnings)
  add_compile_options(
    -Wall
    -Wextra
    -Wno-deprecated-copy
    -Wtype-limits
    -Wnon-virtual-dtor
    -Wvla
    -Wundef
  )
  if (APPLE)
    add_compile_options(-Wno-absolute-value -Wno-inconsistent-missing-override)
  endif()
endif()

if (MSVC)
  # Enable parallel compilation automatically.
  add_compile_options(/MP)
endif()

# Address Sanitizer (GCC/Clang)
option(BUILD_ASAN "Build with Address Sanitizer enabled" OFF)
mark_as_advanced(BUILD_ASAN)

if(BUILD_ASAN AND MSVC)
  message(WARNING "Cannot use BUILD_ASAN option with Microsoft Visual Studio compilers")
  set(BUILD_ASAN OFF)
endif()

if(BUILD_ASAN)
  add_compile_options(-fsanitize=address)
  add_link_options(-fsanitize=address)
endif()

# For valgrind usage (use Debug)
#add_compile_options(-O0)

# C++ header location (keep the trailing '/')
set(INCLUDES 
    ${PROJECT_SOURCE_DIR}/include/)
# C++ source path (prevent using GLOB)
include(src/all_sources.cmake)
set(SOURCES)
foreach(CPP ${SRC})
  set(SOURCES ${SOURCES} ${PROJECT_SOURCE_DIR}/src/${CPP})
endforeach(CPP ${SRC})

# target_include_directories() below is enough to get the code to compile but if
# includes are not explicitly added to SOURCES then Visual Studio doesn't see them
# as part of the project itself. They are GLOB'ed even though it's fragile for
# simplicity (it's just to get the project tree right).
if (CMAKE_GENERATOR MATCHES "Visual Studio")
  file(GLOB_RECURSE INCS RELATIVE ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/include/*.hpp ${PROJECT_SOURCE_DIR}/include/*.h)
  foreach(HPP ${INCS})
    set(SOURCES ${SOURCES} ${PROJECT_SOURCE_DIR}/include/${HPP})
  endforeach(HPP ${INCS})
endif()

# Generation folder (into Release or Debug)
if (NOT IS_MULTI_CONFIG AND NOT WIN32)
  cmake_path(APPEND CMAKE_CURRENT_BINARY_DIR ${CMAKE_BUILD_TYPE} OUTPUT_VARIABLE CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  cmake_path(APPEND CMAKE_CURRENT_BINARY_DIR ${CMAKE_BUILD_TYPE} OUTPUT_VARIABLE CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  cmake_path(APPEND CMAKE_CURRENT_BINARY_DIR ${CMAKE_BUILD_TYPE} OUTPUT_VARIABLE CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
endif()

# Debug find package instruction
#set(CMAKE_FIND_DEBUG_MODE TRUE)

# Look for Boost
#set(Boost_DEBUG 1)
find_package(Boost REQUIRED)
# TODO : If Boost not found, fetch it from the web ?

# Look for OpenMP
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
  add_definitions(-DOPENMP)
endif()

# Look for Eigen
find_package(Eigen3 REQUIRED) 
if(Eigen3_FOUND) 
    message(STATUS "Found Eigen3 version ${Eigen3_VERSION} in ${Eigen3_DIR}")
endif()

# Look for NLOPT
find_package(NLopt REQUIRED)
if (NLopt_FOUND)
    message(STATUS "Found NLopt ${NLOPT_VERSION} from ${NLOPT_CONFIG_FILE}")
else()
    message(FATAL_ERROR "NLopt not found")
endif()

# Look for HDF5
if(USE_HDF5)
  # use MODULE to use CMake's FindHDF5.cmake instead of HDF5 config files
  find_package(HDF5 MODULE REQUIRED COMPONENTS C CXX)
  if(NOT HDF5_FOUND)
    message(FATAL_ERROR "HDF5 not found")
  endif()
  add_definitions(-DHDF5)
endif()

# Shared and Static libraries
add_library(shared                  SHARED ${SOURCES})
set(FLAVORS shared)

############################## Loop on flavors: shared and static
foreach(FLAVOR ${FLAVORS})
  # Convert flavor to uppercase
  string(TOUPPER ${FLAVOR} FLAVOR_UP)

  # Alias target for a better name
  add_library(${PROJECT_NAME}::${FLAVOR} ALIAS ${FLAVOR})

  # Include directories
  # PUBLIC is mandatory for tests and packages (no need to install)
  target_include_directories(${FLAVOR} PUBLIC
    # Add includes path for compiling the library
    "$<BUILD_INTERFACE:${INCLUDES}>"
    # Add binary directory to find generated version.h and export.hpp
    "$<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>"
  )

  # Set some target properties
  set_target_properties(${FLAVOR} PROPERTIES
    # Hide all symbols by default (impose same behavior between Linux and Windows)
    C_VISIBILITY_PRESET hidden
    CXX_VISIBILITY_PRESET hidden
    # Any client who links the library needs -fPIC (static or shared)
    POSITION_INDEPENDENT_CODE 1
  )

  # Rename the output library name
  set_target_properties(${FLAVOR} PROPERTIES OUTPUT_NAME ${PROJECT_NAME})
  # append a 'd' to the output file name of the debug build targets
  set_target_properties(${FLAVOR} PROPERTIES DEBUG_POSTFIX "d")
  
  # Set library version
  set_target_properties(${FLAVOR} PROPERTIES VERSION ${PROJECT_VERSION})

  if(USE_BOOST_SPAN)
    target_compile_definitions(${FLAVOR} PUBLIC USE_BOOST_SPAN)
  endif()

  # Enable OpenMP
  target_link_libraries(${FLAVOR} PRIVATE OpenMP::OpenMP_CXX)

  # Link to csparse and gmtsph
  target_link_libraries(${FLAVOR} PRIVATE csparse gmtsph)
    
  # Link to Eigen
  target_link_libraries(${FLAVOR} PUBLIC Eigen3::Eigen)

  # Link to Boost (use headers)
  # Target for header-only dependencies. (Boost include directory)
  # Currently Boost headers are only used in .cpp so a PRIVATE link minimizes
  # dependencies for projects using gstlearn, except with USE_BOOST_SPAN.
  if(USE_BOOST_SPAN)
    target_link_libraries(${FLAVOR} PUBLIC Boost::boost)
  else()
    target_link_libraries(${FLAVOR} PRIVATE Boost::boost)
  endif()

  # Link to NLopt
  target_link_libraries(${FLAVOR} PRIVATE NLopt::nlopt)

  # Link to HDF5
  if(USE_HDF5)
    target_compile_definitions(${FLAVOR} PRIVATE ${HDF5_DEFINITIONS})
    target_include_directories(${FLAVOR} PRIVATE ${HDF5_INCLUDE_DIRS})
    target_link_libraries(${FLAVOR} PRIVATE ${HDF5_LIBRARIES})
  endif()
  
  # Exclude [L]GPL features from Eigen
  #target_compile_definitions(${FLAVOR} PUBLIC EIGEN_MPL2_ONLY) 

  # Link to specific libraries (only for Microsoft Visual Studio)
  if (MSVC)
    target_link_libraries(${FLAVOR} PUBLIC iphlpapi rpcrt4)
  endif()
  if (MINGW)
    target_link_libraries(${FLAVOR} PUBLIC -liphlpapi -lrpcrt4)
  endif()

  if(CMAKE_COMPILER_IS_GNUCC)
    # Use of std::filesystem needs at least GCC 8.0
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8.0)
      message(SEND_ERROR "GCC>=8.0 is needed to build gstlearn")
    elseif(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
      # GCC 8.0 doesn't link automatically to std::filesystem
      target_link_libraries(${FLAVOR} PUBLIC stdc++fs)
    endif()
  endif()

  # Build a cmake file to be imported by library users
  export(TARGETS ${FLAVOR}
         NAMESPACE ${PROJECT_NAME}::
         FILE ${GSTLEARN_CMAKE_FILE}
         APPEND)
         
endforeach(FLAVOR ${FLAVORS})
############################## End loop on flavors


###################### Shared library specific options

# Generate export header
include(GenerateExportHeader)
set(DISABLE_EXPORT_IF_SWIG "
#ifdef SWIG
#    undef ${PROJECT_NAME_UP}_EXPORT
#    undef ${PROJECT_NAME_UP}_NO_EXPORT
#    undef ${PROJECT_NAME_UP}_DEPRECATED
#    undef ${PROJECT_NAME_UP}_DEPRECATED_EXPORT
#    undef ${PROJECT_NAME_UP}_DEPRECATED_NO_EXPORT
#    define ${PROJECT_NAME_UP}_EXPORT
#    define ${PROJECT_NAME_UP}_NO_EXPORT
#    define ${PROJECT_NAME_UP}_DEPRECATED
#    define ${PROJECT_NAME_UP}_DEPRECATED_EXPORT
#    define ${PROJECT_NAME_UP}_DEPRECATED_NO_EXPORT
#endif
#ifdef ${PROJECT_NAME_UP}_STATIC_DEFINE
#    define ${PROJECT_NAME_UP}_TEMPLATE_EXPORT
#else
#    ifdef shared_EXPORTS
#        define ${PROJECT_NAME_UP}_TEMPLATE_EXPORT
#    else
#        define ${PROJECT_NAME_UP}_TEMPLATE_EXPORT extern
#    endif
#endif
")
generate_export_header(shared
  BASE_NAME ${PROJECT_NAME}
  EXPORT_FILE_NAME ${CMAKE_BINARY_DIR}/${PROJECT_NAME}_export.hpp
  CUSTOM_CONTENT_FROM_VARIABLE DISABLE_EXPORT_IF_SWIG
)

# Set the so version to project major version
set_target_properties(shared PROPERTIES
  SOVERSION ${PROJECT_VERSION_MAJOR}
)
