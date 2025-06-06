# For recent UseSWIG
cmake_minimum_required(VERSION 3.20)

# Define project here
project(gstlearn
        VERSION      1.8.0 # Make it consistent with CHANGES
        DESCRIPTION  "Geostatistics & Machine Learning toolbox"
        HOMEPAGE_URL "https://gstlearn.org"
        LANGUAGES    C CXX) # Enable C language for third party libraries

####################################################
## CONFIGURATION

# By default, activate HDF5 support
option(USE_HDF5 "Using HDF5 support" ON)

# Extract Git commit timestamp, hash
find_package(Git)
execute_process(
  COMMAND ${GIT_EXECUTABLE} show -s --format=%ci
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  OUTPUT_VARIABLE ${PROJECT_NAME}_DATE
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process(
  COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  OUTPUT_VARIABLE ${PROJECT_NAME}_COMMIT
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
# For bib file generation
execute_process(
  COMMAND ${GIT_EXECUTABLE} show -s --format=%cd --date=format:%Y
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  OUTPUT_VARIABLE ${PROJECT_NAME}_YEAR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process(
  COMMAND ${GIT_EXECUTABLE} show -s --format=%cd --date=format:%B
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  OUTPUT_VARIABLE ${PROJECT_NAME}_MONTH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Convert project name to uppercase
string(TOUPPER ${PROJECT_NAME} PROJECT_NAME_UP)

# Fix version.h automatically
# Define a numerical version number (usable in e.g. #if ( VERSION >= ...))
MATH( EXPR gstlearn_VERSION_NUMBER "${gstlearn_VERSION_MAJOR} * 10000 + ${gstlearn_VERSION_MINOR} * 100 + ${gstlearn_VERSION_PATCH}" )
configure_file(version.h.in version.h)

# Detect presence of multi configuration generators
get_property(IS_MULTI_CONFIG GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)

# Activate verbose mode
#set(CMAKE_VERBOSE_MAKEFILE ON)
# Debug find package instruction
#set(CMAKE_FIND_DEBUG_MODE TRUE)

# Options
# By default, deactivate R and python packages
option(BUILD_PYTHON  "Build Python package"                     OFF)
option(BUILD_R       "Build R package"                          OFF)
option(BUILD_TESTING "Build tests"                              OFF)
option(BUILD_DOC     "Build Doxygen and packages documentation" OFF)
if (MINGW)
  set(BUILD_PYTHON OFF)
endif()
if (MSVC)
  set(BUILD_R OFF)
endif()
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  # Clang compiler detected (define the shortcut variable)
  set(CLANG ON)
endif()
message(STATUS "BUILD_PYTHON="  ${BUILD_PYTHON})
message(STATUS "BUILD_R="       ${BUILD_R})
message(STATUS "BUILD_TESTING=" ${BUILD_TESTING})
message(STATUS "BUILD_DOC="     ${BUILD_DOC})
message(STATUS "USE_HDF5="      ${USE_HDF5})

# Create gstlearn cmake file path in the build tree
cmake_path(APPEND GSTLEARN_CMAKE_FILE ${PROJECT_BINARY_DIR} "gstlearn.cmake")

####################################################
## EXTERNAL LIBS
add_subdirectory(3rd-party/gmtsph)
add_subdirectory(3rd-party/csparse)

####################################################
## BUILD & INSTALL

include(cmake/version.cmake)
include(cmake/doc.cmake)
# TODO : Add BUILD_CPP option (not needed for print_version or doxygen targets)
#if (BUILD_CPP)
  include(cmake/cpp.cmake)
  include(cmake/install.cmake)
#endif()
if (BUILD_DOC)
  include(cmake/doxygen.cmake)
endif()

####################################################
## PACKAGES

if (BUILD_PYTHON)
  add_subdirectory(python)
  set(SWIG ON)  
endif()

if (BUILD_R)
  add_subdirectory(r)
  set(SWIG ON)
endif()

if (SWIG)
  add_subdirectory(swig)
endif()

####################################################
## TESTING

# Add non-regression test target
if (BUILD_TESTING)
  include(CTest)
  enable_testing()

  add_subdirectory(tests)
endif()

####################################################
## DEMONSTRATION SCRIPTS

add_subdirectory(doc/demo)
add_subdirectory(doc/courses)
