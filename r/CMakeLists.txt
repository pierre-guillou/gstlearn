if(NOT BUILD_R)
  return()
endif()

# R package name
set(R_PACKAGE_NAME ${PROJECT_NAME})
# Beware : the real package name is the PROJECT_NAME (for better import)

###################################
# DEPENDENCIES

# Look for R
#set(CMAKE_FIND_DEBUG_MODE TRUE)
include(root/FindR.cmake)
if(NOT R_FOUND)
  message(FATAL_ERROR "R not found")
endif()

# Look for SWIG 4.2.0 [customized]
find_package(SWIG 4.2.0 REQUIRED)
message(STATUS "Found SWIG: " ${SWIG_EXECUTABLE} " (found version \"" ${SWIG_VERSION} "\")")
include(${SWIG_USE_FILE})

# Generation folders
if (IS_MULTI_CONFIG)
  set(R_PACKAGE_ROOT_FOLDER        ${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>)
else()
  set(R_PACKAGE_ROOT_FOLDER        ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_BUILD_TYPE})
endif()
set(R_PACKAGE_DESTINATION_FOLDER   ${R_PACKAGE_ROOT_FOLDER}/${R_PACKAGE_NAME})
set(R_PACKAGE_R_FOLDER             ${R_PACKAGE_DESTINATION_FOLDER}/R)
set(R_PACKAGE_DATA_FOLDER          ${R_PACKAGE_DESTINATION_FOLDER}/data)
set(R_PACKAGE_SRC_FOLDER           ${R_PACKAGE_DESTINATION_FOLDER}/src)

######################################  
# DOCUMENTATION
#
# Ensure documentation.i is existing
# TODO : put documentation.i out of the source
#add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/doc/documentation.i
#  COMMAND ${CMAKE_COMMAND} -E touch ${CMAKE_CURRENT_SOURCE_DIR}/doc/documentation.i
#  COMMENT "Ensure that documentation.i exists")
#
# Generate python documentation (before wrapper)
# TODO : Do not regenerate python documentation if doxygen output is up to date
#add_custom_target(python_doc
#  COMMAND ${Python3_EXECUTABLE}
#    ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxy2swig.py -q
#    ${PROJECT_BINARY_DIR}/doxygen/xml/index.xml
#    ${CMAKE_CURRENT_SOURCE_DIR}/doc/documentation.i
#  COMMENT "Generate python documentation from doxygen output"
#)
#
# Tell that we need to generate doxymentation before building python package
#add_dependencies(python_doc doxygen)

# TODO: r_build should call first r_swig (to be created) followed by r_doc (and then r_install)
add_custom_target(r_doc
  COMMAND ${CMAKE_COMMAND} -E echo "Launching create_doc.R"
  COMMAND ${R_EXECUTABLE} CMD BATCH --no-echo --no-timing --no-save --no-restore ${CMAKE_CURRENT_SOURCE_DIR}/doc/create_doc.R ${R_PACKAGE_ROOT_FOLDER}/create_doc.out
  COMMAND ${CMAKE_COMMAND} -E cat ${R_PACKAGE_ROOT_FOLDER}/create_doc.out
  COMMENT "Generating documentation for R package"
  WORKING_DIRECTORY ${R_PACKAGE_DESTINATION_FOLDER}
  VERBATIM
)

####CUSTOM WRAPPER
add_custom_command(OUTPUT ${CMAKE_BINARY_DIR}/swig/generated_r.i
COMMAND ${CMAKE_COMMAND} -E touch ${CMAKE_BINARY_DIR}/swig/generated_r.i
COMMENT "Ensure that generated_r.i exists")
  
add_custom_target(r_wrap_macro
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/../tools/scripts/macrogenerate.py r ${CMAKE_BINARY_DIR}/swig/generated_r.i
  "r"
  COMMENT "Generate macro wrapper package"
  )


# Tell that we need to generate doxymentation before building R package (if required)
if (BUILD_DOC)
  add_dependencies(r_doc doxygen)
endif()

######################################
# WRAPPER GENERATION & BUILD

# R package interface file
set(SRC rgstlearn.i)

# Location for SWIG target language code
set(CMAKE_SWIG_OUTDIR                 ${R_PACKAGE_R_FOLDER})
# Location for SWIG wrapper library
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY    ${R_PACKAGE_SRC_FOLDER}) # Use LIBRARY (gcc)!
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY    ${R_PACKAGE_SRC_FOLDER}) # Use RUNTIME (msvc)!

# SWIG interface source file properties
set_source_files_properties(${SRC} PROPERTIES
  CPLUSPLUS ON
  # Swig wrapper will be regenerated if C++ headers are changed
  # https://stackoverflow.com/a/70375654/3952924 (CMake 3.20 minimum)
  USE_SWIG_DEPENDENCIES TRUE
  # C++ headers will be found by %import swig directive
  USE_TARGET_INCLUDE_DIRECTORIES TRUE
)

# Tell to CMake that ${SRC} depends on some other files
set_property(SOURCE ${SRC} PROPERTY
  DEPENDS
    ${PROJECT_SOURCE_DIR}/swig/swig_exp.i
    ${PROJECT_SOURCE_DIR}/swig/swig_inc.i
    #${CMAKE_CURRENT_SOURCE_DIR}/doc/documentation.i
)

# SWIG flags
#set(CMAKE_SWIG_FLAGS -debug-classes)
#set(CMAKE_SWIG_FLAGS -debug-tmsearch)
set(CMAKE_SWIG_FLAGS)
# Note: Cannot use -doxygen feature because all documentation is in C++ body files
# Note: Use -E option for seeing SWIG preprocessor output

# Forward C++ compile definitions to rgstlearn.i
set_property(
  SOURCE ${SRC}
  PROPERTY COMPILE_DEFINITIONS
           $<TARGET_PROPERTY:${PROJECT_NAME}::shared,COMPILE_DEFINITIONS>)

# SWIG target (generate wrapper & R script then build [so,dll] library)
swig_add_library(
  r_build
  LANGUAGE   r
  TYPE       SHARED
  SOURCES    ${SRC}
)

if (APPLE)
  # To prevent 'import' from crashing under Mac (swig output suffix is .dynlib by default) 
  set_target_properties(r_build PROPERTIES SUFFIX ".so")
endif()

# Reuse shared library objects
target_sources(
  r_build
  PRIVATE
  $<TARGET_OBJECTS:${PROJECT_NAME}::shared>
)
target_include_directories(r_build
  PRIVATE
  $<TARGET_PROPERTY:${PROJECT_NAME}::shared,INCLUDE_DIRECTORIES>
)
target_link_libraries(r_build
  PRIVATE
  $<TARGET_PROPERTY:${PROJECT_NAME}::shared,LINK_LIBRARIES>
)

# Set some properties on SWIG target
set_target_properties(r_build PROPERTIES
  # Modify the generated library name
  OUTPUT_NAME ${PROJECT_NAME}
)

# Transmit the fact that we are linking to static ${PROJECT_NAME} library
target_compile_definitions(r_build PRIVATE -D${PROJECT_NAME_UP}_STATIC_DEFINE)

# Prevent some SWIG C++/R warnings
target_compile_options(r_build PRIVATE
  -Wno-free-nonheap-object
  -Wno-maybe-uninitialized
  -Wno-unused-variable
  -Wno-unused-parameter
  -Wno-uninitialized
  -Wno-free-nonheap-object
  -Wno-cast-function-type
  -Wno-cast-function-type-mismatch)

# Find appropriate R.h
target_include_directories(r_build PRIVATE ${R_INCLUDE_DIRS})

# Link to R's libraries (needed for MSYS and MacOS)
target_link_libraries(r_build PRIVATE ${R_LIBRARIES})

# Tell that we need to generate Rd files before building R package (if required)
if (BUILD_DOC)
  add_dependencies(r_build r_doc)
endif()

if(USE_BOOST_SPAN)
  target_compile_definitions(r_build PUBLIC USE_BOOST_SPAN)
endif()

add_dependencies(r_build r_wrap_macro)

add_dependencies(r_build generate_output)
########################################
# PACKAGING

# Generate NAMESPACE automatically for each configuration
# First step: replace variables (@VAR@)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/NAMESPACE.in 
  ${CMAKE_CURRENT_BINARY_DIR}/NAMESPACE.in
  @ONLY
)
# Second step: replace generator expression ($<>)
file(GENERATE
  OUTPUT ${R_PACKAGE_DESTINATION_FOLDER}/NAMESPACE
  INPUT ${CMAKE_CURRENT_BINARY_DIR}/NAMESPACE.in
)

# Replace DESCRIPTION automatically for each configuration
# First step: replace variables (@VAR@)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/DESCRIPTION.in 
  ${CMAKE_CURRENT_BINARY_DIR}/DESCRIPTION.in
  @ONLY
)
# Second step: replace generator expression ($<>)
file(GENERATE
  OUTPUT ${R_PACKAGE_DESTINATION_FOLDER}/DESCRIPTION
  INPUT ${CMAKE_CURRENT_BINARY_DIR}/DESCRIPTION.in
)

# Copy each modules to package *destination* folder for each configuration
file(GLOB MODULES ${CMAKE_CURRENT_SOURCE_DIR}/modules/*.R)
foreach(MODULE ${MODULES})
  cmake_path(GET MODULE FILENAME MODULE_NAME)
  file(GENERATE
    OUTPUT ${R_PACKAGE_R_FOLDER}/${MODULE_NAME}
    INPUT ${MODULE}
  )
endforeach()

# Copy data folder to package *destination* folder for each configuration
#file(GLOB DATASETS_RDA ${CMAKE_CURRENT_SOURCE_DIR}/data/*.rda)
#foreach(DATASET_RDA ${DATASETS_RDA})
#  cmake_path(GET DATASET_RDA FILENAME DATASET_RDA_NAME)
#  file(COPY
#    OUTPUT ${R_PACKAGE_DATA_FOLDER}/${DATASET_RDA_NAME}
#    INPUT ${DATASET_RDA_NAME}
#  )
#endforeach()

######################################
# INSTALLATION

# Add a custom target for R package installation
# TODO: Do also installation each time README and below files are modified
if (MINGW)
  set(ARCHIVE_FILE_NAME     "${PROJECT_NAME}_${PROJECT_VERSION}.zip")
  set(DST_ARCHIVE_FILE_NAME "${PROJECT_NAME}_${PROJECT_VERSION}.zip")
else()
  if(APPLE)
    set(ARCHIVE_FILE_NAME "${PROJECT_NAME}_${PROJECT_VERSION}.tgz")
    set(DST_ARCHIVE_FILE_NAME "${PROJECT_NAME}_${PROJECT_VERSION}.tgz")
  else()
    set(ARCHIVE_FILE_NAME "${PROJECT_NAME}_${PROJECT_VERSION}_R_x86_64-pc-linux-gnu.tar.gz")
    set(DST_ARCHIVE_FILE_NAME "${PROJECT_NAME}_${PROJECT_VERSION}.tar.gz")
  endif()
endif()
add_custom_target(r_install
  # Empty Makevars (compilation already done by r_build target) - only 64 bits (no multi arch)
  COMMAND ${CMAKE_COMMAND} -E touch ${R_PACKAGE_SRC_FOLDER}/Makevars
  COMMAND ${CMAKE_COMMAND} -E touch ${R_PACKAGE_SRC_FOLDER}/Makevars.win
  COMMAND ${CMAKE_COMMAND} -E touch ${R_PACKAGE_SRC_FOLDER}/Makevars.ucrt
  COMMAND ${R_EXECUTABLE} CMD INSTALL --build --no-multiarch ${R_PACKAGE_DESTINATION_FOLDER}
#  COMMAND ${R_EXECUTABLE} CMD INSTALL --no-multiarch --no-byte-compile ${R_PACKAGE_DESTINATION_FOLDER}
  COMMAND ${CMAKE_COMMAND} -E copy ${ARCHIVE_FILE_NAME} ${R_PACKAGE_ROOT_FOLDER}/${DST_ARCHIVE_FILE_NAME}
  COMMENT "Installing R package"
  VERBATIM
)

# Add r_build target dependency to install target
add_dependencies(r_install r_build)


