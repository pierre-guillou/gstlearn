# Using GLOB could be a bad idea (but use it for the tests)
# https://stackoverflow.com/questions/32411963/why-is-cmake-file-glob-evil

# With GLOB:
file(GLOB TEST_SOURCES *.cpp)
# Without GLOB:
#set(TEST_SOURCES 
#    test_Cholesky.cpp)
#  test_SPDEDrift.cpp)

# Remove test_ppmt (not yet functional)
#list(FILTER TEST_SOURCES EXCLUDE REGEX "test_ppmt.cpp")

# Exclude one test from output comparison
set(EXCLUDE_FROM_DIFF "test_a_template")

# Deactivate HDF5 test if not used
if(NOT USE_HDF5)
  list(FILTER TEST_SOURCES EXCLUDE REGEX "bench_H5.cpp")
endif()

# Generation folder
if (IS_MULTI_CONFIG)
  cmake_path(APPEND CMAKE_CURRENT_BINARY_DIR $<CONFIG>
             OUTPUT_VARIABLE TEST_DST_DIR)
else()
  cmake_path(APPEND CMAKE_CURRENT_BINARY_DIR ${CMAKE_BUILD_TYPE}
             OUTPUT_VARIABLE TEST_DST_DIR)
endif()
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${TEST_DST_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${TEST_DST_DIR})

# Define the prepare target to create the output directory for logs
if (WIN32)
  # Need to copy C++ shared library to tests directory
  if(MSYS)
    add_custom_target(prepare_cpp
      COMMAND ${CMAKE_COMMAND} -E make_directory "${TEST_DST_DIR}"
      COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:${PROJECT_NAME}::shared> ${TEST_DST_DIR}
      # have both lib${PROJECT_NAME}.dll and ${PROJECT_NAME}.dll in ${TEST_DST_DIR}
      # this fixes a bug on MSYS when building R wrappers with shared
      # objects before linking C++ test binaries
      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${TEST_DST_DIR}/lib${PROJECT_NAME}.dll ${TEST_DST_DIR}/${PROJECT_NAME}.dll
    )
  else()
    add_custom_target(prepare_cpp
      COMMAND ${CMAKE_COMMAND} -E make_directory "${TEST_DST_DIR}"
      COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:${PROJECT_NAME}::shared> ${TEST_DST_DIR}
    )
  endif()
else()
  add_custom_target(prepare_cpp
    COMMAND ${CMAKE_COMMAND} -E make_directory "${TEST_DST_DIR}"
  )
endif()

# Compile each test executable
set(TARGETS_EXE "")
foreach(TEST_SOURCE_FILE ${TEST_SOURCES})
    # Retrieve source file name without extension (will become executable name)
    get_filename_component(TEST_NAME ${TEST_SOURCE_FILE} NAME_WE)
    # Add to executable targets list
    list(APPEND TARGETS_EXE ${TEST_NAME})
    # Define sources list for the target executable
    add_executable(${TEST_NAME} ${TEST_SOURCE_FILE})
    # Link each test to shared library
    target_link_libraries(${TEST_NAME} PRIVATE ${PROJECT_NAME}::shared)
    # Pass CMAKE_SOURCE_DIR to executables (for the `getTestData' function)
    target_compile_definitions(${TEST_NAME} PRIVATE -DCMAKE_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
    # Trigger the prepare target each time a test is compiled
    add_dependencies(${TEST_NAME} prepare_cpp)
    # Trigger the build of the test with the target build_tests
    add_dependencies(build_tests ${TEST_NAME})

    # test_Nlopt.cpp explicitly calls NLopt so it must link with it (gstlearn::shared
    # links it as PRIVATE to hide it from projects using it).
    if("${TEST_NAME}" STREQUAL "test_Nlopt")
      target_link_libraries(${TEST_NAME} PRIVATE NLopt::nlopt)
    endif()
endforeach(TEST_SOURCE_FILE ${TEST_SOURCES})

# Display test output in case of failure
set(CTEST_OUTPUT_ON_FAILURE ON)

# Deactivate a particular test (which is under construction)
list(FILTER TARGETS_EXE EXCLUDE REGEX "test_a_template")

# Run each registered executable
foreach(TARGET_EXE ${TARGETS_EXE})
    # Run the test (and generate *.out in ${TEST_DST_DIR})
    add_test(NAME ${TARGET_EXE}
             COMMAND ${TARGET_EXE}
             WORKING_DIRECTORY ${TEST_DST_DIR})
    if (NOT ${TARGET_EXE} STREQUAL ${EXCLUDE_FROM_DIFF})
      # Compare the output result (diff output is shown when CTEST_OUTPUT_ON_FAILURE=1)
      set(REF_FILE ${CMAKE_CURRENT_SOURCE_DIR}/output/${TARGET_EXE}.ref)
      # Use maybe a specific output log for the local computer
      cmake_path(APPEND CMAKE_CURRENT_SOURCE_DIR myoutput ${TARGET_EXE}.ref
                  OUTPUT_VARIABLE REF_FILE_LOCAL)
      if(EXISTS ${REF_FILE_LOCAL})
        set(REF_FILE ${REF_FILE_LOCAL})
      endif()
      # Use maybe a specific output log for MSVC
      if (MSVC)
        set(REF_FILE_MSVC ${CMAKE_CURRENT_SOURCE_DIR}/output/${TARGET_EXE}_msvc.ref)
        if(EXISTS ${REF_FILE_MSVC})
          set(REF_FILE ${REF_FILE_MSVC})
        endif()
      endif()
      # Use maybe a specific output log for MinGW
      if (MINGW)
        set(REF_FILE_MSYS ${CMAKE_CURRENT_SOURCE_DIR}/output/${TARGET_EXE}_msys.ref)
        if(EXISTS ${REF_FILE_MSYS})
          set(REF_FILE ${REF_FILE_MSYS})
        endif()
      endif()
      # Use maybe a specific output log for clang
      if (CLANG)
        set(REF_FILE_CLANG ${CMAKE_CURRENT_SOURCE_DIR}/output/${TARGET_EXE}_clang.ref)
        if(EXISTS ${REF_FILE_CLANG})
          set(REF_FILE ${REF_FILE_CLANG})
        endif()
      endif()
      # Use git diff (cross-platform) with:
      # - no-index (no need to have files under git control)
      # - ignore lines starting with #NO_DIFF#
      # - ignore ALL whitespaces (https://stackoverflow.com/a/65462731)
      add_test(NAME ${TARGET_EXE}_cmp
               COMMAND git diff --no-index -I "#NO_DIFF#" --ignore-blank-lines --ignore-all-space ${REF_FILE} ${TEST_DST_DIR}/${TARGET_EXE}.out)
    endif()
endforeach(TARGET_EXE ${TARGETS_EXE})

# Create the check target to launch the tests
# Look parent CMakeLists for MY_CTEST_COMMAND definition
add_custom_target(check_cpp
  COMMAND ${MY_CTEST_COMMAND} DEPENDS ${TARGETS_EXE}
  JOB_POOL console
)
