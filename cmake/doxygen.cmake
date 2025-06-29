if (NOT BUILD_DOC)
  return()
endif()

# TODO : Do not regenerate doxymentation if nothing has changed in the source code
find_package(Doxygen REQUIRED)

# Configure doxyfile

set(DOXYGEN_OUTPUT_DIRECTORY doxygen)
set(DOXYGEN_PROJECT_BRIEF "Geostatistics & Machine Learning toolbox | <a href=https://gstlearn.org>https://gstlearn.org</a>")
set(DOXYGEN_PROJECT_LOGO ${CMAKE_SOURCE_DIR}/doc/logos/gstlearn_logo_blue_th.png)
set(DOXYGEN_MULTILINE_CPP_IS_BRIEF YES)
set(DOXYGEN_EXTRACT_ALL YES)
set(DOXYGEN_EXTRACT_STATIC YES)
set(DOXYGEN_WARN_NO_PARAMDOC YES)
set(DOXYGEN_USE_MDFILE_AS_MAINPAGE README.md)
set(DOXYGEN_EXCLUDE ${CMAKE_SOURCE_DIR}/include/geoslib_old_f.h
                    ${CMAKE_SOURCE_DIR}/include/geoslib_f_private.h
                    ${CMAKE_SOURCE_DIR}/include/geoslib_d_private.h)
set(DOXYGEN_EXCLUDE_SYMBOLS "FORWARD_METHOD_NON_CONST"
                             "FORWARD_METHOD_CONST")

set(DOXYGEN_VERBATIM_HEADERS NO)
set(DOXYGEN_GENERATE_HTML YES)
set(DOXYGEN_HTML_TIMESTAMP YES)
set(DOXYGEN_HTML_HEADER ${CMAKE_SOURCE_DIR}/cmake/doxygen.header.html)
set(DOXYGEN_LAYOUT_FILE ${CMAKE_SOURCE_DIR}/cmake/DoxygenLayout.xml)
set(DOXYGEN_GENERATE_XML YES)
set(DOXYGEN_GENERATE_TREEVIEW YES)
set(DOXYGEN_MAX_INITIALIZER_LINES 1000) # For very long macros
set(DOXYGEN_MACRO_EXPANSION YES)
set(DOXYGEN_EXPAND_ONLY_PREDEF NO)

set(DOXYGEN_ENABLE_PROCESSING YES)
set(DOXYGEN_MACRO_EXPANSION YES)
set(DOXYGEN_EXPAND_ONLY_PREDEF YES)
set(DOXYGEN_PREDEFINED "protected=private")
set(DOXYGEN_EXTRACT_PRIVATE NO)

set(DOXYGEN_QUIET YES)
set(DOXYGEN_HAVE_DOT NO) # Put NO to reduce generation time (keep YES for UML or better graphs)
# Uncomment if you prefer UML graphs (need DOXYGEN_HAVE_DOT YES)
#set(DOXYGEN_HIDE_UNDOC_RELATIONS NO)
#set(DOXYGEN_UML_LOOK YES)
#set(DOXYGEN_TEMPLATE_RELATIONS YES)

set(DOXYGEN_USE_MATHJAX YES)

# https://stackoverflow.com/questions/25290453/how-do-i-add-a-footnote-in-doxygen
set(DOXYGEN_ALIASES tooltip{1}=\"\\latexonly\\footnote\\{\\1\\}\\endlatexonly\\htmlonly<sup title=\'\\1\'>*</sup>\\endhtmlonly\")

### Following lines create fake hpp to automatically generate documentation for macros
# FORWARD_METHOD_CONST and FORWARD_METHOD_NON_CONST

set(GENERATED_HPP_FILES ${CMAKE_BINARY_DIR}/doxygen/generated_hpp)
set(DOXYGEN_SCRIPT ${CMAKE_SOURCE_DIR}/tools/scripts/macrodoc.py)

add_custom_command(OUTPUT ${GENERATED_HPP_FILES}
  COMMAND ${CMAKE_COMMAND} -E make_directory ${GENERATED_HPP_FILES})

add_custom_target(doc_macro
  COMMAND ${CMAKE_COMMAND} -E make_directory ${GENERATED_HPP_FILES}
  COMMAND ${Python3_EXECUTABLE} ${DOXYGEN_SCRIPT} ${GENERATED_HPP_FILES}
  "python"
  COMMENT "Generate macro wrapper doc"
  )

########

# Add target for generating the doxymentation
doxygen_add_docs(doxygen
                 ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/src 
                 ${CMAKE_SOURCE_DIR}/README.md
                 ${GENERATED_HPP_FILES}
                 COMMENT "Generate doxygen documentation")

add_dependencies(doxygen doc_macro)