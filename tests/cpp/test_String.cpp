/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
/**
 * This test is meant to check the string manipulations
 */

#include "Basic/AStringable.hpp"
#include "Basic/File.hpp"
#include "Basic/String.hpp"
#include "Basic/VectorNumT.hpp"

using namespace gstlrn;

int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  // Testing string conversion Old to New style functions
  char buf_char[] = "abcde";
  message("Buffer (char) = %s\n", buf_char);
  char dst_char[10];
  VectorUChar dst_vec(100);

  // Testing gslStrcpy
  message("Testing gslStrcpy:\n");
  message("- Expected result: abcde\n");

  gslStrcpy(dst_char, buf_char);
  message("- Result = %s\n", dst_char);

  // Testing gslStrcat
  message("Testing gslStrcat:\n");
  message("- Expected result: abcdefgh\n");

  char newbuf_char[] = "fgh";
  gslStrcat(dst_char, newbuf_char);
  message("- Result = %s\n", dst_char);

  // Testing gslSPrintf
  message("Testing gslSPrintf:\n");
  message("- Expected result: A+fgh\n");
  gslSPrintf(dst_char, "A+%s", newbuf_char);
  message("- Result(char) = %s\n", dst_char);

  // Testing gslSPrintf2
  message("Testing gslSPrintf2:\n");
  message("- Expected result: BC+fgh\n");

  gslSPrintf2(dst_vec, "BC+%s", newbuf_char);
  message("- Result(char) = %s\n", dst_vec.data());

  // Testing gslAddSPrintf2
  message("Testing gslAddSPrintf2:\n");
  Id one = 1;
  gslAddSPrintf2(dst_vec, " added %d times", one);
  message("- Result(char) = %s\n", dst_vec.data());

  // Testing gslStrcat2
  message("Testing gslStrcat2:\n");
  gslStrcat2(dst_vec, " (Tested successfully)");
  message("- Result(char) = %s\n", dst_vec.data());

  // Testing gslStrcpy2
  message("Testing gslStrcpy2:\n");
  gslStrcpy2(dst_vec, "After a copy");
  message("- Result(char) = %s\n", dst_vec.data());
  return 0;
}
