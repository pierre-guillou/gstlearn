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
#include "geoslib_old_f.h"

#include "Enum/ECst.hpp"
#include "Enum/ELoadBy.hpp"
#include "Enum/ESpaceType.hpp"

#include "API/SPDE.hpp"

#include "Basic/File.hpp"
#include "Basic/Law.hpp"
#include "Basic/OptCst.hpp"
#include "Basic/OptDbg.hpp"
#include "Covariances/CovContext.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Db/DbStringFormat.hpp"
#include "LinearOp/PrecisionOpMatrix.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Model/Model.hpp"
#include "Space/ASpaceObject.hpp"

using namespace gstlrn;

/****************************************************************************/
/*!
** This test is meant to check the elaboration of the CovGradient class
** and its use in Kriging (or Depth adn Gradient)
**
*****************************************************************************/
int main(int argc, char* argv[])

{
  bool flag_print = false;
  bool flag_save  = true;

  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  /***********************/
  /* 1 - Initializations */
  /***********************/
  int ndim = 2;
  defineDefaultSpace(ESpaceType::RN, ndim);
  ASerializable::setPrefixName("test_Gradient-");

  // Setup constants

  OptDbg::reset();
  OptCst::define(ECst::NTCAR, 10.);
  OptCst::define(ECst::NTDEC, 6.);

  return (0);
}
