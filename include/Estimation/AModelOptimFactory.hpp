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
#pragma once

#include "Estimation/Vecchia.hpp"
#include "Estimation/Likelihood.hpp"
#include "Estimation/AModelOptimNew.hpp"
#include "Model/ModelOptimVMap.hpp"
#include "Model/ModelOptimVario.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

class AModelOptimNew;

class GSTLEARN_EXPORT AModelOptimFactory {
  public:
    /**
     * @brief Instantiate the appropriate AModelOptimNew object based on the provided parameters.
     *
     * @param model ModelGeneric pointer representing the model to be optimized.
     * @param db Db pointer containing experimental data (for standard Likelihood).
     * @param vario Vario pointer containing the variogram (for variogram fitting).
     * @param dbmap DbGrid containing the grid map (for variogram map fitting).
     * @param constraints Constraints (optional)
     * @param mop ModelOptimParam containing fitting options.
     * @param nb_neighVecchia Number of Vecchia neighbors to use (for Vecchia Likelihood).
     * @return AModelOptimNew*
     */
    static AModelOptimNew* create(ModelGeneric* model,
                                  const Db* db,
                                  Vario* vario,
                                  const DbGrid* dbmap,
                                  Constraints* constraints,
                                  const ModelOptimParam& mop,
                                  int nb_neighVecchia = ITEST)
    {
      if (db != nullptr)
      {
        if (nb_neighVecchia != ITEST)
          return Vecchia::createForOptim(model, db, nb_neighVecchia);
        return Likelihood::createForOptim(model, db);
      }
      if (dbmap != nullptr)
        return ModelOptimVMap::createForOptim(model, dbmap, constraints, mop);
      if (vario != nullptr)
        return ModelOptimVario::createForOptim(model, vario, constraints, mop);

      return nullptr;
    }
};
