/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "Basic/Utilities.hpp"
#include "Basic/String.hpp"
#include "Basic/AException.hpp"
#include "LithoRule/RuleShift.hpp"
#include "LithoRule/Rule.hpp"
#include "LithoRule/Node.hpp"
#include "Model/Model.hpp"
#include "geoslib_f.h"
#include "geoslib_enum.h"

/**
 * Definition of the Lithotype RuleShift in the case of Shift
 * @param shift Vector defining the shift
 */
RuleShift::RuleShift(int nfacies, const VectorDouble& shift)
    : Rule(),
      _shDsup(0.),
      _shDown(0.),
      _slope(0.),
      _tgte(TEST),
      _incr(TEST),
      _shift(shift),
      _xyz(),
      _ind1(),
      _ind2()
{
  setModeRule(RULE_SHIFT);
  VectorString nodnames = buildNodNames(nfacies);
  setMainNodeFromNodNames(nodnames);
}

RuleShift::RuleShift(const RuleShift& m)
    : _shDsup(m._shDsup),
      _shDown(m._shDown),
      _slope(m._slope),
      _tgte(m._tgte),
      _incr(m._incr),
      _shift(m._shift)
{
}

RuleShift& RuleShift::operator=(const RuleShift& m)
{
  if (this != &m)
  {
    _shDsup = m._shDsup;
    _shDown = m._shDown;
    _slope = m._slope;
    _tgte = m._tgte;
    _incr = m._incr;
    _shift = m._shift;
  }
  return *this;
}

RuleShift::~RuleShift()
{
}

int RuleShift::deSerializeSpecific()
{
  _shift.resize(3);
  if (_recordRead("Slope for Shadow Rule", "%lf", &_slope)) return 1;
  if (_recordRead("Lower Threshold for Shadow Rule", "%lf", &_shDown)) return 1;
  if (_recordRead("Upper Threshold for Shadow Rule", "%lf", &_shDsup)) return 1;
  if (_recordRead("Shift along first direction", "%lf", &_shift[0]))  return 1;
  if (_recordRead("Shift along second direction", "%lf", &_shift[1])) return 1;
  if (_recordRead("Shift along third direction", "%lf", &_shift[2]))  return 1;
  return 0;
}

void RuleShift::serializeSpecific() const
{
  double slope = (FFFF(_slope)) ? 0. : _slope;
  _recordWrite("%lf", slope);
  double shdown = (FFFF(_shDown)) ? 0. : _shDown;
  _recordWrite("%lf", shdown);
  double shdsup = (FFFF(_shDsup)) ? 0. : _shDsup;
  _recordWrite("%lf", shdsup);
  _recordWrite("#", "Parameters for Shadow option");
  _recordWrite("%lf", _shift[0]);
  _recordWrite("%lf", _shift[1]);
  _recordWrite("%lf", _shift[2]);
  _recordWrite("#", "Parameters for Shift option");
}

String RuleShift::displaySpecific(int flagProp, int flagThresh) const
{
  std::stringstream sstr;
  sstr << toVector("Translation Vector",_shift) << std::endl;
  sstr << "(With the current option, only the first GRF is used)" << std::endl;
  sstr << getMainNode()->nodePrint(flagProp, flagThresh);
  return sstr.str();
}

/****************************************************************************/
/*!
**  Define the particularities of the PGS model
**
** \return  Error return code
**
** \param[in]  db              Db structure
** \param[in]  dbprop          Db structure used for proportions
** \param[in]  model           Model structure (only used for shift option)
** \param[in]  flag_grid_check 1 if grid is compulsory; 0 otherwise
**                             (only for SHIFT)
** \param[in]  flag_stat       1 for stationary; 0 otherwise
**
*****************************************************************************/
int RuleShift::particularities(Db *db,
                               const Db *dbprop,
                               Model *model,
                               int flag_grid_check,
                               int flag_stat)
{
  int ndim = (model != (Model *) NULL) ? model->getDimensionNumber() : 0;
  VectorDouble wxyz(ndim);
  double rhoval;

  /* Dispatch */

  _xyz.resize(ndim);
  double hval = 0.;
  for (int idim = 0; idim < ndim; idim++)
  {
    _xyz[idim] = _shift[idim];
    hval += _xyz[idim] * _xyz[idim];
  }
  hval = sqrt(hval);

  /* Calculate the covariance between the two GRF */

  for (int idim = 0; idim < ndim; idim++)
    wxyz[idim] = _xyz[idim];
  if (model->getVariableNumber() == 1)
    model_evaluate(model, 0, 0, -1, 0, 1, 0, 0, 0, MEMBER_LHS, 1, wxyz, &hval,
                   &rhoval);
  else
    model_evaluate(model, 0, 1, -1, 0, 1, 0, 0, 0, MEMBER_LHS, 1, wxyz, &hval,
                   &rhoval);
  setRho(rhoval);

  /* Translate the shift into grid increments */

  if (_st_shift_on_grid(db, ndim, flag_grid_check)) return (1);
  return (0);
}

int RuleShift::_st_shift_on_grid(Db *db, int ndim, int flag_grid_check)
{
  VectorDouble xyz(ndim);
  _ind1.resize(ndim);

  if (db == (Db *) NULL || ! is_grid(db))
  {
    if (! flag_grid_check) return(0);
    messerr("The shift Rule requires a Grid Db");
    return(1);
  }

  for (int idim=0; idim<ndim; idim++)
    _xyz[idim] = _shift[idim] + db->getX0(idim);

  (void) point_to_grid(db,_xyz.data(),-1,_ind1.data());

  /* Check that the translation is significant */

  int ntot = 0;
  for (int idim=0; idim<ndim; idim++)
    ntot += ABS(_ind1[idim]);
  if (ntot <= 0)
  {
    messerr("The shift of the Lithotype Rule cannot be rendered");
    messerr("using the Output Grid characteristics");
    return(1);
  }
  return(0);
}

bool RuleShift::checkModel(const Model* model, int nvar) const
{
  if (model == (Model *) NULL)
  {
    messerr("No Model is provided");
    return false;
  }
  if (nvar > 0 && model->getVariableNumber() != nvar)
  {
    messerr("The number of variables in the Model (%d) does not match",
            model->getVariableNumber());
    messerr(" the number of variables in the Db (%d)", nvar);
    return false;
  }
  return true;
}

/****************************************************************************/
/*!
**  Combine the underlying GRF into a facies value
**
** \return  Error return code
**
** \param[in]  propdef    Props structure
** \param[in]  dbout      Db output structure
** \param[in]  flag_used  1 if the gaussian is used; 0 otherwise
** \param[in]  ipgs       Rank of the PGS
** \param[in]  isimu      Rank of the simulation
** \param[in]  nbsimu     Number of simulations
**
** \remark Attributes LOC_FACIES and LOC_SIMU are mandatory
**
*****************************************************************************/
int RuleShift::gaus2facResult(PropDef *propdef,
                              Db *dbout,
                              int *flag_used,
                              int ipgs,
                              int isimu,
                              int nbsimu)
{
  int    ndim,iech,jech,idim,igrf,icase;
  double t1min,t1max,t2min,t2max,facies,y[2];

  /* Initializations */

  check_mandatory_attribute("rule_gaus2fac_result",dbout,LOC_FACIES);
  check_mandatory_attribute("rule_gaus2fac_result",dbout,LOC_SIMU);
  ndim   = dbout->getNDim();
  VectorDouble xyz(ndim);
  VectorInt ind1(ndim);
  VectorInt ind2(ndim);

  /* Processing the translation */

  for (iech=0; iech<dbout->getSampleNumber(); iech++)
  {
    if (! dbout->isActive(iech)) continue;

    facies = TEST;
    for (igrf=0; igrf<2; igrf++) y[igrf] = TEST;
    icase = get_rank_from_propdef(propdef, ipgs, 0);
    y[0] = dbout->getSimvar(LOC_SIMU, iech, isimu, 0, icase, nbsimu, 1);
    if (FFFF(y[0])) break;

    if (rule_thresh_define(propdef, dbout, this, ITEST, iech, isimu, nbsimu, 1,
                           &t1min, &t1max, &t2min, &t2max)) return 1;
    db_index_sample_to_grid(dbout, iech, ind2.data());
    for (idim = 0; idim < ndim; idim++)
      ind2[idim] -= ind1[idim];
    jech = db_index_grid_to_sample(dbout, ind2.data());
    if (jech >= 0)
      y[1] = dbout->getSimvar(LOC_SIMU, jech, isimu, 0, icase, nbsimu, 1);
    else
      y[1] = TEST;
    facies = getFaciesFromGaussian(y[0], y[1]);

    /* Combine the underlying GRFs to derive Facies */

    dbout->setSimvar(LOC_FACIES,iech,isimu,0,ipgs,nbsimu,1,facies);
  }
  return 0;
}

/****************************************************************************/
/*!
**  Set the bounds and possibly add replicates
**
** \return  Error return code
**
** \param[in]  propdef    PropDef structure
** \param[in]  dbin       Db structure
** \param[in]  dbout      Db grid structure
** \param[in]  isimu      Rank of the simulation (PROCESS_CONDITIONAL)
** \param[in]  igrf       Rank of the GRF
** \param[in]  ipgs       Rank of the GS
** \param[in]  nbsimu     Number of simulations (PROCESS_CONDITIONAL)
**
*****************************************************************************/
int RuleShift::evaluateBounds(PropDef *propdef,
                              Db *dbin,
                              Db *dbout,
                              int isimu,
                              int igrf,
                              int ipgs,
                              int nbsimu)
{
  int    iech,jech,nadd,nech,idim,facies,nstep;
  double t1min,t1max,t2min,t2max,s1min,s1max,s2min,s2max;

  /* Initializations */

  if (dbin == (Db *) NULL) return(0);
  nadd = nstep = 0;
  nech = dbin->getSampleNumber();

  /* Dispatch */

  if (igrf == 1) return (0);

  /* Loop on the data */
  for (iech = 0; iech < nech; iech++)
  {
    /* Convert the proportions into thresholds for data point */
    if (!dbin->isActive(iech)) continue;
    facies = (int) dbin->getVariable(iech, 0);
    if (rule_thresh_define(propdef, dbin, this, facies, iech, isimu, nbsimu, 1,
                           &t1min, &t1max, &t2min, &t2max)) return (1);
    dbin->setLowerBound(iech, get_rank_from_propdef(propdef, ipgs, igrf),
                        t1min);
    dbin->setUpperBound(iech, get_rank_from_propdef(propdef, ipgs, igrf),
                        t1max);
    if (facies == SHADOW_ISLAND) continue;

    /* Add one replicate */
    jech = dbin->addSamples(1, 0.);
    if (jech < 0) return (1);

    /* Set the coordinates of the replicate */
    for (idim = 0; idim < dbin->getNDim(); idim++)
      dbin->setCoordinate(jech, idim,
                          dbin->getCoordinate(iech, idim) - _shift[idim]);

    /* Can the replicate be added */
    if (replicateInvalid(dbin, dbout, jech))
    {
      dbin->deleteSample(jech);
      return (1);
    }

    /* Convert the proportions into thresholds for replicate */
    if (rule_thresh_define(propdef, dbin, this, facies, jech, isimu, nbsimu, 1,
                           &s1min, &s1max, &s2min, &s2max))
    {
      dbin->deleteSample(jech);
      return (1);
    }

    /* Set the attributes of the replicate */
    if (facies == SHADOW_WATER) dbin->setVariable(jech, 0, SHADOW_WATER);
    if (facies == SHADOW_SHADOW) dbin->setVariable(jech, 0, SHADOW_ISLAND);
    dbin->setLowerBound(jech, get_rank_from_propdef(propdef, ipgs, igrf),
                        s2min);
    dbin->setUpperBound(jech, get_rank_from_propdef(propdef, ipgs, igrf),
                        s2max);
    nadd++;
  }

  if (igrf == 0 && nadd > 0)
  {
    message("Initial count of data = %d\n", nech);
    message("Number of replicates  = %d\n", nadd);
  }
  return (0);
}