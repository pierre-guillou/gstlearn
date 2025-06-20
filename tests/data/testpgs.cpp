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
#include "geoslib_f.h"

#include "Core/Ascii.hpp"
#include "Basic/Law.hpp"
#include "Basic/Limits.hpp"
#include "Basic/File.hpp"
#include "LithoRule/RuleProp.hpp"
#include "LithoRule/Rule.hpp"
#include "LithoRule/RuleStringFormat.hpp"
#include "Space/ASpaceObject.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Db/DbStringFormat.hpp"
#include "Variogram/Vario.hpp"
#include "Model/Model.hpp"
#include "Neigh/NeighUnique.hpp"

/*********************/
/* Program principal */
/*********************/

int main(int argc, char *argv[])
{
  char       filename[BUFFER_LENGTH];
  Db        *dbin;
  DbGrid    *dbout;
  Vario     *vario;
  Model     *model[2][2];
  NeighUnique *neighU;
  Rule      *rule[2];
  Option_VarioFit options;
  RuleStringFormat rulefmt;
  DbStringFormat dbfmt;
  double  total;
  int     i,j,lec,nbsimu,seed,nbtuba,npgs,ntot,nfac[2];
  int     flag_vario,flag_grid,iatt_z,iatt_ind,ifac,nclass;
  VectorDouble props;
  RuleProp* ruleprop;
  static int    niter   = 100;
  static int    nboot   = 10;
  static int    verbose = 0;

  /* Initializations */

  npgs     = ntot = 0;
  dbin     = nullptr;
  dbout    = nullptr;
  vario    = nullptr;
  neighU   = nullptr;
  ruleprop = nullptr;
  for (i=0; i<2; i++)
  {
    rule[i] = nullptr;
    for (j=0; j<2; j++)
      model[i][j] = nullptr;
  }

  /* Standard output redirection to file */

  StdoutRedirect sr("Result.out", argc, argv, 2);

  /* Create the output name (for storage of dump files) */

  VectorString subparts = separateKeywords(argv[1]);
  int nargs             = (int)subparts.size();
  String outname =
    concatenateStrings("", subparts[nargs - 2], subparts[nargs - 1], "-");
  ASerializable::setPrefixName(outname);
  setInternalDebug(true);

  /* Getting the Study name */

  if (argc < 2) messageAbort("Wrong number of arguments");
  ascii_study_define(argv[1]);

  /* Define the environment */

  ascii_filename("Environ",0,0,filename);
  ascii_environ_read(filename,verbose);

  /* Define the data */

  ascii_filename("Data",0,0,filename);
  dbin = Db::createFromNF(filename,verbose);
  if (dbin == nullptr) goto label_end;
  iatt_z = dbin->getUIDByLocator(ELoc::Z,0);
  dbfmt.setFlags(true, false, true, true, true);
  dbin->display(&dbfmt);

  /* Define the Default Space according to the Dimension of the Input Db */

  defineDefaultSpace(ESpaceType::RN,dbin->getNDim());

  /* Define the variogram (optional) */
  
  ascii_filename("Vario",0,0,filename);
  vario = Vario::createFromNF(filename,verbose);
  flag_vario = (vario != nullptr);

  /* Define the output grid file */

  ascii_filename("Grid",0,0,filename);
  dbout = DbGrid::createFromNF(filename,verbose);
  flag_grid = (dbout != nullptr);

  /* Define the rules */

  for (i=lec=0; i<2; i++)
  {

    /* Read the rule */

    ascii_filename("Rule",i,0,filename);
    rule[i] = Rule::createFromNF(filename,verbose);
    if (rule[i] == nullptr) continue;
    npgs++;
    rule[i]->display();
    nfac[i] = rule[i]->getNFacies();

    /* Define the models */
    
    for (j=0; j<2; j++,lec++)
    {
      if (! rule[i]->isYUsed(j)) continue;
      ascii_filename("Model",lec,0,filename);
      model[i][j] = Model::createFromNF(filename,verbose);
      if (model[i][j] == nullptr) goto label_end;
    }

    /* Calculate the experimental variogram of indicators */
    
    if (flag_vario)
    {

      /* Define the indicators */

      nclass   = nfac[i];
      iatt_ind = dbin->getNColumn();
      Limits limits = Limits(nclass);
      limits.toIndicator(dbin);
      dbin->setLocatorsByUID(nclass,iatt_ind,ELoc::Z, 0);
      
      /* Calculate the experimental variograms */
      
      vario->compute(dbin, ECalcVario::VARIOGRAM);
      vario->display();
      ascii_filename("Vario",0,1,filename);
      if (! vario->dumpToNF(filename,verbose))
        messageAbort("ascii_vario_write");
      
      /* Delete the indicator variables */
      
      dbin->clearLocators(ELoc::Z);
      for (ifac=0; ifac<nclass; ifac++)
        dbin->deleteColumnByUID(iatt_ind+ifac);
      dbin->setLocatorByUID(iatt_z,ELoc::Z, 0);
    }
  }

  /* Look for simulations */

  ascii_filename("Simu",0,0,filename);
  ascii_simu_read(filename,verbose,&nbsimu,&nbtuba,&seed);

  /* Create the proportions */

  ntot  = (npgs == 1) ? nfac[0] : nfac[0] * nfac[1];
  props.resize(ntot);
  total = 0.;
  for (i=0; i<ntot; i++)
  {
    props[i] = law_uniform(0.,1.);
    total   += props[i];
  }
  for (i=0; i<ntot; i++) props[i] /= total;

  /* Define the neighborhood */

  neighU = NeighUnique::create();

  /* Perform the Pluri-Gaussian Simulations */

  if (flag_grid)
  {
    if (npgs == 1)
    {
      ruleprop = RuleProp::createFromRule(rule[0],props);
      if (simpgs(dbin,dbout,ruleprop,model[0][0],model[0][1],
                 neighU,nbsimu,seed,0,0,0,0,nbtuba,nboot,niter,1)) goto label_end;
    }
    else
    {
      ruleprop = RuleProp::createFromRules(rule[0],rule[1],props);
      if (simbipgs(dbin,dbout,ruleprop,
                   model[0][0],model[0][1],model[1][0],model[1][1],
                   neighU,nbsimu,seed,0,0,0,0,nbtuba,nboot,niter,1)) goto label_end;
    }
    dbfmt.setFlags(true, false, true, true, true);
    dbout->display(&dbfmt);
  }

  /* Serialization of results (for visual check) */

  (void) dbout->dumpToNF("Result");

  /* Core deallocation */

label_end:
  delete dbin;
  delete dbout;
  for (i=0; i<2; i++)
  {
    rule[i] = rule_free(rule[i]);
    for (j=0; j<2; j++)
      delete model[i][j];
  }
  delete ruleprop;
  delete neighU;

  return(0);
}
