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
#pragma once

#include "gstlearn_export.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Vector.hpp"
#include "Fractures/Description.hpp"
#include "Matrix/MatrixRectangular.hpp"

#define NPART 5
#define NBYFRAC 7
#define NBYWOUT 8

class Environ;
class DbGrid;

class GSTLEARN_EXPORT FracList: public AStringable
{
public:
  FracList(int ndisc = 1000,
           bool flag_check = true,
           double low0 = EPSILON8,
           double low1 = EPSILON6,
           double eps = EPSILON3);
  FracList(const FracList& r);
  FracList& operator=(const FracList& r);
  virtual ~FracList();

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  int getNFracs() const { return (int) _descs.size(); }

  int simulate(const Environ* environ,
               bool flag_sim_layer,
               bool flag_sim_fract,
               int seed,
               bool verbose,
               int nlayers_in,
               const VectorDouble& elevations);
  MatrixRectangular fractureExport() const;
  MatrixRectangular layinfoExport() const { return _layinfo; };
  static FracList* fractureImport(const VectorDouble& frac_segs,
                                  const VectorDouble& layinfo = VectorDouble(),
                                  int nfamilies = 0);

  Description& getDescs(int i) { return _descs[i]; }
  void addDescription(const Description& description = Description()) { _descs.push_back(description); }

  int fractureToBlock(DbGrid *dbgrid,
                      double xmax,
                      VectorDouble& permtab,
                      double perm_mat,
                      double perm_bench,
                      int ndisc);
  VectorDouble fractureToWell(int nval,
                              const VectorDouble& well,
                              double xmax,
                              const VectorDouble& permtab,
                              int *nint,
                              int *ncol);
  int fractureWellToBlock(DbGrid *dbgrid,
                          int col_perm,
                          int col_fluid,
                          int flag_fluid,
                          double val_fluid,
                          const VectorDouble& wellout,
                          int nval,
                          int ndisc,
                          bool verbose);
  VectorDouble fractureExtractLength(int ifam, double cote, double dcote);
  VectorDouble fractureExtractDist(int ifam, double cote, double dcote);

private:
  int getRank(int ifam, int shift) const { return (1 + ifam * NPART + shift); }
  void setMemLayer(int i, double value)             { _layinfo.setValue(i,0,value); }
  void setMemTheta1(int i, int ifam, double value)  { _layinfo.setValue(i,getRank(ifam,0),value); }
  void setMemTheta2(int i, int ifam, double value)  { _layinfo.setValue(i,getRank(ifam,1),value); }
  void setMemPropsur(int i, int ifam, double value) { _layinfo.setValue(i,getRank(ifam,2),value); }
  void setMemFrac(int i, int ifam, double value)    { _layinfo.setValue(i,getRank(ifam,3),value); }
  void setMemTotal(int i, int ifam, double value)   { _layinfo.setValue(i,getRank(ifam,4),value); }
  double getMemLayer(int i)                         { return _layinfo.getValue(i,0); }

  VectorDouble _layersManage(const Environ& environ, double *y0);
  VectorDouble _layersRead(int nlayers_in,
                           const VectorDouble& elevations,
                           double *y0);
  void _manage(int mode);
  int _fracAdd(int ifrac,
               int ifam,
               double xx,
               double cote,
               double thick,
               double orient,
               double* xp);
  void _checkFractureIntersect(double cote, int ifrac0);
  bool _belongToLayer(const Description& desc,
                      double cote,
                      double *xd,
                      double *yd,
                      double *xe,
                      double *ye);
  double _layerIntensity(const Family& family,
                         double thick);
  void _generateDensity(const Environ& environ,
                        const Family& family,
                        int ifam,
                        double cote,
                        VectorDouble& denstab);
  void _correctDensity(const Family& family,
                       int ifam,
                       double cote,
                       VectorDouble& denstab);
  double _deriveIntensity(double theta1,
                          double thetap,
                          double propsur);
  double _extendFractures(const Family& family,
                          int ifam,
                          double cote,
                          double thick,
                          VectorDouble& denstab);
  bool _sameFaultSide(const Environ& environ, int ifault0, double x0);
  double _densityUpdate(const Fault& fault,
                        int side,
                        int ifam,
                        double cote,
                        double xx);
  bool _densityCumulate(const char *title,
                        bool flag_print,
                        const VectorDouble& denstab,
                        double *totdens);
  void _updateRepulsion(double x0, double range, VectorDouble& denstab);
  bool _fractureInterrupt(const Family& family,
                          const Description& desc,
                          double thick);
  double _faultAbscissae(const Fault& fault, double cote);
  double _cubic(double h);
  double _fractureExtension(const Description& desc, double cote, double dcote);
  int _simulateFractures(const Environ& environ,
                         const Family& family,
                         int ifam,
                         double cote,
                         double thick,
                         double theta,
                         VectorDouble& denstab);
  int _getDiscretizedRank(double cumdens, const VectorDouble& denstab);
  int _getEndPointCount() const;
  bool _isValidDisc(int idisc);

  void _plungeSegment(DbGrid *dbgrid,
                      int iptr,
                      double delta,
                      double value,
                      double x1,
                      double y1,
                      double x2,
                      double y2);
  void _welloutAdd(VectorDouble& wellout,
                   double x,
                   double y,
                   int ifrac,
                   int ip,
                   int family,
                   double perm);
  void _trajAdd(VectorDouble& traj, double x, double y);
  void _plungeSegmentGradual(DbGrid *dbgrid,
                             int iptr,
                             double delta,
                             VectorDouble& traj,
                             double perm1,
                             double perm2,
                             double range);

private:
  // Array of fracture descriptions
  std::vector<Description> _descs;
  MatrixRectangular _layinfo;
  int _nlayers;
  // The number of discretization steps used to establish the fracture density
  int _ndisc;
  // The option for checking the fracture intersect or not
  bool _flagCheck;
  // The lower value used for checking the fracture repulsion for the central cell
  double _low0;
  // The lower value used for checking the fracture repulsion for the peripheral cells
  double _low1;
  double _xorigin;
  double _step;
  double _eps;
  bool _verbose;
};