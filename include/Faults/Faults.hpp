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

#include "../Basic/PolyLine2D.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"

class GSTLEARN_EXPORT Faults: public AStringable, public ASerializable
{
public:
  Faults();
  Faults(const Faults& r);
  Faults& operator=(const Faults& r);
  virtual ~Faults();

  /// Interface for AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  static Faults* createFromNF(const String& neutralFilename, bool verbose = true);
  int getNFaults() const { return (int) _faults.size(); }
  void addFault(const PolyLine2D& fault);

protected:
  /// Interface for ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os, bool verbose = false) const override;
  String _getNFName() const override { return "Faults"; }

private:
  std::vector<PolyLine2D> _faults;
};