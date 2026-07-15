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
#include "Basic/ASerializable.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/VectorNumT.hpp"
#include "DataBase/DbData.hpp"

using namespace gstlrn;

/****************************************************************************/
/*!
 ** Main Program
 **
 ** This program is meant to check the manipulation of DbData
 **
 *****************************************************************************/
int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);
  ASerializable::setPrefixName("test_DataBase-");

  DbData data{};
  data.printContents();

  // Checking the different types of Columns that can be added to a DbData
  mestitle(1, "Adding Columns to a DbData");
  data.addColumn("hello", VectorDouble{1., 2., 3.}, RoleID{ERole::X});
  data.addColumn("world", VectorInt{5, 6, 7}, RoleID{ERole::Z});
  data.addColumn(
    "foobar", VectorString{"foo", "bar", "baz"}, RoleID{ERole::Z, 1});
  data.addColumn("gstlearn", VectorBool{true, false, true});
  data.addColumnEmpty<VectorDouble>("Bonjour", 0, 5, RoleID(ERole::F), 3.);
  data.addColumn("MyVar", VH::sequence(15, 3, 2), RoleID(ERole::Z, 2), 5);

  // Checking the use of Neutral Files
  mestitle(1, "Saving and recovering a DbData from a Neutral File");
  data.printContents("Before Saving in a Neutral File");
  std::cout << "C0: " << data.getColumn<VectorDouble>("hello") << std::endl;
  std::cout << "C1: " << data.getColumn<VectorInt>("world") << std::endl;
  std::cout << "C2: " << data.getColumn<VectorString>("foobar") << std::endl;
  std::cout << "C3: " << data.getColumn<VectorBool>("gstlearn") << std::endl;
  std::cout << "C4: " << data.getColumn<VectorDouble>("Bonjour") << std::endl;
  std::cout << "C5: " << data.getColumn<VectorInt>("MyVar") << std::endl;
  data.dumpToNF("test_DataBase.NF");

  auto* data2 = DbData::createFromNF("test_DataBase.NF", true);
  data2->printContents("After recovering from the Neutral File");
  std::cout << "C0: " << data2->getColumn<VectorDouble>("hello") << std::endl;
  std::cout << "C1: " << data2->getColumn<VectorInt>("world") << std::endl;
  std::cout << "C2: " << data2->getColumn<VectorString>("foobar") << std::endl;
  std::cout << "C3: " << data2->getColumn<VectorBool>("gstlearn") << std::endl;
  std::cout << "C4: " << data2->getColumn<VectorDouble>("Bonjour") << std::endl;
  std::cout << "C5: " << data2->getColumn<VectorInt>("MyVar") << std::endl;
  delete data2;

  mestitle(1, "Checking aliases");
  message("Various ways to get value of 'MyVar' for version 0 at sample 2\n");
  auto is = 2;
  auto ind = 2;
  auto iv = 0;
  message("- By Name: %d\n", data.getValue<Id>("MyVar", is).value_or(-1));
  message(
    "- By Name and Version: %d\n",
    data.getValue<Id>({"MyVar", iv}, is).value_or(-1));
  message(
    "- Role: %d\n", data.getValue<Id>(RoleID{ERole::Z, ind}, is).value_or(-1));

  // Checking the different manners to refer to a Column in a DbData
  mestitle(1, "Different manners to refer to a Column");
  message("- by Name: %s\n", data.getName("hello").c_str());
  message("- by Name and Version: %s\n", data.getName({"Bonjour", 1}).c_str());

  message("- by Index: %s\n", data.getName(4).c_str());
  message("- by Index and Version: %s\n", data.getName({4, 1}).c_str());

  auto roleid = RoleID(ERole::F);
  message("- by RoleID: %s\n", data.getName(roleid).c_str());
  message("- by RoleID and Version: %s\n", data.getName({roleid, 0}).c_str());

  message("- by Role: %s\n", data.getName(ERole::Z).c_str());

  mestitle(1, "Retrieving and modifying values in a DbData");
  auto isample = 2;

  message("Initial values\n");
  message("%lf\n", data.getValue<double>(0, isample).value_or(-1.));
  message("%d\n", data.getValue<int>(1, isample).value_or(-1));
  message("%s\n", data.getValue<String>(2, isample).value_or("failed").c_str());

  data.setValue(0, isample, 4.);
  data.setValue(1, isample, 8);
  data.setValue(2, isample, "foobar");

  message("Values after modification\n");
  message("%lf\n", data.getValue<double>(0, isample).value_or(-1.));
  message("%d\n", data.getValue<int>(1, isample).value_or(-1));
  message("%s\n", data.getValue<String>(2, isample).value_or("failed").c_str());

  mestitle(1, "Adding Columns with the same Name and/or Role");
  data.printContents("Initial");
  data.addColumn("hello", VectorDouble{1., 2., 3.}, RoleID{ERole::X, 10});
  data.printContents(
    "\nAfter adding: Same Name ('hello'), Same Role ('X') and different Rank "
    "('10')");
  data.addColumn("world", VectorDouble{1., 2., 3.}, RoleID{ERole::X, 0});
  data.printContents(
    "\nAfter adding: Same Name ('world'), Same Role ('X') and different Rank "
    "('0')");
  data.addColumn("world", VectorDouble{1., 2., 3.}, RoleID{ERole::X, 10});
  data.printContents(
    "\nAfter adding: Same Name ('world'), Same Role ('X') and different Rank "
    "('10')");

  mestitle(1, "Modifying a Column with the wrong type");
  data.getColumn<VectorDouble>("hello").dump("Initial Column");
  data.setValue("hello", isample, 123.);
  data.getColumn<VectorDouble>("hello").dump("After valid change");
  message("Trying to modify a Column with the wrong type\n");
  data.setValue("hello", isample, "String");

  mestitle(1, "Retrieving each Column (in VectorDouble format if possible)");
  data.printContents();
  for (Id icol = 0; icol < data.getNCols(); ++icol)
  {
    data.getColumn<VectorDouble>(icol).dump(
      "Column " + std::to_string(icol + 1), false);
  }

  mestitle(1, "Retrieving series of columns for various criteria");
  data.removeAllColumns();
  data.addColumn("hello", VectorDouble{1., 2., 3.}, RoleID{ERole::X});
  data.addColumn("hellobis", VectorInt{5, 6, 7}, RoleID{ERole::X, 1});
  data.addColumn("helloback", VectorDouble{5, 6, 7}, RoleID{ERole::X, 2});
  data.printContents();

  message("- Matching the name 'hellob*', the following columns are found:\n");
  auto colIDs = data.getColIDs("hellob*");
  for (auto& colID: colIDs)
    message("Matching Column: %s\n", colID.getDescr().c_str());

  message("- Matching the role 'X', the following columns are found:\n");
  colIDs = data.getColIDs(ERole::X);
  for (auto& colID: colIDs)
    message("Matching Column: %s\n", colID.getDescr().c_str());

  return 0;
}
