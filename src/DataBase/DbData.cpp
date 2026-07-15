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
#include "DataBase/DbData.hpp"
#include "Basic/SerializeHDF5.hpp"

namespace gstlrn
{
  DbData* DbData::createFromNF(const String& NFFilename, bool verbose)
  {
    auto* dbdata = new DbData;
    if (dbdata->_fileOpenAndDeserialize(NFFilename, verbose)) return dbdata;
    delete dbdata;
    return nullptr;
  }

#ifdef HDF5
  bool DbData::serializeH5(H5::Group& grp) const
  {
    auto dbG = grp.createGroup("DbData");
    SerializeHDF5::writeValue(dbG, "NColumn", getNCols());

    for (Id i = 0; i < getNCols(); i++)
    {
      auto colG = dbG.createGroup("Column_" + std::to_string(i));

      // Information contained in DbData
      SerializeHDF5::writeValue(colG, "Role", _roleIDs[i].getRole().getValue());
      SerializeHDF5::writeValue(colG, "Rank", _roleIDs[i].getIndex());

      // Information contained in DbCol
      if (!_cols[i].serializeH5(colG)) return false;
    }

    return true;
  }

  bool DbData::deserializeH5(H5::Group& grp)
  {
    auto dbG = grp.openGroup("DbData");

    Id ncols;
    Id roleValue;
    Id rank;

    // Read number of columns
    SerializeHDF5::readValue(dbG, "NColumn", ncols);

    // Clear previous contents of the Data Base
    removeAllColumns();

    for (Id i = 0; i < ncols; i++)
    {
      auto colG = dbG.openGroup("Column_" + std::to_string(i));

      // Read DbData information
      SerializeHDF5::readValue(colG, "Role", roleValue);
      SerializeHDF5::readValue(colG, "Rank", rank);

      // Read DbCol information
      _cols.push_back(DbCol::createEmpty());

      if (!_cols.back().deserializeH5(colG)) return false;

      _roleIDs.emplace_back(ERole::fromValue(roleValue), rank);
    }

    return true;
  }
#endif

  bool DbData::renameColumn(ColID&& colid, const String& newName)
  {
    const auto icol = _getColumnIndex(colid);
    if (!icol) return false;
    this->_cols[*icol].setName(newName);
    return true;
  }

  void DbData::removeColumn(ColID&& colid)
  {
    const auto icol = _getColumnIndex(colid);
    if (!icol) return;

    if (static_cast<Id>(this->_cols.size()) > (*icol))
    {
      this->_cols.erase(this->_cols.begin() + (*icol));
      this->_roleIDs.erase(this->_roleIDs.begin() + (*icol));
    }
  }

  void DbData::removeAllColumns()
  {
    this->_cols.clear();
    this->_roleIDs.clear();
  }

  /**
   * @brief Identify and returns a reference on the Column of interest
   *
   * @param colid Column Indentifier
   * @return std::optional<std::reference_wrapper<DbCol>>
   */
  std::optional<std::reference_wrapper<DbCol>>
    DbData::_identifyColumn(ColID&& colid)
  {
    const auto icol = _getColumnIndex(colid);
    if (!icol) return std::nullopt;
    return {this->_cols[*icol]};
  }

  std::optional<std::reference_wrapper<const DbCol>>
    DbData::_identifyColumn(ColID&& colid) const
  {
    const auto icol = _getColumnIndex(colid);
    if (!icol) return std::nullopt;
    return {this->_cols[*icol]};
  }

  bool DbData::hasColumn(ColID&& colid) const
  {
    const auto icol = _getColumnIndex(colid);
    return static_cast<bool>(icol);
  }

  /**
   * @brief Get the Name object
   *
   * @param colid Column Indentifier
   * @return String
   */
  String DbData::getName(ColID&& colid) const
  {
    const auto icol = _getColumnIndex(colid);
    return icol ? _cols[*icol].getName() : String();
  }

  /**
   * @brief Get the Index of the Column
   *
   * @param colid Column indentifier
   * @return Id
   */
  Id DbData::getICol(ColID&& colid) const
  {
    const auto icol = _getColumnIndex(colid);
    return icol ? *icol : -1;
  }

  /**
   * @brief Get the RoleID of the Column
   *
   * @param colid Column indentifier
   * @return RoleID
   */
  RoleID DbData::getRoleID(ColID&& colid) const
  {
    const auto icol = _getColumnIndex(colid);
    return icol ? _roleIDs[*icol] : RoleID();
  }

  ColID DbData::getColID(const ColID& colid) const
  {
    const auto icol = _getColumnIndex(colid);
    if (!icol) return ColID();

    ColID colIDout;
    colIDout.setName(_cols[*icol].getName());
    colIDout.setICol(*icol);
    colIDout.setRoleID(_roleIDs[*icol]);
    return colIDout;
  }

  /**
   * @brief Return a set of Column Identifiers starting from a Column Name
   *
   * @param name Column Name (with regexp facility)
   * @return std::vector<ColID>
   */
  std::vector<ColID> DbData::getColIDs(const String& name) const
  {
    // Looking for matching names
    VectorString matchNames = expandList(_getNames(), name);

    // Loop to create the list of Column Identifiers
    std::vector<ColID> colIDs;
    for (const auto& matchName: matchNames)
    {
      const auto colID = getColID(ColID(matchName));
      if (colID.getICol() >= 0) colIDs.push_back(colID);
    }
    return colIDs;
  }

  /**
   * @brief Return a set of Column Identifiers starting from a set of Column Names
   *
   * @param names Column Names (with regexp facility)
   * @return std::vector<ColID>
   */
  std::vector<ColID> DbData::getColIDs(const VectorString& names) const
  {
    VectorString matchNames = expandList(_getNames(), names);
    std::vector<ColID> colIDs;
    for (const auto& matchName: matchNames)
    {
      const auto colID = getColID(ColID(matchName));
      if (colID.getICol() >= 0) colIDs.push_back(colID);
    }
    return colIDs;
  }

  /**
   * @brief Returns a set of Column Identifiers matching a given Role
   *
   * @param role Target Role to be searched
   * @return std::vector<ColID>
   */
  std::vector<ColID> DbData::getColIDs(const ERole& role) const
  {
    std::vector<ColID> colIDs;
    for (const auto& id: this->_roleIDs)
    {
      if (id.getRole() == role)
      {
        const auto colID = getColID(ColID(id));
        if (colID.getICol() >= 0) colIDs.push_back(colID);
      }
    }
    return colIDs;
  }

  Id DbData::getNSamples() const
  {
    if (getNCols() <= 0) return 0;
    return _cols[0].getNSamples();
  }

  void DbData::removeRole(ColID&& colid)
  {
    const auto icol = _getColumnIndex(colid);
    if (!icol) return;
    _roleIDs[*icol].removeRole();
  }

  Id DbData::getNVersions(ColID&& colid)
  {
    const auto icol = _getColumnIndex(colid);
    if (!icol) return 0;
    return _cols[*icol].getNVersions();
  }

  /**
   * @brief Erase the Role of all Columns with the given Role
   *
   * @param role Target Role to be erased
   */
  void DbData::clearRole(const ERole& role)
  {
    std::vector<ColID> ids = getColIDs(role);
    for (const auto& id: ids)
    {
      _roleIDs[id.getICol()].removeRole();
    }
  }

  /**
   * @brief Produces a summary of the DbData content
   */
  void DbData::printContents(const String& title) const
  {
    Id ncols = _cols.size();

    if (!title.empty()) std::cout << title << '\n';
    std::cout << "The Data Base contains " << ncols << " columns of "
              << getNSamples() << " samples\n";
    for (Id icol = 0; icol < ncols; ++icol)
    {
      const auto& c = _cols[icol];
      const auto& id = _roleIDs[icol];

      std::cout << "Column " << icol << "/" << ncols;
      std::cout << " : " << c.getDescr();
      std::cout << " {Role: " << id.getDescr() << "}";
      if (c.forbidNA()) std::cout << " [NA forbidden]";
      std::cout << std::endl;
    }
  }

  /**
   * @brief Returns the Column index corresponding to the given Column Indentifier (ColID)
   *
   * @param colid Column Indentifier
   * @return std::optional<Id>
   *
   * @remark: The Column Indentifier (ColID) is searched in the following order:
   * - by Column Name
   * - by Column Role (and Rank)
   * - by Column index
   */
  std::optional<Id> DbData::_getColumnIndex(const ColID& colid) const
  {
    auto ncol = getNCols();

    // Try to identify by Column Name
    const String& localName = colid.getName();
    if (!localName.empty())
    {
      for (Id icol = 0; icol < ncol; ++icol)
      {
        if (_cols[icol].getName() == localName) return icol;
      }
    }

    // Try to identify by Column Role
    if (colid.getRole() != ERole::UNDEFINED)
    {
      const RoleID& roleID = colid.getRoleID();
      for (Id icol = 0; icol < ncol; ++icol)
      {
        if (_roleIDs[icol].match(roleID)) return icol;
      }
    }

    // Try to identify by Column rank
    if (colid.getICol() >= 0) return colid.getICol();

    messerr("Column does not exist.");
    return std::nullopt;
  }

  /**
   * @brief Check if the new name is compatible with existing ones
   *
   * @param name Name of the new column to be added (possibly modified)
   *
   * @remark If the new name matches an already existing one, it is modified to be unique
   */
  void DbData::_updateName(String& name) const
  {
    // Establish the list of already existing names
    VectorString proposedNames = _getNames();
    auto ncol = static_cast<Id>(proposedNames.size());

    // Add the new proposal to the list of already existing names
    proposedNames.push_back(name);

    // Modify the 'ncol' proposal and retrive the modified value
    correctNamesForDuplicates(proposedNames, ncol);
    name = proposedNames[ncol];
  }

  /**
   * @brief Check if the RankID of the new column is compatible with existing ones
   *
   * @param roleID RankID of the new column to be added (possibly modified)
   *
   * @remark If the Role of the new Column is already present in the already defined ones:
   * - if the Rank of the new Column matches the one of the old matching Column:
   *   the Old matching Column is moved to an UNDEFINED Role (and a Rank set to 0).
   *   the New Column keeps its Role and Rank unchanged
   * - if the Rank of the new Column does not match the one of the old matching Column,
   *  this rank is calculated as the largest Rank found in matching Columns incremented by 1.
   */
  void DbData::_updateRoleID(RoleID& roleID)
  {
    const auto& newRole = roleID.getRole();
    const auto newRank = roleID.getIndex();
    Id rankMin = -1;

    // Look for already existing Columns with the same RoleID
    for (auto& id: this->_roleIDs)
    {
      if (id.getRole() == newRole)
      {
        // Same role already exists
        auto oldRank = id.getIndex();
        if (oldRank == newRank)
        {
          // Same role and same rank already exists:
          // Move the old one to UNDEFINED; keep the new one unchanged
          id.setRole(ERole::UNDEFINED);
          id.setIndex(0);
          return;
        }

        // Update the Minimum index
        if (oldRank > rankMin) rankMin = oldRank;
      }
    }

    if (rankMin >= 0)
    {
      // Same role already exist but with different ranks:
      // Set the new rank to the largest rank found + 1
      roleID.setIndex(rankMin + 1);
    }
    else
    {
      // No same role already exists: set the rank to 0 whatever the input rank was
      roleID.setIndex(0);
    }
  }

  // /**
  //  * @brief Check if the new Category Column is Valid when forbidNA is true
  //  *
  //  * @param tab VectorCategory of the new column to be added
  //  */
  // bool DbData::_checkForbidNA(const VectorCategory& tab)
  // {
  //   for (size_t i = 0; i < tab.size(); i++)
  //   {
  //     if (!tab[i].has_value())
  //     {
  //       messerr("Column forbids NA values, but the input tab contains some.");
  //       return false;
  //     }
  //   }
  //   return true;
  // }

  /**
   * @brief Get the Names of all the Columns
   *
   * @return VectorString
   */
  VectorString DbData::_getNames() const
  {
    VectorString names;
    for (const auto& col: _cols)
    {
      names.push_back(col.getName());
    }
    return names;
  }

  void DbData::addSamples(Id nadd, const double valinit)
  {
    if (nadd <= 0) return;
    for (Id icol = 0, ncol = getNCols(); icol < ncol; icol++)
    {
      _cols[icol].addSamples(nadd, valinit);
    }
  }

  void DbData::deleteSample(Id idel)
  {
    for (Id icol = 0, ncol = getNCols(); icol < ncol; icol++)
    {
      _cols[icol].deleteSample(idel);
    }
  }

  void DbData::_checkVersion(Id& nversion)
  {
    if (nversion <= 0)
    {
      messerr(
        "The number of versions (%d) must be strictly positive. It has "
        "been set to 1.",
        nversion);
      nversion = 1;
    }
  }

} // namespace gstlrn
