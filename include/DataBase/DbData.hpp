/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2025) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "Basic/ASerializable.hpp"
#include "DataBase/ColID.hpp"
#include "DataBase/DbCol.hpp"
// #include "DataBase/Dictionary.hpp"
#include "DataBase/RoleID.hpp"
// #include "DataBase/VectorCategory.hpp"
#include "gstlearn_export.hpp"

#include <functional>
#include <optional>

namespace gstlrn
{
  /**
   * @brief The new heterogeneous DbData class is a container for a set of Columns (DbCol)
   *        Each Column can be identified by:
   *        - its name (unique)
   *        - by its RoleID (Role + Rank)
   *        - by its index in the DbData.
   *       The DbData class provides methods to add, remove, and access Columns,
   *        as well as to get and set values in the Columns.
   *
   *       Each Column in turn:
   *       - is characterized by its type: Double, Int, String, etc.
   *       - can have multiple versions
   *
   *       Limitation: all the Columns of the same DbData must have the same number of samples.
   */

  class GSTLEARN_EXPORT DbData: public ASerializable
  {
  public:
#ifndef SWIG
    class ColProxy;
    class ValueProxy;
#endif

    static DbData* createFromNF(const String& NFFilename, bool verbose);

    /// ASerializable interface
    String getNFName() const override { return "DbData"; }
#ifdef HDF5
    bool deserializeH5(H5::Group& grp) override;
    bool serializeH5(H5::Group& grp) const override;
#endif

    /// List of template functions

    /**
     * @brief Initialize a new Column and fill it constantly with default values
     *
     * @tparam VectorType
     * @param name Name of the Column
     * @param nsamples Number of samples (default = 0)
     * @param nversion Number of versions (default = 1)
     * @param roleID RoleID of the Column
     * @param valinit Value to fill the Column with (or std::nullopt for default)
     * @param forbidNA Whether to forbid NA values in the new Column (default = false)
     *
     * @remark: The argument 'nsample' is only needed if the current column is the first one
     * of the current Data Base.
     */
    template<typename VectorType>
    void addColumnEmpty(
      const String& name,
      Id nsamples = 0,
      Id nversion = 1,
      const RoleID& roleID = RoleID(),
      std::optional<typename VectorType::value_type> valinit = std::nullopt,
      bool forbidNA = false)
    // const Dictionary* dict = nullptr)
    {
      if (getNCols() > 0) nsamples = getNSamples();
      if (nsamples <= 0)
      {
        messerr("The number of samples (%d) must be positive.", nsamples);
        return;
      }
      // auto array =
      //   _createEmptyVector<VectorType>(nsamples * nversion, valinit, dict);
      auto array = _createEmptyVector<VectorType>(nsamples * nversion, valinit);
      addColumn(name, std::move(array), roleID, nversion, forbidNA);
    }

    /**
     * @brief Add a new Column
     *
     * @tparam VectorType
     * @param name Name of the Column
     * @param array Array of values to fill the Column with
     * @param roleID RoleID of the Column (optional)
     * @param nversion Number of versions (default = 1)
     * @param forbidNA Whether to forbid NA values in the new Column (default = false)
     */
    template<typename VectorType>
    void addColumn(
      String&& name,
      VectorType&& array,
      const RoleID& roleID = RoleID(),
      Id nversion = 1,
      bool forbidNA = false)
    {
      // Checking the validity of the number of versions
      _checkVersion(nversion);

      // Check the validity of the new name
      auto nameLocal = name;
      for (const auto& col: this->_cols)
      {
        if (col.getName() == name)
        {
          _updateName(nameLocal);
          break;
        }
      }

      // Check the validity of the new Role (if defined)
      auto roleIDLocal = roleID;
      _updateRoleIDAddition(-1, roleIDLocal);

      // Check that input does not contain NA values if forbidNA is true
      if (forbidNA)
      {
        if (!_checkForbidNA(array)) return;
      }

      // Add the new column to the list of columns
      this->_cols.emplace_back(
        std::move(nameLocal), std::forward<VectorType>(array), nversion,
        forbidNA);

      // Add the new RoleID to the list of RoleIDs
      _roleIDs.emplace_back(roleIDLocal);
    }

    template<typename VectorType>
    void addColumn(
      const String& name,
      VectorType&& array,
      const RoleID& roleID = RoleID(),
      Id nversion = 1,
      bool forbidNA = false)
    {
      // Checking the validity of the number of versions
      _checkVersion(nversion);

      // Check the validity of the new name
      auto nameLocal = name;
      for (const auto& col: this->_cols)
      {
        if (col.getName() == name)
        {
          _updateName(nameLocal);
          break;
        }
      }

      // Check the validity of the new Role (if defined)
      auto roleIDLocal = roleID;
      _updateRoleIDAddition(-1, roleIDLocal);

      // Check that input does not contain NA values if forbidNA is true
      if (forbidNA)
      {
        if (!_checkForbidNA(array)) return;
      }

      // Add the new column to the list of columns
      this->_cols.emplace_back(
        std::move(nameLocal), std::forward<VectorType>(array), nversion,
        forbidNA);

      // Add the new RoleID to the list of RoleIDs
      _roleIDs.emplace_back(roleIDLocal);
    }

#ifndef SWIG
    /**
     * @brief Returns the Value of a Column for a given sample and version
     *
     * @tparam T Type of the value
     * @param colid Identification of the Column
     * @param isample Index of the sample
     * @return std::optional<T>
     */
    template<typename T>
    std::optional<T> getValue(ColID&& colid, Id isample) const
    {
      const auto array = this->_identifyColumn(std::move(colid));
      if (!array) return std::nullopt;

      const Id version = colid.getVersion();
      return array->get().getValue<T>(isample, version);
    }

    /**
     * @brief Set the Value of a Column for a given sample
     *
     * @tparam T Type of the value
     * @param colid Identification of the Column
     * @param isample Index of the sample
     * @param value Value to be assigned
     * @return true
     * @return false
     */
    template<typename T>
    bool setValue(ColID&& colid, const Id isample, const T& value)
    {
      const auto array = this->_identifyColumn(std::move(colid));
      if (!array) return false;

      const Id version = colid.getVersion();
      return array->get().setValue<T>(isample, version, value);
    }

    /**
     * @brief returns the whole Target Column
     *
     * @tparam VectorType Type of returned values
     * @param colid Identification of the Column
     * @return VectorType&
     *
     * @remark If a Column contains multiple versions, the returned values correspond to all values of all versions concatenated.
     */
    template<typename VectorType>
    VectorType& getColumn(ColID&& colid)
    {
      static VectorType empty{};

      auto col = this->_identifyColumn(std::move(colid));
      if (!col) return empty;

      auto vec = col->get().template getVector<VectorType>();
      if (!vec) return empty;

      return vec->get();
    }

    template<typename VectorType>
    const VectorType& getColumn(ColID&& colid) const
    {
      static const VectorType empty{};

      auto col = this->_identifyColumn(std::move(colid));
      if (!col) return empty;

      auto vec = col->get().template getVector<VectorType>();
      if (!vec) return empty;

      return vec->get();
    }

    /**
     * @brief returns the Target Column for a given version
     *
     * @tparam VectorType Type of returned values
     * @param colid Identification of the Column
     * @param iversion Version of the Column to get
     * @return std::span<typename VectorType::value_type>
     */
    template<typename VectorType>
    std::span<typename VectorType::value_type>
      getVersion(ColID&& colid, Id iversion = 0)
    {
      auto col = this->_identifyColumn(std::move(colid));
      if (!col) return {};

      auto span = col->get().template getVersion<VectorType>(iversion);
      if (!span) return {};

      return *span;
    }

    template<typename VectorType>
    std::span<const typename VectorType::value_type>
      getVersion(ColID&& colid, Id iversion = 0) const
    {
      auto col = this->_identifyColumn(std::move(colid));
      if (!col) return {};

      auto span = col->get().template getVersion<VectorType>(iversion);
      if (!span) return {};

      return *span;
    }

    /**
     * @brief Set the contents of the Target Column
     *
     * @param colid Identification of the Column
     * @param values Values to set
     * @return bool
     *
     * @remark If a Column contains multiple versions, the values must correspond to all values of all versions concatenated.
     */
    template<typename VectorType>
    bool setColumn(ColID&& colid, const VectorType& values)
    {
      auto col = _identifyColumn(std::move(colid));
      if (!col) return false;

      return col->get().template setVector<VectorType>(values);
    }

    /**
     * @brief Set the contents of the Target Column (for a given version)
     *
     * @param colid Identification of the Column
     * @param values Values to set
     * @param iversion Version of the Column to set
     * @return bool
     */
    template<typename VectorType>
    bool setVersion(ColID&& colid, const VectorType& values, Id iversion = 0)
    {
      auto col = this->_identifyColumn(std::move(colid));
      if (!col) return false;

      return col->get().template setVersion<VectorType>(values, iversion);
    }
#endif

    /*************************************************************************/
    /* Column proxy access.                                                  */
    /* It allows accessing to Columns dedicated to a specific Role and Rank. */
    /* with the following syntax:                                            */
    /*   data.X(ir)[is](iv)                                                  */
    /*   - data is the name of the DbData object                             */
    /*   - X is the Role of the Column (X, Z, W, F)                          */
    /*   - ir is the Rank of the Column (optional, 0-based, default = 0)     */
    /*   - is is the index of the sample (0-based)                           */
    /*   - iv is the index of the version (optional, 0-based, default = 0)   */
    /*************************************************************************/
#ifndef SWIG
    ColProxy X(Id rank = 0);
    ColProxy Z(Id rank = 0);
    ColProxy W(Id rank = 0);
    ColProxy F(Id rank = 0);

    ColProxy col(Id icol);
    ColProxy col(const String& name);
#endif

    /***********************************************************************/
    /* Other public methods                                                */
    /***********************************************************************/

    bool renameColumn(ColID&& colid, const String& newName);

    void deleteColumn(ColID&& colid);

    void deleteAllColumns();

    bool hasColumn(ColID&& colid) const;
    String getName(ColID&& colid) const;
    VectorString getNames() const;

    Id getICol(ColID&& colid) const;

    RoleID getRoleID(ColID&& colid) const;

    const ERole& getRole(ColID&& colid) const;

    ColID getColID(const ColID& colid) const;

    void removeRole(ColID&& colid);

    void removeAllRoles();

    Id getNVersions(ColID&& colid) const;

    Id getNRoles(ColID&& colid) const;

    std::vector<ColID> getColIDs(const String& name) const;
    std::vector<ColID> getColIDs(const VectorString& name) const;
    std::vector<ColID> getColIDs(const ERole& role) const;

    Id getNCols() const { return static_cast<Id>(_cols.size()); }

    Id getNSamples() const;

    void setName(ColID&& colid, const String& newName);
    void setRoleID(ColID&& colid, const RoleID& roleID);

    void printContents(const String& title = "") const;

    void clearRole(const ERole& role);

    void addSamples(Id nadd, const double valinit);
    void deleteSample(Id idel);

    String _summaryRoles(void) const;

  private:
    /***********************************************************************/
    /* Other private methods                                               */
    /***********************************************************************/

    std::optional<std::reference_wrapper<DbCol>> _identifyColumn(ColID&& colid);
    std::optional<std::reference_wrapper<const DbCol>>
      _identifyColumn(ColID&& colid) const;

    std::optional<Id>
      _getColumnIndex(const ColID& colid, bool verbose = true) const;

    template<class VectorType>
    VectorType _createEmptyVector(
      Id n,
      std::optional<typename VectorType::value_type> value)
    // const Dictionary* dict = nullptr)
    {
      // if constexpr (std::is_same_v<VectorType, VectorCategory>)
      // {
      //   if (dict == nullptr)
      //     throw std::invalid_argument("Dictionary is required");

      //   VectorCategory vec(n, *dict);

      //   if (value)
      //   {
      //     for (Id i = 0; i < n; i++) vec[i] = *value;
      //   }
      //   return vec;
      // }
      // else
      {
        const auto actual =
          value.value_or(getNA<typename VectorType::value_type>());

        return VectorType(n, actual);
      }
    }

    void _updateName(String& name) const;

    void _updateRoleIDAddition(Id icol0, RoleID& roleID);

    void _updateRoleIDDeletion(RoleID& roleID);

    template<typename VectorType>
    static bool _checkForbidNA(const VectorType& tab)
    {
      using ValueType = typename VectorType::value_type;

      for (const auto& val: tab)
      {
        if (isNA<ValueType>(val))
        {
          messerr("Column forbids NA values, but the input tab contains some.");
          return false;
        }
      }
      return true;
    }

    // static bool _checkForbidNA(const VectorCategory& tab);

    static void _checkVersion(Id& nversion);
    static void _unknownName(const String& name);
    static void _unknownRoleID(const RoleID& roleID);

  private:
    /***********************************************************************/
    /* DbData members                                                      */
    /***********************************************************************/
    std::vector<DbCol> _cols;
    std::vector<RoleID> _roleIDs;
  };

  /***************************************************************************/
  /*                                                                         */
  /*                     Internal proxy implementation                       */
  /*                                                                         */
  /***************************************************************************/
#ifndef SWIG
  class DbData::ColProxy
  {
  public:
    ColProxy(DbData& db, const ColID& colid)
      : _db(db)
      , _colid(colid)
    {
    }

    ValueProxy operator[](Id isample);

  private:
    DbData& _db;
    ColID _colid;
  };

  class DbData::ValueProxy
  {
  public:
    ValueProxy(DbData& db, const ColID& colid, Id isample)
      : _db(db)
      , _colid(colid)
      , _isample(isample)
    {
    }

    /**
     * @brief Select the version of the value
     *
     * Syntax:
     *   data.X()[isample](iversion) = value;
     */
    ValueProxy operator()(Id version)
    {
      auto copy = _colid;
      copy.setVersion(version);
      return ValueProxy(_db, copy, _isample);
    }

    template<typename T>
    ValueProxy& operator=(const T& value)
    {
      _db.setValue<T>(std::move(_colid), _isample, value);
      return *this;
    }

    template<typename T>
    operator T() const
    {
      auto val = _db.getValue<T>(ColID(_colid), _isample);

      if (val) return *val;

      return getNA<T>();
    }

  private:
    DbData& _db;
    ColID _colid;
    Id _isample;
  };

  inline DbData::ValueProxy DbData::ColProxy::operator[](Id isample)
  {
    return ValueProxy(_db, _colid, isample);
  }
#endif
} // namespace gstlrn
