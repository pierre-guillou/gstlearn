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
#include "DataBase/RoleID.hpp"

namespace gstlrn
{
  RoleID::RoleID(const ERole& role, Id index)
    : _role(role)
    , _index(index)
  {
  }

  bool RoleID::isUnique() const
  {
    const auto& attr = ERoleAttr.at(_role.getKey());
    return attr.isUnique;
  }

  String RoleID::getName() const
  {
    if (!isDefined()) return STRING_NA;

    String name(_role.getKey());

    std::transform(
      name.begin(), name.end(), name.begin(),
      [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

    if (!isUnique()) name += std::to_string(_index + 1);

    return name;
  }

  RoleID roleIDIdentify(const String& name)
  {
    if (name.empty()) return RoleID();

    String lname = name;

    std::transform(
      lname.begin(), lname.end(), lname.begin(),
      [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

    for (Id irank = 0; irank < static_cast<Id>(ERole::getSize()); irank++)
    {
      ERole role = ERole::fromValue(irank);

      String roleName(role.getKey());

      std::transform(
        roleName.begin(), roleName.end(), roleName.begin(),
        [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

      // Test whether this role is unique
      RoleID roleID(role, 0);

      if (roleID.isUnique())
      {
        if (lname == roleName) return roleID;

        continue;
      }

      // Non-unique role: a positive multiplicity is mandatory
      if (lname.size() <= roleName.size()) continue;

      if (lname.compare(0, roleName.size(), roleName) != 0) continue;

      String suffix = lname.substr(roleName.size());

      if (!std::all_of(
            suffix.begin(), suffix.end(),
            [](unsigned char c) { return std::isdigit(c); }))
        continue;

      Id multiplicity = std::stoll(suffix);

      if (multiplicity <= 0) continue;

      return RoleID(role, multiplicity - 1);
    }
    return RoleID();
  }

  /**
   * @brief Check if the current RoleID matches the one provided as argument
   * A match is complete is they have the same Role and the same index
   *
   * @param roleID RoleID to compare with
   * @param checkIndex When True, check the equality of the Index
   * @return true
   * @return false
   */
  bool RoleID::match(const RoleID& roleID, bool checkIndex) const
  {
    if (_role.isDifferent(roleID.getRole())) return false;
    if (!checkIndex) return true;
    return _index == roleID.getIndex();
  }

  /**
   * Given a locator string, create the corresponding RoleID
   * @param string     Locator string
   * @return Error code
   */
  std::optional<RoleID> RoleID::createFromName(String string)
  {
    // Mise en minuscules
    toLower(string);

    // Extracting the Role and the Index
    String roleName;
    String indexString;

    for (char c: string)
    {
      if (std::isdigit(static_cast<unsigned char>(c)))
        indexString += c;
      else
        roleName += c;
    }

    // Recherche du rôle correspondant
    ERole role = ERole::UNDEFINED;

    auto it = ERole::getIterator();
    while (it.hasNext())
    {
      auto current = *it;

      if (current != ERole::UNDEFINED)
      {
        String key = std::string(current.getKey());
        toLower(key);

        if (key == roleName)
        {
          role = current;
          break;
        }
      }

      it.toNext();
    }

    if (role == ERole::UNDEFINED) return RoleID(role, 0);

    // Conversion du rang utilisateur (1-based) en index interne (0-based)
    Id index = 0;
    if (!indexString.empty())
      index = std::max(atoi(indexString.c_str()) - 1, 0);

    // Vérification des rôles uniques
    auto itAttr = ERoleAttr.find(role.getKey());
    if (itAttr != ERoleAttr.end())
    {
      if (itAttr->second.isUnique && index > 0) return std::nullopt;
    }

    return RoleID(role, index);
  }

  void RoleID::removeRole()
  {
    _role = ERole::UNDEFINED;
    _index = 0;
  }

} // namespace gstlrn
