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

  String RoleID::getDescr() const
  {
    if (_role == ERole::UNDEFINED) return STRING_NA;

    String name(std::string(_role.getKey()));

    std::transform(
      name.begin(), name.end(), name.begin(),
      [](unsigned char c) { return std::tolower(c); });

    if (!isUnique()) name += std::to_string(_index + 1);

    return name;
  }

  /**
   * @brief Check if the current RoleID matches the one provided as argument
   * A match is complete is they have the same Role and the same index
   *
   * @param roleID RoleID to compare with
   * @return true
   * @return false
   */
  bool RoleID::match(const RoleID& roleID) const
  {
    return _role == roleID.getRole() && _index == roleID.getIndex();
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

    // Extraction du nom du rôle et du rang éventuel
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

    if (role == ERole::UNDEFINED) return std::nullopt;

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
