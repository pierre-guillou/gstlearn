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

#include "gstlearn_export.hpp"

#include "DataBase/Dictionary.hpp"

#include <optional>
#include <vector>

namespace gstlrn
{
  /**
   * @brief Vector of categorical values associated with a dictionary.
   *
   * A VectorCategory stores one optional category for each sample.
   * Categories are represented by Dictionary::Category and are associated
   * with a Dictionary defining the available category identifiers and
   * their labels.
   *
   * A sample may have no category, in which case its value is represented
   * by an empty optional.
   */
  class GSTLEARN_EXPORT VectorCategory
  {
  public:
    using Category = Dictionary::Category;
    using value_type = Category;

    /**
     * @brief Constructs a categorical vector.
     *
     * @param nsample Number of samples.
     * @param dict Dictionary associated with the categorical values.
     */
    VectorCategory(
      const size_t nsample = 0,
      const Dictionary& dict = Dictionary())
      : _data(nsample)
      , _dict{dict}
    {
    }

    /**
     * @brief Resizes the categorical vector.
     *
     * @param count New number of samples.
     */
    void resize(const size_t count) { this->_data.resize(count); }

    /**
     * @brief Returns the number of samples.
     *
     * @return Number of samples stored in the vector.
     */
    size_t size() const { return this->_data.size(); }

    /**
     * @brief Returns the category associated with a sample.
     *
     * @param isample Sample index.
     *
     * @return The category associated with the sample, or an empty optional
     *         if the sample has no category or if the index is out of range.
     */
    std::optional<Category> getCategory(const size_t isample) const
    {
      if (isample >= this->_data.size()) return {};
      return this->_data[isample];
    }

    /**
     * @brief Sets the category associated with a sample.
     *
     * The operation fails if the sample index is out of range or if the
     * specified category is not defined in the associated dictionary.
     *
     * @param isample Sample index.
     * @param cat Category to assign to the sample.
     *
     * @return @c true if the category was successfully assigned,
     *         @c false otherwise.
     */
    bool setCategory(const size_t isample, const Category& cat);

#ifndef SWIG
    /**
     * @brief Provides access to a sample category.
     *
     * @param isample Sample index.
     *
     * @return Reference to the optional category associated with the sample.
     *
     * @warning No bounds checking is performed.
     */
    std::optional<Category>& operator[](const size_t isample)
    {
      return this->_data[isample];
    }

    /**
     * @brief Provides read-only access to a sample category.
     *
     * @param isample Sample index.
     *
     * @return Constant reference to the optional category associated with
     *         the sample.
     *
     * @warning No bounds checking is performed.
     */
    const std::optional<Category>& operator[](const size_t isample) const
    {
      return this->_data[isample];
    }
#endif

  private:
    std::vector<std::optional<Category>> _data;
    Dictionary _dict;
  };

} // namespace gstlrn
