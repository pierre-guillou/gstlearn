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

#include "geoslib_define.h"

#include "DataBase/VectorCategory.hpp"
#include <optional>

namespace gstlrn
{
  /**
   * @brief Two-dimensional array organized by samples and versions.
   *
   * The array stores @c outer versions of @c inner samples. The number of
   * samples is always the same for all versions.
   *
   * The data are stored contiguously in the underlying vector, with all
   * samples of a given version stored consecutively:
   *
   * @code
   * index(o,i) = o * inner + i
   * @endcode
   *
   * where @c o is the version index and @c i is the sample index.
   *
   * All indices passed to the access and modification methods are expected
   * to be valid. Their validity is checked by the calling code rather than
   * by each individual operation, in order to avoid unnecessary overhead
   * in frequently used operations.
   *
   * @tparam VectorType Type of the vector used to store the data.
   *
   * The class can also be instantiated with VectorCategory, in which case
   * category-specific access is used for getting and setting values.
   *
   * @note This class is temporarily not exported for SWIG.
   */
  template<typename VectorType>
  class Array2D
  {
  public:
    using vector_type = VectorType;
    using value_type = typename VectorType::value_type;

    /**
     * @brief Default constructor.
     *
     * Creates an empty array with zero samples, zero versions and an
     * empty underlying buffer.
     */
    Array2D() = default;

    /**
     * @brief Constructs an array from an existing vector.
     *
     * The vector is moved into the array. Its size must be compatible
     * with the specified number of versions.
     *
     * @param vector Vector containing the array values.
     * @param nversion Number of versions.
     *
     * @note The number of samples is deduced from the vector size:
     *       `inner = vector.size() / nversion`.
     */
    Array2D(VectorType&& vector, Id nversion = 1)
      : _inner{static_cast<Id>(vector.size()) / nversion}
      , _outer{nversion}
      , _buf{std::forward<VectorType>(vector)}
    {
    }

    /**
     * @brief Constructs an array with a specified size and initial value.
     *
     * @param inner Number of samples in each version.
     * @param outer Number of versions.
     * @param val Initial value assigned to all elements.
     */
    Array2D(const Id inner, const Id outer = 1, const value_type val = {})
      : _inner{inner}
      , _outer{outer}
      , _buf(outer * inner, val)
    {
    }

    /**
     * @brief Resizes the array.
     *
     * The underlying buffer is resized to contain @p inner samples for
     * each of the @p outer versions. Newly created elements are initialized
     * with @p val.
     *
     * @param inner New number of samples in each version.
     * @param outer New number of versions.
     * @param val Value used to initialize newly created elements.
     */
    void resize(const Id inner, const Id outer = 1, const value_type val = {})
    {
      this->_inner = inner;
      this->_outer = outer;
      this->_buf.resize(this->_outer * this->_inner, val);
    }

    /**
     * @brief Adds samples to all versions.
     *
     * @param nnewsamp Number of samples to add to each version.
     * @param val Value assigned to the newly added samples.
     *
     * Existing values are preserved.
     */
    void addSamples(const Id nnewsamp, const value_type val)
    {
      VectorType newbuf(this->_outer * (this->_inner + nnewsamp), val);
      for (Id o = 0; o < this->_outer; ++o)
      {
        for (Id i = 0; i < this->_inner; ++i)
        {
          newbuf[(o * this->_inner) + i] = this->_buf[(o * this->_inner) + i];
        }
      }
      this->_inner += nnewsamp;
      std::swap(this->_buf, newbuf);
    }

    /**
     * @brief Deletes one sample from all versions.
     *
     * The sample identified by @p isamp is removed from every version.
     * Existing values of the other samples are preserved.
     *
     * @param isamp Index of the sample to delete.
     */
    void deleteSample(const Id isamp)
    {
      VectorType newbuf(this->_buf);
      newbuf.resize(this->_outer * (this->_inner - 1));

      for (Id o = 0; o < this->_outer; ++o)
      {
        Id a{};
        for (Id i = 0; i < this->_inner; ++i)
        {
          if (i == isamp)
          {
            continue;
          }
          newbuf[(o * (this->_inner - 1)) + a] =
            this->_buf[(o * this->_inner) + i];
          a++;
        }
      }

      this->_inner -= 1;
      std::swap(this->_buf, newbuf);
    }

    /**
     * @brief Returns the value of one element.
     *
     * @param o Version index.
     * @param i Sample index.
     *
     * @return The requested value. An empty optional may be returned when
     *         the underlying vector type represents a missing value.
     */
    std::optional<value_type> getValue(const size_t o, const size_t i) const
    {
      if constexpr (std::is_same_v<VectorType, VectorCategory>)
      {
        return this->_buf.getCategory((o * this->_inner) + i);
      }
      else
      {
        return this->_buf[(o * this->_inner) + i];
      }
    }

    /**
     * @brief Sets the value of one element.
     *
     * @param o Version index.
     * @param i Sample index.
     * @param val New value.
     *
     * @return @c true if the value was successfully set, @c false otherwise.
     */
    bool setValue(const size_t o, const size_t i, const value_type& val)
    {
      if constexpr (std::is_same_v<VectorType, VectorCategory>)
      {
        return this->_buf.setCategory((o * this->_inner) + i, val);
      }
      else
      {
        this->_buf[(o * this->_inner) + i] = val;
        return true;
      }
    }

    /**
     * @brief Returns the number of samples in each version.
     */
    Id inner() const { return this->_inner; }

    /**
     * @brief Returns the number of versions.
     */
    Id outer() const { return this->_outer; }

    /**
     * @brief Returns a pointer to the underlying data.
     */
    value_type* data() { return this->_buf.data(); }

    /**
     * @brief Returns a read-only pointer to the underlying data.
     */
    const value_type* data() const { return this->_buf.data(); }

    /**
     * @brief Returns a reference to the underlying vector.
     *
     * The returned vector can be modified directly.
     */
    VectorType& getBuffer() { return this->_buf; }

    /**
     * @brief Returns a read-only reference to the underlying vector.
     */
    const VectorType& getBuffer() const { return this->_buf; }

  private:
    Id _inner{};
    Id _outer{};
    VectorType _buf;
  };
} // namespace gstlrn
