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
  template<typename VectorType>
  // This class is temporarily not exported (for SWIG)
  class Array2D
  {
  public:
    using vector_type = VectorType;
    using value_type = typename VectorType::value_type;

    Array2D() = default;

    Array2D(VectorType&& vector, Id nversion = 1)
      : _inner{static_cast<Id>(vector.size()) / nversion}
      , _outer{nversion}
      , _buf{std::forward<VectorType>(vector)}
    {
    }

    Array2D(const Id inner, const Id outer = 1, const value_type val = {})
      : _inner{inner}
      , _outer{outer}
      , _buf(outer * inner, val)
    {
    }

    void resize(const Id inner, const Id outer = 1, const value_type val = {})
    {
      this->_inner = inner;
      this->_outer = outer;
      this->_buf.resize(this->_outer * this->_inner, val);
    }

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

    Id inner() const { return this->_inner; }

    Id outer() const { return this->_outer; }

    value_type* data() { return this->_buf.data(); }

    const value_type* data() const { return this->_buf.data(); }

    VectorType& getBuffer() { return this->_buf; }

    const VectorType& getBuffer() const { return this->_buf; }

  private:
    Id _inner{};
    Id _outer{};
    VectorType _buf;
  };
} // namespace gstlrn
