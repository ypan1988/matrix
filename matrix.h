//  matrix.h: this single header matrix lib is a wrapper of std::valarray
//
//  Copyright (C) 2023 Yi Pan <ypan1988@gmail.com>
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.

#ifndef __MATRIX_H
#define __MATRIX_H

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <string>
#include <valarray>

namespace matrix_lib {

//-----------------------------------------------------------------------------

struct _Matrix_error {
  std::string _M_name;
  _Matrix_error(const char* __q) : _M_name(__q) {}
  _Matrix_error(std::string __n) : _M_name(__n) {}
};

//-----------------------------------------------------------------------------

inline void error(const char* __p) { throw _Matrix_error(__p); }

//-----------------------------------------------------------------------------

typedef std::size_t uword;

//-----------------------------------------------------------------------------

// The general Matrix template is simply a prop for its specializations:
template <class _Tp = double, uword _Size = 1>
struct Matrix {
  // multidimensional matrix class
  // ( ) does multidimensional subscripting with range check
  // [ ] does multidimensional subscripting without range check
 private:
  Matrix();  // this should never be compiled
             //	template<class A> Matrix(A);
};

//-----------------------------------------------------------------------------

// Matrix_base represents the common part of the Matrix classes:
template <class _Tp>
struct _Matrix_base {
  typedef _Tp value_type;
  std::valarray<_Tp> _M_elem;

  _Matrix_base() : _M_elem() {}
  _Matrix_base(uword __n) : _M_elem(__n) {}
  _Matrix_base(const value_type& __x, uword __n) : _M_elem(__x, __n) {}
  _Matrix_base(const value_type* __x, uword __n) : _M_elem(__x, __n) {}
  _Matrix_base(const std::valarray<_Tp>& __x) : _M_elem(__x) {}
  _Matrix_base(const std::slice_array<_Tp>& __x) : _M_elem(__x) {}
  _Matrix_base(const std::gslice_array<_Tp>& __x) : _M_elem(__x) {}
  _Matrix_base(const std::mask_array<_Tp>& __x) : _M_elem(__x) {}
  _Matrix_base(const std::indirect_array<_Tp>& __x) : _M_elem(__x) {}
  ~_Matrix_base() {}

  // if necessay, we can get to the raw matrix:
  value_type* data() { return &(_M_elem[0]); }
  const value_type* data() const { return &(_M_elem[0]); }
  void set_elem(const std::valarray<_Tp>& __x) { _M_elem = __x; }
  const std::valarray<_Tp>& get_elem() const { return _M_elem; }
  uword n_elem() const { return _M_elem.size(); }

 public:  // Other member functions.
          // The results are undefined for zero-length arrays
  value_type sum() const { return _M_elem.sum(); }
  value_type min() const { return _M_elem.min(); }
  value_type max() const { return _M_elem.max(); }
};

//-----------------------------------------------------------------------------

template <class _Tp>
struct Matrix<_Tp, 1> : public _Matrix_base<_Tp> {
  uword _M_d1;

 public:
  typedef _Tp value_type;

  Matrix() : _Matrix_base<_Tp>(), _M_d1(0) {}
  Matrix(uword __n1) : _Matrix_base<_Tp>(__n1), _M_d1(__n1) {}
  Matrix(const value_type& __x, uword __n1)
      : _Matrix_base<_Tp>(__x, __n1), _M_d1(__n1) {}
  Matrix(const value_type* __x, uword __n1)
      : _Matrix_base<_Tp>(__x, __n1), _M_d1(__n1) {}
  Matrix(const std::valarray<_Tp>& __x, uword __n1)
      : _Matrix_base<_Tp>(__x), _M_d1(__n1) {
    if (__x.size() != __n1) error("1D Cstor error: dimension mismatch");
  }
  Matrix(const std::slice_array<_Tp>& __x, uword __n1)
      : _Matrix_base<_Tp>(__x), _M_d1(__n1) {}

  uword n_rows() const { return _M_d1; }

  void range_check(uword __n1) const {
    if (_M_d1 <= __n1) error("1D range error: dimension 1");
  }

  // subscripting:
  value_type& operator()(uword __n1) {
    range_check(__n1);
    return this->_M_elem[__n1];
  }
  const value_type& operator()(uword __n1) const {
    range_check(__n1);
    return this->_M_elem[__n1];
  }

 public:
  // clang-format off
  Matrix& operator+=(const value_type& __x) { this->_M_elem += __x; return *this; }
  Matrix& operator-=(const value_type& __x) { this->_M_elem -= __x; return *this; }
  Matrix& operator*=(const value_type& __x) { this->_M_elem *= __x; return *this; }
  Matrix& operator/=(const value_type& __x) { this->_M_elem /= __x; return *this; }
  Matrix& operator%=(const value_type& __x) { this->_M_elem %= __x; return *this; }
  Matrix& operator&=(const value_type& __x) { this->_M_elem &= __x; return *this; }
  Matrix& operator|=(const value_type& __x) { this->_M_elem |= __x; return *this; }
  Matrix& operator^=(const value_type& __x) { this->_M_elem ^= __x; return *this; }
  Matrix& operator<<=(const value_type& __x) { this->_M_elem <<= __x; return *this; }
  Matrix& operator>>=(const value_type& __x) { this->_M_elem >>= __x; return *this; }
  // clang-format on

 public:
  // clang-format off
  Matrix& operator+=(const Matrix& __x) { this->_M_elem += __x._M_elem; return *this; }
  Matrix& operator-=(const Matrix& __x) { this->_M_elem -= __x._M_elem; return *this; }
  Matrix& operator*=(const Matrix& __x) { this->_M_elem *= __x._M_elem; return *this; }
  Matrix& operator/=(const Matrix& __x) { this->_M_elem /= __x._M_elem; return *this; }
  Matrix& operator%=(const Matrix& __x) { this->_M_elem %= __x._M_elem; return *this; }
  Matrix& operator&=(const Matrix& __x) { this->_M_elem &= __x._M_elem; return *this; }
  Matrix& operator|=(const Matrix& __x) { this->_M_elem |= __x._M_elem; return *this; }
  Matrix& operator^=(const Matrix& __x) { this->_M_elem ^= __x._M_elem; return *this; }
  Matrix& operator<<=(const Matrix& __x) { this->_M_elem <<= __x._M_elem; return *this; }
  Matrix& operator>>=(const Matrix& __x) { this->_M_elem >>= __x._M_elem; return *this; }
  // clang-format on
};

//-----------------------------------------------------------------------------

template <class _Tp>
struct Matrix<_Tp, 2> : public _Matrix_base<_Tp> {
  uword _M_d1, _M_d2;

 public:
  typedef _Tp value_type;

  Matrix() : _Matrix_base<_Tp>(), _M_d1(0), _M_d2(0) {}
  Matrix(uword __n1, uword __n2)
      : _Matrix_base<_Tp>(__n1 * __n2), _M_d1(__n1), _M_d2(__n2) {}
  Matrix(const value_type& __x, uword __n1, uword __n2)
      : _Matrix_base<_Tp>(__x, __n1 * __n2), _M_d1(__n1), _M_d2(__n2) {}
  Matrix(const value_type* __x, uword __n1, uword __n2)
      : _Matrix_base<_Tp>(__x, __n1 * __n2), _M_d1(__n1), _M_d2(__n2) {}
  Matrix(const std::valarray<_Tp>& __x, uword __n1, uword __n2)
      : _Matrix_base<_Tp>(__x), _M_d1(__n1), _M_d2(__n2) {
    if (__x.size() != __n1 * __n2) error("2D Cstor error: dimension mismatch");
  }
  Matrix(const std::slice_array<_Tp>& __x, uword __n1, uword __n2)
      : _Matrix_base<_Tp>(__x), _M_d1(__n1), _M_d2(__n2) {}
  uword n_rows() const { return _M_d1; }
  uword n_cols() const { return _M_d2; }

  void range_check(uword __n1, uword __n2) const {
    if (_M_d1 <= __n1) error("2D range error: dimension 1");
    if (_M_d2 <= __n2) error("2D range error: dimension 2");
  }

  // subscripting:
  value_type& operator()(uword __n1, uword __n2) {
    range_check(__n1, __n2);
    return this->_M_elem[__n1 + __n2 * _M_d1];
  }
  const value_type& operator()(uword __n1, uword __n2) const {
    range_check(__n1, __n2);
    return this->_M_elem[__n1 + __n2 * _M_d1];
  }

  std::slice_array<_Tp> diag() {
    return this->_M_elem[std::slice(0, std::min(_M_d1, _M_d2), _M_d1 + 1)];
  }
  std::slice_array<_Tp> row(uword __r) {
    range_check(__r, 0);
    return this->_M_elem[std::slice(__r, _M_d2, _M_d1)];
  }
  std::slice_array<_Tp> col(uword __c) {
    range_check(0, __c);
    return this->_M_elem[std::slice(__c * _M_d1, _M_d1, 1)];
  }

 public:
  // clang-format off
  Matrix& operator+=(const value_type& __x) { this->_M_elem += __x; return *this; }
  Matrix& operator-=(const value_type& __x) { this->_M_elem -= __x; return *this; }
  Matrix& operator*=(const value_type& __x) { this->_M_elem *= __x; return *this; }
  Matrix& operator/=(const value_type& __x) { this->_M_elem /= __x; return *this; }
  Matrix& operator%=(const value_type& __x) { this->_M_elem %= __x; return *this; }
  Matrix& operator&=(const value_type& __x) { this->_M_elem &= __x; return *this; }
  Matrix& operator|=(const value_type& __x) { this->_M_elem |= __x; return *this; }
  Matrix& operator^=(const value_type& __x) { this->_M_elem ^= __x; return *this; }
  Matrix& operator<<=(const value_type& __x) { this->_M_elem <<= __x; return *this; }
  Matrix& operator>>=(const value_type& __x) { this->_M_elem >>= __x; return *this; }
  // clang-format on

 public:
  // clang-format off
  Matrix& operator+=(const Matrix& __x) { this->_M_elem += __x._M_elem; return *this; }
  Matrix& operator-=(const Matrix& __x) { this->_M_elem -= __x._M_elem; return *this; }
  Matrix& operator*=(const Matrix& __x) { this->_M_elem *= __x._M_elem; return *this; }
  Matrix& operator/=(const Matrix& __x) { this->_M_elem /= __x._M_elem; return *this; }
  Matrix& operator%=(const Matrix& __x) { this->_M_elem %= __x._M_elem; return *this; }
  Matrix& operator&=(const Matrix& __x) { this->_M_elem &= __x._M_elem; return *this; }
  Matrix& operator|=(const Matrix& __x) { this->_M_elem |= __x._M_elem; return *this; }
  Matrix& operator^=(const Matrix& __x) { this->_M_elem ^= __x._M_elem; return *this; }
  Matrix& operator<<=(const Matrix& __x) { this->_M_elem <<= __x._M_elem; return *this; }
  Matrix& operator>>=(const Matrix& __x) { this->_M_elem >>= __x._M_elem; return *this; }
  // clang-format on
};

//-----------------------------------------------------------------------------

template <class _Tp>
struct Matrix<_Tp, 3> : public _Matrix_base<_Tp> {
  uword _M_d1, _M_d2, _M_d3;

 private:
  uword _M_d1xd2;

 public:
  typedef _Tp value_type;

  Matrix() : _Matrix_base<_Tp>(), _M_d1(0), _M_d2(0), _M_d3(0) {}
  Matrix(uword __n1, uword __n2, uword __n3)
      : _Matrix_base<_Tp>(__n1 * __n2 * __n3),
        _M_d1(__n1),
        _M_d2(__n2),
        _M_d3(__n3),
        _M_d1xd2(__n1 * __n2) {}
  Matrix(const value_type& __x, uword __n1, uword __n2, uword __n3)
      : _Matrix_base<_Tp>(__x, __n1 * __n2 * __n3),
        _M_d1(__n1),
        _M_d2(__n2),
        _M_d3(__n3),
        _M_d1xd2(__n1 * __n2) {}
  Matrix(const value_type* __x, uword __n1, uword __n2, uword __n3)
      : _Matrix_base<_Tp>(__x, __n1 * __n2 * __n3),
        _M_d1(__n1),
        _M_d2(__n2),
        _M_d3(__n3),
        _M_d1xd2(__n1 * __n2) {}
  Matrix(const std::valarray<_Tp>& __x, uword __n1, uword __n2, uword __n3)
      : _Matrix_base<_Tp>(__x),
        _M_d1(__n1),
        _M_d2(__n2),
        _M_d3(__n3),
        _M_d1xd2(__n1 * __n2) {
    if (__x.size() != __n1 * __n2 * __n3)
      error("3D Cstor error: dimension mismatch");
  }
  Matrix(const std::slice_array<_Tp>& __x, uword __n1, uword __n2, uword __n3)
      : _Matrix_base<_Tp>(__x),
        _M_d1(__n1),
        _M_d2(__n2),
        _M_d3(__n3),
        _M_d1xd2(__n1 * __n2) {}

  uword n_rows() const { return _M_d1; }
  uword n_cols() const { return _M_d2; }
  uword n_slices() const { return _M_d3; }

  void range_check(uword __n1, uword __n2, uword __n3) const {
    if (_M_d1 <= __n1) error("3D range error: dimension 1");
    if (_M_d2 <= __n2) error("3D range error: dimension 2");
    if (_M_d3 <= __n3) error("3D range error: dimension 3");
  }

  // subscripting:
  value_type& operator()(uword __n1, uword __n2, uword __n3) {
    range_check(__n1, __n2, __n3);
    return this->_M_elem[__n1 + __n2 * _M_d1 + __n3 * _M_d1xd2];
  }
  const value_type& operator()(uword __n1, uword __n2, uword __n3) const {
    range_check(__n1, __n2, __n3);
    return this->_M_elem[__n1 + __n2 * _M_d1 + __n3 * _M_d1xd2];
  }

 public:
  // clang-format off
  Matrix& operator+=(const value_type& __x) { this->_M_elem += __x; return *this; }
  Matrix& operator-=(const value_type& __x) { this->_M_elem -= __x; return *this; }
  Matrix& operator*=(const value_type& __x) { this->_M_elem *= __x; return *this; }
  Matrix& operator/=(const value_type& __x) { this->_M_elem /= __x; return *this; }
  Matrix& operator%=(const value_type& __x) { this->_M_elem %= __x; return *this; }
  Matrix& operator&=(const value_type& __x) { this->_M_elem &= __x; return *this; }
  Matrix& operator|=(const value_type& __x) { this->_M_elem |= __x; return *this; }
  Matrix& operator^=(const value_type& __x) { this->_M_elem ^= __x; return *this; }
  Matrix& operator<<=(const value_type& __x) { this->_M_elem <<= __x; return *this; }
  Matrix& operator>>=(const value_type& __x) { this->_M_elem >>= __x; return *this; }
  // clang-format on

 public:
  // clang-format off
  Matrix& operator+=(const Matrix& __x) { this->_M_elem += __x._M_elem; return *this; }
  Matrix& operator-=(const Matrix& __x) { this->_M_elem -= __x._M_elem; return *this; }
  Matrix& operator*=(const Matrix& __x) { this->_M_elem *= __x._M_elem; return *this; }
  Matrix& operator/=(const Matrix& __x) { this->_M_elem /= __x._M_elem; return *this; }
  Matrix& operator%=(const Matrix& __x) { this->_M_elem %= __x._M_elem; return *this; }
  Matrix& operator&=(const Matrix& __x) { this->_M_elem &= __x._M_elem; return *this; }
  Matrix& operator|=(const Matrix& __x) { this->_M_elem |= __x._M_elem; return *this; }
  Matrix& operator^=(const Matrix& __x) { this->_M_elem ^= __x._M_elem; return *this; }
  Matrix& operator<<=(const Matrix& __x) { this->_M_elem <<= __x._M_elem; return *this; }
  Matrix& operator>>=(const Matrix& __x) { this->_M_elem >>= __x._M_elem; return *this; }
  // clang-format on
};

//-----------------------------------------------------------------------------

}  // namespace matrix_lib

#endif