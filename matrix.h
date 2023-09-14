//  matrix.h: this single header matrix lib is a wrapper of std::valarray
//
//  Copyright (C) 2022-2023 Yi Pan <ypan1988@gmail.com>
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

#if __cplusplus >= 201103L || (defined(_MSC_VER) && _MSC_VER >= 1900)
#define __MATRIX_LIB_USE_CPP11
#endif

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <string>
#include <valarray>

#if defined(__MATRIX_LIB_USE_CPP11)
#include <initializer_list>
#endif

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

template <class _Tp = double>
struct MatrixRef;

//-----------------------------------------------------------------------------

// Matrix_base represents the common part of the Matrix classes:
template <class _Tp>
struct _Matrix_base {
 protected:
  std::valarray<_Tp> _M_elem;

 public:
  typedef _Tp value_type;

  // construct/destroy:
  _Matrix_base() : _M_elem() {}
  _Matrix_base(uword __n) : _M_elem(__n) {}
  _Matrix_base(const value_type& __x, uword __n) : _M_elem(__x, __n) {}
  _Matrix_base(const value_type* __x, uword __n) : _M_elem(__x, __n) {}
  _Matrix_base(const std::valarray<_Tp>& __x) : _M_elem(__x) {}
#if defined(__MATRIX_LIB_USE_CPP11)
  _Matrix_base(const _Matrix_base& __x) = default;
  _Matrix_base(_Matrix_base&& __x) = default;
  _Matrix_base(std::initializer_list<_Tp> __x) : _M_elem(__x) {}
#endif
  _Matrix_base(const std::slice_array<_Tp>& __x) : _M_elem(__x) {}
  _Matrix_base(const std::gslice_array<_Tp>& __x) : _M_elem(__x) {}
  _Matrix_base(const std::mask_array<_Tp>& __x) : _M_elem(__x) {}
  _Matrix_base(const std::indirect_array<_Tp>& __x) : _M_elem(__x) {}
  virtual ~_Matrix_base() {}

  // assignment
#if defined(__MATRIX_LIB_USE_CPP11)
  _Matrix_base& operator=(const _Matrix_base& __x) = default;
  _Matrix_base& operator=(_Matrix_base&& __x) = default;
#endif

  // clang-format off
  value_type  operator[](uword __n) const { return _M_elem[__n]; }
  value_type& operator[](uword __n) { return _M_elem[__n]; }
  std::valarray<_Tp> operator[](std::slice __x) const { return _M_elem[__x]; }
  std::slice_array<_Tp> operator[](std::slice __x) { return _M_elem[__x]; }
  std::valarray<_Tp> operator[](std::gslice __x) const { return _M_elem[__x]; }
  std::gslice_array<_Tp> operator[](std::gslice __x) { return _M_elem[__x]; }
  std::valarray<_Tp> operator[](const std::valarray<bool> &__x) const { return _M_elem[__x]; }
  std::mask_array<_Tp> operator[](const std::valarray<bool> &__x) { return _M_elem[__x]; }
  std::valarray<_Tp> operator[](const std::valarray<uword> &__x) const { return _M_elem[__x]; }
  std::indirect_array<_Tp> operator[](const std::valarray<uword> &__x) { return _M_elem[__x]; }
  // clang-format on

  // if necessay, we can get to the raw matrix:
  value_type* data() { return &(_M_elem[0]); }
  const value_type* data() const { return &(_M_elem[0]); }
  void set_elem(const std::valarray<_Tp>& __x) { _M_elem = __x; }
  std::valarray<_Tp> get_elem() const { return _M_elem; }
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
 public:
  typedef _Tp value_type;

  // clang-format off
  // construct/destroy:
  Matrix() : _Matrix_base<_Tp>() { _M_init(); }
  Matrix(const Matrix& __x) : _Matrix_base<_Tp>(__x._M_elem) { _M_init(); }
#if defined(__MATRIX_LIB_USE_CPP11)
  Matrix(Matrix&& __x) = default;
  Matrix(std::initializer_list<_Tp> __x) : _Matrix_base<_Tp>(__x) { _M_init(); }
#endif
  explicit Matrix(uword __n1) : _Matrix_base<_Tp>(__n1) { _M_init(); }
  Matrix(const value_type& __x, uword __n1) : _Matrix_base<_Tp>(__x, __n1) { _M_init(); }
  Matrix(const value_type* __x, uword __n1) : _Matrix_base<_Tp>(__x, __n1) { _M_init(); }
  Matrix(const std::valarray<_Tp>& __x) : _Matrix_base<_Tp>(__x) { _M_init(); }
  Matrix(const std::slice_array<_Tp>& __x) : _Matrix_base<_Tp>(__x) { _M_init(); }
  Matrix(const std::gslice_array<_Tp>& __x) : _Matrix_base<_Tp>(__x) { _M_init(); }
  Matrix(const std::mask_array<_Tp>& __x) : _Matrix_base<_Tp>(__x) { _M_init(); }
  Matrix(const std::indirect_array<_Tp>& __x) : _Matrix_base<_Tp>(__x) { _M_init(); }
  Matrix(const std::valarray<_Tp>& __x, const uword __dims[1]) : _Matrix_base<_Tp>(__x) { _M_init(__dims); }
  Matrix(const std::valarray<_Tp>& __x, const std::valarray<uword>& __dims) : _Matrix_base<_Tp>(__x) { _M_init(__dims); }
  ~Matrix() {}

  // assignment
  Matrix& operator=(const Matrix& __x) { this->_M_elem = __x._M_elem; _M_init(); return *this; }
#if defined(__MATRIX_LIB_USE_CPP11)
  Matrix& operator=(Matrix&& __x) = default;
#endif
  // clang-format on

  std::valarray<uword> get_dims() const {
    return std::valarray<uword>(_M_dims, 1);
  }
  uword n_rows() const { return _M_dims[0]; }

  void range_check(uword __n1) const {
    if (_M_dims[0] <= __n1) error("1D range error: dimension 1");
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

  std::valarray<_Tp> subvec(std::size_t __i, std::size_t __j) const {
    if (__i > __j || __j >= _M_dims[0]) error("1D subscription error");
    return this->_M_elem[std::slice(__i, __j - __i + 1, 1)];
  }
  std::slice_array<_Tp> subvec(std::size_t __i, std::size_t __j) {
    if (__i > __j || __j >= _M_dims[0]) error("1D subscription error");
    return this->_M_elem[std::slice(__i, __j - __i + 1, 1)];
  }

  // clang-format off
 public:
  Matrix operator+() const { return *this; }
  Matrix operator-() const { return Matrix(-this->_M_elem, _M_dims); }
  Matrix operator~() const { return Matrix(~this->_M_elem, _M_dims); }
  Matrix<bool, 1> operator!() const { return Matrix<bool, 1>(!this->_M_elem, _M_dims); }

 public:
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

 public:
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

 private:
  uword _M_dims[1];

  void _M_init() { _M_dims[0] = this->n_elem(); }
  void _M_init(const uword __dims[1]) {
    if (this->n_elem() != __dims[0])
      error("1D Cstor error: dimension mismatch");
    _M_dims[0] = __dims[0];
  }
  void _M_init(const std::valarray<uword>& __dims) {
    if (this->n_elem() != __dims[0])
      error("1D Cstor error: dimension mismatch");
    _M_dims[0] = __dims[0];
  }
};

//-----------------------------------------------------------------------------

template <class _Tp>
struct Matrix<_Tp, 2> : public _Matrix_base<_Tp> {
 public:
  typedef _Tp value_type;

  // clang-format off
  // construct/destroy:
  Matrix() : _Matrix_base<_Tp>() { _M_init(0, 0); }
  Matrix(const Matrix& __x) : _Matrix_base<_Tp>(__x._M_elem) { _M_init(__x._M_dims); }
#if defined(__MATRIX_LIB_USE_CPP11)
  Matrix(Matrix&& __x) = default;
  Matrix(std::initializer_list<_Tp> __x, uword __n1, uword __n2) : _Matrix_base<_Tp>(__x) { _M_init(__n1, __n2); }
#endif
  Matrix(uword __n1, uword __n2) : _Matrix_base<_Tp>(__n1 * __n2) { _M_init(__n1, __n2); }
  Matrix(const value_type& __x, uword __n1, uword __n2) : _Matrix_base<_Tp>(__x, __n1 * __n2) { _M_init(__n1, __n2); }
  Matrix(const value_type* __x, uword __n1, uword __n2) : _Matrix_base<_Tp>(__x, __n1 * __n2) { _M_init(__n1, __n2); }
  Matrix(const std::valarray<_Tp>& __x, uword __n1, uword __n2) : _Matrix_base<_Tp>(__x) { _M_init(__n1, __n2); }
  Matrix(const std::slice_array<_Tp>& __x, uword __n1, uword __n2) : _Matrix_base<_Tp>(__x) { _M_init(__n1, __n2); }
  Matrix(const std::gslice_array<_Tp>& __x, uword __n1, uword __n2) : _Matrix_base<_Tp>(__x) { _M_init(__n1, __n2); }
  Matrix(const std::mask_array<_Tp>& __x, uword __n1, uword __n2) : _Matrix_base<_Tp>(__x) { _M_init(__n1, __n2); }
  Matrix(const std::indirect_array<_Tp>& __x, uword __n1, uword __n2) : _Matrix_base<_Tp>(__x) { _M_init(__n1, __n2); }
  Matrix(const std::valarray<_Tp>& __x, const uword __dims[2]) : _Matrix_base<_Tp>(__x) { _M_init(__dims); }
  Matrix(const std::valarray<_Tp>& __x, const std::valarray<uword>& __dims) : _Matrix_base<_Tp>(__x) { _M_init(__dims); }
  ~Matrix() {}

  // assignment
  Matrix& operator=(const Matrix& __x) { this->_M_elem = __x._M_elem; _M_init(__x._M_dims); return *this; }
#if defined(__MATRIX_LIB_USE_CPP11)
  Matrix& operator=(Matrix&& __x) = default;
#endif
  // clang-format on

  std::valarray<uword> get_dims() const {
    return std::valarray<uword>(_M_dims, 2);
  }

  uword n_rows() const { return _M_dims[0]; }
  uword n_cols() const { return _M_dims[1]; }

  void range_check(uword __n1, uword __n2) const {
    if (_M_dims[0] <= __n1) error("2D range error: dimension 1");
    if (_M_dims[1] <= __n2) error("2D range error: dimension 2");
  }

  // subscripting:
  value_type& operator()(uword __n1, uword __n2) {
    range_check(__n1, __n2);
    return this->_M_elem[__n1 + __n2 * _M_dims[0]];
  }
  const value_type& operator()(uword __n1, uword __n2) const {
    range_check(__n1, __n2);
    return this->_M_elem[__n1 + __n2 * _M_dims[0]];
  }

  std::valarray<_Tp> diag() const {
    return this->_M_elem[std::slice(0, std::min(_M_dims[0], _M_dims[1]),
                                    _M_dims[0] + 1)];
  }
  std::slice_array<_Tp> diag() {
    return this->_M_elem[std::slice(0, std::min(_M_dims[0], _M_dims[1]),
                                    _M_dims[0] + 1)];
  }
  std::valarray<_Tp> row(uword __r) const {
    range_check(__r, 0);
    return this->_M_elem[std::slice(__r, _M_dims[1], _M_dims[0])];
  }
  std::slice_array<_Tp> row(uword __r) {
    range_check(__r, 0);
    return this->_M_elem[std::slice(__r, _M_dims[1], _M_dims[0])];
  }
  std::valarray<_Tp> col(uword __c) const {
    range_check(0, __c);
    return this->_M_elem[std::slice(__c * _M_dims[0], _M_dims[0], 1)];
  }
  std::slice_array<_Tp> col(uword __c) {
    range_check(0, __c);
    return this->_M_elem[std::slice(__c * _M_dims[0], _M_dims[0], 1)];
  }
  Matrix<_Tp, 2> submat(uword __first_row, uword __first_col, uword __last_row,
                        uword __last_col) const;
  MatrixRef<_Tp> submat(uword __first_row, uword __first_col, uword __last_row,
                        uword __last_col);

  // clang-format off
 public:
  Matrix operator+() const { return *this; }
  Matrix operator-() const { return Matrix(-this->_M_elem, _M_dims); }
  Matrix operator~() const { return Matrix(~this->_M_elem, _M_dims); }
  Matrix<bool, 2> operator!() const { return Matrix<bool, 2>(!this->_M_elem, _M_dims); }

 public:
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

 public:
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

 private:
  uword _M_dims[2];

  void _M_init(uword __n1, uword __n2) {
    if (this->n_elem() != __n1 * __n2)
      error("2D Cstor error: dimension mismatch");
    _M_dims[0] = __n1, _M_dims[1] = __n2;
  }
  void _M_init(const uword __dims[2]) {
    _M_dims[0] = __dims[0], _M_dims[1] = __dims[1];
  }
  void _M_init(const std::valarray<uword>& __dims) {
    _M_dims[0] = __dims[0], _M_dims[1] = __dims[1];
  }
};

//-----------------------------------------------------------------------------

template <class _Tp>
struct Matrix<_Tp, 3> : public _Matrix_base<_Tp> {
 public:
  typedef _Tp value_type;

  // clang-format off
  // construct/destroy:
  Matrix() : _Matrix_base<_Tp>() { _M_init(0, 0, 0); }
  Matrix(const Matrix& __x) : _Matrix_base<_Tp>(__x._M_elem) { _M_init(__x._M_dims); }
#if defined(__MATRIX_LIB_USE_CPP11)
  Matrix(Matrix&& __x) = default;
  Matrix(std::initializer_list<_Tp> __x, uword __n1, uword __n2, uword __n3) : _Matrix_base<_Tp>(__x) { _M_init(__n1, __n2, __n3); }
#endif
  Matrix(uword __n1, uword __n2, uword __n3) : _Matrix_base<_Tp>(__n1 * __n2 * __n3) { _M_init(__n1, __n2, __n3); }
  Matrix(const value_type& __x, uword __n1, uword __n2, uword __n3) : _Matrix_base<_Tp>(__x, __n1 * __n2 * __n3) { _M_init(__n1, __n2, __n3); }
  Matrix(const value_type* __x, uword __n1, uword __n2, uword __n3) : _Matrix_base<_Tp>(__x, __n1 * __n2 * __n3) { _M_init(__n1, __n2, __n3); }
  Matrix(const std::valarray<_Tp>& __x, uword __n1, uword __n2, uword __n3) : _Matrix_base<_Tp>(__x) { _M_init(__n1, __n2, __n3); }
  Matrix(const std::slice_array<_Tp>& __x, uword __n1, uword __n2, uword __n3) : _Matrix_base<_Tp>(__x) { _M_init(__n1, __n2, __n3); }
  Matrix(const std::gslice_array<_Tp>& __x, uword __n1, uword __n2, uword __n3) : _Matrix_base<_Tp>(__x) { _M_init(__n1, __n2, __n3); }
  Matrix(const std::mask_array<_Tp>& __x, uword __n1, uword __n2, uword __n3) : _Matrix_base<_Tp>(__x) { _M_init(__n1, __n2, __n3); }
  Matrix(const std::indirect_array<_Tp>& __x, uword __n1, uword __n2, uword __n3) : _Matrix_base<_Tp>(__x) { _M_init(__n1, __n2, __n3); }
  Matrix(const std::valarray<_Tp>& __x, const uword __dims[3]) : _Matrix_base<_Tp>(__x) { _M_init(__dims); }
  Matrix(const std::valarray<_Tp>& __x, const std::valarray<uword>& __dims) : _Matrix_base<_Tp>(__x) { _M_init(__dims); }
  ~Matrix() {}

  // assignment
  Matrix& operator=(const Matrix& __x) { this->_M_elem = __x._M_elem; _M_init(__x._M_dims); return *this; }
#if defined(__MATRIX_LIB_USE_CPP11)
  Matrix& operator=(Matrix&& __x) = default;
#endif
  // clang-format on

  std::valarray<uword> get_dims() const {
    return std::valarray<uword>(_M_dims, 3);
  }

  uword n_rows() const { return _M_dims[0]; }
  uword n_cols() const { return _M_dims[1]; }
  uword n_slices() const { return _M_dims[2]; }

  void range_check(uword __n1, uword __n2, uword __n3) const {
    if (_M_dims[0] <= __n1) error("3D range error: dimension 1");
    if (_M_dims[1] <= __n2) error("3D range error: dimension 2");
    if (_M_dims[2] <= __n3) error("3D range error: dimension 3");
  }

  // subscripting:
  value_type& operator()(uword __n1, uword __n2, uword __n3) {
    range_check(__n1, __n2, __n3);
    return this->_M_elem[__n1 + __n2 * _M_dims[0] + __n3 * _M_d1xd2];
  }
  const value_type& operator()(uword __n1, uword __n2, uword __n3) const {
    range_check(__n1, __n2, __n3);
    return this->_M_elem[__n1 + __n2 * _M_dims[0] + __n3 * _M_d1xd2];
  }

  // clang-format off
 public:
  Matrix operator+() const { return *this; }
  Matrix operator-() const { return Matrix(-this->_M_elem, _M_dims); }
  Matrix operator~() const { return Matrix(~this->_M_elem, _M_dims); }
  Matrix<bool, 3> operator!() const { return Matrix<bool, 3>(!this->_M_elem, _M_dims); }

 public:
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

 public:
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

 private:
  uword _M_dims[3], _M_d1xd2;

  void _M_init(uword __n1, uword __n2, uword __n3) {
    if (this->n_elem() != __n1 * __n2 * __n3)
      error("3D Cstor error: dimension mismatch");
    _M_dims[0] = __n1, _M_dims[1] = __n2, _M_dims[2] = __n3;
    _M_d1xd2 = __n1 * __n2;
  }
  void _M_init(const uword __dims[3]) {
    _M_dims[0] = __dims[0], _M_dims[1] = __dims[1], _M_dims[2] = __dims[2];
    _M_d1xd2 = __dims[0] * __dims[1];
  }
  void _M_init(const std::valarray<uword>& __dims) {
    _M_dims[0] = __dims[0], _M_dims[1] = __dims[1], _M_dims[2] = __dims[2];
    _M_d1xd2 = __dims[0] * __dims[1];
  }
};

//-----------------------------------------------------------------------------

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator+(const Matrix<_Tp, _Size>& __x,
                                    const Matrix<_Tp, _Size>& __y) {
  Matrix<_Tp, _Size> __tmp(__x);
  return __tmp += __y;
}

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator-(const Matrix<_Tp, _Size>& __x,
                                    const Matrix<_Tp, _Size>& __y) {
  Matrix<_Tp, _Size> __tmp(__x);
  return __tmp -= __y;
}

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator*(const Matrix<_Tp, _Size>& __x,
                                    const Matrix<_Tp, _Size>& __y) {
  Matrix<_Tp, _Size> __tmp(__x);
  return __tmp *= __y;
}

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator/(const Matrix<_Tp, _Size>& __x,
                                    const Matrix<_Tp, _Size>& __y) {
  Matrix<_Tp, _Size> __tmp(__x);
  return __tmp /= __y;
}

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator%(const Matrix<_Tp, _Size>& __x,
                                    const Matrix<_Tp, _Size>& __y) {
  Matrix<_Tp, _Size> __tmp(__x);
  return __tmp %= __y;
}

//-----------------------------------------------------------------------------

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator+(const Matrix<_Tp, _Size>& __x,
                                    const _Tp& __c) {
  Matrix<_Tp, _Size> __tmp(__x);
  return __tmp += __c;
}

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator-(const Matrix<_Tp, _Size>& __x,
                                    const _Tp& __c) {
  Matrix<_Tp, _Size> __tmp(__x);
  return __tmp -= __c;
}

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator*(const Matrix<_Tp, _Size>& __x,
                                    const _Tp& __c) {
  Matrix<_Tp, _Size> __tmp(__x);
  return __tmp *= __c;
}

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator/(const Matrix<_Tp, _Size>& __x,
                                    const _Tp& __c) {
  Matrix<_Tp, _Size> __tmp(__x);
  return __tmp /= __c;
}

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator%(const Matrix<_Tp, _Size>& __x,
                                    const _Tp& __c) {
  Matrix<_Tp, _Size> __tmp(__x);
  return __tmp %= __c;
}

//-----------------------------------------------------------------------------

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator+(const _Tp& __c,
                                    const Matrix<_Tp, _Size>& __x) {
  Matrix<_Tp, _Size> __tmp(__c + __x.get_elem(), __x.get_elem());
  return __tmp;
}

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator-(const _Tp& __c,
                                    const Matrix<_Tp, _Size>& __x) {
  Matrix<_Tp, _Size> __tmp(__c - __x.get_elem(), __x.get_elem());
  return __tmp;
}

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator*(const _Tp& __c,
                                    const Matrix<_Tp, _Size>& __x) {
  Matrix<_Tp, _Size> __tmp(__c * __x.get_elem(), __x.get_elem());
  return __tmp;
}

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator/(const _Tp& __c,
                                    const Matrix<_Tp, _Size>& __x) {
  Matrix<_Tp, _Size> __tmp(__c / __x.get_elem(), __x.get_elem());
  return __tmp;
}

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator%(const _Tp& __c,
                                    const Matrix<_Tp, _Size>& __x) {
  Matrix<_Tp, _Size> __tmp(__c % __x.get_elem(), __x.get_elem());
  return __tmp;
}

//-----------------------------------------------------------------------------

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> abs(const Matrix<_Tp, _Size>& __x) {
  Matrix<_Tp, _Size> __tmp(std::abs(__x.get_elem()), __x.get_dims());
  return __tmp;
}

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> exp(const Matrix<_Tp, _Size>& __x) {
  Matrix<_Tp, _Size> __tmp(std::exp(__x.get_elem()), __x.get_dims());
  return __tmp;
}

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> log(const Matrix<_Tp, _Size>& __x) {
  Matrix<_Tp, _Size> __tmp(std::log(__x.get_elem()), __x.get_dims());
  return __tmp;
}

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> log10(const Matrix<_Tp, _Size>& __x) {
  Matrix<_Tp, _Size> __tmp(std::log10(__x.get_elem()), __x.get_dims());
  return __tmp;
}

//-----------------------------------------------------------------------------

template <class _Tp>
struct MatrixRef {
  std::gslice_array<_Tp> _M_elem;
  std::valarray<uword> _M_dims;
  MatrixRef(std::valarray<_Tp>& __va, uword __start,
            const std::valarray<std::size_t>& __size,
            const std::valarray<std::size_t>& __stride)
      : _M_elem(__va[std::gslice(__start, __size, __stride)]),
        _M_dims(__size.size()) {
    uword n = __size.size();
    for (uword idx = 0; idx < n; ++idx) {
      _M_dims[idx] = __size[n - idx - 1];
    }
  }
  void operator=(const _Tp& __value) { _M_elem = __value; }
};

template <class _Tp>
inline Matrix<_Tp, 2> Matrix<_Tp, 2>::submat(uword __first_row,
                                             uword __first_col,
                                             uword __last_row,
                                             uword __last_col) const {
  const uword __sz[2] = {__last_col - __first_col + 1,
                         __last_row - __first_row + 1};
  const uword __st[2] = {n_rows(), 1};

  const uword __start = n_rows() * __first_col + __first_row;
  const std::valarray<uword> __size(__sz, 2);
  const std::valarray<uword> __stride(__st, 2);

  return Matrix<_Tp, 2>(this->_M_elem[std::gslice(__start, __size, __stride)],
                        __size[1], __size[0]);
}

template <class _Tp>
inline MatrixRef<_Tp> Matrix<_Tp, 2>::submat(uword __first_row,
                                             uword __first_col,
                                             uword __last_row,
                                             uword __last_col) {
  const uword __sz[2] = {__last_col - __first_col + 1,
                         __last_row - __first_row + 1};
  const uword __st[2] = {n_rows(), 1};

  const uword __start = n_rows() * __first_col + __first_row;
  const std::valarray<uword> __size(__sz, 2);
  const std::valarray<uword> __stride(__st, 2);

  return MatrixRef<_Tp>(this->_M_elem, __start, __size, __stride);
}

}  // namespace matrix_lib

#endif