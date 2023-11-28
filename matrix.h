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
#include <sstream>
#include <string>
#include <valarray>

#if defined(__MATRIX_LIB_USE_CPP11)
#include <initializer_list>
#endif

#if defined(__MATRIX_LIB_USE_R)
#include <R_ext/BLAS.h>
#include <RcppCommon.h>
#else
#include <cstdlib>  // std::abort()
#endif

namespace matrix_lib {

//-----------------------------------------------------------------------------

#if defined(__MATRIX_LIB_USE_R)
inline void error(const char* __p) { Rcpp::stop(__p); }
#else
inline void error(const char* __p) {
  std::cerr << __p << std::endl;
  std::abort();
}
#endif

#ifndef NDEBUG
#define __assert(expr, msg)                                             \
  if (!(expr)) {                                                        \
    std::stringstream ss;                                               \
    ss << std::endl                                                     \
       << __FILE__ << ":" << __LINE__ << " assertion failed: " << #expr \
       << "\n"                                                          \
       << msg << std::endl;                                             \
    error(ss.str().c_str());                                            \
  }
#else
#define __assert(expr)
#endif

//-----------------------------------------------------------------------------

typedef std::size_t uword;
typedef std::valarray<bool> bool_array;
typedef std::valarray<std::size_t> index_array;

//-----------------------------------------------------------------------------

// The general Matrix template is simply a prop for its specializations:
template <class _Tp = double, uword _Size = 1>
struct Matrix {
  // multidimensional matrix class
  // ( ) does multidimensional subscripting with range check
  // [ ] retrieves single element or portions of the matrix
 private:
  Matrix();  // this should never be compiled
             //	template<class A> Matrix(A);
};

template <class _Tp = double>
struct SliceMatrix;
template <class _Tp = double, uword _Size = 1>
struct GsliceMatrix;
template <class _Tp = double, uword _Size = 1>
struct IndirectMatrix;
template <class _Tp = double>
struct MaskMatrix;

//-----------------------------------------------------------------------------

// Matrix_base represents the common part of the Matrix classes:
template <class _Tp>
struct _Matrix_base {
 protected:
  std::valarray<_Tp> _M_elem;

 public:
  typedef _Tp elem_type;

  // clang-format off
  // construct/destroy:
  _Matrix_base() : _M_elem() {}
  _Matrix_base(uword __n) : _M_elem(__n) {}
  _Matrix_base(const elem_type& __x, uword __n) : _M_elem(__x, __n) {}
  _Matrix_base(const elem_type* __x, uword __n) : _M_elem(__x, __n) {}
  _Matrix_base(const std::valarray<_Tp>& __x) : _M_elem(__x) {}
#if defined(__MATRIX_LIB_USE_CPP11)
  _Matrix_base(const _Matrix_base& __x) = default;
  _Matrix_base(_Matrix_base&& __x) = default;
  _Matrix_base(std::initializer_list<_Tp> __x) : _M_elem(__x) {}
#endif
  _Matrix_base(const std::slice_array<_Tp>   & __x) : _M_elem(__x) {}
  _Matrix_base(const std::gslice_array<_Tp>  & __x) : _M_elem(__x) {}
  _Matrix_base(const std::mask_array<_Tp>    & __x) : _M_elem(__x) {}
  _Matrix_base(const std::indirect_array<_Tp>& __x) : _M_elem(__x) {}
#if defined(__MATRIX_LIB_USE_R)
  template <typename _U = _Tp,
  typename std::enable_if<std::is_same<_U, double>::value, int>::type=0>
    _Matrix_base(SEXP __x) : _M_elem(REAL(__x), Rf_length(__x)) {}
#endif

  virtual ~_Matrix_base() {}

  // assignments
#if defined(__MATRIX_LIB_USE_CPP11)
  _Matrix_base& operator=(const _Matrix_base& __x) = default;
  _Matrix_base& operator=(_Matrix_base&& __x) = default;
#endif

  // elements accessor and mutator functions
  const std::valarray<_Tp>& get_elem() const { return _M_elem; }
  std::valarray<_Tp> get_elem(const std::slice & __x) const { return _M_elem[__x]; }
  std::valarray<_Tp> get_elem(const std::gslice& __x) const { return _M_elem[__x]; }
  std::valarray<_Tp> get_elem(const bool_array & __x) const { return _M_elem[__x]; }
  std::valarray<_Tp> get_elem(const index_array& __x) const { return _M_elem[__x]; }

  void set_elem(const std::valarray<_Tp>      & __x) { _M_elem.resize(__x.size()); _M_elem = __x; }
  void set_elem(const std::slice_array<_Tp>   & __x) { _M_elem = __x; }
  void set_elem(const std::gslice_array<_Tp>  & __x) { _M_elem = __x; }
  void set_elem(const std::mask_array<_Tp>    & __x) { _M_elem = __x; }
  void set_elem(const std::indirect_array<_Tp>& __x) { _M_elem = __x; }
#if defined(__MATRIX_LIB_USE_CPP11)
  void set_elem(std::initializer_list<_Tp>    & __x) { _M_elem = __x; }
#endif
  uword n_elem() const { return _M_elem.size(); }

  // element access
  elem_type                operator[](uword              __n) const { return _M_elem[__n]; }
  elem_type&               operator[](uword              __n)       { return _M_elem[__n]; }

  // subsetting operations with auxiliary type
  std::valarray<_Tp>       operator[](std::slice         __x) const { return _M_elem[__x]; }
  std::slice_array<_Tp>    operator[](std::slice         __x)       { return _M_elem[__x]; }
  std::valarray<_Tp>       operator[](std::gslice        __x) const { return _M_elem[__x]; }
  std::gslice_array<_Tp>   operator[](std::gslice        __x)       { return _M_elem[__x]; }
  std::valarray<_Tp>       operator[](const bool_array & __x) const { return _M_elem[__x]; }
  std::mask_array<_Tp>     operator[](const bool_array & __x)       { return _M_elem[__x]; }
  std::valarray<_Tp>       operator[](const index_array& __x) const { return _M_elem[__x]; }
  std::indirect_array<_Tp> operator[](const index_array& __x)       { return _M_elem[__x]; }

  // if necessay, we can get to the raw matrix:
        elem_type* data()       { return &(_M_elem[0]); }
  const elem_type* data() const { return &(_M_elem[0]); }
  // clang-format on

 public:  // Other member functions.
          // The results are undefined for zero-length arrays
  elem_type sum() const { return _M_elem.sum(); }
  elem_type min() const { return _M_elem.min(); }
  elem_type max() const { return _M_elem.max(); }

  template <typename F>
  void for_each(F f) {
    for (uword i = 0; i < n_elem(); ++i) f(_M_elem[i]);
  }
};

//-----------------------------------------------------------------------------

template <class _Tp>
struct Matrix<_Tp, 1> : public _Matrix_base<_Tp> {
 public:
  typedef _Tp elem_type;
  bool is_column_vector;
  friend struct Matrix<_Tp, 2>;
  friend struct Matrix<_Tp, 3>;

  // clang-format off
  // construct/destroy:
  Matrix() : _Matrix_base<_Tp>() { _M_init(); }
  Matrix(const         Matrix        & __x) : _Matrix_base<_Tp>(__x._M_elem) { _M_init(); }
  Matrix(const    SliceMatrix<_Tp   >& __x) : _Matrix_base<_Tp>(__x._M_elem[__x._M_desc]) { _M_init(__x.is_column_vector); }
  Matrix(const   GsliceMatrix<_Tp, 1>& __x) : _Matrix_base<_Tp>(__x._M_elem[__x._M_desc]) { _M_init(); }
  Matrix(const IndirectMatrix<_Tp, 1>& __x) : _Matrix_base<_Tp>(__x._M_elem[__x._M_desc]) { _M_init(); }
  Matrix(const     MaskMatrix<_Tp   >& __x) : _Matrix_base<_Tp>(__x._M_elem[__x._M_desc]) { _M_init(); }
#if defined(__MATRIX_LIB_USE_CPP11)
  Matrix(Matrix&& __x) = default;
  Matrix(std::initializer_list<_Tp> __x) : _Matrix_base<_Tp>(__x) { _M_init(); }
  Matrix(Matrix<_Tp, 2>&& __x);
#endif
  explicit Matrix(uword __n1) : _Matrix_base<_Tp>(__n1) { _M_init(); }
  Matrix(const elem_type& __x, uword __n1) : _Matrix_base<_Tp>(__x, __n1) { _M_init(); }
  Matrix(const elem_type* __x, uword __n1) : _Matrix_base<_Tp>(__x, __n1) { _M_init(); }
  Matrix(const std::valarray<_Tp>      & __x) : _Matrix_base<_Tp>(__x) { _M_init(); }
  Matrix(const std::slice_array<_Tp>   & __x) : _Matrix_base<_Tp>(__x) { _M_init(); }
  Matrix(const std::gslice_array<_Tp>  & __x) : _Matrix_base<_Tp>(__x) { _M_init(); }
  Matrix(const std::mask_array<_Tp>    & __x) : _Matrix_base<_Tp>(__x) { _M_init(); }
  Matrix(const std::indirect_array<_Tp>& __x) : _Matrix_base<_Tp>(__x) { _M_init(); }
  Matrix(const std::valarray<_Tp>& __x, const uword      __dims[1]) : _Matrix_base<_Tp>(__x) { _M_init(__dims); }
  Matrix(const std::valarray<_Tp>& __x, const index_array& __dims ) : _Matrix_base<_Tp>(__x) { _M_init(__dims); }
  Matrix(const Matrix<_Tp, 2>& __x);
#if defined(__MATRIX_LIB_USE_R)
  Matrix(SEXP __x) : _Matrix_base<_Tp>(__x) { _M_init();  }
#endif
  ~Matrix() {}

  // assignment
  Matrix& operator=(const Matrix              & __x) { this->set_elem(__x._M_elem); _M_init(); return *this; }
  Matrix& operator=(const GsliceMatrix<_Tp, 1>& __x) { this->set_elem(__x._M_elem); _M_init(); return *this; }
#if defined(__MATRIX_LIB_USE_CPP11)
  Matrix& operator=(Matrix&& __x) = default;
  Matrix& operator=(std::initializer_list<_Tp>       __x) { this->_M_elem =  __x; _M_init(); return *this; }
  Matrix& operator=(Matrix<_Tp, 2>&& __x);
#endif
  Matrix& operator=(const elem_type               & __x ) { this->_M_elem =  __x;            return *this; }
  Matrix& operator=(const std::slice_array<_Tp>   & __sa) { this->_M_elem = __sa; _M_init(); return *this; }
  Matrix& operator=(const std::gslice_array<_Tp>  & __ga) { this->_M_elem = __ga; _M_init(); return *this; }
  Matrix& operator=(const std::mask_array<_Tp>    & __ma) { this->_M_elem = __ma; _M_init(); return *this; }
  Matrix& operator=(const std::indirect_array<_Tp>& __ia) { this->_M_elem = __ia; _M_init(); return *this; }
  Matrix& operator=(const Matrix<_Tp, 2>& __x);

#if defined(__MATRIX_LIB_USE_R)
  // this operator enables implicit Rcpp::wrap
  operator SEXP()       { return __export_matrix_to_sexp(); }
  operator SEXP() const { return __export_matrix_to_sexp(); }
  SEXP __export_matrix_to_sexp() const;
#endif

  index_array get_dims() const { return index_array(_M_dims, 1); }
  uword         n_rows() const { return  is_column_vector ? _M_dims[0] : (_M_dims[0] ? 1 : 0); }
  uword         n_cols() const { return !is_column_vector ? _M_dims[0] : (_M_dims[0] ? 1 : 0); }

  void _M_range_check(uword __n1) const {
    __assert(__n1 < this->n_elem(),
            "Matrix<T, 1>::_M_range_check:\n" <<
            "  Index " << __n1 << " is out of bound for axis 0 with size " << n_rows());
  }

  // subscripting:
        elem_type& operator()(uword __n1)       { _M_range_check(__n1); return this->_M_elem[__n1]; }
  const elem_type& operator()(uword __n1) const { _M_range_check(__n1); return this->_M_elem[__n1]; }

  // GsliceMatrix related member functions
          Matrix<_Tp, 1> subvec(uword __first_index, uword __last_index) const;
    GsliceMatrix<_Tp, 1> subvec(uword __first_index, uword __last_index);
          Matrix<_Tp, 1> operator()(std::slice __s) const;
    GsliceMatrix<_Tp, 1> operator()(std::slice __s);

  // IndirectMatrix related member functions
          Matrix<_Tp, 1> elem(const index_array& __idx_arr) const;
  IndirectMatrix<_Tp, 1> elem(const index_array& __idx_arr);
          Matrix<_Tp, 1> operator()(const index_array& __idx_arr) const;
  IndirectMatrix<_Tp, 1> operator()(const index_array& __idx_arr);

  // MaskMatrix related member functions
          Matrix<_Tp, 1> operator()(const bool_array& __bool_arr) const;
      MaskMatrix<_Tp   > operator()(const bool_array& __bool_arr);

  Matrix<_Tp, 2> t() const;

 public:
  Matrix operator+() const { return *this; }
  Matrix operator-() const { return Matrix(-this->_M_elem, _M_dims); }
  Matrix operator~() const { return Matrix(~this->_M_elem, _M_dims); }
  Matrix<bool, 1> operator!() const { return Matrix<bool, 1>(!this->_M_elem, _M_dims); }

 public:
  Matrix& operator+=(const elem_type& __x) { this->_M_elem += __x; return *this; }
  Matrix& operator-=(const elem_type& __x) { this->_M_elem -= __x; return *this; }
  Matrix& operator*=(const elem_type& __x) { this->_M_elem *= __x; return *this; }
  Matrix& operator/=(const elem_type& __x) { this->_M_elem /= __x; return *this; }
  Matrix& operator%=(const elem_type& __x) { this->_M_elem %= __x; return *this; }
  Matrix& operator&=(const elem_type& __x) { this->_M_elem &= __x; return *this; }
  Matrix& operator|=(const elem_type& __x) { this->_M_elem |= __x; return *this; }
  Matrix& operator^=(const elem_type& __x) { this->_M_elem ^= __x; return *this; }
  Matrix& operator<<=(const elem_type& __x) { this->_M_elem <<= __x; return *this; }
  Matrix& operator>>=(const elem_type& __x) { this->_M_elem >>= __x; return *this; }

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

  void _M_init(bool __is_colvec = true) {
    _M_dims[0] = this->n_elem();
    is_column_vector = __is_colvec;
  }
  void _M_init(const uword __dims[1], bool __is_colvec = true) {
    if (this->n_elem() != __dims[0])
      error("1D Cstor error: dimension mismatch");
    _M_dims[0] = __dims[0];
    is_column_vector = __is_colvec;
  }
  void _M_init(const index_array& __dims, bool __is_colvec = true) {
    if (this->n_elem() != __dims[0])
      error("1D Cstor error: dimension mismatch");
    _M_dims[0] = __dims[0];
    is_column_vector = __is_colvec;
  }
  Matrix(const std::valarray<_Tp>& __va, uword __start, const uword __size,
         const uword __stride, bool __is_colvec = true)
      : _Matrix_base<_Tp>(__va[std::slice(__start, __size, __stride)]) {
    _M_init(__is_colvec);
  }
  Matrix(const std::valarray<_Tp>& __va, uword __start, const uword __size[1],
         const uword __stride[1], bool __is_colvec = true)
      : _Matrix_base<_Tp>(__va[std::gslice(__start, index_array(__size, 1),
                                           index_array(__stride, 1))]) {
    _M_init(__is_colvec);
  }
  Matrix(const std::valarray<_Tp>& __va, const index_array& __idx_arr,
         const uword __dims[1])
      : _Matrix_base<_Tp>(__va[__idx_arr]) {
    _M_init();
  }
  Matrix(const std::valarray<_Tp>& __va, const bool_array& __bool_arr)
      : _Matrix_base<_Tp>(__va[__bool_arr]) {
    _M_init();
  }
  uword sub2ind(uword __i) const { return __i; }
};

//-----------------------------------------------------------------------------

template <class _Tp>
struct Matrix<_Tp, 2> : public _Matrix_base<_Tp> {
 public:
  typedef _Tp elem_type;
  friend struct Matrix<_Tp, 1>;
  friend struct Matrix<_Tp, 3>;

  // clang-format off
  // construct/destroy:
  Matrix() : _Matrix_base<_Tp>() { _M_init(0, 0); }
  Matrix(const         Matrix        & __x) : _Matrix_base<_Tp>(__x._M_elem) { _M_init(__x._M_dims); }
  Matrix(const    SliceMatrix<_Tp   >& __x) : _Matrix_base<_Tp>(__x._M_elem[__x._M_desc]) { __x.is_column_vector? _M_init(this->n_elem(), 1) : _M_init(1, this->n_elem()); }
  Matrix(const   GsliceMatrix<_Tp, 2>& __x) : _Matrix_base<_Tp>(__x._M_elem[__x._M_desc]) { _M_init(__x._M_dims); }
  Matrix(const IndirectMatrix<_Tp, 2>& __x) : _Matrix_base<_Tp>(__x._M_elem[__x._M_desc]) { _M_init(__x._M_dims); }

#if defined(__MATRIX_LIB_USE_CPP11)
  Matrix(Matrix&& __x) = default;
  Matrix(std::initializer_list<_Tp>      __x, uword __n1, uword __n2) : _Matrix_base<_Tp>(__x) { _M_init(__n1, __n2); }
  Matrix(Matrix<_Tp, 1>&& __x);
#endif
  Matrix(uword __n1, uword __n2) : _Matrix_base<_Tp>(__n1 * __n2) { _M_init(__n1, __n2); }
  Matrix(const elem_type& __x, uword __n1, uword __n2) : _Matrix_base<_Tp>(__x, __n1 * __n2)   { _M_init(__n1, __n2); }
  Matrix(const elem_type* __x, uword __n1, uword __n2) : _Matrix_base<_Tp>(__x, __n1 * __n2)   { _M_init(__n1, __n2); }
  Matrix(const std::valarray<_Tp>      & __x, uword __n1, uword __n2) : _Matrix_base<_Tp>(__x) { _M_init(__n1, __n2); }
  Matrix(const std::slice_array<_Tp>   & __x, uword __n1, uword __n2) : _Matrix_base<_Tp>(__x) { _M_init(__n1, __n2); }
  Matrix(const std::gslice_array<_Tp>  & __x, uword __n1, uword __n2) : _Matrix_base<_Tp>(__x) { _M_init(__n1, __n2); }
  Matrix(const std::mask_array<_Tp>    & __x, uword __n1, uword __n2) : _Matrix_base<_Tp>(__x) { _M_init(__n1, __n2); }
  Matrix(const std::indirect_array<_Tp>& __x, uword __n1, uword __n2) : _Matrix_base<_Tp>(__x) { _M_init(__n1, __n2); }
  Matrix(const std::valarray<_Tp>& __x, const uword        __dims[2]) : _Matrix_base<_Tp>(__x) { _M_init(__dims); }
  Matrix(const std::valarray<_Tp>& __x, const index_array&   __dims ) : _Matrix_base<_Tp>(__x) { _M_init(__dims); }
  Matrix(const Matrix<_Tp, 1>& __x);
#if defined(__MATRIX_LIB_USE_R)
  Matrix(SEXP __x) : _Matrix_base<_Tp>(__x) {
    SEXP __dims = Rf_getAttrib(__x, R_DimSymbol);
    _M_init(INTEGER(__dims)[0], INTEGER(__dims)[1]);
  }
#endif
  ~Matrix() {}

  // assignment
  Matrix& operator=(const Matrix              & __x) { this->set_elem(__x._M_elem); _M_init(__x._M_dims); return *this; }
  Matrix& operator=(const GsliceMatrix<_Tp, 2>& __x) { this->set_elem(__x._M_elem); _M_init(__x._M_dims); return *this; }
#if defined(__MATRIX_LIB_USE_CPP11)
  Matrix& operator=(Matrix&& __x) = default;
  Matrix& operator=(Matrix<_Tp, 1>&& __x);
#endif
  Matrix& operator=(const elem_type& __x) { this->_M_elem = __x; return *this; }
  Matrix& operator=(const Matrix<_Tp, 1>& __x);

#if defined(__MATRIX_LIB_USE_R)
  // this operator enables implicit Rcpp::wrap
  operator SEXP()       { return __export_matrix_to_sexp(); }
  operator SEXP() const { return __export_matrix_to_sexp(); }
  SEXP __export_matrix_to_sexp() const;
#endif

  index_array get_dims() const { return index_array(_M_dims, 2); }
  uword         n_rows() const { return _M_dims[0]; }
  uword         n_cols() const { return _M_dims[1]; }

  void _M_range_check(uword __n1, uword __n2) const {
    __assert(__n1 < n_rows(), "Matrix<T, 2>::_M_range_check: index is out of bound for dimension 1");
    __assert(__n2 < n_cols(), "Matrix<T, 2>::_M_range_check: index is out of bound for dimension 2");
  }

  // subscripting:
        elem_type& operator()(uword __n1, uword __n2)
  { _M_range_check(__n1, __n2); return this->_M_elem[sub2ind(__n1, __n2)]; }
  const elem_type& operator()(uword __n1, uword __n2) const
  { _M_range_check(__n1, __n2); return this->_M_elem[sub2ind(__n1, __n2)]; }

  std::valarray<_Tp> diag() const {
    return this->_M_elem[std::slice(0, std::min(_M_dims[0], _M_dims[1]),
                                    _M_dims[0] + 1)];
  }
  std::slice_array<_Tp> diag() {
    return this->_M_elem[std::slice(0, std::min(_M_dims[0], _M_dims[1]),
                                    _M_dims[0] + 1)];
  }

  // SliceMatrix related member functions
  Matrix<_Tp, 1> row(uword __row_number) const;
  SliceMatrix<_Tp> row(uword __row_number);
  Matrix<_Tp, 1> col(uword __col_number) const;
  SliceMatrix<_Tp> col(uword __col_number);

  // GsliceMatrix related member functions
          Matrix<_Tp, 2> rows(uword __first_row, uword __last_row) const;
    GsliceMatrix<_Tp, 2> rows(uword __first_row, uword __last_row);
          Matrix<_Tp, 2> cols(uword __first_col, uword __last_col) const;
    GsliceMatrix<_Tp, 2> cols(uword __first_col, uword __last_col);
          Matrix<_Tp, 2> submat(uword __first_row, uword __first_col, uword __last_row, uword __last_col) const;
    GsliceMatrix<_Tp, 2> submat(uword __first_row, uword __first_col, uword __last_row, uword __last_col);
          Matrix<_Tp, 2> operator()(std::slice __s1, std::slice __s2) const;
    GsliceMatrix<_Tp, 2> operator()(std::slice __s1, std::slice __s2);

  // IndirectMatrix related member functions
          Matrix<_Tp, 1> elem(const index_array& __idx_arr) const;
  IndirectMatrix<_Tp, 1> elem(const index_array& __idx_arr);
          Matrix<_Tp, 2> operator()(const index_array& __idx_arr1, const index_array& __idx_arr2) const;
  IndirectMatrix<_Tp, 2> operator()(const index_array& __idx_arr1, const index_array& __idx_arr2);

  // MaskMatrix related member functions
          Matrix<_Tp, 1> operator()(const bool_array& __bool_arr) const;
      MaskMatrix<_Tp   > operator()(const bool_array& __bool_arr);

  Matrix<_Tp, 2> t() const;

 public:
  Matrix operator+() const { return *this; }
  Matrix operator-() const { return Matrix(-this->_M_elem, _M_dims); }
  Matrix operator~() const { return Matrix(~this->_M_elem, _M_dims); }
  Matrix<bool, 2> operator!() const { return Matrix<bool, 2>(!this->_M_elem, _M_dims); }

 public:
  Matrix& operator+=(const elem_type& __x) { this->_M_elem += __x; return *this; }
  Matrix& operator-=(const elem_type& __x) { this->_M_elem -= __x; return *this; }
  Matrix& operator*=(const elem_type& __x) { this->_M_elem *= __x; return *this; }
  Matrix& operator/=(const elem_type& __x) { this->_M_elem /= __x; return *this; }
  Matrix& operator%=(const elem_type& __x) { this->_M_elem %= __x; return *this; }
  Matrix& operator&=(const elem_type& __x) { this->_M_elem &= __x; return *this; }
  Matrix& operator|=(const elem_type& __x) { this->_M_elem |= __x; return *this; }
  Matrix& operator^=(const elem_type& __x) { this->_M_elem ^= __x; return *this; }
  Matrix& operator<<=(const elem_type& __x) { this->_M_elem <<= __x; return *this; }
  Matrix& operator>>=(const elem_type& __x) { this->_M_elem >>= __x; return *this; }

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
  void _M_init(const index_array& __dims) {
    _M_dims[0] = __dims[0], _M_dims[1] = __dims[1];
  }
  Matrix(const std::valarray<_Tp>& __va, uword __start, const uword __size[2],
         const uword __stride[2])
      : _Matrix_base<_Tp>(__va[std::gslice(__start, index_array(__size, 2),
                                           index_array(__stride, 2))]) {
    _M_init(__size[1], __size[0]);
  }
  Matrix(const std::valarray<_Tp>& __va, const index_array& __idx_arr,
         const uword __dims[2])
      : _Matrix_base<_Tp>(__va[__idx_arr]) {
    _M_init(__dims[0], __dims[1]);
  }
  uword sub2ind(uword __i, uword __j) const { return __i + __j * _M_dims[0]; }
};

//-----------------------------------------------------------------------------

template <class _Tp>
struct Matrix<_Tp, 3> : public _Matrix_base<_Tp> {
 public:
  typedef _Tp elem_type;

  // clang-format off
  // construct/destroy:
  Matrix() : _Matrix_base<_Tp>() { _M_init(0, 0, 0); }
  Matrix(const         Matrix        & __x) : _Matrix_base<_Tp>(__x._M_elem) { _M_init(__x._M_dims); }
  Matrix(const   GsliceMatrix<_Tp, 3>& __x) : _Matrix_base<_Tp>(__x._M_elem[__x._M_desc]) { _M_init(__x._M_dims); }
  Matrix(const IndirectMatrix<_Tp, 3>& __x) : _Matrix_base<_Tp>(__x._M_elem[__x._M_desc]) { _M_init(__x._M_dims); }

#if defined(__MATRIX_LIB_USE_CPP11)
  Matrix(Matrix&& __x) = default;
  Matrix(std::initializer_list<_Tp> __x,      uword __n1, uword __n2, uword __n3) : _Matrix_base<_Tp>(__x) { _M_init(__n1, __n2, __n3); }
#endif
  Matrix(uword __n1, uword __n2, uword __n3) : _Matrix_base<_Tp>(__n1 * __n2 * __n3) { _M_init(__n1, __n2, __n3); }
  Matrix(const elem_type& __x, uword __n1, uword __n2, uword __n3) : _Matrix_base<_Tp>(__x, __n1 * __n2 * __n3) { _M_init(__n1, __n2, __n3); }
  Matrix(const elem_type* __x, uword __n1, uword __n2, uword __n3) : _Matrix_base<_Tp>(__x, __n1 * __n2 * __n3) { _M_init(__n1, __n2, __n3); }
  Matrix(const std::valarray<_Tp>      & __x, uword __n1, uword __n2, uword __n3) : _Matrix_base<_Tp>(__x) { _M_init(__n1, __n2, __n3); }
  Matrix(const std::slice_array<_Tp>   & __x, uword __n1, uword __n2, uword __n3) : _Matrix_base<_Tp>(__x) { _M_init(__n1, __n2, __n3); }
  Matrix(const std::gslice_array<_Tp>  & __x, uword __n1, uword __n2, uword __n3) : _Matrix_base<_Tp>(__x) { _M_init(__n1, __n2, __n3); }
  Matrix(const std::mask_array<_Tp>    & __x, uword __n1, uword __n2, uword __n3) : _Matrix_base<_Tp>(__x) { _M_init(__n1, __n2, __n3); }
  Matrix(const std::indirect_array<_Tp>& __x, uword __n1, uword __n2, uword __n3) : _Matrix_base<_Tp>(__x) { _M_init(__n1, __n2, __n3); }
  Matrix(const std::valarray<_Tp>& __x, const uword      __dims[3]) : _Matrix_base<_Tp>(__x) { _M_init(__dims); }
  Matrix(const std::valarray<_Tp>& __x, const index_array& __dims ) : _Matrix_base<_Tp>(__x) { _M_init(__dims); }
#if defined(__MATRIX_LIB_USE_R)
  Matrix(SEXP __x) : _Matrix_base<_Tp>(__x) {
    SEXP __dims = Rf_getAttrib(__x, R_DimSymbol);
    _M_init(INTEGER(__dims)[0], INTEGER(__dims)[1], INTEGER(__dims)[2]);
  }
#endif
  ~Matrix() {}

  // assignment
  Matrix& operator=(const Matrix              & __x) { this->set_elem(__x._M_elem); _M_init(__x._M_dims); return *this; }
  Matrix& operator=(const GsliceMatrix<_Tp, 3>& __x) { this->set_elem(__x._M_elem); _M_init(__x._M_dims); return *this; }
#if defined(__MATRIX_LIB_USE_CPP11)
  Matrix& operator=(Matrix&& __x) = default;
#endif
  Matrix& operator=(const elem_type& __x) { this->_M_elem = __x; return *this; }

#if defined(__MATRIX_LIB_USE_R)
  // this operator enables implicit Rcpp::wrap
  operator SEXP()       { return __export_matrix_to_sexp(); }
  operator SEXP() const { return __export_matrix_to_sexp(); }
  SEXP __export_matrix_to_sexp() const;
#endif

  index_array get_dims() const { return index_array(_M_dims, 3); }
  uword         n_rows() const { return _M_dims[0]; }
  uword         n_cols() const { return _M_dims[1]; }
  uword       n_slices() const { return _M_dims[2]; }

  void _M_range_check(uword __n1, uword __n2, uword __n3) const {
    __assert(__n1 < n_rows()  , "Matrix<T, 3>::_M_range_check: index is out of bound for dimension 1");
    __assert(__n2 < n_cols()  , "Matrix<T, 3>::_M_range_check: index is out of bound for dimension 2");
    __assert(__n3 < n_slices(), "Matrix<T, 3>::_M_range_check: index is out of bound for dimension 3");
  }

  // subscripting:
        elem_type& operator()(uword __n1, uword __n2, uword __n3)
  { _M_range_check(__n1, __n2, __n3); return this->_M_elem[sub2ind(__n1, __n2, __n3)]; }
  const elem_type& operator()(uword __n1, uword __n2, uword __n3) const
  { _M_range_check(__n1, __n2, __n3); return this->_M_elem[sub2ind(__n1, __n2, __n3)]; }

  // GsliceMatrix related member functions
          Matrix<_Tp, 2> slice(uword __slice_number) const;
    GsliceMatrix<_Tp, 2> slice(uword __slice_number);
          Matrix<_Tp, 3> slices(uword __first_slice, uword __last_slice) const;
    GsliceMatrix<_Tp, 3> slices(uword __first_slice, uword __last_slice);
          Matrix<_Tp, 3> subcube(uword __first_row, uword __first_col, uword __first_slice,
                                 uword __last_row,  uword __last_col,  uword __last_slice) const;
    GsliceMatrix<_Tp, 3> subcube(uword __first_row, uword __first_col, uword __first_slice,
                                 uword __last_row,  uword __last_col,  uword __last_slice);
          Matrix<_Tp, 3> operator()(std::slice __s1, std::slice __s2, std::slice __s3) const;
    GsliceMatrix<_Tp, 3> operator()(std::slice __s1, std::slice __s2, std::slice __s3);

  // IndirectMatrix related member functions
          Matrix<_Tp, 1> elem(const index_array& __idx_arr) const;
  IndirectMatrix<_Tp, 1> elem(const index_array& __idx_arr);
          Matrix<_Tp, 3> operator()(const index_array&, const index_array&, const index_array&) const;
  IndirectMatrix<_Tp, 3> operator()(const index_array&, const index_array&, const index_array&);

  // MaskMatrix related member functions
          Matrix<_Tp, 1> operator()(const bool_array& __bool_arr) const;
      MaskMatrix<_Tp   > operator()(const bool_array& __bool_arr);

 public:
  Matrix operator+() const { return *this; }
  Matrix operator-() const { return Matrix(-this->_M_elem, _M_dims); }
  Matrix operator~() const { return Matrix(~this->_M_elem, _M_dims); }
  Matrix<bool, 3> operator!() const { return Matrix<bool, 3>(!this->_M_elem, _M_dims); }

 public:
  Matrix& operator+=(const elem_type& __x) { this->_M_elem += __x; return *this; }
  Matrix& operator-=(const elem_type& __x) { this->_M_elem -= __x; return *this; }
  Matrix& operator*=(const elem_type& __x) { this->_M_elem *= __x; return *this; }
  Matrix& operator/=(const elem_type& __x) { this->_M_elem /= __x; return *this; }
  Matrix& operator%=(const elem_type& __x) { this->_M_elem %= __x; return *this; }
  Matrix& operator&=(const elem_type& __x) { this->_M_elem &= __x; return *this; }
  Matrix& operator|=(const elem_type& __x) { this->_M_elem |= __x; return *this; }
  Matrix& operator^=(const elem_type& __x) { this->_M_elem ^= __x; return *this; }
  Matrix& operator<<=(const elem_type& __x) { this->_M_elem <<= __x; return *this; }
  Matrix& operator>>=(const elem_type& __x) { this->_M_elem >>= __x; return *this; }

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
  void _M_init(const index_array& __dims) {
    _M_dims[0] = __dims[0], _M_dims[1] = __dims[1], _M_dims[2] = __dims[2];
    _M_d1xd2 = __dims[0] * __dims[1];
  }
  Matrix(const std::valarray<_Tp>& __va, uword __start, const uword __size[3],
         const uword __stride[3])
      : _Matrix_base<_Tp>(__va[std::gslice(__start, index_array(__size, 3),
                                           index_array(__stride, 3))]) {
    _M_init(__size[2], __size[1], __size[0]);
  }
  Matrix(const std::valarray<_Tp>& __va, const index_array& __idx_arr,
         const uword __dims[3])
      : _Matrix_base<_Tp>(__va[__idx_arr]) {
    _M_init(__dims[0], __dims[1], __dims[2]);
  }
  uword sub2ind(uword __i, uword __j, uword __k) const {
    return __i + __j * _M_dims[0] + __k * _M_d1xd2;
  }
};

//-----------------------------------------------------------------------------

template <class _Tp>
Matrix<_Tp, 1>::Matrix(const Matrix<_Tp, 2>& __x)
    : _Matrix_base<_Tp>(__x.get_elem()) {
  if (__x.n_cols() != 1) error("__x is not a n x 1 matrix");
  _M_init();
}

template <class _Tp>
Matrix<_Tp, 2>::Matrix(const Matrix<_Tp, 1>& __x)
    : _Matrix_base<_Tp>(__x.get_elem()) {
  __x.is_column_vector ? _M_init(__x.n_elem(), 1) : _M_init(1, __x.n_elem());
}

template <class _Tp>
Matrix<_Tp, 1>& Matrix<_Tp, 1>::operator=(const Matrix<_Tp, 2>& __x) {
  if (__x.n_cols() != 1) error("__x is not a n x 1 matrix");
  if (__x.n_rows() != n_rows()) this->_M_elem.resize(__x.n_rows());
  this->_M_elem = __x.get_elem();
  _M_init();
  return *this;
}

template <class _Tp>
Matrix<_Tp, 2>& Matrix<_Tp, 2>::operator=(const Matrix<_Tp, 1>& __x) {
  if (__x.n_rows() != n_rows()) this->_M_elem.resize(__x.n_rows());
  this->_M_elem = __x.get_elem();
  _M_init(__x.n_elem(), 1);
  return *this;
}

#if defined(__MATRIX_LIB_USE_CPP11)
template <class _Tp>
Matrix<_Tp, 1>::Matrix(Matrix<_Tp, 2>&& __x)
    : _Matrix_base<_Tp>(std::move(__x._M_elem)) {
  if (__x.n_cols() != 1) error("__x is not a n x 1 matrix");
  _M_init();
}

template <class _Tp>
Matrix<_Tp, 2>::Matrix(Matrix<_Tp, 1>&& __x)
    : _Matrix_base<_Tp>(std::move(__x._M_elem)) {
  __x.is_column_vector ? _M_init(__x.n_elem(), 1) : _M_init(1, __x.n_elem());
}

template <class _Tp>
Matrix<_Tp, 1>& Matrix<_Tp, 1>::operator=(Matrix<_Tp, 2>&& __x) {
  if (__x.n_cols() != 1) error("__x is not a n x 1 matrix");
  this->_M_elem = std::move(__x._M_elem);
  _M_init();
  return *this;
}

template <class _Tp>
Matrix<_Tp, 2>& Matrix<_Tp, 2>::operator=(Matrix<_Tp, 1>&& __x) {
  this->_M_elem = std::move(__x._M_elem);
  _M_init(this->n_elem(), 1);
  return *this;
}
#endif

// Matrix Rcpp-related member function

#if defined(__MATRIX_LIB_USE_R)
template <class _Tp>
SEXP Matrix<_Tp, 1>::__export_matrix_to_sexp() const {
  SEXP __x = PROTECT(Rcpp::wrap(this->data(), this->data() + this->n_elem()));
  SEXP __dims = PROTECT(Rf_allocVector(INTSXP, 1));
  INTEGER(__dims)[0] = n_rows();
  Rf_setAttrib(__x, R_DimSymbol, __dims);
  UNPROTECT(2);
  return __x;
}

template <class _Tp>
SEXP Matrix<_Tp, 2>::__export_matrix_to_sexp() const {
  SEXP __x = PROTECT(Rcpp::wrap(this->data(), this->data() + this->n_elem()));
  SEXP __dims = PROTECT(Rf_allocVector(INTSXP, 2));
  INTEGER(__dims)[0] = n_rows();
  INTEGER(__dims)[1] = n_cols();
  Rf_setAttrib(__x, R_DimSymbol, __dims);
  UNPROTECT(2);
  return __x;
}

template <class _Tp>
SEXP Matrix<_Tp, 3>::__export_matrix_to_sexp() const {
  SEXP __x = PROTECT(Rcpp::wrap(this->data(), this->data() + this->n_elem()));
  SEXP __dims = PROTECT(Rf_allocVector(INTSXP, 3));
  INTEGER(__dims)[0] = n_rows();
  INTEGER(__dims)[1] = n_cols();
  INTEGER(__dims)[2] = n_slices();
  Rf_setAttrib(__x, R_DimSymbol, __dims);
  UNPROTECT(2);
  return __x;
}
#endif

//----------------------------------------------------------------------
// Matrix non-member functions.

bool __same_dims(const std::valarray<uword>& __d1,
                 const std::valarray<uword>& __d2) {
  uword __n = __d1.size();
  for (uword __i = 0; __i != __n; ++__i)
    if (__d1[__i] != __d2[__i]) return false;
  return true;
}

// Binary arithmetic operations between two Matrix.

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator+(const Matrix<_Tp, _Size>& __x,
                                    const Matrix<_Tp, _Size>& __y) {
  __assert(__same_dims(__x.get_dims(), __y.get_dims()), "dimension mismatch");
  Matrix<_Tp, _Size> __tmp(__x);
  return __tmp += __y;
}

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator-(const Matrix<_Tp, _Size>& __x,
                                    const Matrix<_Tp, _Size>& __y) {
  __assert(__same_dims(__x.get_dims(), __y.get_dims()), "dimension mismatch");
  Matrix<_Tp, _Size> __tmp(__x);
  return __tmp -= __y;
}

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator*(const Matrix<_Tp, _Size>& __x,
                                    const Matrix<_Tp, _Size>& __y) {
  __assert(__same_dims(__x.get_dims(), __y.get_dims()), "dimension mismatch");
  Matrix<_Tp, _Size> __tmp(__x);
  return __tmp *= __y;
}

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator/(const Matrix<_Tp, _Size>& __x,
                                    const Matrix<_Tp, _Size>& __y) {
  __assert(__same_dims(__x.get_dims(), __y.get_dims()), "dimension mismatch");
  Matrix<_Tp, _Size> __tmp(__x);
  return __tmp /= __y;
}

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator%(const Matrix<_Tp, _Size>& __x,
                                    const Matrix<_Tp, _Size>& __y) {
  __assert(__same_dims(__x.get_dims(), __y.get_dims()), "dimension mismatch");
  Matrix<_Tp, _Size> __tmp(__x);
  return __tmp %= __y;
}

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator&(const Matrix<_Tp, _Size>& __x,
                                    const Matrix<_Tp, _Size>& __y) {
  __assert(__same_dims(__x.get_dims(), __y.get_dims()), "dimension mismatch");
  Matrix<_Tp, _Size> __tmp(__x);
  return __tmp &= __y;
}

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator|(const Matrix<_Tp, _Size>& __x,
                                    const Matrix<_Tp, _Size>& __y) {
  __assert(__same_dims(__x.get_dims(), __y.get_dims()), "dimension mismatch");
  Matrix<_Tp, _Size> __tmp(__x);
  return __tmp |= __y;
}

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator^(const Matrix<_Tp, _Size>& __x,
                                    const Matrix<_Tp, _Size>& __y) {
  __assert(__same_dims(__x.get_dims(), __y.get_dims()), "dimension mismatch");
  Matrix<_Tp, _Size> __tmp(__x);
  return __tmp ^= __y;
}

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator<<(const Matrix<_Tp, _Size>& __x,
                                     const Matrix<_Tp, _Size>& __y) {
  __assert(__same_dims(__x.get_dims(), __y.get_dims()), "dimension mismatch");
  Matrix<_Tp, _Size> __tmp(__x);
  return __tmp <<= __y;
}

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator>>(const Matrix<_Tp, _Size>& __x,
                                     const Matrix<_Tp, _Size>& __y) {
  __assert(__same_dims(__x.get_dims(), __y.get_dims()), "dimension mismatch");
  Matrix<_Tp, _Size> __tmp(__x);
  return __tmp >>= __y;
}

template <class _Tp, uword _Size>
inline bool_array operator&&(const Matrix<_Tp, _Size>& __x,
                             const Matrix<_Tp, _Size>& __y) {
  __assert(__same_dims(__x.get_dims(), __y.get_dims()), "dimension mismatch");
  return __x.get_elem() && __y.get_elem();
}

template <class _Tp, uword _Size>
inline bool_array operator||(const Matrix<_Tp, _Size>& __x,
                             const Matrix<_Tp, _Size>& __y) {
  __assert(__same_dims(__x.get_dims(), __y.get_dims()), "dimension mismatch");
  return __x.get_elem() || __y.get_elem();
}

// Binary arithmetic operations between an array and a scalar.

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator+(const Matrix<_Tp, _Size>& __x,
                                    const _Tp& __c) {
  Matrix<_Tp, _Size> __tmp(__x);
  return __tmp += __c;
}

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator+(const _Tp& __c,
                                    const Matrix<_Tp, _Size>& __x) {
  return Matrix<_Tp, _Size>(__c + __x.get_elem(), __x.get_dims());
}

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator-(const Matrix<_Tp, _Size>& __x,
                                    const _Tp& __c) {
  Matrix<_Tp, _Size> __tmp(__x);
  return __tmp -= __c;
}

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator-(const _Tp& __c,
                                    const Matrix<_Tp, _Size>& __x) {
  return Matrix<_Tp, _Size>(__c - __x.get_elem(), __x.get_dims());
}

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator*(const Matrix<_Tp, _Size>& __x,
                                    const _Tp& __c) {
  Matrix<_Tp, _Size> __tmp(__x);
  return __tmp *= __c;
}

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator*(const _Tp& __c,
                                    const Matrix<_Tp, _Size>& __x) {
  return Matrix<_Tp, _Size>(__c * __x.get_elem(), __x.get_dims());
}

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator/(const Matrix<_Tp, _Size>& __x,
                                    const _Tp& __c) {
  Matrix<_Tp, _Size> __tmp(__x);
  return __tmp /= __c;
}

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator/(const _Tp& __c,
                                    const Matrix<_Tp, _Size>& __x) {
  Matrix<_Tp, _Size> __tmp(__c / __x.get_elem(), __x.get_dims());
  return __tmp;
}

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator%(const Matrix<_Tp, _Size>& __x,
                                    const _Tp& __c) {
  Matrix<_Tp, _Size> __tmp(__x);
  return __tmp %= __c;
}

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator%(const _Tp& __c,
                                    const Matrix<_Tp, _Size>& __x) {
  Matrix<_Tp, _Size> __tmp(__c % __x.get_elem(), __x.get_dims());
  return __tmp;
}

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator&(const Matrix<_Tp, _Size>& __x,
                                    const _Tp& __c) {
  Matrix<_Tp, _Size> __tmp(__x);
  return __tmp &= __c;
}

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator&(const _Tp& __c,
                                    const Matrix<_Tp, _Size>& __x) {
  Matrix<_Tp, _Size> __tmp(__c & __x.get_elem(), __x.get_dims());
  return __tmp;
}

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator|(const Matrix<_Tp, _Size>& __x,
                                    const _Tp& __c) {
  Matrix<_Tp, _Size> __tmp(__x);
  return __tmp |= __c;
}

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator|(const _Tp& __c,
                                    const Matrix<_Tp, _Size>& __x) {
  Matrix<_Tp, _Size> __tmp(__c | __x.get_elem(), __x.get_dims());
  return __tmp;
}

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator^(const Matrix<_Tp, _Size>& __x,
                                    const _Tp& __c) {
  Matrix<_Tp, _Size> __tmp(__x);
  return __tmp ^= __c;
}

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator^(const _Tp& __c,
                                    const Matrix<_Tp, _Size>& __x) {
  Matrix<_Tp, _Size> __tmp(__c ^ __x.get_elem(), __x.get_dims());
  return __tmp;
}

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator<<(const Matrix<_Tp, _Size>& __x,
                                     const _Tp& __c) {
  Matrix<_Tp, _Size> __tmp(__x);
  return __tmp <<= __c;
}

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator<<(const _Tp& __c,
                                     const Matrix<_Tp, _Size>& __x) {
  Matrix<_Tp, _Size> __tmp(__c << __x.get_elem(), __x.get_dims());
  return __tmp;
}

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator>>(const Matrix<_Tp, _Size>& __x,
                                     const _Tp& __c) {
  Matrix<_Tp, _Size> __tmp(__x);
  return __tmp >>= __c;
}

template <class _Tp, uword _Size>
inline Matrix<_Tp, _Size> operator>>(const _Tp& __c,
                                     const Matrix<_Tp, _Size>& __x) {
  Matrix<_Tp, _Size> __tmp(__c >> __x.get_elem(), __x.get_dims());
  return __tmp;
}

template <class _Tp, uword _Size>
inline bool_array operator&&(const Matrix<_Tp, _Size>& __x, const _Tp& __c) {
  return __x.get_elem() && __c;
}

template <class _Tp, uword _Size>
inline bool_array operator&&(const _Tp& __c, const Matrix<_Tp, _Size>& __x) {
  return __c && __x.get_elem();
}

template <class _Tp, uword _Size>
inline bool_array operator||(const Matrix<_Tp, _Size>& __x, const _Tp& __c) {
  return __x.get_elem() || __c;
}

template <class _Tp, uword _Size>
inline bool_array operator||(const _Tp& __c, const Matrix<_Tp, _Size>& __x) {
  return __c || __x.get_elem();
}

// Binary logical operations between two Matrices.

template <class _Tp, uword _Size>
inline bool_array operator==(const Matrix<_Tp, _Size>& __x,
                             const Matrix<_Tp, _Size>& __y) {
  __assert(__same_dims(__x.get_dims(), __y.get_dims()), "dimension mismatch");
  return __x.get_elem() == __y.get_elem();
}

template <class _Tp, uword _Size>
inline bool_array operator!=(const Matrix<_Tp, _Size>& __x,
                             const Matrix<_Tp, _Size>& __y) {
  __assert(__same_dims(__x.get_dims(), __y.get_dims()), "dimension mismatch");
  return __x.get_elem() != __y.get_elem();
}

template <class _Tp, uword _Size>
inline bool_array operator<(const Matrix<_Tp, _Size>& __x,
                            const Matrix<_Tp, _Size>& __y) {
  __assert(__same_dims(__x.get_dims(), __y.get_dims()), "dimension mismatch");
  return __x.get_elem() < __y.get_elem();
}

template <class _Tp, uword _Size>
inline bool_array operator<=(const Matrix<_Tp, _Size>& __x,
                             const Matrix<_Tp, _Size>& __y) {
  __assert(__same_dims(__x.get_dims(), __y.get_dims()), "dimension mismatch");
  return __x.get_elem() <= __y.get_elem();
}

template <class _Tp, uword _Size>
inline bool_array operator>(const Matrix<_Tp, _Size>& __x,
                            const Matrix<_Tp, _Size>& __y) {
  __assert(__same_dims(__x.get_dims(), __y.get_dims()), "dimension mismatch");
  return __x.get_elem() > __y.get_elem();
}

template <class _Tp, uword _Size>
inline bool_array operator>=(const Matrix<_Tp, _Size>& __x,
                             const Matrix<_Tp, _Size>& __y) {
  __assert(__same_dims(__x.get_dims(), __y.get_dims()), "dimension mismatch");
  return __x.get_elem() >= __y.get_elem();
}

// Logical operations between a Matrix and a scalar.

template <class _Tp, uword _Size>
inline bool_array operator==(const Matrix<_Tp, _Size>& __x, const _Tp& __c) {
  return __x.get_elem() == __c;
}

template <class _Tp, uword _Size>
inline bool_array operator==(const _Tp& __c, const Matrix<_Tp, _Size>& __x) {
  return __c == __x.get_elem();
}

template <class _Tp, uword _Size>
inline bool_array operator!=(const Matrix<_Tp, _Size>& __x, const _Tp& __c) {
  return __x.get_elem() != __c;
}

template <class _Tp, uword _Size>
inline bool_array operator!=(const _Tp& __c, const Matrix<_Tp, _Size>& __x) {
  return __c != __x.get_elem();
}

template <class _Tp, uword _Size>
inline bool_array operator<(const Matrix<_Tp, _Size>& __x, const _Tp& __c) {
  return __x.get_elem() < __c;
}

template <class _Tp, uword _Size>
inline bool_array operator<(const _Tp& __c, const Matrix<_Tp, _Size>& __x) {
  return __c < __x.get_elem();
}

template <class _Tp, uword _Size>
inline bool_array operator<=(const Matrix<_Tp, _Size>& __x, const _Tp& __c) {
  return __x.get_elem() <= __c;
}

template <class _Tp, uword _Size>
inline bool_array operator<=(const _Tp& __c, const Matrix<_Tp, _Size>& __x) {
  return __c <= __x.get_elem();
}

template <class _Tp, uword _Size>
inline bool_array operator>(const Matrix<_Tp, _Size>& __x, const _Tp& __c) {
  return __x.get_elem() > __c;
}

template <class _Tp, uword _Size>
inline bool_array operator>(const _Tp& __c, const Matrix<_Tp, _Size>& __x) {
  return __c > __x.get_elem();
}

template <class _Tp, uword _Size>
inline bool_array operator>=(const Matrix<_Tp, _Size>& __x, const _Tp& __c) {
  return __x.get_elem() >= __c;
}

template <class _Tp, uword _Size>
inline bool_array operator>=(const _Tp& __c, const Matrix<_Tp, _Size>& __x) {
  return __c >= __x.get_elem();
}

// Matrix "transcendentals" (the list includes abs and sqrt, which,
// of course, are not transcendental).

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

template <class _Tp>
_Tp dot(const Matrix<_Tp, 1>& __x, const Matrix<_Tp, 1>& __y) {
  return (__x.get_elem() * __y.get_elem()).sum();
}

//----------------------------------------------------------------------
// SliceMatrix

template <class _Tp>
struct SliceMatrix {
 public:
  typedef _Tp elem_type;
  std::valarray<_Tp>& _M_elem;
  std::slice _M_desc;
  bool is_column_vector;
  SliceMatrix(std::valarray<_Tp>& __va, uword __start, uword __size,
              uword __stride, bool __is_colvec = true)
      : _M_elem(__va),
        _M_desc(__start, __size, __stride),
        is_column_vector(__is_colvec) {}
  void operator=(const _Tp& __value) { _M_elem[_M_desc] = __value; }

  std::valarray<_Tp> get_elem() const {
    return std::valarray<_Tp>(_M_elem[_M_desc]);
  }

  elem_type sum() const { return get_elem().sum(); }
  elem_type min() const { return get_elem().min(); }
  elem_type max() const { return get_elem().max(); }
};

//----------------------------------------------------------------------
// GsliceMatrix

template <class _Tp, uword _Size>
struct GsliceMatrix {
 public:
  typedef _Tp elem_type;
  std::valarray<_Tp>& _M_elem;
  std::gslice _M_desc;
  uword _M_dims[_Size];
  GsliceMatrix(std::valarray<_Tp>& __va, uword __start,
               const index_array& __size, const index_array& __stride)
      : _M_elem(__va), _M_desc(__start, __size, __stride) {
    uword n = __size.size();
    for (uword idx = 0; idx < n; ++idx) {
      _M_dims[idx] = __size[n - idx - 1];
    }
  }
  GsliceMatrix(std::valarray<_Tp>& __va, uword __start,
               const uword __size[_Size], const uword __stride[_Size])
      : _M_elem(__va),
        _M_desc(__start, index_array(__size, _Size),
                index_array(__stride, _Size)) {
    uword n = _Size;
    for (uword idx = 0; idx < n; ++idx) {
      _M_dims[idx] = __size[n - idx - 1];
    }
  }
  void operator=(const _Tp& __value) { _M_elem[_M_desc] = __value; }

  std::valarray<_Tp> get_elem() const {
    return std::valarray<_Tp>(_M_elem[_M_desc]);
  }

  elem_type sum() const { return get_elem().sum(); }
  elem_type min() const { return get_elem().min(); }
  elem_type max() const { return get_elem().max(); }
};

//----------------------------------------------------------------------
// IndirectMatrix

template <class _Tp, uword _Size>
struct IndirectMatrix {
 public:
  typedef _Tp elem_type;
  std::valarray<_Tp>& _M_elem;
  index_array _M_desc;
  uword _M_dims[_Size];
  IndirectMatrix(std::valarray<_Tp>& __va, const index_array& __ind_arr,
                 const uword __dims[_Size])
      : _M_elem(__va), _M_desc(__ind_arr) {
    for (uword idx = 0; idx < _Size; ++idx) {
      _M_dims[idx] = __dims[idx];
    }
  }
  void operator=(const _Tp& __value) { _M_elem[_M_desc] = __value; }

  std::valarray<_Tp> get_elem() const {
    return std::valarray<_Tp>(_M_elem[_M_desc]);
  }

  elem_type sum() const { return get_elem().sum(); }
  elem_type min() const { return get_elem().min(); }
  elem_type max() const { return get_elem().max(); }
};

//----------------------------------------------------------------------
// MaskMatrix

template <class _Tp>
struct MaskMatrix {
 public:
  typedef _Tp elem_type;
  std::valarray<_Tp>& _M_elem;
  bool_array _M_desc;
  MaskMatrix(std::valarray<_Tp>& __va, const bool_array& __boolarr)
      : _M_elem(__va), _M_desc(__boolarr) {}
  void operator=(const _Tp& __value) { _M_elem[_M_desc] = __value; }

  std::valarray<_Tp> get_elem() const {
    return std::valarray<_Tp>(_M_elem[_M_desc]);
  }

  elem_type sum() const { return get_elem().sum(); }
  elem_type min() const { return get_elem().min(); }
  elem_type max() const { return get_elem().max(); }
};

// Matrix member functions dealing with SliceMatrix

template <class _Tp>
inline Matrix<_Tp, 1> Matrix<_Tp, 2>::row(uword __r) const {
  _M_range_check(__r, 0);
  const uword __start = __r;
  const uword __size = _M_dims[1];
  const uword __stride = _M_dims[0];
  return Matrix<_Tp, 1>(this->_M_elem, __start, __size, __stride, false);
}

template <class _Tp>
inline SliceMatrix<_Tp> Matrix<_Tp, 2>::row(uword __r) {
  _M_range_check(__r, 0);
  const uword __start = __r;
  const uword __size = _M_dims[1];
  const uword __stride = _M_dims[0];
  return SliceMatrix<_Tp>(this->_M_elem, __start, __size, __stride, false);
}

template <class _Tp>
inline Matrix<_Tp, 1> Matrix<_Tp, 2>::col(uword __c) const {
  _M_range_check(0, __c);
  const uword __start = __c * _M_dims[0];
  const uword __size = _M_dims[0];
  const uword __stride = 1;
  return Matrix<_Tp, 1>(this->_M_elem, __start, __size, __stride);
}

template <class _Tp>
inline SliceMatrix<_Tp> Matrix<_Tp, 2>::col(uword __c) {
  _M_range_check(0, __c);
  const uword __start = __c * _M_dims[0];
  const uword __size = _M_dims[0];
  const uword __stride = 1;
  return SliceMatrix<_Tp>(this->_M_elem, __start, __size, __stride);
}

// Matrix member functions dealing with GsliceMatrix

template <class _Tp>
inline Matrix<_Tp, 2> Matrix<_Tp, 3>::slice(uword __s) const {
  const uword __start = __s * _M_d1xd2;
  const uword __size[2] = {n_cols(), n_rows()};
  const uword __stride[2] = {n_rows(), 1};
  return Matrix<_Tp, 2>(this->_M_elem, __start, __size, __stride);
}

template <class _Tp>
inline GsliceMatrix<_Tp, 2> Matrix<_Tp, 3>::slice(uword __s) {
  const uword __start = __s * _M_d1xd2;
  const uword __size[2] = {n_cols(), n_rows()};
  const uword __stride[2] = {n_rows(), 1};
  return GsliceMatrix<_Tp, 2>(this->_M_elem, __start, __size, __stride);
}

template <class _Tp>
inline Matrix<_Tp, 2> Matrix<_Tp, 2>::rows(uword __fr, uword __lr) const {
  return submat(__fr, 0, __lr, n_cols() - 1);
}

template <class _Tp>
inline GsliceMatrix<_Tp, 2> Matrix<_Tp, 2>::rows(uword __fr, uword __lr) {
  return submat(__fr, 0, __lr, n_cols() - 1);
}

template <class _Tp>
inline Matrix<_Tp, 2> Matrix<_Tp, 2>::cols(uword __fc, uword __lc) const {
  return submat(0, __fc, n_rows() - 1, __lc);
}

template <class _Tp>
inline GsliceMatrix<_Tp, 2> Matrix<_Tp, 2>::cols(uword __fc, uword __lc) {
  return submat(0, __fc, n_rows() - 1, __lc);
}

template <class _Tp>
inline Matrix<_Tp, 3> Matrix<_Tp, 3>::slices(uword __fs, uword __ls) const {
  return subcube(0, 0, __fs, n_rows() - 1, n_cols() - 1, __ls);
}

template <class _Tp>
inline GsliceMatrix<_Tp, 3> Matrix<_Tp, 3>::slices(uword __fs, uword __ls) {
  return subcube(0, 0, __fs, n_rows() - 1, n_cols() - 1, __ls);
}

template <class _Tp>
inline Matrix<_Tp, 1> Matrix<_Tp, 1>::subvec(uword __i, uword __j) const {
  if (__i > __j || __j >= _M_dims[0]) error("1D subscription error");
  const uword __start = __i;
  const uword __size[1] = {__j - __i + 1};
  const uword __stride[1] = {1};
  return Matrix<_Tp, 1>(this->_M_elem, __start, __size, __stride);
}

template <class _Tp>
inline GsliceMatrix<_Tp, 1> Matrix<_Tp, 1>::subvec(uword __i, uword __j) {
  if (__i > __j || __j >= _M_dims[0]) error("1D subscription error");
  const uword __start = __i;
  const uword __size[1] = {__j - __i + 1};
  const uword __stride[1] = {1};
  return GsliceMatrix<_Tp, 1>(this->_M_elem, __start, __size, __stride);
}

template <class _Tp>
inline Matrix<_Tp, 2> Matrix<_Tp, 2>::submat(uword __fr, uword __fc, uword __lr,
                                             uword __lc) const {
  const uword __start = n_rows() * __fc + __fr;
  const uword __size[2] = {__lc - __fc + 1, __lr - __fr + 1};
  const uword __stride[2] = {n_rows(), 1};
  return Matrix<_Tp, 2>(this->_M_elem, __start, __size, __stride);
}

template <class _Tp>
inline GsliceMatrix<_Tp, 2> Matrix<_Tp, 2>::submat(uword __fr, uword __fc,
                                                   uword __lr, uword __lc) {
  const uword __start = n_rows() * __fc + __fr;
  const uword __size[2] = {__lc - __fc + 1, __lr - __fr + 1};
  const uword __stride[2] = {n_rows(), 1};
  return GsliceMatrix<_Tp, 2>(this->_M_elem, __start, __size, __stride);
}

template <class _Tp>
inline Matrix<_Tp, 3> Matrix<_Tp, 3>::subcube(uword __fr, uword __fc,
                                              uword __fs, uword __lr,
                                              uword __lc, uword __ls) const {
  const uword __start = sub2ind(__fr, __fc, __fs);
  const uword __size[3] = {__ls - __fs + 1, __lc - __fc + 1, __lr - __fr + 1};
  const uword __stride[3] = {_M_d1xd2, n_rows(), 1};
  return Matrix<_Tp, 3>(this->_M_elem, __start, __size, __stride);
}

template <class _Tp>
inline GsliceMatrix<_Tp, 3> Matrix<_Tp, 3>::subcube(uword __fr, uword __fc,
                                                    uword __fs, uword __lr,
                                                    uword __lc, uword __ls) {
  const uword __start = sub2ind(__fr, __fc, __fs);
  const uword __size[3] = {__ls - __fs + 1, __lc - __fc + 1, __lr - __fr + 1};
  const uword __stride[3] = {_M_d1xd2, n_rows(), 1};
  return GsliceMatrix<_Tp, 3>(this->_M_elem, __start, __size, __stride);
}

template <class _Tp>
inline Matrix<_Tp, 1> Matrix<_Tp, 1>::operator()(std::slice __s) const {
  const uword __start = __s.start();
  const uword __size[1] = {__s.size()};
  const uword __stride[1] = {__s.stride()};
  return Matrix<_Tp, 1>(this->_M_elem, __start, __size, __stride);
}

template <class _Tp>
inline GsliceMatrix<_Tp, 1> Matrix<_Tp, 1>::operator()(std::slice __s) {
  const uword __start = __s.start();
  const uword __size[1] = {__s.size()};
  const uword __stride[1] = {__s.stride()};
  return GsliceMatrix<_Tp, 1>(this->_M_elem, __start, __size, __stride);
}

template <class _Tp>
inline Matrix<_Tp, 2> Matrix<_Tp, 2>::operator()(std::slice __s1,
                                                 std::slice __s2) const {
  const uword __start = sub2ind(__s1.start(), __s2.start());
  const uword __size[2] = {__s2.size(), __s1.size()};
  const uword __stride[2] = {n_rows() * __s2.stride(), 1 * __s1.stride()};
  return Matrix<_Tp, 2>(this->_M_elem, __start, __size, __stride);
}

template <class _Tp>
inline GsliceMatrix<_Tp, 2> Matrix<_Tp, 2>::operator()(std::slice __s1,
                                                       std::slice __s2) {
  const uword __start = sub2ind(__s1.start(), __s2.start());
  const uword __size[2] = {__s2.size(), __s1.size()};
  const uword __stride[2] = {n_rows() * __s2.stride(), 1 * __s1.stride()};
  return GsliceMatrix<_Tp, 2>(this->_M_elem, __start, __size, __stride);
}

template <class _Tp>
inline Matrix<_Tp, 3> Matrix<_Tp, 3>::operator()(std::slice __s1,
                                                 std::slice __s2,
                                                 std::slice __s3) const {
  const uword __start = sub2ind(__s1.start(), __s2.start(), __s3.start());
  const uword __size[3] = {__s3.size(), __s2.size(), __s1.size()};
  const uword __stride[3] = {_M_d1xd2 * __s3.stride(), n_rows() * __s2.stride(),
                             1 * __s1.stride()};
  return Matrix<_Tp, 3>(this->_M_elem, __start, __size, __stride);
}

template <class _Tp>
inline GsliceMatrix<_Tp, 3> Matrix<_Tp, 3>::operator()(std::slice __s1,
                                                       std::slice __s2,
                                                       std::slice __s3) {
  const uword __start = sub2ind(__s1.start(), __s2.start(), __s3.start());
  const uword __size[3] = {__s3.size(), __s2.size(), __s1.size()};
  const uword __stride[3] = {_M_d1xd2 * __s3.stride(), n_rows() * __s2.stride(),
                             1 * __s1.stride()};
  return GsliceMatrix<_Tp, 3>(this->_M_elem, __start, __size, __stride);
}

// Matrix member functions dealing with index_array and IndirectMatrix

template <class _Tp>
inline Matrix<_Tp, 1> Matrix<_Tp, 1>::elem(const index_array& __idx_arr) const {
  uword dims[1] = {__idx_arr.size()};
  return Matrix<_Tp, 1>(this->_M_elem, __idx_arr, dims);
}

template <class _Tp>
inline IndirectMatrix<_Tp, 1> Matrix<_Tp, 1>::elem(
    const index_array& __idx_arr) {
  uword dims[1] = {__idx_arr.size()};
  return IndirectMatrix<_Tp, 1>(this->_M_elem, __idx_arr, dims);
}
template <class _Tp>
inline Matrix<_Tp, 1> Matrix<_Tp, 2>::elem(const index_array& __idx_arr) const {
  uword dims[1] = {__idx_arr.size()};
  return Matrix<_Tp, 1>(this->_M_elem, __idx_arr, dims);
}

template <class _Tp>
inline IndirectMatrix<_Tp, 1> Matrix<_Tp, 2>::elem(
    const index_array& __idx_arr) {
  uword dims[1] = {__idx_arr.size()};
  return IndirectMatrix<_Tp, 1>(this->_M_elem, __idx_arr, dims);
}

template <class _Tp>
inline Matrix<_Tp, 1> Matrix<_Tp, 3>::elem(const index_array& __idx_arr) const {
  uword dims[1] = {__idx_arr.size()};
  return Matrix<_Tp, 1>(this->_M_elem, __idx_arr, dims);
}

template <class _Tp>
inline IndirectMatrix<_Tp, 1> Matrix<_Tp, 3>::elem(
    const index_array& __idx_arr) {
  uword dims[1] = {__idx_arr.size()};
  return IndirectMatrix<_Tp, 1>(this->_M_elem, __idx_arr, dims);
}

template <class _Tp>
inline Matrix<_Tp, 1> Matrix<_Tp, 1>::operator()(
    const index_array& __idx_arr) const {
  uword __dims[1] = {__idx_arr.size()};
  return Matrix(this->_M_elem, __idx_arr, __dims);
}

template <class _Tp>
inline IndirectMatrix<_Tp, 1> Matrix<_Tp, 1>::operator()(
    const index_array& __idx_arr) {
  uword __dims[1] = {__idx_arr.size()};
  return IndirectMatrix<_Tp, 1>(this->_M_elem, __idx_arr, __dims);
}

template <class _Tp>
inline Matrix<_Tp, 2> Matrix<_Tp, 2>::operator()(
    const index_array& __idx_arr1, const index_array& __idx_arr2) const {
  uword __dims[2] = {__idx_arr1.size(), __idx_arr2.size()};
  index_array __idx_arr(__dims[0] * __dims[1]);
  uword __idx = 0;
  for (uword __j = 0; __j < __idx_arr2.size(); ++__j) {
    for (uword __i = 0; __i < __idx_arr1.size(); ++__i) {
      __idx_arr[__idx++] = sub2ind(__idx_arr1[__i], __idx_arr2[__j]);
    }
  }
  return Matrix<_Tp, 2>(this->_M_elem, __idx_arr, __dims);
}

template <class _Tp>
inline IndirectMatrix<_Tp, 2> Matrix<_Tp, 2>::operator()(
    const index_array& __idx_arr1, const index_array& __idx_arr2) {
  uword __dims[2] = {__idx_arr1.size(), __idx_arr2.size()};
  index_array __idx_arr(__dims[0] * __dims[1]);
  uword __idx = 0;
  for (uword __j = 0; __j < __idx_arr2.size(); ++__j) {
    for (uword __i = 0; __i < __idx_arr1.size(); ++__i) {
      __idx_arr[__idx++] = sub2ind(__idx_arr1[__i], __idx_arr2[__j]);
    }
  }
  return IndirectMatrix<_Tp, 2>(this->_M_elem, __idx_arr, __dims);
}

template <class _Tp>
inline Matrix<_Tp, 3> Matrix<_Tp, 3>::operator()(
    const index_array& __idx_arr1, const index_array& __idx_arr2,
    const index_array& __idx_arr3) const {
  uword __dims[3] = {__idx_arr1.size(), __idx_arr2.size(), __idx_arr3.size()};
  index_array __idx_arr(__dims[0] * __dims[1] * __dims[2]);
  uword __idx = 0;
  for (uword __k = 0; __k < __idx_arr3.size(); ++__k) {
    for (uword __j = 0; __j < __idx_arr2.size(); ++__j) {
      for (uword __i = 0; __i < __idx_arr1.size(); ++__i) {
        __idx_arr[__idx++] =
            sub2ind(__idx_arr1[__i], __idx_arr2[__j], __idx_arr3[__k]);
      }
    }
  }
  return Matrix<_Tp, 3>(this->_M_elem, __idx_arr, __dims);
}

template <class _Tp>
inline IndirectMatrix<_Tp, 3> Matrix<_Tp, 3>::operator()(
    const index_array& __idx_arr1, const index_array& __idx_arr2,
    const index_array& __idx_arr3) {
  uword __dims[3] = {__idx_arr1.size(), __idx_arr2.size(), __idx_arr3.size()};
  index_array __idx_arr(__dims[0] * __dims[1] * __dims[2]);
  uword __idx = 0;
  for (uword __k = 0; __k < __idx_arr3.size(); ++__k) {
    for (uword __j = 0; __j < __idx_arr2.size(); ++__j) {
      for (uword __i = 0; __i < __idx_arr1.size(); ++__i) {
        __idx_arr[__idx++] =
            sub2ind(__idx_arr1[__i], __idx_arr2[__j], __idx_arr3[__k]);
      }
    }
  }
  return IndirectMatrix<_Tp, 3>(this->_M_elem, __idx_arr, __dims);
}

// Matrix member functions dealing with bool_array and MaskMatrix

template <class _Tp>
inline Matrix<_Tp, 1> Matrix<_Tp, 1>::operator()(
    const bool_array& __bool_arr) const {
  return Matrix<_Tp>(this->_M_elem, __bool_arr);
}

template <class _Tp>
inline MaskMatrix<_Tp> Matrix<_Tp, 1>::operator()(
    const bool_array& __bool_arr) {
  return MaskMatrix<_Tp>(this->_M_elem, __bool_arr);
}

template <class _Tp>
inline Matrix<_Tp, 1> Matrix<_Tp, 2>::operator()(
    const bool_array& __bool_arr) const {
  return Matrix<_Tp>(this->_M_elem, __bool_arr);
}

template <class _Tp>
inline MaskMatrix<_Tp> Matrix<_Tp, 2>::operator()(
    const bool_array& __bool_arr) {
  return MaskMatrix<_Tp>(this->_M_elem, __bool_arr);
}

template <class _Tp>
inline Matrix<_Tp, 1> Matrix<_Tp, 3>::operator()(
    const bool_array& __bool_arr) const {
  return Matrix<_Tp>(this->_M_elem, __bool_arr);
}

template <class _Tp>
inline MaskMatrix<_Tp> Matrix<_Tp, 3>::operator()(
    const bool_array& __bool_arr) {
  return MaskMatrix<_Tp>(this->_M_elem, __bool_arr);
}

//----------------------------------------------------------------------
// Other Matrix noninline member functions

template <class _Tp>
Matrix<_Tp, 2> Matrix<_Tp, 1>::t() const {
  return Matrix<_Tp, 2>(this->_M_elem, 1, n_rows());
}

template <class _Tp>
Matrix<_Tp, 2> Matrix<_Tp, 2>::t() const {
  uword n = n_rows(), m = n_cols();
  Matrix<_Tp, 2> __res(m, n);
  for (uword idx = 0; idx < this->n_elem(); ++idx) {
    uword i = idx / m;
    uword j = idx % m;
    __res._M_elem[idx] = (*this)[n * j + i];
  }

  return __res;
}

//----------------------------------------------------------------------
// Matrix functions

// matrix multiplication

template <typename _Tp>
Matrix<_Tp, 2> matmul(const Matrix<_Tp, 1>& __x, const Matrix<_Tp, 2>& __y) {
  __assert(__y.n_rows() == 1, "matmul(x, y): non-conformable arguments");
  const uword __n = __x.n_rows();
  const uword __m = __y.n_cols();
  Matrix<_Tp, 2> __res(__n, __m);
  for (uword __i = 0; __i != __n; ++__i)
    for (uword __j = 0; __j != __m; ++__j)
      __res(__i, __j) = __x(__i) * __y(0, __j);
  return __res;
}

template <typename _Tp>
Matrix<_Tp, 1> matmul(const Matrix<_Tp, 2>& __x, const Matrix<_Tp, 1>& __y) {
  __assert(__x.n_cols() == __y.n_rows(),
           "matmul(x, y): non-conformable arguments");
  const uword __nr = __x.n_rows();
  const uword __nc = __x.n_cols();
  Matrix<_Tp, 1> __res(__nr);
  for (uword __i = 0; __i != __nr; ++__i)
    for (uword __j = 0; __j != __nc; ++__j)
      __res(__i) += __x(__i, __j) * __y(__j);
  return __res;
}

template <typename _Tp>
Matrix<_Tp, 2> matmul(const Matrix<_Tp, 2>& __x, const Matrix<_Tp, 2>& __y) {
  __assert(__x.n_cols() == __y.n_rows(),
           "matmul(x, y): non-conformable arguments");
  const uword __nr = __x.n_rows();
  const uword __nc = __x.n_cols();
  const uword __p = __y.n_cols();
  Matrix<_Tp, 2> __res(__nr, __p);
  for (uword __i = 0; __i != __nr; ++__i)
    for (uword __j = 0; __j != __p; ++__j)
      for (uword __k = 0; __k != __nc; ++__k)
        __res(__i, __j) += __x(__i, __k) * __y(__k, __j);
  return __res;
}

//----------------------------------------------------------------------
// A list of typedefs

typedef Matrix<double, 1> vec;
typedef Matrix<double, 2> mat;
typedef Matrix<double, 3> cube;

typedef Matrix<double, 1> dvec;
typedef Matrix<double, 2> dmat;
typedef Matrix<double, 3> dcube;

typedef Matrix<float, 1> fvec;
typedef Matrix<float, 2> fmat;
typedef Matrix<float, 3> fcube;

typedef SliceMatrix<double> slice_vec;
typedef SliceMatrix<double> slice_dvec;
typedef SliceMatrix<float> slice_fvec;

typedef GsliceMatrix<double, 1> gslice_vec;
typedef GsliceMatrix<double, 2> gslice_mat;
typedef GsliceMatrix<double, 3> gslice_cube;

typedef GsliceMatrix<double, 1> gslice_dvec;
typedef GsliceMatrix<double, 2> gslice_dmat;
typedef GsliceMatrix<double, 3> gslice_dcube;

typedef GsliceMatrix<float, 1> gslice_fvec;
typedef GsliceMatrix<float, 2> gslice_fmat;
typedef GsliceMatrix<float, 3> gslice_fcube;

typedef IndirectMatrix<double, 1> indirect_vec;
typedef IndirectMatrix<double, 2> indirect_mat;
typedef IndirectMatrix<double, 3> indirect_cube;

typedef IndirectMatrix<double, 1> indirect_dvec;
typedef IndirectMatrix<double, 2> indirect_dmat;
typedef IndirectMatrix<double, 3> indirect_dcube;

typedef IndirectMatrix<float, 1> indirect_fvec;
typedef IndirectMatrix<float, 2> indirect_fmat;
typedef IndirectMatrix<float, 3> indirect_fcube;

typedef MaskMatrix<double> mask_vec;
typedef MaskMatrix<double> mask_dvec;
typedef MaskMatrix<float> mask_fvec;

}  // namespace matrix_lib

#endif
