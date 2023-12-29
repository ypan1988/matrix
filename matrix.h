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

#ifndef MATRIX_H
#define MATRIX_H

#if __cplusplus >= 201103L || (defined(_MSC_VER) && _MSC_VER >= 1900)
#define MATRIX_LIB_USE_CPP11
#endif

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <sstream>
#include <string>
#include <valarray>

#if defined(MATRIX_LIB_USE_CPP11)
#include <initializer_list>
#endif

#if defined(MATRIX_LIB_USE_R)
#include <R_ext/BLAS.h>
#include <RcppCommon.h>
#else
#include <cstdlib>  // std::abort()
#endif

namespace matrix_lib {

//-----------------------------------------------------------------------------

#if defined(MATRIX_LIB_USE_R)
inline void error(const char* p) { Rcpp::stop(p); }
#else
inline void error(const char* p) {
  std::cerr << p << std::endl;
  std::abort();
}
#endif

#ifndef NDEBUG
#define matrix_assert(expr, msg)                                          \
  do {                                                                    \
    if (!(expr)) {                                                        \
      std::stringstream ss;                                               \
      ss << std::endl                                                     \
         << __FILE__ << ":" << __LINE__ << " assertion failed: " << #expr \
         << "\n"                                                          \
         << msg << std::endl;                                             \
      error(ss.str().c_str());                                            \
    }                                                                     \
  } while (false)
#else
#define matrix_assert(expr)
#endif

#ifndef MATRIX_LIB_USE_CPP11
#define INIT_ARR1E(VAR_NAME, X1) \
  index_array VAR_NAME(1);       \
  VAR_NAME[0] = X1
#define INIT_ARR2E(VAR_NAME, X1, X2) \
  index_array VAR_NAME(2);           \
  VAR_NAME[0] = X1, VAR_NAME[1] = X2
#define INIT_ARR3E(VAR_NAME, X1, X2, X3) \
  index_array VAR_NAME(3);               \
  VAR_NAME[0] = X1, VAR_NAME[1] = X2, VAR_NAME[2] = X3
#else
#define INIT_ARR1E(VAR_NAME, X1) index_array VAR_NAME = {X1}
#define INIT_ARR2E(VAR_NAME, X1, X2) index_array VAR_NAME = {X1, X2}
#define INIT_ARR3E(VAR_NAME, X1, X2, X3) index_array VAR_NAME = {X1, X2, X3}
#endif

//-----------------------------------------------------------------------------

typedef std::size_t uword;
typedef std::valarray<bool> bool_array;
typedef std::valarray<std::size_t> index_array;

//-----------------------------------------------------------------------------

// The general Matrix template is simply a prop for its specializations:
template <class Tp = double, uword Size = 1>
struct Matrix {
  // multidimensional matrix class
  // ( ) does multidimensional subscripting with range check
  // [ ] retrieves single element or portions of the matrix
#if defined(MATRIX_LIB_USE_CPP11)
  Matrix() = delete;
#else
 private:
  Matrix();  // this should never be compiled
             //	template<class A> Matrix(A);
#endif
};

template <class Tp = double>
struct SliceMatrix;
template <class Tp = double>
struct GsliceMatrix;
template <class Tp = double>
struct MaskMatrix;
template <class Tp = double>
struct IndirectMatrix;

//-----------------------------------------------------------------------------

// Matrix_base represents the common part of the Matrix classes:
template <class Tp>
struct Matrix_base {
 protected:
  std::valarray<Tp> M_elem;

 public:
  typedef Tp elem_type;

  // clang-format off
  // construct/destroy:
  Matrix_base() : M_elem() {}
  Matrix_base(uword n) : M_elem(n) {}
  Matrix_base(const elem_type& x, uword n) : M_elem(x, n) {}
  Matrix_base(const elem_type* x, uword n) : M_elem(x, n) {}
  Matrix_base(const std::valarray<Tp>& x) : M_elem(x) {}
#if defined(MATRIX_LIB_USE_CPP11)
  Matrix_base(const Matrix_base& x) = default;
  Matrix_base(Matrix_base&& x) = default;
  Matrix_base(std::initializer_list<Tp> x) : M_elem(x) {}
#endif
  Matrix_base(const std::slice_array<Tp>   & x) : M_elem(x) {}
  Matrix_base(const std::gslice_array<Tp>  & x) : M_elem(x) {}
  Matrix_base(const std::mask_array<Tp>    & x) : M_elem(x) {}
  Matrix_base(const std::indirect_array<Tp>& x) : M_elem(x) {}
#if defined(MATRIX_LIB_USE_R)
  template <typename _U = Tp,
  typename std::enable_if<std::is_same<_U, double>::value, int>::type=0>
    Matrix_base(SEXP x) : M_elem(REAL(x), Rf_length(x)) {}
#endif

  virtual ~Matrix_base() {}

  // assignments
#if defined(MATRIX_LIB_USE_CPP11)
  Matrix_base& operator=(const Matrix_base& x) = default;
  Matrix_base& operator=(Matrix_base&& x) = default;
#endif

  // elements accessor and mutator functions
  const std::valarray<Tp>& get_elem() const { return M_elem; }
  std::valarray<Tp> get_elem(const std::slice & x) const { return M_elem[x]; }
  std::valarray<Tp> get_elem(const std::gslice& x) const { return M_elem[x]; }
  std::valarray<Tp> get_elem(const bool_array & x) const { return M_elem[x]; }
  std::valarray<Tp> get_elem(const index_array& x) const { return M_elem[x]; }

  void set_elem(const std::valarray<Tp>      & x) { M_elem.resize(x.size()); M_elem = x; }
  void set_elem(const std::slice_array<Tp>   & x, uword sz) { M_elem.resize(sz); M_elem = x; }
  void set_elem(const std::gslice_array<Tp>  & x, uword sz) { M_elem.resize(sz); M_elem = x; }
  void set_elem(const std::mask_array<Tp>    & x, uword sz) { M_elem.resize(sz); M_elem = x; }
  void set_elem(const std::indirect_array<Tp>& x, uword sz) { M_elem.resize(sz); M_elem = x; }
#if defined(MATRIX_LIB_USE_CPP11)
  void set_elem(std::initializer_list<Tp>    & x) { M_elem = x; }
  std::valarray<Tp> make_2d_array(std::initializer_list<std::initializer_list<Tp>> x, uword nr, uword nc);
  std::valarray<Tp> make_3d_array(std::initializer_list<std::initializer_list<std::initializer_list<Tp>>> x, uword nr, uword nc, uword ns);
#endif
  uword n_elem() const { return M_elem.size(); }

  // element access
  elem_type               operator[](uword              n) const { return M_elem[n]; }
  elem_type&              operator[](uword              n)       { return M_elem[n]; }

  // subsetting operations with auxiliary type
  std::valarray<Tp>       operator[](std::slice         x) const { return M_elem[x]; }
  std::slice_array<Tp>    operator[](std::slice         x)       { return M_elem[x]; }
  std::valarray<Tp>       operator[](std::gslice        x) const { return M_elem[x]; }
  std::gslice_array<Tp>   operator[](std::gslice        x)       { return M_elem[x]; }
  std::valarray<Tp>       operator[](const bool_array & x) const { return M_elem[x]; }
  std::mask_array<Tp>     operator[](const bool_array & x)       { return M_elem[x]; }
  std::valarray<Tp>       operator[](const index_array& x) const { return M_elem[x]; }
  std::indirect_array<Tp> operator[](const index_array& x)       { return M_elem[x]; }

  // if necessary, we can get to the raw matrix:
        elem_type* data()       { return &(M_elem[0]); }
  const elem_type* data() const { return &(M_elem[0]); }
  // clang-format on

 public:  // Other member functions.
          // The results are undefined for zero-length arrays
  elem_type sum() const { return M_elem.sum(); }
  elem_type min() const { return M_elem.min(); }
  elem_type max() const { return M_elem.max(); }

  template <typename F>
  void for_each(F f) {
    for (uword i = 0; i < n_elem(); ++i) f(M_elem[i]);
  }
};

#if defined(MATRIX_LIB_USE_CPP11)
template <class Tp>
std::valarray<Tp> Matrix_base<Tp>::make_2d_array(
    std::initializer_list<std::initializer_list<Tp>> x, uword nr, uword nc) {
  std::valarray<Tp> res(nr * nc);
  uword r = 0;
  for (auto riter = x.begin(); riter != x.end(); ++riter, ++r) {
    uword c = 0;
    for (auto citer = riter->begin(); citer != riter->end(); ++citer, ++c) {
      res[r + c * nr] = *citer;
    }
  }
  return res;
}
template <class Tp>
std::valarray<Tp> Matrix_base<Tp>::make_3d_array(
    std::initializer_list<std::initializer_list<std::initializer_list<Tp>>> x,
    uword nr, uword nc, uword ns) {
  std::valarray<Tp> res(nr * nc * ns);
  uword s = 0;
  for (auto siter = x.begin(); siter != x.end(); ++siter, ++s) {
    uword r = 0;
    for (auto riter = siter->begin(); riter != siter->end(); ++riter, ++r) {
      uword c = 0;
      for (auto citer = riter->begin(); citer != riter->end(); ++citer, ++c) {
        res[r + c * nr + s * nr * nc] = *citer;
      }
    }
  }
  return res;
}

#endif

//-----------------------------------------------------------------------------

template <class Tp>
struct Matrix<Tp, 1> : public Matrix_base<Tp> {
 public:
  typedef Tp elem_type;
  friend struct Matrix<Tp, 2>;
  friend struct Matrix<Tp, 3>;

  // clang-format off
  // construct/destroy:
  Matrix()                             : Matrix_base<Tp>(                  ) { M_init(); }    // (1)
  explicit Matrix(uword n1)            : Matrix_base<Tp>(                n1) { M_init(); }    // (2)
  Matrix(const elem_type& x, uword n1) : Matrix_base<Tp>(x,              n1) { M_init(); }    // (3)
  Matrix(const elem_type* x, uword n1) : Matrix_base<Tp>(x,              n1) { M_init(); }    // (4)
  Matrix(const         Matrix    &  x) : Matrix_base<Tp>(x.M_elem          ) { M_init(); }    // (5)
  Matrix(const    SliceMatrix<Tp>&  x) : Matrix_base<Tp>(x.M_elem[x.M_desc]) { M_init(x.is_column_vector); }    // (7)
  Matrix(const   GsliceMatrix<Tp>&  x) : Matrix_base<Tp>(x.M_elem[x.M_desc]) { M_init(); }    // (8)
  Matrix(const     MaskMatrix<Tp>&  x) : Matrix_base<Tp>(x.M_elem[x.M_desc]) { M_init(); }    // (9)
  Matrix(const IndirectMatrix<Tp>&  x) : Matrix_base<Tp>(x.M_elem[x.M_desc]) { M_init(); }    // (10)
#if defined(MATRIX_LIB_USE_CPP11)
  Matrix(Matrix&& x) = default;                                                               // (6)
  Matrix(std::initializer_list<Tp>  x) : Matrix_base<Tp>(x                 ) { M_init(); }    // (11)
  Matrix(Matrix<Tp, 2>&& x);
#endif
  Matrix(const std::valarray<Tp>  & x) : Matrix_base<Tp>(x                 ) { M_init(); }    // (12)

#if defined(MATRIX_LIB_USE_R)
  Matrix(SEXP                         x) : Matrix_base<Tp>(x               ) { M_init(); }
#endif
  Matrix(const std::valarray<Tp>& x, const uword        dims[1]) : Matrix_base<Tp>(x) { M_init(dims); }
  Matrix(const std::valarray<Tp>& x, const index_array& dims   ) : Matrix_base<Tp>(x) { M_init(dims); }
  Matrix(const Matrix<Tp, 2>& x);
  ~Matrix() {}

  // assignment:
  Matrix& operator=(const         Matrix    & x) { this->set_elem(x.M_elem                      ); M_init();     return *this; }    // (1)
  Matrix& operator=(const       elem_type   & x) { this->M_elem = x                              ;               return *this; }    // (3)
  Matrix& operator=(const    SliceMatrix<Tp>& x) { this->set_elem(x.M_elem[x.M_desc], x.n_elem()); M_init(x.is_column_vector); return *this; }    // (4)
  Matrix& operator=(const   GsliceMatrix<Tp>& x) { this->set_elem(x.M_elem[x.M_desc], x.n_elem()); M_init();     return *this; }    // (5)
  Matrix& operator=(const     MaskMatrix<Tp>& x) { this->set_elem(x.M_elem[x.M_desc], x.n_elem()); M_init(true); return *this; }    // (6)
  Matrix& operator=(const IndirectMatrix<Tp>& x) { this->set_elem(x.M_elem[x.M_desc], x.n_elem()); M_init();     return *this; }    // (7)
#if defined(MATRIX_LIB_USE_CPP11)
  Matrix& operator=(Matrix&& x) = default;                                                                                          // (2)
  Matrix& operator=(std::initializer_list<Tp> x) { this->M_elem = x;                               M_init();     return *this; }    // (8)
  Matrix& operator=(Matrix<Tp, 2>&& x);                                                                                             // (10)
#endif
  Matrix& operator=(const Matrix<Tp, 2>     & x);                                                                                   // (9)

#if defined(MATRIX_LIB_USE_R)
  // this operator enables implicit Rcpp::wrap
  operator SEXP()       { return export_matrix_to_sexp(); }
  operator SEXP() const { return export_matrix_to_sexp(); }
  SEXP export_matrix_to_sexp() const;
#endif

  index_array get_dims() const { return index_array(M_dims, 1); }
  uword         n_rows() const { return  is_column_vector ? M_dims[0] : (M_dims[0] ? 1 : 0); }
  uword         n_cols() const { return !is_column_vector ? M_dims[0] : (M_dims[0] ? 1 : 0); }

  void _M_range_check(uword n1) const {
    matrix_assert(n1 < this->n_elem(),
            "Matrix<T, 1>::_M_range_check:\n" <<
            "  Index " << n1 << " is out of bound for axis 0 with size " << n_rows());
  }

  // subscripting:
        elem_type& operator()(uword n1)       { _M_range_check(n1); return this->M_elem[n1]; }
  const elem_type& operator()(uword n1) const { _M_range_check(n1); return this->M_elem[n1]; }

  // SliceMatrix related member functions
          Matrix<Tp, 1> operator()(std::slice s) const;
     SliceMatrix<Tp   > operator()(std::slice s);
          Matrix<Tp, 1> subvec(uword first_index, uword last_index) const;
     SliceMatrix<Tp   > subvec(uword first_index, uword last_index);

  // MaskMatrix related member functions
          Matrix<Tp, 1> operator()(const bool_array& bool_arr) const;
      MaskMatrix<Tp   > operator()(const bool_array& bool_arr);

  // IndirectMatrix related member functions
          Matrix<Tp, 1> elem(const index_array& idx_arr) const;
  IndirectMatrix<Tp   > elem(const index_array& idx_arr);
          Matrix<Tp, 1> operator()(const index_array& idx_arr) const;
  IndirectMatrix<Tp   > operator()(const index_array& idx_arr);

  Matrix<Tp, 2> t() const;

 public:
  Matrix operator+() const { return *this; }
  Matrix operator-() const { return Matrix(-this->M_elem, M_dims); }
  Matrix operator~() const { return Matrix(~this->M_elem, M_dims); }
  Matrix<bool, 1> operator!() const { return Matrix<bool, 1>(!this->M_elem, M_dims); }

 public:
  Matrix& operator+=(const elem_type& x) { this->M_elem += x; return *this; }
  Matrix& operator-=(const elem_type& x) { this->M_elem -= x; return *this; }
  Matrix& operator*=(const elem_type& x) { this->M_elem *= x; return *this; }
  Matrix& operator/=(const elem_type& x) { this->M_elem /= x; return *this; }
  Matrix& operator%=(const elem_type& x) { this->M_elem %= x; return *this; }
  Matrix& operator&=(const elem_type& x) { this->M_elem &= x; return *this; }
  Matrix& operator|=(const elem_type& x) { this->M_elem |= x; return *this; }
  Matrix& operator^=(const elem_type& x) { this->M_elem ^= x; return *this; }
  Matrix& operator<<=(const elem_type& x) { this->M_elem <<= x; return *this; }
  Matrix& operator>>=(const elem_type& x) { this->M_elem >>= x; return *this; }

 public:
  Matrix& operator+=(const Matrix& x) { this->M_elem += x.M_elem; return *this; }
  Matrix& operator-=(const Matrix& x) { this->M_elem -= x.M_elem; return *this; }
  Matrix& operator*=(const Matrix& x) { this->M_elem *= x.M_elem; return *this; }
  Matrix& operator/=(const Matrix& x) { this->M_elem /= x.M_elem; return *this; }
  Matrix& operator%=(const Matrix& x) { this->M_elem %= x.M_elem; return *this; }
  Matrix& operator&=(const Matrix& x) { this->M_elem &= x.M_elem; return *this; }
  Matrix& operator|=(const Matrix& x) { this->M_elem |= x.M_elem; return *this; }
  Matrix& operator^=(const Matrix& x) { this->M_elem ^= x.M_elem; return *this; }
  Matrix& operator<<=(const Matrix& x) { this->M_elem <<= x.M_elem; return *this; }
  Matrix& operator>>=(const Matrix& x) { this->M_elem >>= x.M_elem; return *this; }
  // clang-format on

#if defined(MATRIX_LIB_USE_CPP11)
 public:
  bool is_column_vector = true;

 private:
  uword M_dims[1] = {0};
#else
 public:
  bool is_column_vector;

 private:
  uword M_dims[1];
#endif
 private:
  void M_init(bool is_colvec = true) {
    M_dims[0] = this->n_elem();
    is_column_vector = is_colvec;
  }
  void M_init(const uword dims[1], bool is_colvec = true) {
    if (this->n_elem() != dims[0]) error("1D Cstor error: dimension mismatch");
    M_dims[0] = dims[0];
    is_column_vector = is_colvec;
  }
  void M_init(const index_array& dims, bool is_colvec = true) {
    if (this->n_elem() != dims[0]) error("1D Cstor error: dimension mismatch");
    M_dims[0] = dims[0];
    is_column_vector = is_colvec;
  }
  Matrix(const std::valarray<Tp>& va, uword start, const uword size,
         const uword stride, bool is_colvec = true)
      : Matrix_base<Tp>(va[std::slice(start, size, stride)]) {
    M_init(is_colvec);
  }
  Matrix(const std::valarray<Tp>& va, uword start, const index_array& size,
         const index_array& stride, bool is_colvec = true)
      : Matrix_base<Tp>(va[std::gslice(start, size, stride)]) {
    M_init(is_colvec);
  }
  Matrix(const std::valarray<Tp>& va, const bool_array& bool_arr)
      : Matrix_base<Tp>(va[bool_arr]) {
    M_init();
  }
  Matrix(const std::valarray<Tp>& va, const index_array& idx_arr,
         const index_array& dims)
      : Matrix_base<Tp>(va[idx_arr]) {
    M_init();
  }
  uword sub2ind(uword i) const { return i; }
};

//-----------------------------------------------------------------------------

template <class Tp>
struct Matrix<Tp, 2> : public Matrix_base<Tp> {
 public:
  typedef Tp elem_type;
  friend struct Matrix<Tp, 1>;
  friend struct Matrix<Tp, 3>;

  // clang-format off
  // construct/destroy:
  Matrix()                                       : Matrix_base<Tp>(                  ) { M_init( 0,  0  ); }    // (1)
  Matrix(uword n1, uword n2)                     : Matrix_base<Tp>(           n1 * n2) { M_init(n1, n2  ); }    // (2)
  Matrix(const elem_type& x, uword n1, uword n2) : Matrix_base<Tp>(x,         n1 * n2) { M_init(n1, n2  ); }    // (3)
  Matrix(const elem_type* x, uword n1, uword n2) : Matrix_base<Tp>(x,         n1 * n2) { M_init(n1, n2  ); }    // (4)
  Matrix(const                    Matrix    & x) : Matrix_base<Tp>(x.M_elem          ) { M_init(x.M_dims); }    // (5)
  Matrix(const               SliceMatrix<Tp>& x) : Matrix_base<Tp>(x.M_elem[x.M_desc]) { M_init(x.is_column_vector); }    // (7)
  Matrix(const              GsliceMatrix<Tp>& x) : Matrix_base<Tp>(x.M_elem[x.M_desc]) { M_init(x.M_dims); }    // (8)
  Matrix(const                MaskMatrix<Tp>& x) : Matrix_base<Tp>(x.M_elem[x.M_desc]) { M_init(true    ); }    // (9)
  Matrix(const            IndirectMatrix<Tp>& x) : Matrix_base<Tp>(x.M_elem[x.M_desc]) { M_init(x.M_dims); }    // (10)
#if defined(MATRIX_LIB_USE_CPP11)
  Matrix(Matrix&& x) = default;                                                                                 // (6)
  Matrix(std::initializer_list<std::initializer_list<Tp>> x);                                                   // (11)
  Matrix(Matrix<Tp, 1>&& x);
#endif
  Matrix(const std::valarray<Tp>&  x, uword n1, uword n2) : Matrix_base<Tp>(x)         { M_init(n1, n2  ); }    // (12)
#if defined(MATRIX_LIB_USE_R)
  Matrix(SEXP x) : Matrix_base<Tp>(x) { SEXP dims = Rf_getAttrib(x, R_DimSymbol); M_init(INTEGER(dims)[0], INTEGER(dims)[1]); }
#endif
  Matrix(const std::valarray<Tp>& x, const uword        dims[2]) : Matrix_base<Tp>(x)  { M_init(dims)    ; }
  Matrix(const std::valarray<Tp>& x, const index_array& dims   ) : Matrix_base<Tp>(x)  { M_init(dims)    ; }
  Matrix(const Matrix<Tp, 1>& x);
  ~Matrix() {}

  // assignment:
  Matrix& operator=(const         Matrix    & x) { this->set_elem(x.M_elem                      ); M_init(x.M_dims); return *this; }    // (1)
  Matrix& operator=(const       elem_type   & x) { this->M_elem = x                              ;                   return *this; }    // (3)
  Matrix& operator=(const    SliceMatrix<Tp>& x) { this->set_elem(x.M_elem[x.M_desc], x.n_elem()); M_init(x.is_column_vector); return *this; }    // (4)
  Matrix& operator=(const   GsliceMatrix<Tp>& x) { this->set_elem(x.M_elem[x.M_desc], x.n_elem()); M_init(x.M_dims); return *this; }    // (5)
  Matrix& operator=(const     MaskMatrix<Tp>& x) { this->set_elem(x.M_elem[x.M_desc], x.n_elem()); M_init(true    ); return *this; }    // (6)
  Matrix& operator=(const IndirectMatrix<Tp>& x) { this->set_elem(x.M_elem[x.M_desc], x.n_elem()); M_init(x.M_dims); return *this; }    // (7)
#if defined(MATRIX_LIB_USE_CPP11)
  Matrix& operator=(Matrix&& x) = default;                                                                                              // (2)
  Matrix& operator=(std::initializer_list<std::initializer_list<Tp>> x);                                                                // (8)
  Matrix& operator=(Matrix<Tp, 1>          && x);                                                                                       // (10)
#endif
  Matrix& operator=(const      Matrix<Tp, 1>& x);                                                                                       // (9)

#if defined(MATRIX_LIB_USE_R)
  // this operator enables implicit Rcpp::wrap
  operator SEXP()       { return export_matrix_to_sexp(); }
  operator SEXP() const { return export_matrix_to_sexp(); }
  SEXP export_matrix_to_sexp() const;
#endif

  index_array get_dims() const { return index_array(M_dims, 2); }
  uword         n_rows() const { return M_dims[0]; }
  uword         n_cols() const { return M_dims[1]; }

  void _M_range_check(uword n1, uword n2) const {
    matrix_assert(n1 < n_rows(), "Matrix<T, 2>::_M_range_check: index is out of bound for dimension 1");
    matrix_assert(n2 < n_cols(), "Matrix<T, 2>::_M_range_check: index is out of bound for dimension 2");
  }

  // subscripting:
        elem_type& operator()(uword n1, uword n2)
  { _M_range_check(n1, n2); return this->M_elem[sub2ind(n1, n2)]; }
  const elem_type& operator()(uword n1, uword n2) const
  { _M_range_check(n1, n2); return this->M_elem[sub2ind(n1, n2)]; }

  // SliceMatrix related member functions
          Matrix<Tp, 1> row(uword row_number) const;
     SliceMatrix<Tp   > row(uword row_number);
          Matrix<Tp, 1> col(uword col_number) const;
     SliceMatrix<Tp   > col(uword col_number);
          Matrix<Tp, 1> diag(int k = 0) const;
     SliceMatrix<Tp   > diag(int k = 0);

  // GsliceMatrix related member functions
          Matrix<Tp, 2> operator()(std::slice s1, std::slice s2) const;
    GsliceMatrix<Tp   > operator()(std::slice s1, std::slice s2);
          Matrix<Tp, 2> submat(uword first_row, uword first_col, uword last_row, uword last_col) const;
    GsliceMatrix<Tp   > submat(uword first_row, uword first_col, uword last_row, uword last_col);
          Matrix<Tp, 2> rows(uword fr, uword lr) const { return submat(fr, 0, lr, n_cols() - 1); }
    GsliceMatrix<Tp   > rows(uword fr, uword lr)       { return submat(fr, 0, lr, n_cols() - 1); }
          Matrix<Tp, 2> cols(uword fc, uword lc) const { return submat(0, fc, n_rows() - 1, lc); }
    GsliceMatrix<Tp   > cols(uword fc, uword lc)       { return submat(0, fc, n_rows() - 1, lc); }

  // MaskMatrix related member functions
          Matrix<Tp, 1> operator()(const bool_array& bool_arr) const;
      MaskMatrix<Tp   > operator()(const bool_array& bool_arr);

  // IndirectMatrix related member functions
          Matrix<Tp, 1> elem(const index_array& idx_arr) const;
  IndirectMatrix<Tp   > elem(const index_array& idx_arr);
          Matrix<Tp, 2> operator()(const index_array& idx_arr1, const index_array& idx_arr2) const;
  IndirectMatrix<Tp   > operator()(const index_array& idx_arr1, const index_array& idx_arr2);

  Matrix<Tp, 2> t() const;

 public:
  Matrix operator+() const { return *this; }
  Matrix operator-() const { return Matrix(-this->M_elem, M_dims); }
  Matrix operator~() const { return Matrix(~this->M_elem, M_dims); }
  Matrix<bool, 2> operator!() const { return Matrix<bool, 2>(!this->M_elem, M_dims); }

 public:
  Matrix& operator+=(const elem_type& x) { this->M_elem += x; return *this; }
  Matrix& operator-=(const elem_type& x) { this->M_elem -= x; return *this; }
  Matrix& operator*=(const elem_type& x) { this->M_elem *= x; return *this; }
  Matrix& operator/=(const elem_type& x) { this->M_elem /= x; return *this; }
  Matrix& operator%=(const elem_type& x) { this->M_elem %= x; return *this; }
  Matrix& operator&=(const elem_type& x) { this->M_elem &= x; return *this; }
  Matrix& operator|=(const elem_type& x) { this->M_elem |= x; return *this; }
  Matrix& operator^=(const elem_type& x) { this->M_elem ^= x; return *this; }
  Matrix& operator<<=(const elem_type& x) { this->M_elem <<= x; return *this; }
  Matrix& operator>>=(const elem_type& x) { this->M_elem >>= x; return *this; }

 public:
  Matrix& operator+=(const Matrix& x) { this->M_elem += x.M_elem; return *this; }
  Matrix& operator-=(const Matrix& x) { this->M_elem -= x.M_elem; return *this; }
  Matrix& operator*=(const Matrix& x) { this->M_elem *= x.M_elem; return *this; }
  Matrix& operator/=(const Matrix& x) { this->M_elem /= x.M_elem; return *this; }
  Matrix& operator%=(const Matrix& x) { this->M_elem %= x.M_elem; return *this; }
  Matrix& operator&=(const Matrix& x) { this->M_elem &= x.M_elem; return *this; }
  Matrix& operator|=(const Matrix& x) { this->M_elem |= x.M_elem; return *this; }
  Matrix& operator^=(const Matrix& x) { this->M_elem ^= x.M_elem; return *this; }
  Matrix& operator<<=(const Matrix& x) { this->M_elem <<= x.M_elem; return *this; }
  Matrix& operator>>=(const Matrix& x) { this->M_elem >>= x.M_elem; return *this; }
  // clang-format on

 private:
#if defined(MATRIX_LIB_USE_CPP11)
  uword M_dims[2] = {0, 0};
#else
  uword M_dims[2];
#endif

  void M_init(bool is_colvec = true) {
    is_colvec ? M_init(this->n_elem(), 1) : M_init(1, this->n_elem());
  }
  void M_init(uword n1, uword n2) {
    if (this->n_elem() != n1 * n2) error("2D Cstor error: dimension mismatch");
    M_dims[0] = n1, M_dims[1] = n2;
  }
  void M_init(const uword dims[2]) { M_dims[0] = dims[0], M_dims[1] = dims[1]; }
  void M_init(const index_array& dims) {
    M_dims[0] = dims[0], M_dims[1] = dims[1];
  }
  Matrix(const std::valarray<Tp>& va, uword start, const index_array& size,
         const index_array& stride)
      : Matrix_base<Tp>(va[std::gslice(start, size, stride)]) {
    M_init(size[1], size[0]);
  }
  Matrix(const std::valarray<Tp>& va, const index_array& idx_arr,
         const index_array& dims)
      : Matrix_base<Tp>(va[idx_arr]) {
    M_init(dims);
  }
  uword sub2ind(uword i, uword j) const { return i + j * M_dims[0]; }
};

//-----------------------------------------------------------------------------

template <class Tp>
struct Matrix<Tp, 3> : public Matrix_base<Tp> {
 public:
  typedef Tp elem_type;

  // clang-format off
  // construct/destroy:
  Matrix()                                                 : Matrix_base<Tp>(                  ) { M_init( 0,  0,  0); }    // (1)
  Matrix(                    uword n1, uword n2, uword n3) : Matrix_base<Tp>(      n1 * n2 * n3) { M_init(n1, n2, n3); }    // (2)
  Matrix(const elem_type& x, uword n1, uword n2, uword n3) : Matrix_base<Tp>(x,    n1 * n2 * n3) { M_init(n1, n2, n3); }    // (3)
  Matrix(const elem_type* x, uword n1, uword n2, uword n3) : Matrix_base<Tp>(x,    n1 * n2 * n3) { M_init(n1, n2, n3); }    // (4)
  Matrix(const         Matrix    & x                     ) : Matrix_base<Tp>(x.M_elem          ) { M_init(x.M_dims  ); }    // (5)
  Matrix(const   GsliceMatrix<Tp>& x                     ) : Matrix_base<Tp>(x.M_elem[x.M_desc]) { M_init(x.M_dims  ); }    // (8)
  Matrix(const IndirectMatrix<Tp>& x                     ) : Matrix_base<Tp>(x.M_elem[x.M_desc]) { M_init(x.M_dims  ); }    // (10)
#if defined(MATRIX_LIB_USE_CPP11)
  Matrix(Matrix&& x) = default;                                                                                             // (6)
  Matrix(std::initializer_list<std::initializer_list<std::initializer_list<Tp>>> x);                                        // (11)

#endif
  Matrix(const std::valarray<Tp>&  x, uword n1, uword n2, uword n3) : Matrix_base<Tp>(x)         { M_init(n1, n2, n3); }    // (12)
  Matrix(const std::valarray<Tp>&  x, const uword          dims[3]) : Matrix_base<Tp>(x)         { M_init(dims      ); }
  Matrix(const std::valarray<Tp>&  x, const index_array&   dims   ) : Matrix_base<Tp>(x)         { M_init(dims      ); }
#if defined(MATRIX_LIB_USE_R)
  Matrix(SEXP x) : Matrix_base<Tp>(x) { SEXP dims = Rf_getAttrib(x, R_DimSymbol); M_init(INTEGER(dims)[0], INTEGER(dims)[1], INTEGER(dims)[2]); }
#endif
  ~Matrix() {}

  // assignment:
  Matrix& operator=(const         Matrix    & x) { this->set_elem(x.M_elem                      ); M_init(x.M_dims); return *this; }    // (1)
  Matrix& operator=(const       elem_type   & x) { this->M_elem = x                              ;                   return *this; }    // (3)
  Matrix& operator=(const   GsliceMatrix<Tp>& x) { this->set_elem(x.M_elem[x.M_desc], x.n_elem()); M_init(x.M_dims); return *this; }    // (5)
  Matrix& operator=(const IndirectMatrix<Tp>& x) { this->set_elem(x.M_elem[x.M_desc], x.n_elem()); M_init(x.M_dims); return *this; }    // (7)
#if defined(MATRIX_LIB_USE_CPP11)
  Matrix& operator=(Matrix&& x) = default;                                                                                              // (2)
  Matrix& operator=(std::initializer_list<std::initializer_list<std::initializer_list<Tp>>> x);                                         // (8)
#endif

#if defined(MATRIX_LIB_USE_R)
  // this operator enables implicit Rcpp::wrap
  operator SEXP()       { return export_matrix_to_sexp(); }
  operator SEXP() const { return export_matrix_to_sexp(); }
  SEXP export_matrix_to_sexp() const;
#endif

  index_array get_dims() const { return index_array(M_dims, 3); }
  uword         n_rows() const { return M_dims[0]; }
  uword         n_cols() const { return M_dims[1]; }
  uword       n_slices() const { return M_dims[2]; }

  void _M_range_check(uword n1, uword n2, uword n3) const {
    matrix_assert(n1 < n_rows()  , "Matrix<T, 3>::_M_range_check: index is out of bound for dimension 1");
    matrix_assert(n2 < n_cols()  , "Matrix<T, 3>::_M_range_check: index is out of bound for dimension 2");
    matrix_assert(n3 < n_slices(), "Matrix<T, 3>::_M_range_check: index is out of bound for dimension 3");
  }

  // subscripting:
        elem_type& operator()(uword n1, uword n2, uword n3)
  { _M_range_check(n1, n2, n3); return this->M_elem[sub2ind(n1, n2, n3)]; }
  const elem_type& operator()(uword n1, uword n2, uword n3) const
  { _M_range_check(n1, n2, n3); return this->M_elem[sub2ind(n1, n2, n3)]; }

  // GsliceMatrix related member functions
          Matrix<Tp, 3> operator()(std::slice s1, std::slice s2, std::slice s3) const;
    GsliceMatrix<Tp   > operator()(std::slice s1, std::slice s2, std::slice s3);
          Matrix<Tp, 3> subcube(uword first_row, uword first_col, uword first_slice,
                                uword  last_row, uword  last_col, uword  last_slice) const;
    GsliceMatrix<Tp   > subcube(uword first_row, uword first_col, uword first_slice,
                      uword  last_row, uword  last_col, uword  last_slice);
          Matrix<Tp, 2> row(uword slice_number) const;
    GsliceMatrix<Tp   > row(uword slice_number);
          Matrix<Tp, 2> col(uword slice_number) const;
    GsliceMatrix<Tp   > col(uword slice_number);
          Matrix<Tp, 2> slice(uword slice_number) const;
    GsliceMatrix<Tp   > slice(uword slice_number);
          Matrix<Tp, 3> slices(uword fs, uword ls) const { return subcube(0, 0, fs, n_rows() - 1, n_cols() - 1, ls); }
    GsliceMatrix<Tp   > slices(uword fs, uword ls)       { return subcube(0, 0, fs, n_rows() - 1, n_cols() - 1, ls); }

  // MaskMatrix related member functions
          Matrix<Tp, 1> operator()(const bool_array& bool_arr) const;
      MaskMatrix<Tp   > operator()(const bool_array& bool_arr);

  // IndirectMatrix related member functions
          Matrix<Tp, 1> elem(const index_array& idx_arr) const;
  IndirectMatrix<Tp   > elem(const index_array& idx_arr);
          Matrix<Tp, 3> operator()(const index_array&, const index_array&, const index_array&) const;
  IndirectMatrix<Tp   > operator()(const index_array&, const index_array&, const index_array&);

 public:
  Matrix operator+() const { return *this; }
  Matrix operator-() const { return Matrix(-this->M_elem, M_dims); }
  Matrix operator~() const { return Matrix(~this->M_elem, M_dims); }
  Matrix<bool, 3> operator!() const { return Matrix<bool, 3>(!this->M_elem, M_dims); }

 public:
  Matrix& operator+=(const elem_type& x) { this->M_elem += x; return *this; }
  Matrix& operator-=(const elem_type& x) { this->M_elem -= x; return *this; }
  Matrix& operator*=(const elem_type& x) { this->M_elem *= x; return *this; }
  Matrix& operator/=(const elem_type& x) { this->M_elem /= x; return *this; }
  Matrix& operator%=(const elem_type& x) { this->M_elem %= x; return *this; }
  Matrix& operator&=(const elem_type& x) { this->M_elem &= x; return *this; }
  Matrix& operator|=(const elem_type& x) { this->M_elem |= x; return *this; }
  Matrix& operator^=(const elem_type& x) { this->M_elem ^= x; return *this; }
  Matrix& operator<<=(const elem_type& x) { this->M_elem <<= x; return *this; }
  Matrix& operator>>=(const elem_type& x) { this->M_elem >>= x; return *this; }

 public:
  Matrix& operator+=(const Matrix& x) { this->M_elem += x.M_elem; return *this; }
  Matrix& operator-=(const Matrix& x) { this->M_elem -= x.M_elem; return *this; }
  Matrix& operator*=(const Matrix& x) { this->M_elem *= x.M_elem; return *this; }
  Matrix& operator/=(const Matrix& x) { this->M_elem /= x.M_elem; return *this; }
  Matrix& operator%=(const Matrix& x) { this->M_elem %= x.M_elem; return *this; }
  Matrix& operator&=(const Matrix& x) { this->M_elem &= x.M_elem; return *this; }
  Matrix& operator|=(const Matrix& x) { this->M_elem |= x.M_elem; return *this; }
  Matrix& operator^=(const Matrix& x) { this->M_elem ^= x.M_elem; return *this; }
  Matrix& operator<<=(const Matrix& x) { this->M_elem <<= x.M_elem; return *this; }
  Matrix& operator>>=(const Matrix& x) { this->M_elem >>= x.M_elem; return *this; }
  // clang-format on

 private:
#if defined(MATRIX_LIB_USE_CPP11)
  uword M_dims[3] = {0, 0, 0}, _M_d1xd2 = 0;
#else
  uword M_dims[3], _M_d1xd2;
#endif
  void M_init(uword n1, uword n2, uword n3) {
    if (this->n_elem() != n1 * n2 * n3)
      error("3D Cstor error: dimension mismatch");
    M_dims[0] = n1, M_dims[1] = n2, M_dims[2] = n3;
    _M_d1xd2 = n1 * n2;
  }
  void M_init(const uword dims[3]) {
    M_dims[0] = dims[0], M_dims[1] = dims[1], M_dims[2] = dims[2];
    _M_d1xd2 = dims[0] * dims[1];
  }
  void M_init(const index_array& dims) {
    M_dims[0] = dims[0], M_dims[1] = dims[1], M_dims[2] = dims[2];
    _M_d1xd2 = dims[0] * dims[1];
  }
  Matrix(const std::valarray<Tp>& va, uword start, const index_array& size,
         const index_array& stride)
      : Matrix_base<Tp>(va[std::gslice(start, size, stride)]) {
    M_init(size[2], size[1], size[0]);
  }
  Matrix(const std::valarray<Tp>& va, const index_array& idx_arr,
         const index_array& dims)
      : Matrix_base<Tp>(va[idx_arr]) {
    M_init(dims);
  }
  uword sub2ind(uword i, uword j, uword k) const {
    return i + j * M_dims[0] + k * _M_d1xd2;
  }
};

//-----------------------------------------------------------------------------

#if defined(MATRIX_LIB_USE_CPP11)
template <class Tp>
Matrix<Tp, 2>::Matrix(std::initializer_list<std::initializer_list<Tp>> x)
    : Matrix_base<Tp>() {
  // init number of rows and columns
  uword nr = x.size();
  uword nc = x.begin()->size();

  // check dimensions
  for (auto iter = x.begin(); iter != x.end(); ++iter)
    if (iter->size() != nc) error("dimension mismatch");

  // init elements and dimensions
  this->M_elem = this->make_2d_array(x, nr, nc);
  M_init(nr, nc);
}

template <class Tp>
Matrix<Tp, 2>& Matrix<Tp, 2>::operator=(
    std::initializer_list<std::initializer_list<Tp>> x) {
  // init number of rows and columns
  uword nr = x.size();
  uword nc = x.begin()->size();

  // check dimensions
  for (auto iter = x.begin(); iter != x.end(); ++iter)
    if (iter->size() != nc) error("dimension mismatch");

  // init elements and dimensions
  this->M_elem = this->make_2d_array(x, nr, nc);
  M_init(nr, nc);

  return *this;
}

template <class Tp>
Matrix<Tp, 3>::Matrix(
    std::initializer_list<std::initializer_list<std::initializer_list<Tp>>> x)
    : Matrix_base<Tp>() {
  // init number of rows, columns and slides
  uword nr = x.begin()->size();
  uword nc = x.begin()->begin()->size();
  uword ns = x.size();

  // check dimensions
  for (auto siter = x.begin(); siter != x.end(); ++siter) {
    if (siter->size() != nr) error("dimension mismatch");
    for (auto riter = siter->begin(); riter != siter->end(); ++riter) {
      if (riter->size() != nc) error("dimension mismatch");
    }
  }

  // init elements and dimensions
  this->M_elem = this->make_3d_array(x, nr, nc, ns);
  M_init(nr, nc, ns);
}

template <class Tp>
Matrix<Tp, 3>& Matrix<Tp, 3>::operator=(
    std::initializer_list<std::initializer_list<std::initializer_list<Tp>>> x) {
  // init number of rows, columns and slides
  uword nr = x.begin()->size();
  uword nc = x.begin()->begin()->size();
  uword ns = x.size();

  // check dimensions
  for (auto siter = x.begin(); siter != x.end(); ++siter) {
    if (siter->size() != nr) error("dimension mismatch");
    for (auto riter = siter->begin(); riter != siter->end(); ++riter) {
      if (riter->size() != nc) error("dimension mismatch");
    }
  }

  // init elements and dimensions
  this->M_elem = this->make_3d_array(x, nr, nc, ns);
  M_init(nr, nc, ns);

  return *this;
}
#endif

template <class Tp>
Matrix<Tp, 1>::Matrix(const Matrix<Tp, 2>& x) : Matrix_base<Tp>(x.get_elem()) {
  if (x.n_cols() != 1 && x.n_rows() != 1) error("x is not a col/row vector");
  M_init(x.n_cols() == 1);
}

template <class Tp>
Matrix<Tp, 2>::Matrix(const Matrix<Tp, 1>& x) : Matrix_base<Tp>(x.get_elem()) {
  M_init(x.is_column_vector);
}

template <class Tp>
Matrix<Tp, 1>& Matrix<Tp, 1>::operator=(const Matrix<Tp, 2>& x) {
  if (x.n_cols() != 1 && x.n_rows() != 1) error("x is not a col/row vector");
  this->set_elem(x.get_elem());
  M_init(x.n_cols() == 1);
  return *this;
}

template <class Tp>
Matrix<Tp, 2>& Matrix<Tp, 2>::operator=(const Matrix<Tp, 1>& x) {
  this->set_elem(x.get_elem());
  M_init(x.is_column_vector);
  return *this;
}

#if defined(MATRIX_LIB_USE_CPP11)
template <class Tp>
Matrix<Tp, 1>::Matrix(Matrix<Tp, 2>&& x)
    : Matrix_base<Tp>(std::move(x.M_elem)) {
  if (x.n_cols() != 1 && x.n_rows() != 1) error("x is not a col/row matrix");
  M_init(x.n_cols() == 1);
}

template <class Tp>
Matrix<Tp, 2>::Matrix(Matrix<Tp, 1>&& x)
    : Matrix_base<Tp>(std::move(x.M_elem)) {
  M_init(x.is_column_vector);
}

template <class Tp>
Matrix<Tp, 1>& Matrix<Tp, 1>::operator=(Matrix<Tp, 2>&& x) {
  if (x.n_cols() != 1 && x.n_rows() != 1) error("x is not a col/row matrix");
  this->M_elem = std::move(x.M_elem);
  M_init(x.n_cols() == 1);
  return *this;
}

template <class Tp>
Matrix<Tp, 2>& Matrix<Tp, 2>::operator=(Matrix<Tp, 1>&& x) {
  this->M_elem = std::move(x.M_elem);
  M_init(x.is_column_vector);
  return *this;
}
#endif

// Matrix Rcpp-related member function

#if defined(MATRIX_LIB_USE_R)
template <class Tp>
SEXP Matrix<Tp, 1>::export_matrix_to_sexp() const {
  SEXP x = PROTECT(Rcpp::wrap(this->data(), this->data() + this->n_elem()));
  SEXP dims = PROTECT(Rf_allocVector(INTSXP, 1));
  INTEGER(dims)[0] = n_elem();
  Rf_setAttrib(x, R_DimSymbol, dims);
  UNPROTECT(2);
  return x;
}

template <class Tp>
SEXP Matrix<Tp, 2>::export_matrix_to_sexp() const {
  SEXP x = PROTECT(Rcpp::wrap(this->data(), this->data() + this->n_elem()));
  SEXP dims = PROTECT(Rf_allocVector(INTSXP, 2));
  INTEGER(dims)[0] = n_rows();
  INTEGER(dims)[1] = n_cols();
  Rf_setAttrib(x, R_DimSymbol, dims);
  UNPROTECT(2);
  return x;
}

template <class Tp>
SEXP Matrix<Tp, 3>::export_matrix_to_sexp() const {
  SEXP x = PROTECT(Rcpp::wrap(this->data(), this->data() + this->n_elem()));
  SEXP dims = PROTECT(Rf_allocVector(INTSXP, 3));
  INTEGER(dims)[0] = n_rows();
  INTEGER(dims)[1] = n_cols();
  INTEGER(dims)[2] = n_slices();
  Rf_setAttrib(x, R_DimSymbol, dims);
  UNPROTECT(2);
  return x;
}
#endif

//----------------------------------------------------------------------
// Matrix non-member functions.

bool same_dims(const std::valarray<uword>& d1, const std::valarray<uword>& d2) {
  uword n = d1.size();
  for (uword i = 0; i != n; ++i)
    if (d1[i] != d2[i]) return false;
  return true;
}

// Binary arithmetic operations between two Matrix.

template <class Tp, uword Size>
inline Matrix<Tp, Size> operator+(const Matrix<Tp, Size>& x,
                                  const Matrix<Tp, Size>& y) {
  matrix_assert(same_dims(x.get_dims(), y.get_dims()), "dimension mismatch");
  Matrix<Tp, Size> tmp(x);
  return tmp += y;
}

template <class Tp, uword Size>
inline Matrix<Tp, Size> operator-(const Matrix<Tp, Size>& x,
                                  const Matrix<Tp, Size>& y) {
  matrix_assert(same_dims(x.get_dims(), y.get_dims()), "dimension mismatch");
  Matrix<Tp, Size> tmp(x);
  return tmp -= y;
}

template <class Tp, uword Size>
inline Matrix<Tp, Size> operator*(const Matrix<Tp, Size>& x,
                                  const Matrix<Tp, Size>& y) {
  matrix_assert(same_dims(x.get_dims(), y.get_dims()), "dimension mismatch");
  Matrix<Tp, Size> tmp(x);
  return tmp *= y;
}

template <class Tp, uword Size>
inline Matrix<Tp, Size> operator/(const Matrix<Tp, Size>& x,
                                  const Matrix<Tp, Size>& y) {
  matrix_assert(same_dims(x.get_dims(), y.get_dims()), "dimension mismatch");
  Matrix<Tp, Size> tmp(x);
  return tmp /= y;
}

template <class Tp, uword Size>
inline Matrix<Tp, Size> operator%(const Matrix<Tp, Size>& x,
                                  const Matrix<Tp, Size>& y) {
  matrix_assert(same_dims(x.get_dims(), y.get_dims()), "dimension mismatch");
  Matrix<Tp, Size> tmp(x);
  return tmp %= y;
}

template <class Tp, uword Size>
inline Matrix<Tp, Size> operator&(const Matrix<Tp, Size>& x,
                                  const Matrix<Tp, Size>& y) {
  matrix_assert(same_dims(x.get_dims(), y.get_dims()), "dimension mismatch");
  Matrix<Tp, Size> tmp(x);
  return tmp &= y;
}

template <class Tp, uword Size>
inline Matrix<Tp, Size> operator|(const Matrix<Tp, Size>& x,
                                  const Matrix<Tp, Size>& y) {
  matrix_assert(same_dims(x.get_dims(), y.get_dims()), "dimension mismatch");
  Matrix<Tp, Size> tmp(x);
  return tmp |= y;
}

template <class Tp, uword Size>
inline Matrix<Tp, Size> operator^(const Matrix<Tp, Size>& x,
                                  const Matrix<Tp, Size>& y) {
  matrix_assert(same_dims(x.get_dims(), y.get_dims()), "dimension mismatch");
  Matrix<Tp, Size> tmp(x);
  return tmp ^= y;
}

template <class Tp, uword Size>
inline Matrix<Tp, Size> operator<<(const Matrix<Tp, Size>& x,
                                   const Matrix<Tp, Size>& y) {
  matrix_assert(same_dims(x.get_dims(), y.get_dims()), "dimension mismatch");
  Matrix<Tp, Size> tmp(x);
  return tmp <<= y;
}

template <class Tp, uword Size>
inline Matrix<Tp, Size> operator>>(const Matrix<Tp, Size>& x,
                                   const Matrix<Tp, Size>& y) {
  matrix_assert(same_dims(x.get_dims(), y.get_dims()), "dimension mismatch");
  Matrix<Tp, Size> tmp(x);
  return tmp >>= y;
}

template <class Tp, uword Size>
inline bool_array operator&&(const Matrix<Tp, Size>& x,
                             const Matrix<Tp, Size>& y) {
  matrix_assert(same_dims(x.get_dims(), y.get_dims()), "dimension mismatch");
  return x.get_elem() && y.get_elem();
}

template <class Tp, uword Size>
inline bool_array operator||(const Matrix<Tp, Size>& x,
                             const Matrix<Tp, Size>& y) {
  matrix_assert(same_dims(x.get_dims(), y.get_dims()), "dimension mismatch");
  return x.get_elem() || y.get_elem();
}

// Binary arithmetic operations between an array and a scalar.

template <class Tp, uword Size>
inline Matrix<Tp, Size> operator+(const Matrix<Tp, Size>& x, const Tp& c) {
  Matrix<Tp, Size> tmp(x);
  return tmp += c;
}

template <class Tp, uword Size>
inline Matrix<Tp, Size> operator+(const Tp& c, const Matrix<Tp, Size>& x) {
  return Matrix<Tp, Size>(c + x.get_elem(), x.get_dims());
}

template <class Tp, uword Size>
inline Matrix<Tp, Size> operator-(const Matrix<Tp, Size>& x, const Tp& c) {
  Matrix<Tp, Size> tmp(x);
  return tmp -= c;
}

template <class Tp, uword Size>
inline Matrix<Tp, Size> operator-(const Tp& c, const Matrix<Tp, Size>& x) {
  return Matrix<Tp, Size>(c - x.get_elem(), x.get_dims());
}

template <class Tp, uword Size>
inline Matrix<Tp, Size> operator*(const Matrix<Tp, Size>& x, const Tp& c) {
  Matrix<Tp, Size> tmp(x);
  return tmp *= c;
}

template <class Tp, uword Size>
inline Matrix<Tp, Size> operator*(const Tp& c, const Matrix<Tp, Size>& x) {
  return Matrix<Tp, Size>(c * x.get_elem(), x.get_dims());
}

template <class Tp, uword Size>
inline Matrix<Tp, Size> operator/(const Matrix<Tp, Size>& x, const Tp& c) {
  Matrix<Tp, Size> tmp(x);
  return tmp /= c;
}

template <class Tp, uword Size>
inline Matrix<Tp, Size> operator/(const Tp& c, const Matrix<Tp, Size>& x) {
  Matrix<Tp, Size> tmp(c / x.get_elem(), x.get_dims());
  return tmp;
}

template <class Tp, uword Size>
inline Matrix<Tp, Size> operator%(const Matrix<Tp, Size>& x, const Tp& c) {
  Matrix<Tp, Size> tmp(x);
  return tmp %= c;
}

template <class Tp, uword Size>
inline Matrix<Tp, Size> operator%(const Tp& c, const Matrix<Tp, Size>& x) {
  Matrix<Tp, Size> tmp(c % x.get_elem(), x.get_dims());
  return tmp;
}

template <class Tp, uword Size>
inline Matrix<Tp, Size> operator&(const Matrix<Tp, Size>& x, const Tp& c) {
  Matrix<Tp, Size> tmp(x);
  return tmp &= c;
}

template <class Tp, uword Size>
inline Matrix<Tp, Size> operator&(const Tp& c, const Matrix<Tp, Size>& x) {
  Matrix<Tp, Size> tmp(c & x.get_elem(), x.get_dims());
  return tmp;
}

template <class Tp, uword Size>
inline Matrix<Tp, Size> operator|(const Matrix<Tp, Size>& x, const Tp& c) {
  Matrix<Tp, Size> tmp(x);
  return tmp |= c;
}

template <class Tp, uword Size>
inline Matrix<Tp, Size> operator|(const Tp& c, const Matrix<Tp, Size>& x) {
  Matrix<Tp, Size> tmp(c | x.get_elem(), x.get_dims());
  return tmp;
}

template <class Tp, uword Size>
inline Matrix<Tp, Size> operator^(const Matrix<Tp, Size>& x, const Tp& c) {
  Matrix<Tp, Size> tmp(x);
  return tmp ^= c;
}

template <class Tp, uword Size>
inline Matrix<Tp, Size> operator^(const Tp& c, const Matrix<Tp, Size>& x) {
  Matrix<Tp, Size> tmp(c ^ x.get_elem(), x.get_dims());
  return tmp;
}

template <class Tp, uword Size>
inline Matrix<Tp, Size> operator<<(const Matrix<Tp, Size>& x, const Tp& c) {
  Matrix<Tp, Size> tmp(x);
  return tmp <<= c;
}

template <class Tp, uword Size>
inline Matrix<Tp, Size> operator<<(const Tp& c, const Matrix<Tp, Size>& x) {
  Matrix<Tp, Size> tmp(c << x.get_elem(), x.get_dims());
  return tmp;
}

template <class Tp, uword Size>
inline Matrix<Tp, Size> operator>>(const Matrix<Tp, Size>& x, const Tp& c) {
  Matrix<Tp, Size> tmp(x);
  return tmp >>= c;
}

template <class Tp, uword Size>
inline Matrix<Tp, Size> operator>>(const Tp& c, const Matrix<Tp, Size>& x) {
  Matrix<Tp, Size> tmp(c >> x.get_elem(), x.get_dims());
  return tmp;
}

template <class Tp, uword Size>
inline bool_array operator&&(const Matrix<Tp, Size>& x, const Tp& c) {
  return x.get_elem() && c;
}

template <class Tp, uword Size>
inline bool_array operator&&(const Tp& c, const Matrix<Tp, Size>& x) {
  return c && x.get_elem();
}

template <class Tp, uword Size>
inline bool_array operator||(const Matrix<Tp, Size>& x, const Tp& c) {
  return x.get_elem() || c;
}

template <class Tp, uword Size>
inline bool_array operator||(const Tp& c, const Matrix<Tp, Size>& x) {
  return c || x.get_elem();
}

// Binary logical operations between two Matrices.

template <class Tp, uword Size>
inline bool_array operator==(const Matrix<Tp, Size>& x,
                             const Matrix<Tp, Size>& y) {
  matrix_assert(same_dims(x.get_dims(), y.get_dims()), "dimension mismatch");
  return x.get_elem() == y.get_elem();
}

template <class Tp, uword Size>
inline bool_array operator!=(const Matrix<Tp, Size>& x,
                             const Matrix<Tp, Size>& y) {
  matrix_assert(same_dims(x.get_dims(), y.get_dims()), "dimension mismatch");
  return x.get_elem() != y.get_elem();
}

template <class Tp, uword Size>
inline bool_array operator<(const Matrix<Tp, Size>& x,
                            const Matrix<Tp, Size>& y) {
  matrix_assert(same_dims(x.get_dims(), y.get_dims()), "dimension mismatch");
  return x.get_elem() < y.get_elem();
}

template <class Tp, uword Size>
inline bool_array operator<=(const Matrix<Tp, Size>& x,
                             const Matrix<Tp, Size>& y) {
  matrix_assert(same_dims(x.get_dims(), y.get_dims()), "dimension mismatch");
  return x.get_elem() <= y.get_elem();
}

template <class Tp, uword Size>
inline bool_array operator>(const Matrix<Tp, Size>& x,
                            const Matrix<Tp, Size>& y) {
  matrix_assert(same_dims(x.get_dims(), y.get_dims()), "dimension mismatch");
  return x.get_elem() > y.get_elem();
}

template <class Tp, uword Size>
inline bool_array operator>=(const Matrix<Tp, Size>& x,
                             const Matrix<Tp, Size>& y) {
  matrix_assert(same_dims(x.get_dims(), y.get_dims()), "dimension mismatch");
  return x.get_elem() >= y.get_elem();
}

// Logical operations between a Matrix and a scalar.

template <class Tp, uword Size>
inline bool_array operator==(const Matrix<Tp, Size>& x, const Tp& c) {
  return x.get_elem() == c;
}

template <class Tp, uword Size>
inline bool_array operator==(const Tp& c, const Matrix<Tp, Size>& x) {
  return c == x.get_elem();
}

template <class Tp, uword Size>
inline bool_array operator!=(const Matrix<Tp, Size>& x, const Tp& c) {
  return x.get_elem() != c;
}

template <class Tp, uword Size>
inline bool_array operator!=(const Tp& c, const Matrix<Tp, Size>& x) {
  return c != x.get_elem();
}

template <class Tp, uword Size>
inline bool_array operator<(const Matrix<Tp, Size>& x, const Tp& c) {
  return x.get_elem() < c;
}

template <class Tp, uword Size>
inline bool_array operator<(const Tp& c, const Matrix<Tp, Size>& x) {
  return c < x.get_elem();
}

template <class Tp, uword Size>
inline bool_array operator<=(const Matrix<Tp, Size>& x, const Tp& c) {
  return x.get_elem() <= c;
}

template <class Tp, uword Size>
inline bool_array operator<=(const Tp& c, const Matrix<Tp, Size>& x) {
  return c <= x.get_elem();
}

template <class Tp, uword Size>
inline bool_array operator>(const Matrix<Tp, Size>& x, const Tp& c) {
  return x.get_elem() > c;
}

template <class Tp, uword Size>
inline bool_array operator>(const Tp& c, const Matrix<Tp, Size>& x) {
  return c > x.get_elem();
}

template <class Tp, uword Size>
inline bool_array operator>=(const Matrix<Tp, Size>& x, const Tp& c) {
  return x.get_elem() >= c;
}

template <class Tp, uword Size>
inline bool_array operator>=(const Tp& c, const Matrix<Tp, Size>& x) {
  return c >= x.get_elem();
}

// Matrix "transcendentals" (the list includes abs and sqrt, which,
// of course, are not transcendental).

template <class Tp, uword Size>
inline Matrix<Tp, Size> abs(const Matrix<Tp, Size>& x) {
  Matrix<Tp, Size> tmp(std::abs(x.get_elem()), x.get_dims());
  return tmp;
}

template <class Tp, uword Size>
inline Matrix<Tp, Size> exp(const Matrix<Tp, Size>& x) {
  Matrix<Tp, Size> tmp(std::exp(x.get_elem()), x.get_dims());
  return tmp;
}

template <class Tp, uword Size>
inline Matrix<Tp, Size> log(const Matrix<Tp, Size>& x) {
  Matrix<Tp, Size> tmp(std::log(x.get_elem()), x.get_dims());
  return tmp;
}

template <class Tp, uword Size>
inline Matrix<Tp, Size> log10(const Matrix<Tp, Size>& x) {
  Matrix<Tp, Size> tmp(std::log10(x.get_elem()), x.get_dims());
  return tmp;
}

template <class Tp>
Tp dot(const Matrix<Tp, 1>& x, const Matrix<Tp, 1>& y) {
  return (x.get_elem() * y.get_elem()).sum();
}

template <class Tp, class Tpx>
struct SubMatrix_base {
 public:
  typedef Tp elem_type;
  std::valarray<Tp>& M_elem;
  SubMatrix_base(std::valarray<Tp>& va) : M_elem(va) {}

  virtual uword n_elem() const = 0;
  virtual Tpx& desc() = 0;
  virtual const Tpx& desc() const = 0;

  std::valarray<Tp> get_elem() const {
    return std::valarray<Tp>(M_elem[this->desc()]);
  }

  elem_type sum() const { return get_elem().sum(); }
  elem_type min() const { return get_elem().min(); }
  elem_type max() const { return get_elem().max(); }
};

//----------------------------------------------------------------------
// SliceMatrix

template <class Tp>
struct SliceMatrix : public SubMatrix_base<Tp, std::slice> {
 public:
  std::slice M_desc;
  bool is_column_vector;
  SliceMatrix(std::valarray<Tp>& va, uword start, uword size, uword stride,
              bool is_colvec = true)
      : SubMatrix_base<Tp, std::slice>(va),
        M_desc(start, size, stride),
        is_column_vector(is_colvec) {}
  uword n_elem() const { return M_desc.size(); }
  std::slice& desc() { return M_desc; }
  const std::slice& desc() const { return M_desc; }

  void operator=(const Tp& value) { this->M_elem[M_desc] = value; }
};

//----------------------------------------------------------------------
// GsliceMatrix

template <class Tp>
struct GsliceMatrix : public SubMatrix_base<Tp, std::gslice> {
 public:
  std::gslice M_desc;
  uword M_order;
  index_array M_dims;
  uword M_size;
  GsliceMatrix(std::valarray<Tp>& va, uword start, const index_array& size,
               const index_array& stride)
      : SubMatrix_base<Tp, std::gslice>(va),
        M_desc(start, size, stride),
        M_order(size.size()),
        M_dims(M_order),
        M_size(1) {
    for (uword idx = 0; idx < M_order; ++idx) {
      M_dims[idx] = size[M_order - idx - 1];
      M_size *= M_dims[idx];
    }
  }
  uword n_elem() const { return M_size; }
  std::gslice& desc() { return M_desc; }
  const std::gslice& desc() const { return M_desc; }

  void operator=(const Tp& value) { this->M_elem[M_desc] = value; }
};

//----------------------------------------------------------------------
// MaskMatrix

template <class Tp>
struct MaskMatrix : public SubMatrix_base<Tp, bool_array> {
 public:
  bool_array M_desc;
  uword M_size;
  MaskMatrix(std::valarray<Tp>& va, const bool_array& boolarr)
      : SubMatrix_base<Tp, bool_array>(va), M_desc(boolarr), M_size(0) {
    std::valarray<uword> tmp((uword)0, boolarr.size());
    tmp[boolarr] = 1;
    M_size = tmp.sum();
  }
  uword n_elem() const { return M_size; }
  bool_array& desc() { return M_desc; }
  const bool_array& desc() const { return M_desc; }

  void operator=(const Tp& value) { this->M_elem[M_desc] = value; }
};

//----------------------------------------------------------------------
// IndirectMatrix

template <class Tp>
struct IndirectMatrix : public SubMatrix_base<Tp, index_array> {
 public:
  index_array M_desc;
  index_array M_dims;
  uword M_order;
  IndirectMatrix(std::valarray<Tp>& va, const index_array& ind_arr,
                 const index_array& dims)
      : SubMatrix_base<Tp, index_array>(va),
        M_desc(ind_arr),
        M_dims(dims),
        M_order(dims.size()) {}
  uword n_elem() const { return M_desc.size(); }
  index_array& desc() { return M_desc; }
  const index_array& desc() const { return M_desc; }

  void operator=(const Tp& value) { this->M_elem[M_desc] = value; }
};

// Matrix member functions dealing with SliceMatrix

template <class Tp>
inline Matrix<Tp, 1> Matrix<Tp, 1>::operator()(std::slice s) const {
  const uword start = s.start();
  const uword size = s.size();
  const uword stride = s.stride();
  return Matrix<Tp, 1>(this->M_elem, start, size, stride, is_column_vector);
}

template <class Tp>
inline SliceMatrix<Tp> Matrix<Tp, 1>::operator()(std::slice s) {
  const uword start = s.start();
  const uword size = s.size();
  const uword stride = s.stride();
  return SliceMatrix<Tp>(this->M_elem, start, size, stride, is_column_vector);
}

template <class Tp>
inline Matrix<Tp, 1> Matrix<Tp, 1>::subvec(uword i, uword j) const {
  if (i > j || j >= M_dims[0]) error("1D subscription error");
  const uword start = i;
  const uword size = j - i + 1;
  const uword stride = 1;
  return Matrix<Tp, 1>(this->M_elem, start, size, stride, is_column_vector);
}

template <class Tp>
inline SliceMatrix<Tp> Matrix<Tp, 1>::subvec(uword i, uword j) {
  if (i > j || j >= M_dims[0]) error("1D subscription error");
  const uword start = i;
  const uword size = j - i + 1;
  const uword stride = 1;
  return SliceMatrix<Tp>(this->M_elem, start, size, stride, is_column_vector);
}

template <class Tp>
inline Matrix<Tp, 1> Matrix<Tp, 2>::row(uword r) const {
  _M_range_check(r, 0);
  const uword start = r;
  const uword size = M_dims[1];
  const uword stride = M_dims[0];
  return Matrix<Tp, 1>(this->M_elem, start, size, stride, false);
}

template <class Tp>
inline SliceMatrix<Tp> Matrix<Tp, 2>::row(uword r) {
  _M_range_check(r, 0);
  const uword start = r;
  const uword size = M_dims[1];
  const uword stride = M_dims[0];
  return SliceMatrix<Tp>(this->M_elem, start, size, stride, false);
}

template <class Tp>
inline Matrix<Tp, 1> Matrix<Tp, 2>::col(uword c) const {
  _M_range_check(0, c);
  const uword start = c * M_dims[0];
  const uword size = M_dims[0];
  const uword stride = 1;
  return Matrix<Tp, 1>(this->M_elem, start, size, stride);
}

template <class Tp>
inline SliceMatrix<Tp> Matrix<Tp, 2>::col(uword c) {
  _M_range_check(0, c);
  const uword start = c * M_dims[0];
  const uword size = M_dims[0];
  const uword stride = 1;
  return SliceMatrix<Tp>(this->M_elem, start, size, stride);
}

template <class Tp>
inline Matrix<Tp, 1> Matrix<Tp, 2>::diag(int k) const {
  uword start = 0, size = 0;
  if (k == 0) {
    start = 0, size = std::min(M_dims[0], M_dims[1]);
  } else if (k < 0) {
    start = -k, size = std::min(M_dims[0] + k, M_dims[1]);
  } else {
    start = k * M_dims[0], size = std::min(M_dims[0], M_dims[1] - k);
  }
  const uword stride = M_dims[0] + 1;
  return Matrix<Tp, 1>(this->M_elem, start, size, stride);
}

template <class Tp>
inline SliceMatrix<Tp> Matrix<Tp, 2>::diag(int k) {
  uword start = 0, size = 0;
  if (k == 0) {
    start = 0, size = std::min(M_dims[0], M_dims[1]);
  } else if (k < 0) {
    start = -k, size = std::min(M_dims[0] + k, M_dims[1]);
  } else {
    start = k * M_dims[0], size = std::min(M_dims[0], M_dims[1] - k);
  }
  const uword stride = M_dims[0] + 1;
  return SliceMatrix<Tp>(this->M_elem, start, size, stride);
}

// Matrix member functions dealing with GsliceMatrix
template <class Tp>
inline Matrix<Tp, 2> Matrix<Tp, 2>::operator()(std::slice s1,
                                               std::slice s2) const {
  const uword start = sub2ind(s1.start(), s2.start());
  INIT_ARR2E(size, s2.size(), s1.size());
  INIT_ARR2E(stride, n_rows() * s2.stride(), 1 * s1.stride());
  return Matrix<Tp, 2>(this->M_elem, start, size, stride);
}

template <class Tp>
inline GsliceMatrix<Tp> Matrix<Tp, 2>::operator()(std::slice s1,
                                                  std::slice s2) {
  const uword start = sub2ind(s1.start(), s2.start());
  INIT_ARR2E(size, s2.size(), s1.size());
  INIT_ARR2E(stride, n_rows() * s2.stride(), 1 * s1.stride());
  return GsliceMatrix<Tp>(this->M_elem, start, size, stride);
}

template <class Tp>
inline Matrix<Tp, 3> Matrix<Tp, 3>::operator()(std::slice s1, std::slice s2,
                                               std::slice s3) const {
  const uword start = sub2ind(s1.start(), s2.start(), s3.start());
  INIT_ARR3E(size, s3.size(), s2.size(), s1.size());
  INIT_ARR3E(stride, _M_d1xd2 * s3.stride(), n_rows() * s2.stride(),
             s1.stride());
  return Matrix<Tp, 3>(this->M_elem, start, size, stride);
}

template <class Tp>
inline GsliceMatrix<Tp> Matrix<Tp, 3>::operator()(std::slice s1, std::slice s2,
                                                  std::slice s3) {
  const uword start = sub2ind(s1.start(), s2.start(), s3.start());
  INIT_ARR3E(size, s3.size(), s2.size(), s1.size());
  INIT_ARR3E(stride, _M_d1xd2 * s3.stride(), n_rows() * s2.stride(),
             s1.stride());
  return GsliceMatrix<Tp>(this->M_elem, start, size, stride);
}

template <class Tp>
inline Matrix<Tp, 2> Matrix<Tp, 2>::submat(uword fr, uword fc, uword lr,
                                           uword lc) const {
  const uword start = n_rows() * fc + fr;
  INIT_ARR2E(size, lc - fc + 1, lr - fr + 1);
  INIT_ARR2E(stride, n_rows(), 1);
  return Matrix<Tp, 2>(this->M_elem, start, size, stride);
}

template <class Tp>
inline GsliceMatrix<Tp> Matrix<Tp, 2>::submat(uword fr, uword fc, uword lr,
                                              uword lc) {
  const uword start = n_rows() * fc + fr;
  INIT_ARR2E(size, lc - fc + 1, lr - fr + 1);
  INIT_ARR2E(stride, n_rows(), 1);
  return GsliceMatrix<Tp>(this->M_elem, start, size, stride);
}

template <class Tp>
inline Matrix<Tp, 3> Matrix<Tp, 3>::subcube(uword fr, uword fc, uword fs,
                                            uword lr, uword lc,
                                            uword ls) const {
  const uword start = sub2ind(fr, fc, fs);
  INIT_ARR3E(size, ls - fs + 1, lc - fc + 1, lr - fr + 1);
  INIT_ARR3E(stride, _M_d1xd2, n_rows(), 1);
  return Matrix<Tp, 3>(this->M_elem, start, size, stride);
}

template <class Tp>
inline GsliceMatrix<Tp> Matrix<Tp, 3>::subcube(uword fr, uword fc, uword fs,
                                               uword lr, uword lc, uword ls) {
  const uword start = sub2ind(fr, fc, fs);
  INIT_ARR3E(size, ls - fs + 1, lc - fc + 1, lr - fr + 1);
  INIT_ARR3E(stride, _M_d1xd2, n_rows(), 1);
  return GsliceMatrix<Tp>(this->M_elem, start, size, stride);
}

template <class Tp>
inline Matrix<Tp, 2> Matrix<Tp, 3>::row(uword i) const {
  const uword start = i;
  INIT_ARR2E(size, n_slices(), n_cols());
  INIT_ARR2E(stride, _M_d1xd2, n_rows());
  return Matrix<Tp, 2>(this->M_elem, start, size, stride);
}

template <class Tp>
inline GsliceMatrix<Tp> Matrix<Tp, 3>::row(uword i) {
  const uword start = i;
  INIT_ARR2E(size, n_slices(), n_cols());
  INIT_ARR2E(stride, _M_d1xd2, n_rows());
  return GsliceMatrix<Tp>(this->M_elem, start, size, stride);
}

template <class Tp>
inline Matrix<Tp, 2> Matrix<Tp, 3>::col(uword j) const {
  const uword start = j * n_rows();
  INIT_ARR2E(size, n_slices(), n_rows());
  INIT_ARR2E(stride, _M_d1xd2, 1);
  return Matrix<Tp, 2>(this->M_elem, start, size, stride);
}

template <class Tp>
inline GsliceMatrix<Tp> Matrix<Tp, 3>::col(uword j) {
  const uword start = j * n_rows();
  INIT_ARR2E(size, n_slices(), n_rows());
  INIT_ARR2E(stride, _M_d1xd2, 1);
  return GsliceMatrix<Tp>(this->M_elem, start, size, stride);
}

template <class Tp>
inline Matrix<Tp, 2> Matrix<Tp, 3>::slice(uword k) const {
  const uword start = k * _M_d1xd2;
  INIT_ARR2E(size, n_cols(), n_rows());
  INIT_ARR2E(stride, n_rows(), 1);
  return Matrix<Tp, 2>(this->M_elem, start, size, stride);
}

template <class Tp>
inline GsliceMatrix<Tp> Matrix<Tp, 3>::slice(uword k) {
  const uword start = k * _M_d1xd2;
  INIT_ARR2E(size, n_cols(), n_rows());
  INIT_ARR2E(stride, n_rows(), 1);
  return GsliceMatrix<Tp>(this->M_elem, start, size, stride);
}

// Matrix member functions dealing with bool_array and MaskMatrix

template <class Tp>
inline Matrix<Tp, 1> Matrix<Tp, 1>::operator()(
    const bool_array& bool_arr) const {
  return Matrix<Tp, 1>(this->M_elem, bool_arr);
}

template <class Tp>
inline MaskMatrix<Tp> Matrix<Tp, 1>::operator()(const bool_array& bool_arr) {
  return MaskMatrix<Tp>(this->M_elem, bool_arr);
}

template <class Tp>
inline Matrix<Tp, 1> Matrix<Tp, 2>::operator()(
    const bool_array& bool_arr) const {
  return Matrix<Tp, 1>(this->M_elem, bool_arr);
}

template <class Tp>
inline MaskMatrix<Tp> Matrix<Tp, 2>::operator()(const bool_array& bool_arr) {
  return MaskMatrix<Tp>(this->M_elem, bool_arr);
}

template <class Tp>
inline Matrix<Tp, 1> Matrix<Tp, 3>::operator()(
    const bool_array& bool_arr) const {
  return Matrix<Tp, 1>(this->M_elem, bool_arr);
}

template <class Tp>
inline MaskMatrix<Tp> Matrix<Tp, 3>::operator()(const bool_array& bool_arr) {
  return MaskMatrix<Tp>(this->M_elem, bool_arr);
}

// Matrix member functions dealing with index_array and IndirectMatrix

template <class Tp>
inline Matrix<Tp, 1> Matrix<Tp, 1>::elem(const index_array& idx_arr) const {
  INIT_ARR1E(dims, idx_arr.size());
  return Matrix<Tp, 1>(this->M_elem, idx_arr, dims);
}

template <class Tp>
inline IndirectMatrix<Tp> Matrix<Tp, 1>::elem(const index_array& idx_arr) {
  INIT_ARR1E(dims, idx_arr.size());
  return IndirectMatrix<Tp>(this->M_elem, idx_arr, dims);
}
template <class Tp>
inline Matrix<Tp, 1> Matrix<Tp, 2>::elem(const index_array& idx_arr) const {
  INIT_ARR1E(dims, idx_arr.size());
  return Matrix<Tp, 1>(this->M_elem, idx_arr, dims);
}

template <class Tp>
inline IndirectMatrix<Tp> Matrix<Tp, 2>::elem(const index_array& idx_arr) {
  INIT_ARR1E(dims, idx_arr.size());
  return IndirectMatrix<Tp>(this->M_elem, idx_arr, dims);
}

template <class Tp>
inline Matrix<Tp, 1> Matrix<Tp, 3>::elem(const index_array& idx_arr) const {
  INIT_ARR1E(dims, idx_arr.size());
  return Matrix<Tp, 1>(this->M_elem, idx_arr, dims);
}

template <class Tp>
inline IndirectMatrix<Tp> Matrix<Tp, 3>::elem(const index_array& idx_arr) {
  INIT_ARR1E(dims, idx_arr.size());
  return IndirectMatrix<Tp>(this->M_elem, idx_arr, dims);
}

template <class Tp>
inline Matrix<Tp, 1> Matrix<Tp, 1>::operator()(
    const index_array& idx_arr) const {
  INIT_ARR1E(dims, idx_arr.size());
  return Matrix(this->M_elem, idx_arr, dims);
}

template <class Tp>
inline IndirectMatrix<Tp> Matrix<Tp, 1>::operator()(
    const index_array& idx_arr) {
  INIT_ARR1E(dims, idx_arr.size());
  return IndirectMatrix<Tp>(this->M_elem, idx_arr, dims);
}

template <class Tp>
inline Matrix<Tp, 2> Matrix<Tp, 2>::operator()(
    const index_array& idx_arr1, const index_array& idx_arr2) const {
  INIT_ARR2E(dims, idx_arr1.size(), idx_arr2.size());
  index_array idx_arr(dims[0] * dims[1]);
  uword idx = 0;
  for (uword j = 0; j < idx_arr2.size(); ++j) {
    for (uword i = 0; i < idx_arr1.size(); ++i) {
      idx_arr[idx++] = sub2ind(idx_arr1[i], idx_arr2[j]);
    }
  }
  return Matrix<Tp, 2>(this->M_elem, idx_arr, dims);
}

template <class Tp>
inline IndirectMatrix<Tp> Matrix<Tp, 2>::operator()(
    const index_array& idx_arr1, const index_array& idx_arr2) {
  INIT_ARR2E(dims, idx_arr1.size(), idx_arr2.size());
  index_array idx_arr(dims[0] * dims[1]);
  uword idx = 0;
  for (uword j = 0; j < idx_arr2.size(); ++j) {
    for (uword i = 0; i < idx_arr1.size(); ++i) {
      idx_arr[idx++] = sub2ind(idx_arr1[i], idx_arr2[j]);
    }
  }
  return IndirectMatrix<Tp>(this->M_elem, idx_arr, dims);
}

template <class Tp>
inline Matrix<Tp, 3> Matrix<Tp, 3>::operator()(
    const index_array& idx_arr1, const index_array& idx_arr2,
    const index_array& idx_arr3) const {
  INIT_ARR3E(dims, idx_arr1.size(), idx_arr2.size(), idx_arr3.size());
  index_array idx_arr(dims[0] * dims[1] * dims[2]);
  uword idx = 0;
  for (uword k = 0; k < idx_arr3.size(); ++k) {
    for (uword j = 0; j < idx_arr2.size(); ++j) {
      for (uword i = 0; i < idx_arr1.size(); ++i) {
        idx_arr[idx++] = sub2ind(idx_arr1[i], idx_arr2[j], idx_arr3[k]);
      }
    }
  }
  return Matrix<Tp, 3>(this->M_elem, idx_arr, dims);
}

template <class Tp>
inline IndirectMatrix<Tp> Matrix<Tp, 3>::operator()(
    const index_array& idx_arr1, const index_array& idx_arr2,
    const index_array& idx_arr3) {
  INIT_ARR3E(dims, idx_arr1.size(), idx_arr2.size(), idx_arr3.size());
  index_array idx_arr(dims[0] * dims[1] * dims[2]);
  uword idx = 0;
  for (uword k = 0; k < idx_arr3.size(); ++k) {
    for (uword j = 0; j < idx_arr2.size(); ++j) {
      for (uword i = 0; i < idx_arr1.size(); ++i) {
        idx_arr[idx++] = sub2ind(idx_arr1[i], idx_arr2[j], idx_arr3[k]);
      }
    }
  }
  return IndirectMatrix<Tp>(this->M_elem, idx_arr, dims);
}

//----------------------------------------------------------------------
// Other Matrix noninline member functions

template <class Tp>
Matrix<Tp, 2> Matrix<Tp, 1>::t() const {
  return Matrix<Tp, 2>(this->M_elem, 1, n_rows());
}

template <class Tp>
Matrix<Tp, 2> Matrix<Tp, 2>::t() const {
  uword n = n_rows(), m = n_cols();
  Matrix<Tp, 2> res(m, n);
  for (uword idx = 0; idx < this->n_elem(); ++idx) {
    uword i = idx / m;
    uword j = idx % m;
    res.M_elem[idx] = (*this)[n * j + i];
  }

  return res;
}

//----------------------------------------------------------------------
// Matrix functions

// matrix multiplication

template <typename Tp>
Matrix<Tp, 2> matmul(const Matrix<Tp, 1>& x, const Matrix<Tp, 2>& y) {
  matrix_assert(y.n_rows() == 1, "matmul(x, y): non-conformable arguments");
  const uword n = x.n_rows();
  const uword m = y.n_cols();
  Matrix<Tp, 2> res(n, m);
  for (uword i = 0; i != n; ++i)
    for (uword j = 0; j != m; ++j) res(i, j) = x(i) * y(0, j);
  return res;
}

template <typename Tp>
Matrix<Tp, 1> matmul(const Matrix<Tp, 2>& x, const Matrix<Tp, 1>& y) {
  matrix_assert(x.n_cols() == y.n_rows(),
                "matmul(x, y): non-conformable arguments");
  const uword nr = x.n_rows();
  const uword nc = x.n_cols();
  Matrix<Tp, 1> res(nr);
  for (uword i = 0; i != nr; ++i)
    for (uword j = 0; j != nc; ++j) res(i) += x(i, j) * y(j);
  return res;
}

template <typename Tp>
Matrix<Tp, 2> matmul(const Matrix<Tp, 2>& x, const Matrix<Tp, 2>& y) {
  matrix_assert(x.n_cols() == y.n_rows(),
                "matmul(x, y): non-conformable arguments");
  const uword nr = x.n_rows();
  const uword nc = x.n_cols();
  const uword p = y.n_cols();
  Matrix<Tp, 2> res(nr, p);
  for (uword i = 0; i != nr; ++i)
    for (uword j = 0; j != p; ++j)
      for (uword k = 0; k != nc; ++k) res(i, j) += x(i, k) * y(k, j);
  return res;
}

//----------------------------------------------------------------------
// A list of typedefs

typedef Matrix<double, 1> vec;
typedef Matrix<double, 2> mat;
typedef Matrix<double, 3> cube;

}  // namespace matrix_lib

#endif
