//  matrix.h: this single header matrix lib is a wrapper of std::valarray
//
//  Copyright (C) 2023-2024 Yi Pan <ypan1988@gmail.com>
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

#ifndef MATRIX_LIB_USE_CPP11
bool all(const bool_array& ba) {
  for (uword i = 0; i != ba.size(); ++i) {
    if (!ba[i]) return false;
  }
  return true;
}
bool any(const bool_array& ba) {
  for (uword i = 0; i != ba.size(); ++i) {
    if (ba[i]) return true;
  }
  return false;
}
#else
bool all(const bool_array& ba) {
  if (std::end(ba) == std::find(std::begin(ba), std::end(ba), false))
    return true;
  else
    return false;
}
bool any(const bool_array& ba) {
  if (std::begin(ba) == std::find(std::begin(ba), std::end(ba), true))
    return false;
  else
    return true;
}
#endif

template <class Tp>
bool all_equal(const std::valarray<Tp>& va1, const std::valarray<Tp>& va2) {
  return all(va1 == va2);
}

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

template <class Tp>
struct SliceMatrix;
template <class Tp, uword Size>
struct GsliceMatrix;
template <class Tp>
struct MaskMatrix;
template <class Tp, uword Size>
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
  Matrix_base(const std::valarray<Tp>&  x) : M_elem(x   ) {}
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
  const std::valarray<Tp>& e()    const { return M_elem; }
  const std::valarray<Tp>& elem() const { return M_elem; }

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
  Matrix()                                : Matrix_base<Tp>(            ) { M_init(); }    // (1)
  explicit Matrix(              uword n1) : Matrix_base<Tp>(          n1) { M_init(); }    // (2)
  Matrix(const elem_type& x,    uword n1) : Matrix_base<Tp>(x,        n1) { M_init(); }    // (3)
  Matrix(const elem_type* x,    uword n1) : Matrix_base<Tp>(x,        n1) { M_init(); }    // (4)
  Matrix(const         Matrix       &  x) : Matrix_base<Tp>(x.M_elem    ) { M_init(); }    // (5)
  Matrix(const    SliceMatrix<Tp   >&  x) : Matrix_base<Tp>(x.sub_elem()) { M_init(x.is_col()); }    // (7)
  Matrix(const   GsliceMatrix<Tp, 1>&  x) : Matrix_base<Tp>(x.sub_elem()) { M_init(); }    // (8)
  Matrix(const     MaskMatrix<Tp   >&  x) : Matrix_base<Tp>(x.sub_elem()) { M_init(); }    // (9)
  Matrix(const IndirectMatrix<Tp, 1>&  x) : Matrix_base<Tp>(x.sub_elem()) { M_init(); }    // (10)
#if defined(MATRIX_LIB_USE_CPP11)
  Matrix(Matrix&& x) = default;                                                            // (6)
  Matrix(std::initializer_list<Tp>     x) : Matrix_base<Tp>(x           ) { M_init(); }    // (11)
  Matrix(Matrix<Tp, 2>&& x);
#endif
  Matrix(const std::valarray<Tp>     & x) : Matrix_base<Tp>(x           ) { M_init(); }    // (12)

#if defined(MATRIX_LIB_USE_R)
  Matrix(SEXP                          x) : Matrix_base<Tp>(x           ) { M_init(); }
#endif
  Matrix(const std::valarray<Tp>& x, const uword        dims[1]) : Matrix_base<Tp>(x) { M_init(dims); }
  Matrix(const std::valarray<Tp>& x, const index_array& dims   ) : Matrix_base<Tp>(x) { M_init(dims); }
  Matrix(const Matrix<Tp, 2>& x);
  ~Matrix() {}

  // assignment:
  Matrix& operator=(const         Matrix       & x) { this->set_elem(x.M_elem                ); M_init();     return *this; }    // (1)
  Matrix& operator=(const       elem_type      & x) { this->M_elem = x                        ;               return *this; }    // (3)
  Matrix& operator=(const    SliceMatrix<Tp   >& x) { this->set_elem(x.sub_elem(), x.n_elem()); M_init(x.is_col()); return *this; }    // (4)
  Matrix& operator=(const   GsliceMatrix<Tp, 1>& x) { this->set_elem(x.sub_elem(), x.n_elem()); M_init();     return *this; }    // (5)
  Matrix& operator=(const     MaskMatrix<Tp   >& x) { this->set_elem(x.sub_elem(), x.n_elem()); M_init(true); return *this; }    // (6)
  Matrix& operator=(const IndirectMatrix<Tp, 1>& x) { this->set_elem(x.sub_elem(), x.n_elem()); M_init();     return *this; }    // (7)
#if defined(MATRIX_LIB_USE_CPP11)
  Matrix& operator=(Matrix&& x) = default;                                                                                       // (2)
  Matrix& operator=(std::initializer_list<Tp>    x) { this->M_elem = x                        ; M_init();     return *this; }    // (8)
  Matrix& operator=(Matrix<Tp, 2>&& x);                                                                                          // (10)
#endif
  Matrix& operator=(const Matrix<Tp, 2>        & x);                                                                             // (9)

#if defined(MATRIX_LIB_USE_R)
  // this operator enables implicit Rcpp::wrap
  operator SEXP()       { return export_matrix_to_sexp(); }
  operator SEXP() const { return export_matrix_to_sexp(); }
  SEXP export_matrix_to_sexp() const;
#endif

  index_array size() const {
    if (is_column_vector) return index_array(M_dims, 1);
    INIT_ARR2E(res, 1, this->n_elem());
    return res;
  }
  uword n_rows() const { return  is_column_vector ? M_dims[0] : (M_dims[0] ? 1 : 0); }
  uword n_cols() const { return !is_column_vector ? M_dims[0] : (M_dims[0] ? 1 : 0); }

  void M_range_check(uword n1) const {
    matrix_assert(n1 < this->n_elem(),
            "Matrix<T, 1>::M_range_check:\n" <<
            "  Index " << n1 << " is out of bound for axis 0 with size " << n_rows());
  }

  // subscripting:
        elem_type& operator()(uword n1)       { M_range_check(n1); return this->M_elem[n1]; }
  const elem_type& operator()(uword n1) const { M_range_check(n1); return this->M_elem[n1]; }

  // SliceMatrix related member functions
          Matrix<Tp, 1> operator()(std::slice s) const;
     SliceMatrix<Tp   > operator()(std::slice s);
          Matrix<Tp, 1> subvec(uword fi, uword li) const;
     SliceMatrix<Tp   > subvec(uword fi, uword li);

  // MaskMatrix related member functions
          Matrix<Tp, 1> operator()(const bool_array& ba) const;
      MaskMatrix<Tp   > operator()(const bool_array& ba);

  // IndirectMatrix related member functions
          Matrix<Tp, 1> operator()(const index_array& ia) const;
  IndirectMatrix<Tp, 1> operator()(const index_array& ia);

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
    M_init(is_colvec);
  }
  void M_init(const index_array& dims) {
    bool is_rowvec = (dims.size() == 2) && (dims[0] == 1);
    M_init(!is_rowvec);
  }
  Matrix(const std::valarray<Tp>& va, uword start, const uword size,
         const uword stride, bool is_colvec = true)
      : Matrix_base<Tp>(va[std::slice(start, size, stride)]) {
    M_init(is_colvec);
  }
  Matrix(const std::valarray<Tp>& va, uword start, const index_array& size,
         const index_array& stride)
      : Matrix_base<Tp>(va[std::gslice(start, size, stride)]) {
    M_init();
  }
  Matrix(const std::valarray<Tp>& va, const bool_array& ba)
      : Matrix_base<Tp>(va[ba]) {
    M_init();
  }
  Matrix(const std::valarray<Tp>& va, const index_array& ia,
         const index_array& dims)
      : Matrix_base<Tp>(va[ia]) {
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
  Matrix()                                       : Matrix_base<Tp>(            ) { M_init( 0,  0    ); }    // (1)
  Matrix(                    uword n1, uword n2) : Matrix_base<Tp>(     n1 * n2) { M_init(n1, n2    ); }    // (2)
  Matrix(const elem_type& x, uword n1, uword n2) : Matrix_base<Tp>(x,   n1 * n2) { M_init(n1, n2    ); }    // (3)
  Matrix(const elem_type* x, uword n1, uword n2) : Matrix_base<Tp>(x,   n1 * n2) { M_init(n1, n2    ); }    // (4)
  Matrix(const         Matrix       & x        ) : Matrix_base<Tp>(x.M_elem    ) { M_init(x.M_dims  ); }    // (5)
  Matrix(const    SliceMatrix<Tp   >& x        ) : Matrix_base<Tp>(x.sub_elem()) { M_init(x.is_col()); }    // (7)
  Matrix(const   GsliceMatrix<Tp, 2>& x        ) : Matrix_base<Tp>(x.sub_elem()) { M_init(x.size()  ); }    // (8)
  Matrix(const     MaskMatrix<Tp   >& x        ) : Matrix_base<Tp>(x.sub_elem()) { M_init(true      ); }    // (9)
  Matrix(const IndirectMatrix<Tp, 2>& x        ) : Matrix_base<Tp>(x.sub_elem()) { M_init(x.size()  ); }    // (10)
#if defined(MATRIX_LIB_USE_CPP11)
  Matrix(Matrix&& x) = default;                                                                             // (6)
  Matrix(std::initializer_list<std::initializer_list<Tp>> x);                                               // (11)
  Matrix(Matrix<Tp, 1>&& x);
#endif
  Matrix(const std::valarray<Tp>&  x, uword n1, uword n2) : Matrix_base<Tp>(x)   { M_init(n1, n2    ); }    // (12)
#if defined(MATRIX_LIB_USE_R)
  Matrix(SEXP x) : Matrix_base<Tp>(x) { SEXP dims = Rf_getAttrib(x, R_DimSymbol); M_init(INTEGER(dims)[0], INTEGER(dims)[1]); }
#endif
  Matrix(const std::valarray<Tp>& x, const uword        dims[2]) : Matrix_base<Tp>(x)  { M_init(dims)    ; }
  Matrix(const std::valarray<Tp>& x, const index_array& dims   ) : Matrix_base<Tp>(x)  { M_init(dims)    ; }
  Matrix(const Matrix<Tp, 1>& x);
  ~Matrix() {}

  // assignment:
  Matrix& operator=(const         Matrix       & x) { this->set_elem(x.M_elem                ); M_init(x.M_dims  ); return *this; }    // (1)
  Matrix& operator=(const       elem_type      & x) { this->M_elem = x                        ;                     return *this; }    // (3)
  Matrix& operator=(const    SliceMatrix<Tp   >& x) { this->set_elem(x.sub_elem(), x.n_elem()); M_init(x.is_col()); return *this; }    // (4)
  Matrix& operator=(const   GsliceMatrix<Tp, 2>& x) { this->set_elem(x.sub_elem(), x.n_elem()); M_init(x.size()  ); return *this; }    // (5)
  Matrix& operator=(const     MaskMatrix<Tp   >& x) { this->set_elem(x.sub_elem(), x.n_elem()); M_init(true      ); return *this; }    // (6)
  Matrix& operator=(const IndirectMatrix<Tp, 2>& x) { this->set_elem(x.sub_elem(), x.n_elem()); M_init(x.size()  ); return *this; }    // (7)
#if defined(MATRIX_LIB_USE_CPP11)
  Matrix& operator=(Matrix&& x) = default;                                                                                             // (2)
  Matrix& operator=(std::initializer_list<std::initializer_list<Tp>> x);                                                               // (8)
  Matrix& operator=(Matrix<Tp, 1>            && x);                                                                                    // (10)
#endif
  Matrix& operator=(const        Matrix<Tp, 1>& x);                                                                                    // (9)

#if defined(MATRIX_LIB_USE_R)
  // this operator enables implicit Rcpp::wrap
  operator SEXP()       { return export_matrix_to_sexp(); }
  operator SEXP() const { return export_matrix_to_sexp(); }
  SEXP export_matrix_to_sexp() const;
#endif

  index_array size() const { return index_array(M_dims, 2); }
  uword     n_rows() const { return M_dims[0]; }
  uword     n_cols() const { return M_dims[1]; }

  void M_range_check(uword n1, uword n2) const {
    matrix_assert(n1 < n_rows(), "Matrix<T, 2>::M_range_check: index is out of bound for dimension 1");
    matrix_assert(n2 < n_cols(), "Matrix<T, 2>::M_range_check: index is out of bound for dimension 2");
  }

  // subscripting:
        elem_type& operator()(uword n1, uword n2)
  { M_range_check(n1, n2); return this->M_elem[sub2ind(n1, n2)]; }
  const elem_type& operator()(uword n1, uword n2) const
  { M_range_check(n1, n2); return this->M_elem[sub2ind(n1, n2)]; }

  // SliceMatrix related member functions
          Matrix<Tp, 1> row(uword i) const;
     SliceMatrix<Tp   > row(uword i);
          Matrix<Tp, 1> col(uword j) const;
     SliceMatrix<Tp   > col(uword j);
          Matrix<Tp, 1> diag(int k = 0) const;
     SliceMatrix<Tp   > diag(int k = 0);

  // GsliceMatrix related member functions
          Matrix<Tp, 2> operator()(std::slice s1, std::slice s2) const;
    GsliceMatrix<Tp, 2> operator()(std::slice s1, std::slice s2);
          Matrix<Tp, 2> submat(uword fr, uword fc, uword lr, uword lc) const;
    GsliceMatrix<Tp, 2> submat(uword fr, uword fc, uword lr, uword lc);
          Matrix<Tp, 2> rows(uword fr, uword lr) const { return submat(fr, 0, lr, n_cols() - 1); }
    GsliceMatrix<Tp, 2> rows(uword fr, uword lr)       { return submat(fr, 0, lr, n_cols() - 1); }
          Matrix<Tp, 2> cols(uword fc, uword lc) const { return submat(0, fc, n_rows() - 1, lc); }
    GsliceMatrix<Tp, 2> cols(uword fc, uword lc)       { return submat(0, fc, n_rows() - 1, lc); }

  // MaskMatrix related member functions
          Matrix<Tp, 1> operator()(const bool_array& ba) const;
      MaskMatrix<Tp   > operator()(const bool_array& ba);

  // IndirectMatrix related member functions
          Matrix<Tp, 1> operator()(const index_array& ia) const;
  IndirectMatrix<Tp, 1> operator()(const index_array& ia);
          Matrix<Tp, 2> operator()(const index_array& ia1, const index_array& ia2) const;
  IndirectMatrix<Tp, 2> operator()(const index_array& ia1, const index_array& ia2);

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
  Matrix(const std::valarray<Tp>& va, const index_array& ia,
         const index_array& dims)
      : Matrix_base<Tp>(va[ia]) {
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
  Matrix()                                                 : Matrix_base<Tp>(            ) { M_init( 0,  0,  0); }    // (1)
  Matrix(                    uword n1, uword n2, uword n3) : Matrix_base<Tp>(    n1*n2*n3) { M_init(n1, n2, n3); }    // (2)
  Matrix(const elem_type& x, uword n1, uword n2, uword n3) : Matrix_base<Tp>(x,  n1*n2*n3) { M_init(n1, n2, n3); }    // (3)
  Matrix(const elem_type* x, uword n1, uword n2, uword n3) : Matrix_base<Tp>(x,  n1*n2*n3) { M_init(n1, n2, n3); }    // (4)
  Matrix(const         Matrix       & x                  ) : Matrix_base<Tp>(x.M_elem    ) { M_init(x.M_dims  ); }    // (5)
  Matrix(const   GsliceMatrix<Tp, 3>& x                  ) : Matrix_base<Tp>(x.sub_elem()) { M_init(x.size()  ); }    // (8)
  Matrix(const IndirectMatrix<Tp, 3>& x                  ) : Matrix_base<Tp>(x.sub_elem()) { M_init(x.size()  ); }    // (10)
#if defined(MATRIX_LIB_USE_CPP11)
  Matrix(Matrix&& x) = default;                                                                                       // (6)
  Matrix(std::initializer_list<std::initializer_list<std::initializer_list<Tp>>> x);                                  // (11)

#endif
  Matrix(const std::valarray<Tp>&  x, uword n1, uword n2, uword n3) : Matrix_base<Tp>(x)   { M_init(n1, n2, n3); }    // (12)
  Matrix(const std::valarray<Tp>&  x, const uword          dims[3]) : Matrix_base<Tp>(x)   { M_init(dims      ); }
  Matrix(const std::valarray<Tp>&  x, const index_array&   dims   ) : Matrix_base<Tp>(x)   { M_init(dims      ); }
#if defined(MATRIX_LIB_USE_R)
  Matrix(SEXP x) : Matrix_base<Tp>(x) { SEXP dims = Rf_getAttrib(x, R_DimSymbol); M_init(INTEGER(dims)[0], INTEGER(dims)[1], INTEGER(dims)[2]); }
#endif
  ~Matrix() {}

  // assignment:
  Matrix& operator=(const         Matrix       & x) { this->set_elem(x.M_elem                ); M_init(x.M_dims); return *this; }    // (1)
  Matrix& operator=(const       elem_type      & x) { this->M_elem = x                        ;                   return *this; }    // (3)
  Matrix& operator=(const   GsliceMatrix<Tp, 3>& x) { this->set_elem(x.sub_elem(), x.n_elem()); M_init(x.size()); return *this; }    // (5)
  Matrix& operator=(const IndirectMatrix<Tp, 3>& x) { this->set_elem(x.sub_elem(), x.n_elem()); M_init(x.size()); return *this; }    // (7)
#if defined(MATRIX_LIB_USE_CPP11)
  Matrix& operator=(Matrix&& x) = default;                                                                                           // (2)
  Matrix& operator=(std::initializer_list<std::initializer_list<std::initializer_list<Tp>>> x);                                      // (8)
#endif

#if defined(MATRIX_LIB_USE_R)
  // this operator enables implicit Rcpp::wrap
  operator SEXP()       { return export_matrix_to_sexp(); }
  operator SEXP() const { return export_matrix_to_sexp(); }
  SEXP export_matrix_to_sexp() const;
#endif

  index_array size() const { return index_array(M_dims, 3); }
  uword     n_rows() const { return M_dims[0]; }
  uword     n_cols() const { return M_dims[1]; }
  uword   n_slices() const { return M_dims[2]; }

  void M_range_check(uword n1, uword n2, uword n3) const {
    matrix_assert(n1 < n_rows()  , "Matrix<T, 3>::M_range_check: index is out of bound for dimension 1");
    matrix_assert(n2 < n_cols()  , "Matrix<T, 3>::M_range_check: index is out of bound for dimension 2");
    matrix_assert(n3 < n_slices(), "Matrix<T, 3>::M_range_check: index is out of bound for dimension 3");
  }

  // subscripting:
        elem_type& operator()(uword n1, uword n2, uword n3)
  { M_range_check(n1, n2, n3); return this->M_elem[sub2ind(n1, n2, n3)]; }
  const elem_type& operator()(uword n1, uword n2, uword n3) const
  { M_range_check(n1, n2, n3); return this->M_elem[sub2ind(n1, n2, n3)]; }

  // GsliceMatrix related member functions
          Matrix<Tp, 3> operator()(std::slice s1, std::slice s2, std::slice s3) const;
    GsliceMatrix<Tp, 3> operator()(std::slice s1, std::slice s2, std::slice s3);
          Matrix<Tp, 3> subcube(uword fr, uword fc, uword fs, uword lr, uword lc, uword ls) const;
    GsliceMatrix<Tp, 3> subcube(uword fr, uword fc, uword fs, uword lr, uword lc, uword ls);
          Matrix<Tp, 2>    row(uword i) const;
    GsliceMatrix<Tp, 2>    row(uword i);
          Matrix<Tp, 2>    col(uword j) const;
    GsliceMatrix<Tp, 2>    col(uword j);
          Matrix<Tp, 2>  slice(uword k) const;
    GsliceMatrix<Tp, 2>  slice(uword k);
          Matrix<Tp, 3>   rows(uword fr, uword lr) const { return subcube(fr,  0,  0,           lr, n_cols() - 1, n_slices() - 1); }
    GsliceMatrix<Tp, 3>   rows(uword fr, uword lr)       { return subcube(fr,  0,  0,           lr, n_cols() - 1, n_slices() - 1); }
          Matrix<Tp, 3>   cols(uword fc, uword lc) const { return subcube( 0, fc,  0, n_rows() - 1,           lc, n_slices() - 1); }
    GsliceMatrix<Tp, 3>   cols(uword fc, uword lc)       { return subcube( 0, fc,  0, n_rows() - 1,           lc, n_slices() - 1); }
          Matrix<Tp, 3> slices(uword fs, uword ls) const { return subcube( 0,  0, fs, n_rows() - 1, n_cols() - 1,             ls); }
    GsliceMatrix<Tp, 3> slices(uword fs, uword ls)       { return subcube( 0,  0, fs, n_rows() - 1, n_cols() - 1,             ls); }

  // MaskMatrix related member functions
          Matrix<Tp, 1> operator()(const bool_array& ba) const;
      MaskMatrix<Tp   > operator()(const bool_array& ba);

  // IndirectMatrix related member functions
          Matrix<Tp, 1> operator()(const index_array& ia) const;
  IndirectMatrix<Tp, 1> operator()(const index_array& ia);
          Matrix<Tp, 3> operator()(const index_array& ia1, const index_array& ia2, const index_array& ia3) const;
  IndirectMatrix<Tp, 3> operator()(const index_array& ia1, const index_array& ia2, const index_array& ia3);

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
  uword M_dims[3] = {0, 0, 0}, M_nrxnc = 0;
#else
  uword M_dims[3], M_nrxnc;
#endif
  void M_init(uword n1, uword n2, uword n3) {
    if (this->n_elem() != n1 * n2 * n3)
      error("3D Cstor error: dimension mismatch");
    M_dims[0] = n1, M_dims[1] = n2, M_dims[2] = n3;
    M_nrxnc = n1 * n2;
  }
  void M_init(const uword dims[3]) {
    M_dims[0] = dims[0], M_dims[1] = dims[1], M_dims[2] = dims[2];
    M_nrxnc = dims[0] * dims[1];
  }
  void M_init(const index_array& dims) {
    M_dims[0] = dims[0], M_dims[1] = dims[1], M_dims[2] = dims[2];
    M_nrxnc = dims[0] * dims[1];
  }
  Matrix(const std::valarray<Tp>& va, uword start, const index_array& size,
         const index_array& stride)
      : Matrix_base<Tp>(va[std::gslice(start, size, stride)]) {
    M_init(size[2], size[1], size[0]);
  }
  Matrix(const std::valarray<Tp>& va, const index_array& ia,
         const index_array& dims)
      : Matrix_base<Tp>(va[ia]) {
    M_init(dims);
  }
  uword sub2ind(uword i, uword j, uword k) const {
    return i + j * M_dims[0] + k * M_nrxnc;
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
Matrix<Tp, 1>::Matrix(const Matrix<Tp, 2>& x) : Matrix_base<Tp>(x.elem()) {
  if (x.n_cols() != 1 && x.n_rows() != 1) error("x is not a col/row vector");
  M_init(x.n_cols() == 1);
}

template <class Tp>
Matrix<Tp, 2>::Matrix(const Matrix<Tp, 1>& x) : Matrix_base<Tp>(x.elem()) {
  M_init(x.is_column_vector);
}

template <class Tp>
Matrix<Tp, 1>& Matrix<Tp, 1>::operator=(const Matrix<Tp, 2>& x) {
  if (x.n_cols() != 1 && x.n_rows() != 1) error("x is not a col/row vector");
  this->set_elem(x.elem());
  M_init(x.n_cols() == 1);
  return *this;
}

template <class Tp>
Matrix<Tp, 2>& Matrix<Tp, 2>::operator=(const Matrix<Tp, 1>& x) {
  this->set_elem(x.elem());
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

template <class Tp, uword Size>
struct SubMatrix {
 public:
  typedef Tp elem_type;

  SubMatrix(std::valarray<Tp>& va) : M_elem(va), M_dims(Size), M_n_elem(0) {}
  virtual ~SubMatrix() {}

  virtual std::valarray<Tp> elem() const { return std::valarray<Tp>(); }
  virtual index_array size() const { return this->M_dims; }
  uword n_elem() const { return M_n_elem; }

  elem_type sum() const { return elem().sum(); }
  elem_type min() const { return elem().min(); }
  elem_type max() const { return elem().max(); }

 protected:
  std::valarray<Tp>& M_elem;
  index_array M_dims;
  uword M_n_elem;
};

//----------------------------------------------------------------------
// SliceMatrix

template <class Tp>
struct SliceMatrix : public SubMatrix<Tp, 1> {
 public:
  typedef Tp elem_type;

  SliceMatrix(std::valarray<Tp>& va, uword start, uword size, uword stride,
              bool is_colvec = true)
      : SubMatrix<Tp, 1>(va),
        M_desc(start, size, stride),
        is_column_vector(is_colvec) {
    this->M_dims[0] = size;
    this->M_n_elem = size;
  }
  std::valarray<Tp> elem() const {
    return std::valarray<Tp>(this->M_elem[M_desc]);
  }
  std::slice_array<Tp> sub_elem() const { return this->M_elem[M_desc]; }
  bool is_col() const { return is_column_vector; }
  index_array size() const {
    if (is_column_vector) return this->M_dims;
    INIT_ARR2E(res, 1, this->n_elem());
    return res;
  }

  // clang-format off
  SliceMatrix& operator= (const elem_type& x) { this->M_elem[M_desc]  =                   x                 ; return *this; }
  SliceMatrix& operator+=(const elem_type& x) { this->M_elem[M_desc] += std::valarray<Tp>(x, this->n_elem()); return *this; }
  SliceMatrix& operator-=(const elem_type& x) { this->M_elem[M_desc] -= std::valarray<Tp>(x, this->n_elem()); return *this; }
  SliceMatrix& operator*=(const elem_type& x) { this->M_elem[M_desc] *= std::valarray<Tp>(x, this->n_elem()); return *this; }
  SliceMatrix& operator/=(const elem_type& x) { this->M_elem[M_desc] /= std::valarray<Tp>(x, this->n_elem()); return *this; }
  SliceMatrix& operator%=(const elem_type& x) { this->M_elem[M_desc] %= std::valarray<Tp>(x, this->n_elem()); return *this; }
  // clang-format on

 private:
  std::slice M_desc;
  bool is_column_vector;
};

//----------------------------------------------------------------------
// GsliceMatrix

template <class Tp, uword Size>
struct GsliceMatrix : public SubMatrix<Tp, Size> {
 public:
  typedef Tp elem_type;

  GsliceMatrix(std::valarray<Tp>& va, uword start, const index_array& size,
               const index_array& stride)
      : SubMatrix<Tp, Size>(va), M_desc(start, size, stride) {
    this->M_n_elem = 1;
    for (uword idx = 0; idx < Size; ++idx) {
      this->M_dims[idx] = size[Size - idx - 1];
      this->M_n_elem *= this->M_dims[idx];
    }
  }
  std::valarray<Tp> elem() const {
    return std::valarray<Tp>(this->M_elem[M_desc]);
  }
  std::gslice_array<Tp> sub_elem() const { return this->M_elem[M_desc]; }

  // clang-format off
  GsliceMatrix& operator= (const elem_type& x) { this->M_elem[M_desc]  =                   x                 ; return *this; }
  GsliceMatrix& operator+=(const elem_type& x) { this->M_elem[M_desc] += std::valarray<Tp>(x, this->n_elem()); return *this; }
  GsliceMatrix& operator-=(const elem_type& x) { this->M_elem[M_desc] -= std::valarray<Tp>(x, this->n_elem()); return *this; }
  GsliceMatrix& operator*=(const elem_type& x) { this->M_elem[M_desc] *= std::valarray<Tp>(x, this->n_elem()); return *this; }
  GsliceMatrix& operator/=(const elem_type& x) { this->M_elem[M_desc] /= std::valarray<Tp>(x, this->n_elem()); return *this; }
  GsliceMatrix& operator%=(const elem_type& x) { this->M_elem[M_desc] %= std::valarray<Tp>(x, this->n_elem()); return *this; }
  // clang-format on

 private:
  std::gslice M_desc;
};

//----------------------------------------------------------------------
// MaskMatrix

template <class Tp>
struct MaskMatrix : public SubMatrix<Tp, 1> {
 public:
  typedef Tp elem_type;

  MaskMatrix(std::valarray<Tp>& va, const bool_array& ba)
      : SubMatrix<Tp, 1>(va), M_desc(ba) {
    std::valarray<uword> tmp((uword)0, ba.size());
    tmp[ba] = 1;
    this->M_dims[0] = this->M_n_elem = tmp.sum();
  }
  std::valarray<Tp> elem() const {
    return std::valarray<Tp>(this->M_elem[M_desc]);
  }
  std::mask_array<Tp> sub_elem() const { return this->M_elem[M_desc]; }

  // clang-format off
  MaskMatrix& operator= (const elem_type& x) { this->M_elem[M_desc]  =                   x                 ; return *this; }
  MaskMatrix& operator+=(const elem_type& x) { this->M_elem[M_desc] += std::valarray<Tp>(x, this->n_elem()); return *this; }
  MaskMatrix& operator-=(const elem_type& x) { this->M_elem[M_desc] -= std::valarray<Tp>(x, this->n_elem()); return *this; }
  MaskMatrix& operator*=(const elem_type& x) { this->M_elem[M_desc] *= std::valarray<Tp>(x, this->n_elem()); return *this; }
  MaskMatrix& operator/=(const elem_type& x) { this->M_elem[M_desc] /= std::valarray<Tp>(x, this->n_elem()); return *this; }
  MaskMatrix& operator%=(const elem_type& x) { this->M_elem[M_desc] %= std::valarray<Tp>(x, this->n_elem()); return *this; }
  // clang-format on

 private:
  bool_array M_desc;
};

//----------------------------------------------------------------------
// IndirectMatrix

template <class Tp, uword Size>
struct IndirectMatrix : public SubMatrix<Tp, Size> {
 public:
  typedef Tp elem_type;

  IndirectMatrix(std::valarray<Tp>& va, const index_array& ia,
                 const index_array& dims)
      : SubMatrix<Tp, Size>(va), M_desc(ia) {
    this->M_dims = dims;
    this->M_n_elem = M_desc.size();
  }
  std::valarray<Tp> elem() const {
    return std::valarray<Tp>(this->M_elem[M_desc]);
  }
  std::indirect_array<Tp> sub_elem() const { return this->M_elem[M_desc]; }

  // clang-format off
  IndirectMatrix& operator= (const elem_type& x) { this->M_elem[M_desc]  =                   x                 ; return *this; }
  IndirectMatrix& operator+=(const elem_type& x) { this->M_elem[M_desc] += std::valarray<Tp>(x, this->n_elem()); return *this; }
  IndirectMatrix& operator-=(const elem_type& x) { this->M_elem[M_desc] -= std::valarray<Tp>(x, this->n_elem()); return *this; }
  IndirectMatrix& operator*=(const elem_type& x) { this->M_elem[M_desc] *= std::valarray<Tp>(x, this->n_elem()); return *this; }
  IndirectMatrix& operator/=(const elem_type& x) { this->M_elem[M_desc] /= std::valarray<Tp>(x, this->n_elem()); return *this; }
  IndirectMatrix& operator%=(const elem_type& x) { this->M_elem[M_desc] %= std::valarray<Tp>(x, this->n_elem()); return *this; }
  // clang-format on

 private:
  index_array M_desc;
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
inline Matrix<Tp, 1> Matrix<Tp, 1>::subvec(uword fi, uword li) const {
  if (fi > li || li >= M_dims[0]) error("1D subscription error");
  const uword start = fi;
  const uword size = li - fi + 1;
  const uword stride = 1;
  return Matrix<Tp, 1>(this->M_elem, start, size, stride, is_column_vector);
}

template <class Tp>
inline SliceMatrix<Tp> Matrix<Tp, 1>::subvec(uword fi, uword li) {
  if (fi > li || li >= M_dims[0]) error("1D subscription error");
  const uword start = fi;
  const uword size = li - fi + 1;
  const uword stride = 1;
  return SliceMatrix<Tp>(this->M_elem, start, size, stride, is_column_vector);
}

template <class Tp>
inline Matrix<Tp, 1> Matrix<Tp, 2>::row(uword i) const {
  M_range_check(i, 0);
  const uword start = i;
  const uword size = M_dims[1];
  const uword stride = M_dims[0];
  return Matrix<Tp, 1>(this->M_elem, start, size, stride, false);
}

template <class Tp>
inline SliceMatrix<Tp> Matrix<Tp, 2>::row(uword i) {
  M_range_check(i, 0);
  const uword start = i;
  const uword size = M_dims[1];
  const uword stride = M_dims[0];
  return SliceMatrix<Tp>(this->M_elem, start, size, stride, false);
}

template <class Tp>
inline Matrix<Tp, 1> Matrix<Tp, 2>::col(uword j) const {
  M_range_check(0, j);
  const uword start = j * M_dims[0];
  const uword size = M_dims[0];
  const uword stride = 1;
  return Matrix<Tp, 1>(this->M_elem, start, size, stride);
}

template <class Tp>
inline SliceMatrix<Tp> Matrix<Tp, 2>::col(uword j) {
  M_range_check(0, j);
  const uword start = j * M_dims[0];
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
inline GsliceMatrix<Tp, 2> Matrix<Tp, 2>::operator()(std::slice s1,
                                                     std::slice s2) {
  const uword start = sub2ind(s1.start(), s2.start());
  INIT_ARR2E(size, s2.size(), s1.size());
  INIT_ARR2E(stride, n_rows() * s2.stride(), 1 * s1.stride());
  return GsliceMatrix<Tp, 2>(this->M_elem, start, size, stride);
}

template <class Tp>
inline Matrix<Tp, 3> Matrix<Tp, 3>::operator()(std::slice s1, std::slice s2,
                                               std::slice s3) const {
  const uword start = sub2ind(s1.start(), s2.start(), s3.start());
  INIT_ARR3E(size, s3.size(), s2.size(), s1.size());
  INIT_ARR3E(stride, M_nrxnc * s3.stride(), n_rows() * s2.stride(),
             s1.stride());
  return Matrix<Tp, 3>(this->M_elem, start, size, stride);
}

template <class Tp>
inline GsliceMatrix<Tp, 3> Matrix<Tp, 3>::operator()(std::slice s1,
                                                     std::slice s2,
                                                     std::slice s3) {
  const uword start = sub2ind(s1.start(), s2.start(), s3.start());
  INIT_ARR3E(size, s3.size(), s2.size(), s1.size());
  INIT_ARR3E(stride, M_nrxnc * s3.stride(), n_rows() * s2.stride(),
             s1.stride());
  return GsliceMatrix<Tp, 3>(this->M_elem, start, size, stride);
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
inline GsliceMatrix<Tp, 2> Matrix<Tp, 2>::submat(uword fr, uword fc, uword lr,
                                                 uword lc) {
  const uword start = n_rows() * fc + fr;
  INIT_ARR2E(size, lc - fc + 1, lr - fr + 1);
  INIT_ARR2E(stride, n_rows(), 1);
  return GsliceMatrix<Tp, 2>(this->M_elem, start, size, stride);
}

template <class Tp>
inline Matrix<Tp, 3> Matrix<Tp, 3>::subcube(uword fr, uword fc, uword fs,
                                            uword lr, uword lc,
                                            uword ls) const {
  const uword start = sub2ind(fr, fc, fs);
  INIT_ARR3E(size, ls - fs + 1, lc - fc + 1, lr - fr + 1);
  INIT_ARR3E(stride, M_nrxnc, n_rows(), 1);
  return Matrix<Tp, 3>(this->M_elem, start, size, stride);
}

template <class Tp>
inline GsliceMatrix<Tp, 3> Matrix<Tp, 3>::subcube(uword fr, uword fc, uword fs,
                                                  uword lr, uword lc,
                                                  uword ls) {
  const uword start = sub2ind(fr, fc, fs);
  INIT_ARR3E(size, ls - fs + 1, lc - fc + 1, lr - fr + 1);
  INIT_ARR3E(stride, M_nrxnc, n_rows(), 1);
  return GsliceMatrix<Tp, 3>(this->M_elem, start, size, stride);
}

template <class Tp>
inline Matrix<Tp, 2> Matrix<Tp, 3>::row(uword i) const {
  const uword start = i;
  INIT_ARR2E(size, n_slices(), n_cols());
  INIT_ARR2E(stride, M_nrxnc, n_rows());
  return Matrix<Tp, 2>(this->M_elem, start, size, stride);
}

template <class Tp>
inline GsliceMatrix<Tp, 2> Matrix<Tp, 3>::row(uword i) {
  const uword start = i;
  INIT_ARR2E(size, n_slices(), n_cols());
  INIT_ARR2E(stride, M_nrxnc, n_rows());
  return GsliceMatrix<Tp, 2>(this->M_elem, start, size, stride);
}

template <class Tp>
inline Matrix<Tp, 2> Matrix<Tp, 3>::col(uword j) const {
  const uword start = j * n_rows();
  INIT_ARR2E(size, n_slices(), n_rows());
  INIT_ARR2E(stride, M_nrxnc, 1);
  return Matrix<Tp, 2>(this->M_elem, start, size, stride);
}

template <class Tp>
inline GsliceMatrix<Tp, 2> Matrix<Tp, 3>::col(uword j) {
  const uword start = j * n_rows();
  INIT_ARR2E(size, n_slices(), n_rows());
  INIT_ARR2E(stride, M_nrxnc, 1);
  return GsliceMatrix<Tp, 2>(this->M_elem, start, size, stride);
}

template <class Tp>
inline Matrix<Tp, 2> Matrix<Tp, 3>::slice(uword k) const {
  const uword start = k * M_nrxnc;
  INIT_ARR2E(size, n_cols(), n_rows());
  INIT_ARR2E(stride, n_rows(), 1);
  return Matrix<Tp, 2>(this->M_elem, start, size, stride);
}

template <class Tp>
inline GsliceMatrix<Tp, 2> Matrix<Tp, 3>::slice(uword k) {
  const uword start = k * M_nrxnc;
  INIT_ARR2E(size, n_cols(), n_rows());
  INIT_ARR2E(stride, n_rows(), 1);
  return GsliceMatrix<Tp, 2>(this->M_elem, start, size, stride);
}

// Matrix member functions dealing with bool_array and MaskMatrix

template <class Tp>
inline Matrix<Tp, 1> Matrix<Tp, 1>::operator()(const bool_array& ba) const {
  return Matrix<Tp, 1>(this->M_elem, ba);
}

template <class Tp>
inline MaskMatrix<Tp> Matrix<Tp, 1>::operator()(const bool_array& ba) {
  return MaskMatrix<Tp>(this->M_elem, ba);
}

template <class Tp>
inline Matrix<Tp, 1> Matrix<Tp, 2>::operator()(const bool_array& ba) const {
  return Matrix<Tp, 1>(this->M_elem, ba);
}

template <class Tp>
inline MaskMatrix<Tp> Matrix<Tp, 2>::operator()(const bool_array& ba) {
  return MaskMatrix<Tp>(this->M_elem, ba);
}

template <class Tp>
inline Matrix<Tp, 1> Matrix<Tp, 3>::operator()(const bool_array& ba) const {
  return Matrix<Tp, 1>(this->M_elem, ba);
}

template <class Tp>
inline MaskMatrix<Tp> Matrix<Tp, 3>::operator()(const bool_array& ba) {
  return MaskMatrix<Tp>(this->M_elem, ba);
}

// Matrix member functions dealing with index_array and IndirectMatrix

template <class Tp>
inline Matrix<Tp, 1> Matrix<Tp, 1>::operator()(const index_array& ia) const {
  INIT_ARR1E(dims, ia.size());
  return Matrix<Tp, 1>(this->M_elem, ia, dims);
}

template <class Tp>
inline IndirectMatrix<Tp, 1> Matrix<Tp, 1>::operator()(const index_array& ia) {
  INIT_ARR1E(dims, ia.size());
  return IndirectMatrix<Tp, 1>(this->M_elem, ia, dims);
}
template <class Tp>
inline Matrix<Tp, 1> Matrix<Tp, 2>::operator()(const index_array& ia) const {
  INIT_ARR1E(dims, ia.size());
  return Matrix<Tp, 1>(this->M_elem, ia, dims);
}

template <class Tp>
inline IndirectMatrix<Tp, 1> Matrix<Tp, 2>::operator()(const index_array& ia) {
  INIT_ARR1E(dims, ia.size());
  return IndirectMatrix<Tp, 1>(this->M_elem, ia, dims);
}

template <class Tp>
inline Matrix<Tp, 1> Matrix<Tp, 3>::operator()(const index_array& ia) const {
  INIT_ARR1E(dims, ia.size());
  return Matrix<Tp, 1>(this->M_elem, ia, dims);
}

template <class Tp>
inline IndirectMatrix<Tp, 1> Matrix<Tp, 3>::operator()(const index_array& ia) {
  INIT_ARR1E(dims, ia.size());
  return IndirectMatrix<Tp, 1>(this->M_elem, ia, dims);
}

template <class Tp>
inline Matrix<Tp, 2> Matrix<Tp, 2>::operator()(const index_array& ia1,
                                               const index_array& ia2) const {
  INIT_ARR2E(dims, ia1.size(), ia2.size());
  index_array ia(dims[0] * dims[1]);
  uword idx = 0;
  for (uword j = 0; j < ia2.size(); ++j) {
    for (uword i = 0; i < ia1.size(); ++i) {
      ia[idx++] = sub2ind(ia1[i], ia2[j]);
    }
  }
  return Matrix<Tp, 2>(this->M_elem, ia, dims);
}

template <class Tp>
inline IndirectMatrix<Tp, 2> Matrix<Tp, 2>::operator()(const index_array& ia1,
                                                       const index_array& ia2) {
  INIT_ARR2E(dims, ia1.size(), ia2.size());
  index_array ia(dims[0] * dims[1]);
  uword idx = 0;
  for (uword j = 0; j < ia2.size(); ++j) {
    for (uword i = 0; i < ia1.size(); ++i) {
      ia[idx++] = sub2ind(ia1[i], ia2[j]);
    }
  }
  return IndirectMatrix<Tp, 2>(this->M_elem, ia, dims);
}

template <class Tp>
inline Matrix<Tp, 3> Matrix<Tp, 3>::operator()(const index_array& ia1,
                                               const index_array& ia2,
                                               const index_array& ia3) const {
  INIT_ARR3E(dims, ia1.size(), ia2.size(), ia3.size());
  index_array ia(dims[0] * dims[1] * dims[2]);
  uword idx = 0;
  for (uword k = 0; k < ia3.size(); ++k) {
    for (uword j = 0; j < ia2.size(); ++j) {
      for (uword i = 0; i < ia1.size(); ++i) {
        ia[idx++] = sub2ind(ia1[i], ia2[j], ia3[k]);
      }
    }
  }
  return Matrix<Tp, 3>(this->M_elem, ia, dims);
}

template <class Tp>
inline IndirectMatrix<Tp, 3> Matrix<Tp, 3>::operator()(const index_array& ia1,
                                                       const index_array& ia2,
                                                       const index_array& ia3) {
  INIT_ARR3E(dims, ia1.size(), ia2.size(), ia3.size());
  index_array ia(dims[0] * dims[1] * dims[2]);
  uword idx = 0;
  for (uword k = 0; k < ia3.size(); ++k) {
    for (uword j = 0; j < ia2.size(); ++j) {
      for (uword i = 0; i < ia1.size(); ++i) {
        ia[idx++] = sub2ind(ia1[i], ia2[j], ia3[k]);
      }
    }
  }
  return IndirectMatrix<Tp, 3>(this->M_elem, ia, dims);
}

//----------------------------------------------------------------------
// Matrix / Sub-Matrix non-member functions.

// Binary arithmetic operations between two Matrix.

// clang-format off
template <class Tp, uword Size> void check_size(const    Matrix<Tp, Size>& x, const    Matrix<Tp, Size>& y) { matrix_assert(all_equal(x.size(), y.size()), "dimension mismatch"); }
template <class Tp, uword Size> void check_size(const    Matrix<Tp, Size>& x, const SubMatrix<Tp, Size>& y) { matrix_assert(all_equal(x.size(), y.size()), "dimension mismatch"); }
template <class Tp, uword Size> void check_size(const SubMatrix<Tp, Size>& x, const    Matrix<Tp, Size>& y) { matrix_assert(all_equal(x.size(), y.size()), "dimension mismatch"); }
template <class Tp, uword Size> void check_size(const SubMatrix<Tp, Size>& x, const SubMatrix<Tp, Size>& y) { matrix_assert(all_equal(x.size(), y.size()), "dimension mismatch"); }

template <class Tp, uword Size> inline Matrix<Tp, Size> operator+(const    Matrix<Tp, Size>& x, const    Matrix<Tp, Size>& y) { check_size(x, y); return Matrix<Tp, Size>(x.elem() + y.elem(), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator+(const SubMatrix<Tp, Size>& x, const    Matrix<Tp, Size>& y) { check_size(x, y); return Matrix<Tp, Size>(x.elem() + y.elem(), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator+(const    Matrix<Tp, Size>& x, const SubMatrix<Tp, Size>& y) { check_size(x, y); return Matrix<Tp, Size>(x.elem() + y.elem(), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator+(const SubMatrix<Tp, Size>& x, const SubMatrix<Tp, Size>& y) { check_size(x, y); return Matrix<Tp, Size>(x.elem() + y.elem(), x.size()); }

template <class Tp, uword Size> inline Matrix<Tp, Size> operator-(const    Matrix<Tp, Size>& x, const    Matrix<Tp, Size>& y) { check_size(x, y); return Matrix<Tp, Size>(x.elem() - y.elem(), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator-(const SubMatrix<Tp, Size>& x, const    Matrix<Tp, Size>& y) { check_size(x, y); return Matrix<Tp, Size>(x.elem() - y.elem(), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator-(const    Matrix<Tp, Size>& x, const SubMatrix<Tp, Size>& y) { check_size(x, y); return Matrix<Tp, Size>(x.elem() - y.elem(), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator-(const SubMatrix<Tp, Size>& x, const SubMatrix<Tp, Size>& y) { check_size(x, y); return Matrix<Tp, Size>(x.elem() - y.elem(), x.size()); }

template <class Tp, uword Size> inline Matrix<Tp, Size> operator*(const    Matrix<Tp, Size>& x, const    Matrix<Tp, Size>& y) { check_size(x, y); return Matrix<Tp, Size>(x.elem() * y.elem(), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator*(const SubMatrix<Tp, Size>& x, const    Matrix<Tp, Size>& y) { check_size(x, y); return Matrix<Tp, Size>(x.elem() * y.elem(), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator*(const    Matrix<Tp, Size>& x, const SubMatrix<Tp, Size>& y) { check_size(x, y); return Matrix<Tp, Size>(x.elem() * y.elem(), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator*(const SubMatrix<Tp, Size>& x, const SubMatrix<Tp, Size>& y) { check_size(x, y); return Matrix<Tp, Size>(x.elem() * y.elem(), x.size()); }

template <class Tp, uword Size> inline Matrix<Tp, Size> operator/(const    Matrix<Tp, Size>& x, const    Matrix<Tp, Size>& y) { check_size(x, y); return Matrix<Tp, Size>(x.elem() / y.elem(), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator/(const SubMatrix<Tp, Size>& x, const    Matrix<Tp, Size>& y) { check_size(x, y); return Matrix<Tp, Size>(x.elem() / y.elem(), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator/(const    Matrix<Tp, Size>& x, const SubMatrix<Tp, Size>& y) { check_size(x, y); return Matrix<Tp, Size>(x.elem() / y.elem(), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator/(const SubMatrix<Tp, Size>& x, const SubMatrix<Tp, Size>& y) { check_size(x, y); return Matrix<Tp, Size>(x.elem() / y.elem(), x.size()); }

template <class Tp, uword Size> inline Matrix<Tp, Size> operator%(const    Matrix<Tp, Size>& x, const    Matrix<Tp, Size>& y) { check_size(x, y); return Matrix<Tp, Size>(x.elem() % y.elem(), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator%(const SubMatrix<Tp, Size>& x, const    Matrix<Tp, Size>& y) { check_size(x, y); return Matrix<Tp, Size>(x.elem() % y.elem(), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator%(const    Matrix<Tp, Size>& x, const SubMatrix<Tp, Size>& y) { check_size(x, y); return Matrix<Tp, Size>(x.elem() % y.elem(), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator%(const SubMatrix<Tp, Size>& x, const SubMatrix<Tp, Size>& y) { check_size(x, y); return Matrix<Tp, Size>(x.elem() % y.elem(), x.size()); }

template <class Tp, uword Size> inline Matrix<Tp, Size> operator&(const    Matrix<Tp, Size>& x, const    Matrix<Tp, Size>& y) { check_size(x, y); return Matrix<Tp, Size>(x.elem() & y.elem(), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator&(const SubMatrix<Tp, Size>& x, const    Matrix<Tp, Size>& y) { check_size(x, y); return Matrix<Tp, Size>(x.elem() & y.elem(), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator&(const    Matrix<Tp, Size>& x, const SubMatrix<Tp, Size>& y) { check_size(x, y); return Matrix<Tp, Size>(x.elem() & y.elem(), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator&(const SubMatrix<Tp, Size>& x, const SubMatrix<Tp, Size>& y) { check_size(x, y); return Matrix<Tp, Size>(x.elem() & y.elem(), x.size()); }

template <class Tp, uword Size> inline Matrix<Tp, Size> operator|(const    Matrix<Tp, Size>& x, const    Matrix<Tp, Size>& y) { check_size(x, y); return Matrix<Tp, Size>(x.elem() | y.elem(), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator|(const SubMatrix<Tp, Size>& x, const    Matrix<Tp, Size>& y) { check_size(x, y); return Matrix<Tp, Size>(x.elem() | y.elem(), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator|(const    Matrix<Tp, Size>& x, const SubMatrix<Tp, Size>& y) { check_size(x, y); return Matrix<Tp, Size>(x.elem() | y.elem(), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator|(const SubMatrix<Tp, Size>& x, const SubMatrix<Tp, Size>& y) { check_size(x, y); return Matrix<Tp, Size>(x.elem() | y.elem(), x.size()); }

template <class Tp, uword Size> inline Matrix<Tp, Size> operator^(const    Matrix<Tp, Size>& x, const    Matrix<Tp, Size>& y) { check_size(x, y); return Matrix<Tp, Size>(x.elem() ^ y.elem(), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator^(const SubMatrix<Tp, Size>& x, const    Matrix<Tp, Size>& y) { check_size(x, y); return Matrix<Tp, Size>(x.elem() ^ y.elem(), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator^(const    Matrix<Tp, Size>& x, const SubMatrix<Tp, Size>& y) { check_size(x, y); return Matrix<Tp, Size>(x.elem() ^ y.elem(), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator^(const SubMatrix<Tp, Size>& x, const SubMatrix<Tp, Size>& y) { check_size(x, y); return Matrix<Tp, Size>(x.elem() ^ y.elem(), x.size()); }

template <class Tp, uword Size> inline Matrix<Tp, Size> operator<<(const    Matrix<Tp, Size>& x, const    Matrix<Tp, Size>& y) { check_size(x, y); return Matrix<Tp, Size>(x.elem() << y.elem(), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator<<(const SubMatrix<Tp, Size>& x, const    Matrix<Tp, Size>& y) { check_size(x, y); return Matrix<Tp, Size>(x.elem() << y.elem(), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator<<(const    Matrix<Tp, Size>& x, const SubMatrix<Tp, Size>& y) { check_size(x, y); return Matrix<Tp, Size>(x.elem() << y.elem(), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator<<(const SubMatrix<Tp, Size>& x, const SubMatrix<Tp, Size>& y) { check_size(x, y); return Matrix<Tp, Size>(x.elem() << y.elem(), x.size()); }

template <class Tp, uword Size> inline Matrix<Tp, Size> operator>>(const    Matrix<Tp, Size>& x, const    Matrix<Tp, Size>& y) { check_size(x, y); return Matrix<Tp, Size>(x.elem() >> y.elem(), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator>>(const SubMatrix<Tp, Size>& x, const    Matrix<Tp, Size>& y) { check_size(x, y); return Matrix<Tp, Size>(x.elem() >> y.elem(), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator>>(const    Matrix<Tp, Size>& x, const SubMatrix<Tp, Size>& y) { check_size(x, y); return Matrix<Tp, Size>(x.elem() >> y.elem(), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator>>(const SubMatrix<Tp, Size>& x, const SubMatrix<Tp, Size>& y) { check_size(x, y); return Matrix<Tp, Size>(x.elem() >> y.elem(), x.size()); }

template <class Tp, uword Size> inline bool_array operator&&(const    Matrix<Tp, Size>& x, const    Matrix<Tp, Size>& y) { check_size(x, y); return x.elem() && y.elem(); }
template <class Tp, uword Size> inline bool_array operator&&(const SubMatrix<Tp, Size>& x, const    Matrix<Tp, Size>& y) { check_size(x, y); return x.elem() && y.elem(); }
template <class Tp, uword Size> inline bool_array operator&&(const    Matrix<Tp, Size>& x, const SubMatrix<Tp, Size>& y) { check_size(x, y); return x.elem() && y.elem(); }
template <class Tp, uword Size> inline bool_array operator&&(const SubMatrix<Tp, Size>& x, const SubMatrix<Tp, Size>& y) { check_size(x, y); return x.elem() && y.elem(); }

template <class Tp, uword Size> inline bool_array operator||(const    Matrix<Tp, Size>& x, const    Matrix<Tp, Size>& y) { check_size(x, y); return x.elem() || y.elem(); }
template <class Tp, uword Size> inline bool_array operator||(const SubMatrix<Tp, Size>& x, const    Matrix<Tp, Size>& y) { check_size(x, y); return x.elem() || y.elem(); }
template <class Tp, uword Size> inline bool_array operator||(const    Matrix<Tp, Size>& x, const SubMatrix<Tp, Size>& y) { check_size(x, y); return x.elem() || y.elem(); }
template <class Tp, uword Size> inline bool_array operator||(const SubMatrix<Tp, Size>& x, const SubMatrix<Tp, Size>& y) { check_size(x, y); return x.elem() || y.elem(); }
// clang-format on

// Binary arithmetic operations between an array and a scalar.

// clang-format off
template <class Tp, uword Size> inline Matrix<Tp, Size> operator+(const    Matrix<Tp, Size>& x, const Tp& c) { return Matrix<Tp, Size>(x.elem() + c, x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator+(const SubMatrix<Tp, Size>& x, const Tp& c) { return Matrix<Tp, Size>(x.elem() + c, x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator+(const Tp& c, const    Matrix<Tp, Size>& x) { return Matrix<Tp, Size>(c + x.elem(), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator+(const Tp& c, const SubMatrix<Tp, Size>& x) { return Matrix<Tp, Size>(c + x.elem(), x.size());}

template <class Tp, uword Size> inline Matrix<Tp, Size> operator-(const    Matrix<Tp, Size>& x, const Tp& c) { return Matrix<Tp, Size>(x.elem() - c, x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator-(const SubMatrix<Tp, Size>& x, const Tp& c) { return Matrix<Tp, Size>(x.elem() - c, x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator-(const Tp& c, const    Matrix<Tp, Size>& x) { return Matrix<Tp, Size>(c - x.elem(), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator-(const Tp& c, const SubMatrix<Tp, Size>& x) { return Matrix<Tp, Size>(c - x.elem(), x.size()); }

template <class Tp, uword Size> inline Matrix<Tp, Size> operator*(const    Matrix<Tp, Size>& x, const Tp& c) { return Matrix<Tp, Size>(x.elem() * c, x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator*(const SubMatrix<Tp, Size>& x, const Tp& c) { return Matrix<Tp, Size>(x.elem() * c, x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator*(const Tp& c, const    Matrix<Tp, Size>& x) { return Matrix<Tp, Size>(c * x.elem(), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator*(const Tp& c, const SubMatrix<Tp, Size>& x) { return Matrix<Tp, Size>(c * x.elem(), x.size()); }

template <class Tp, uword Size> inline Matrix<Tp, Size> operator/(const    Matrix<Tp, Size>& x, const Tp& c) { return Matrix<Tp, Size>(x.elem() / c, x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator/(const SubMatrix<Tp, Size>& x, const Tp& c) { return Matrix<Tp, Size>(x.elem() / c, x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator/(const Tp& c, const    Matrix<Tp, Size>& x) { return Matrix<Tp, Size>(c / x.elem(), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator/(const Tp& c, const SubMatrix<Tp, Size>& x) { return Matrix<Tp, Size>(c / x.elem(), x.size()); }

template <class Tp, uword Size> inline Matrix<Tp, Size> operator%(const    Matrix<Tp, Size>& x, const Tp& c) { return Matrix<Tp, Size>(x.elem() % c, x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator%(const SubMatrix<Tp, Size>& x, const Tp& c) { return Matrix<Tp, Size>(x.elem() % c, x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator%(const Tp& c, const    Matrix<Tp, Size>& x) { return Matrix<Tp, Size>(c % x.elem(), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator%(const Tp& c, const SubMatrix<Tp, Size>& x) { return Matrix<Tp, Size>(c % x.elem(), x.size()); }

template <class Tp, uword Size> inline Matrix<Tp, Size> operator&(const    Matrix<Tp, Size>& x, const Tp& c) { return Matrix<Tp, Size>(x.elem() & c, x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator&(const SubMatrix<Tp, Size>& x, const Tp& c) { return Matrix<Tp, Size>(x.elem() & c, x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator&(const Tp& c, const    Matrix<Tp, Size>& x) { return Matrix<Tp, Size>(c & x.elem(), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator&(const Tp& c, const SubMatrix<Tp, Size>& x) { return Matrix<Tp, Size>(c & x.elem(), x.size()); }

template <class Tp, uword Size> inline Matrix<Tp, Size> operator|(const    Matrix<Tp, Size>& x, const Tp& c) { return Matrix<Tp, Size>(x.elem() | c, x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator|(const SubMatrix<Tp, Size>& x, const Tp& c) { return Matrix<Tp, Size>(x.elem() | c, x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator|(const Tp& c, const    Matrix<Tp, Size>& x) { return Matrix<Tp, Size>(c | x.elem(), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator|(const Tp& c, const SubMatrix<Tp, Size>& x) { return Matrix<Tp, Size>(c | x.elem(), x.size()); }

template <class Tp, uword Size> inline Matrix<Tp, Size> operator^(const    Matrix<Tp, Size>& x, const Tp& c) { return Matrix<Tp, Size>(x.elem() ^ c, x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator^(const SubMatrix<Tp, Size>& x, const Tp& c) { return Matrix<Tp, Size>(x.elem() ^ c, x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator^(const Tp& c, const    Matrix<Tp, Size>& x) { return Matrix<Tp, Size>(c ^ x.elem(), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator^(const Tp& c, const SubMatrix<Tp, Size>& x) { return Matrix<Tp, Size>(c ^ x.elem(), x.size()); }

template <class Tp, uword Size> inline Matrix<Tp, Size> operator<<(const    Matrix<Tp, Size>& x, const Tp& c) { return Matrix<Tp, Size>(x.elem() << c, x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator<<(const SubMatrix<Tp, Size>& x, const Tp& c) { return Matrix<Tp, Size>(x.elem() << c, x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator<<(const Tp& c, const    Matrix<Tp, Size>& x) { return Matrix<Tp, Size>(c << x.elem(), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator<<(const Tp& c, const SubMatrix<Tp, Size>& x) { return Matrix<Tp, Size>(c << x.elem(), x.size()); }

template <class Tp, uword Size> inline Matrix<Tp, Size> operator>>(const    Matrix<Tp, Size>& x, const Tp& c) { return Matrix<Tp, Size>(x.elem() >> c, x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator>>(const SubMatrix<Tp, Size>& x, const Tp& c) { return Matrix<Tp, Size>(x.elem() >> c, x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator>>(const Tp& c, const    Matrix<Tp, Size>& x) { return Matrix<Tp, Size>(c >> x.elem(), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> operator>>(const Tp& c, const SubMatrix<Tp, Size>& x) { return Matrix<Tp, Size>(c >> x.elem(), x.size()); }

template <class Tp, uword Size> inline bool_array operator&&(const    Matrix<Tp, Size>& x, const Tp& c) { return Matrix<Tp, Size>(x.elem() && c, x.size()); }
template <class Tp, uword Size> inline bool_array operator&&(const SubMatrix<Tp, Size>& x, const Tp& c) { return Matrix<Tp, Size>(x.elem() && c, x.size()); }
template <class Tp, uword Size> inline bool_array operator&&(const Tp& c, const    Matrix<Tp, Size>& x) { return Matrix<Tp, Size>(c && x.elem(), x.size()); }
template <class Tp, uword Size> inline bool_array operator&&(const Tp& c, const SubMatrix<Tp, Size>& x) { return Matrix<Tp, Size>(c && x.elem(), x.size()); }

template <class Tp, uword Size> inline bool_array operator||(const    Matrix<Tp, Size>& x, const Tp& c) { return Matrix<Tp, Size>(x.elem() || c, x.size()); }
template <class Tp, uword Size> inline bool_array operator||(const SubMatrix<Tp, Size>& x, const Tp& c) { return Matrix<Tp, Size>(x.elem() || c, x.size()); }
template <class Tp, uword Size> inline bool_array operator||(const Tp& c, const    Matrix<Tp, Size>& x) { return Matrix<Tp, Size>(c || x.elem(), x.size()); }
template <class Tp, uword Size> inline bool_array operator||(const Tp& c, const SubMatrix<Tp, Size>& x) { return Matrix<Tp, Size>(c || x.elem(), x.size()); }
// clang-format on

// Binary logical operations between two Matrices.

// clang-format off
template <class Tp, uword Size> inline bool_array operator==(const    Matrix<Tp, Size>& x, const    Matrix<Tp, Size>& y) { check_size(x, y); return x.elem() == y.elem(); }
template <class Tp, uword Size> inline bool_array operator==(const SubMatrix<Tp, Size>& x, const    Matrix<Tp, Size>& y) { check_size(x, y); return x.elem() == y.elem(); }
template <class Tp, uword Size> inline bool_array operator==(const    Matrix<Tp, Size>& x, const SubMatrix<Tp, Size>& y) { check_size(x, y); return x.elem() == y.elem(); }
template <class Tp, uword Size> inline bool_array operator==(const SubMatrix<Tp, Size>& x, const SubMatrix<Tp, Size>& y) { check_size(x, y); return x.elem() == y.elem(); }

template <class Tp, uword Size> inline bool_array operator!=(const    Matrix<Tp, Size>& x, const    Matrix<Tp, Size>& y) { check_size(x, y); return x.elem() != y.elem(); }
template <class Tp, uword Size> inline bool_array operator!=(const SubMatrix<Tp, Size>& x, const    Matrix<Tp, Size>& y) { check_size(x, y); return x.elem() != y.elem(); }
template <class Tp, uword Size> inline bool_array operator!=(const    Matrix<Tp, Size>& x, const SubMatrix<Tp, Size>& y) { check_size(x, y); return x.elem() != y.elem(); }
template <class Tp, uword Size> inline bool_array operator!=(const SubMatrix<Tp, Size>& x, const SubMatrix<Tp, Size>& y) { check_size(x, y); return x.elem() != y.elem(); }

template <class Tp, uword Size> inline bool_array operator< (const    Matrix<Tp, Size>& x, const    Matrix<Tp, Size>& y) { check_size(x, y); return x.elem() <  y.elem(); }
template <class Tp, uword Size> inline bool_array operator< (const SubMatrix<Tp, Size>& x, const    Matrix<Tp, Size>& y) { check_size(x, y); return x.elem() <  y.elem(); }
template <class Tp, uword Size> inline bool_array operator< (const    Matrix<Tp, Size>& x, const SubMatrix<Tp, Size>& y) { check_size(x, y); return x.elem() <  y.elem(); }
template <class Tp, uword Size> inline bool_array operator< (const SubMatrix<Tp, Size>& x, const SubMatrix<Tp, Size>& y) { check_size(x, y); return x.elem() <  y.elem(); }

template <class Tp, uword Size> inline bool_array operator<=(const    Matrix<Tp, Size>& x, const    Matrix<Tp, Size>& y) { check_size(x, y); return x.elem() <= y.elem(); }
template <class Tp, uword Size> inline bool_array operator<=(const SubMatrix<Tp, Size>& x, const    Matrix<Tp, Size>& y) { check_size(x, y); return x.elem() <= y.elem(); }
template <class Tp, uword Size> inline bool_array operator<=(const    Matrix<Tp, Size>& x, const SubMatrix<Tp, Size>& y) { check_size(x, y); return x.elem() <= y.elem(); }
template <class Tp, uword Size> inline bool_array operator<=(const SubMatrix<Tp, Size>& x, const SubMatrix<Tp, Size>& y) { check_size(x, y); return x.elem() <= y.elem(); }

template <class Tp, uword Size> inline bool_array operator> (const    Matrix<Tp, Size>& x, const    Matrix<Tp, Size>& y) { check_size(x, y); return x.elem() >  y.elem(); }
template <class Tp, uword Size> inline bool_array operator> (const SubMatrix<Tp, Size>& x, const    Matrix<Tp, Size>& y) { check_size(x, y); return x.elem() >  y.elem(); }
template <class Tp, uword Size> inline bool_array operator> (const    Matrix<Tp, Size>& x, const SubMatrix<Tp, Size>& y) { check_size(x, y); return x.elem() >  y.elem(); }
template <class Tp, uword Size> inline bool_array operator> (const SubMatrix<Tp, Size>& x, const SubMatrix<Tp, Size>& y) { check_size(x, y); return x.elem() >  y.elem(); }

template <class Tp, uword Size> inline bool_array operator>=(const    Matrix<Tp, Size>& x, const    Matrix<Tp, Size>& y) { check_size(x, y); return x.elem() >= y.elem(); }
template <class Tp, uword Size> inline bool_array operator>=(const SubMatrix<Tp, Size>& x, const    Matrix<Tp, Size>& y) { check_size(x, y); return x.elem() >= y.elem(); }
template <class Tp, uword Size> inline bool_array operator>=(const    Matrix<Tp, Size>& x, const SubMatrix<Tp, Size>& y) { check_size(x, y); return x.elem() >= y.elem(); }
template <class Tp, uword Size> inline bool_array operator>=(const SubMatrix<Tp, Size>& x, const SubMatrix<Tp, Size>& y) { check_size(x, y); return x.elem() >= y.elem(); }
// clang-format on

// Logical operations between a Matrix and a scalar.
// clang-format off
template <class Tp, uword Size> inline bool_array operator==(const    Matrix<Tp, Size>& x, const Tp& c) { return x.elem() == c; }
template <class Tp, uword Size> inline bool_array operator==(const SubMatrix<Tp, Size>& x, const Tp& c) { return x.elem() == c; }
template <class Tp, uword Size> inline bool_array operator==(const Tp& c, const    Matrix<Tp, Size>& x) { return c == x.elem(); }
template <class Tp, uword Size> inline bool_array operator==(const Tp& c, const SubMatrix<Tp, Size>& x) { return c == x.elem(); }

template <class Tp, uword Size> inline bool_array operator!=(const    Matrix<Tp, Size>& x, const Tp& c) { return x.elem() != c; }
template <class Tp, uword Size> inline bool_array operator!=(const SubMatrix<Tp, Size>& x, const Tp& c) { return x.elem() != c; }
template <class Tp, uword Size> inline bool_array operator!=(const Tp& c, const    Matrix<Tp, Size>& x) { return c != x.elem(); }
template <class Tp, uword Size> inline bool_array operator!=(const Tp& c, const SubMatrix<Tp, Size>& x) { return c != x.elem(); }

template <class Tp, uword Size> inline bool_array operator< (const    Matrix<Tp, Size>& x, const Tp& c) { return x.elem() < c; }
template <class Tp, uword Size> inline bool_array operator< (const SubMatrix<Tp, Size>& x, const Tp& c) { return x.elem() < c; }
template <class Tp, uword Size> inline bool_array operator< (const Tp& c, const    Matrix<Tp, Size>& x) { return c < x.elem(); }
template <class Tp, uword Size> inline bool_array operator< (const Tp& c, const SubMatrix<Tp, Size>& x) { return c < x.elem(); }

template <class Tp, uword Size> inline bool_array operator<=(const    Matrix<Tp, Size>& x, const Tp& c) { return x.elem() <= c; }
template <class Tp, uword Size> inline bool_array operator<=(const SubMatrix<Tp, Size>& x, const Tp& c) { return x.elem() <= c; }
template <class Tp, uword Size> inline bool_array operator<=(const Tp& c, const    Matrix<Tp, Size>& x) { return c <= x.elem(); }
template <class Tp, uword Size> inline bool_array operator<=(const Tp& c, const SubMatrix<Tp, Size>& x) { return c <= x.elem(); }

template <class Tp, uword Size> inline bool_array operator> (const    Matrix<Tp, Size>& x, const Tp& c) { return x.elem() > c; }
template <class Tp, uword Size> inline bool_array operator> (const SubMatrix<Tp, Size>& x, const Tp& c) { return x.elem() > c; }
template <class Tp, uword Size> inline bool_array operator> (const Tp& c, const    Matrix<Tp, Size>& x) { return c > x.elem(); }
template <class Tp, uword Size> inline bool_array operator> (const Tp& c, const SubMatrix<Tp, Size>& x) { return c > x.elem(); }

template <class Tp, uword Size> inline bool_array operator>=(const    Matrix<Tp, Size>& x, const Tp& c) { return x.elem() >= c; }
template <class Tp, uword Size> inline bool_array operator>=(const SubMatrix<Tp, Size>& x, const Tp& c) { return x.elem() >= c; }
template <class Tp, uword Size> inline bool_array operator>=(const Tp& c, const    Matrix<Tp, Size>& x) { return c >= x.elem(); }
template <class Tp, uword Size> inline bool_array operator>=(const Tp& c, const SubMatrix<Tp, Size>& x) { return c >= x.elem(); }
// clang-format on

// Matrix "transcendentals" (the list includes abs and sqrt, which,
// of course, are not transcendental).

// clang-format off
template <class Tp, uword Size> inline Matrix<Tp, Size> abs(const    Matrix<Tp, Size>& x) { return Matrix<Tp, Size>(std::abs(x.elem()), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> abs(const SubMatrix<Tp, Size>& x) { return Matrix<Tp, Size>(std::abs(x.elem()), x.size()); }

template <class Tp, uword Size> inline Matrix<Tp, Size> exp(const    Matrix<Tp, Size>& x) { return Matrix<Tp, Size>(std::exp(x.elem()), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> exp(const SubMatrix<Tp, Size>& x) { return Matrix<Tp, Size>(std::exp(x.elem()), x.size()); }

template <class Tp, uword Size> inline Matrix<Tp, Size> log(const    Matrix<Tp, Size>& x) { return Matrix<Tp, Size>(std::log(x.elem()), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> log(const SubMatrix<Tp, Size>& x) { return Matrix<Tp, Size>(std::log(x.elem()), x.size()); }

template <class Tp, uword Size> inline Matrix<Tp, Size> log10(const    Matrix<Tp, Size>& x) { return Matrix<Tp, Size>(std::log10(x.elem()), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> log10(const SubMatrix<Tp, Size>& x) { return Matrix<Tp, Size>(std::log10(x.elem()), x.size()); }

template <class Tp, uword Size> inline Matrix<Tp, Size> pow(const    Matrix<Tp, Size>&  x, const    Matrix<Tp, Size>&  y) { check_size(x, y); return Matrix<Tp, Size>(std::pow( x.elem(),  y.elem()), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> pow(const SubMatrix<Tp, Size>&  x, const    Matrix<Tp, Size>&  y) { check_size(x, y); return Matrix<Tp, Size>(std::pow( x.elem(),  y.elem()), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> pow(const    Matrix<Tp, Size>&  x, const SubMatrix<Tp, Size>&  y) { check_size(x, y); return Matrix<Tp, Size>(std::pow( x.elem(),  y.elem()), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> pow(const SubMatrix<Tp, Size>&  x, const SubMatrix<Tp, Size>&  y) { check_size(x, y); return Matrix<Tp, Size>(std::pow( x.elem(),  y.elem()), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> pow(const    Matrix<Tp, Size>&  x, const           Tp       & vy) {                   return Matrix<Tp, Size>(std::pow( x.elem(), vy       ), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> pow(const SubMatrix<Tp, Size>&  x, const           Tp       & vy) {                   return Matrix<Tp, Size>(std::pow( x.elem(), vy       ), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> pow(const           Tp       & vx, const    Matrix<Tp, Size>&  y) {                   return Matrix<Tp, Size>(std::pow(vx       ,  y.elem()), y.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> pow(const           Tp       & vx, const SubMatrix<Tp, Size>&  y) {                   return Matrix<Tp, Size>(std::pow(vx       ,  y.elem()), y.size()); }

template <class Tp, uword Size> inline Matrix<Tp, Size> sqrt(const    Matrix<Tp, Size>& x) { return Matrix<Tp, Size>(std::sqrt(x.elem()), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> sqrt(const SubMatrix<Tp, Size>& x) { return Matrix<Tp, Size>(std::sqrt(x.elem()), x.size()); }

template <class Tp, uword Size> inline Matrix<Tp, Size> sin(const    Matrix<Tp, Size>& x) { return Matrix<Tp, Size>(std::sin(x.elem()), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> sin(const SubMatrix<Tp, Size>& x) { return Matrix<Tp, Size>(std::sin(x.elem()), x.size()); }

template <class Tp, uword Size> inline Matrix<Tp, Size> cos(const    Matrix<Tp, Size>& x) { return Matrix<Tp, Size>(std::cos(x.elem()), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> cos(const SubMatrix<Tp, Size>& x) { return Matrix<Tp, Size>(std::cos(x.elem()), x.size()); }

template <class Tp, uword Size> inline Matrix<Tp, Size> tan(const    Matrix<Tp, Size>& x) { return Matrix<Tp, Size>(std::tan(x.elem()), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> tan(const SubMatrix<Tp, Size>& x) { return Matrix<Tp, Size>(std::tan(x.elem()), x.size()); }

template <class Tp, uword Size> inline Matrix<Tp, Size> asin(const    Matrix<Tp, Size>& x) { return Matrix<Tp, Size>(std::asin(x.elem()), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> asin(const SubMatrix<Tp, Size>& x) { return Matrix<Tp, Size>(std::asin(x.elem()), x.size()); }

template <class Tp, uword Size> inline Matrix<Tp, Size> acos(const    Matrix<Tp, Size>& x) { return Matrix<Tp, Size>(std::acos(x.elem()), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> acos(const SubMatrix<Tp, Size>& x) { return Matrix<Tp, Size>(std::acos(x.elem()), x.size()); }

template <class Tp, uword Size> inline Matrix<Tp, Size> atan(const    Matrix<Tp, Size>& x) { return Matrix<Tp, Size>(std::atan(x.elem()), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> atan(const SubMatrix<Tp, Size>& x) { return Matrix<Tp, Size>(std::atan(x.elem()), x.size()); }

template <class Tp, uword Size> inline Matrix<Tp, Size> atan2(const    Matrix<Tp, Size>&  y, const    Matrix<Tp, Size>&  x) { check_size(x, y); return Matrix<Tp, Size>(std::atan2( y.elem(),  x.elem()), y.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> atan2(const SubMatrix<Tp, Size>&  y, const    Matrix<Tp, Size>&  x) { check_size(x, y); return Matrix<Tp, Size>(std::atan2( y.elem(),  x.elem()), y.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> atan2(const    Matrix<Tp, Size>&  y, const SubMatrix<Tp, Size>&  x) { check_size(x, y); return Matrix<Tp, Size>(std::atan2( y.elem(),  x.elem()), y.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> atan2(const SubMatrix<Tp, Size>&  y, const SubMatrix<Tp, Size>&  x) { check_size(x, y); return Matrix<Tp, Size>(std::atan2( y.elem(),  x.elem()), y.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> atan2(const    Matrix<Tp, Size>&  y, const           Tp       & vx) {                   return Matrix<Tp, Size>(std::atan2( y.elem(), vx       ), y.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> atan2(const SubMatrix<Tp, Size>&  y, const           Tp       & vx) {                   return Matrix<Tp, Size>(std::atan2( y.elem(), vx       ), y.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> atan2(const           Tp       & vy, const    Matrix<Tp, Size>&  x) {                   return Matrix<Tp, Size>(std::atan2(vy       ,  x.elem()), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> atan2(const           Tp       & vy, const SubMatrix<Tp, Size>&  x) {                   return Matrix<Tp, Size>(std::atan2(vy       ,  x.elem()), x.size()); }

template <class Tp, uword Size> inline Matrix<Tp, Size> sinh(const    Matrix<Tp, Size>& x) { return Matrix<Tp, Size>(std::sinh(x.elem()), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> sinh(const SubMatrix<Tp, Size>& x) { return Matrix<Tp, Size>(std::sinh(x.elem()), x.size()); }

template <class Tp, uword Size> inline Matrix<Tp, Size> cosh(const    Matrix<Tp, Size>& x) { return Matrix<Tp, Size>(std::cosh(x.elem()), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> cosh(const SubMatrix<Tp, Size>& x) { return Matrix<Tp, Size>(std::cosh(x.elem()), x.size()); }

template <class Tp, uword Size> inline Matrix<Tp, Size> tanh(const    Matrix<Tp, Size>& x) { return Matrix<Tp, Size>(std::tanh(x.elem()), x.size()); }
template <class Tp, uword Size> inline Matrix<Tp, Size> tanh(const SubMatrix<Tp, Size>& x) { return Matrix<Tp, Size>(std::tanh(x.elem()), x.size()); }
// clang-format on

template <class Tp>
Tp dot(const Matrix<Tp, 1>& x, const Matrix<Tp, 1>& y) {
  return (x.elem() * y.elem()).sum();
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
