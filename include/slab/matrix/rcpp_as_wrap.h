//
// Copyright 2018-2019 The Statslabs Authors.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//

/// @file rcpp_as_wrap.h
/// @brief Data interchange between R and C++ for the Matrix template

inline void Rf_copyCubeByRow(SEXP target, SEXP source) {
  SEXP dims = Rf_getAttrib(source, R_DimSymbol);

  int dim0 = INTEGER(dims)[0];
  int dim1 = INTEGER(dims)[1];
  int dim2 = INTEGER(dims)[2];

  int s_idx = 0, t_idx = 0;
  for (int i = 0; i != dim0; ++i) {
    for (int j = 0; j != dim1; ++j) {
      for (int k = 0; k != dim2; ++k, ++t_idx) {
        s_idx = k * (dim0 * dim1) + j * dim0 + i;

        if (Rf_isReal(source))
          REAL(target)[t_idx] = REAL(source)[s_idx];
        else if (Rf_isInteger(source))
          INTEGER(target)[t_idx] = INTEGER(source)[s_idx];
      }
    }
  }
}

template <typename T, std::size_t N>
Matrix<T, N>::Matrix(SEXP s) {
  SEXP dims = Rf_getAttrib(s, R_DimSymbol);
  _SLAB_ASSERT((Rf_isNull(dims) && this->order() == 1) ||
                   Rf_length(dims) == this->order(),
               "Matrix(SEXP): unmatched dimensions.");

  int num_elem = Rf_length(s);
  int num_dims = 0;
  if (Rf_isNull(dims))
    num_dims = 1;
  else
    num_dims = Rf_length(dims);

  SEXP s2 = PROTECT(Rf_allocVector(TYPEOF(s), (R_xlen_t)num_elem));
  SEXP dims2 = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t)num_dims));

  std::array<std::size_t, N> exts;
  if (num_dims == 1) {
    exts[0] = INTEGER(dims2)[0] = num_elem;
  } else if (num_dims == 2) {
    exts[0] = INTEGER(dims2)[1] = INTEGER(dims)[0];
    exts[1] = INTEGER(dims2)[0] = INTEGER(dims)[1];
  } else {
    for (int i = 0; i <= N; ++i)
      exts[i] = INTEGER(dims2)[i] = INTEGER(dims)[N - i - 1];
  }
  this->desc_ = MatrixSlice<N>(exts);

  Rf_setAttrib(s2, R_DimSymbol, dims2);
  if (num_dims <= 2)
    Rf_copyMatrix(s2, s, TRUE);
  else if (num_dims == 3)
    Rf_copyCubeByRow(s2, s);

  elems_.reserve(num_elem);
  for (int i = 0; i != num_elem; ++i) {
    if (Rf_isReal(s))
      elems_.push_back(REAL(s2)[i]);
    else if (Rf_isInteger(s))
      elems_.push_back(INTEGER(s2)[i]);
    else
      _SLAB_ERROR("Matrix(SEXP): unsupported SEXP type");
  }

  UNPROTECT(2);
}

template <typename T, std::size_t N>
Matrix<T, N>::operator SEXP() {
  int num_elems = this->size();
  SEXP res = PROTECT(Rcpp::wrap(this->data(), this->data() + num_elems));
  SEXP res2 = PROTECT(Rf_allocVector(Rcpp::traits::r_sexptype_traits<T>::rtype,
                                     (R_xlen_t)num_elems));

  int num_dims = this->order();
  SEXP dim = PROTECT(Rf_allocVector(INTSXP, num_dims));
  int *idim = INTEGER(dim);
  if (num_dims <= 2) {
    for (int i = 0; i != num_dims; ++i) idim[i] = this->desc_.extents[i];
  } else {
    idim[0] = this->desc_.extents[num_dims - 2];
    idim[1] = this->desc_.extents[num_dims - 1];
    for (int i = 2; i != num_dims; ++i)
      idim[i] = this->desc_.extents[num_dims - i - 1];
  }
  Rf_setAttrib(res, R_DimSymbol, dim);
  Rf_setAttrib(res2, R_DimSymbol, dim);

  if (num_dims <= 2) Rf_copyMatrix(res2, res, TRUE);
  // TODO: else CopyMatrixByRow(res2, res); // for num_dims > 2

  UNPROTECT(3);

  return res2;
}

template <typename T, std::size_t N>
Matrix<T, N>::operator SEXP() const {
  int num_elems = this->size();
  SEXP res = PROTECT(Rcpp::wrap(this->data(), this->data() + num_elems));
  SEXP res2 = PROTECT(Rf_allocVector(Rcpp::traits::r_sexptype_traits<T>::rtype,
                                     (R_xlen_t)num_elems));

  int num_dims = this->order();
  SEXP dim = PROTECT(Rf_allocVector(INTSXP, num_dims));
  int *idim = INTEGER(dim);
  if (num_dims <= 2) {
    for (int i = 0; i != num_dims; ++i) idim[i] = this->desc_.extents[i];
  } else {
    idim[0] = this->desc_.extents[num_dims - 2];
    idim[1] = this->desc_.extents[num_dims - 1];
    for (int i = 2; i != num_dims; ++i)
      idim[i] = this->desc_.extents[num_dims - i - 1];
  }
  Rf_setAttrib(res, R_DimSymbol, dim);
  Rf_setAttrib(res2, R_DimSymbol, dim);

  if (num_dims <= 2) Rf_copyMatrix(res2, res, TRUE);
  // TODO: else CopyMatrixByRow(res2, res); // for num_dims > 2

  UNPROTECT(3);

  return res2;
}
