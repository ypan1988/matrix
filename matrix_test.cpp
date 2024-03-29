#include "matrix.h"

#include <cassert>
#include <iostream>
#include <string>
#include <valarray>

using namespace matrix_lib;

const double arr_1d[] = {1, 2, 3, 4};
const double arr_2d[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
const double arr_3d[] = {1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12,
                         13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24};
const double arr_3d_2[] = {1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12,
                           13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,
                           25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36};

const bool bool_arr_1d[4] = {false, true, false, true};
const bool bool_arr_2d[12] = {false, true,  false, true, false, false,
                              false, false, false, true, false, true};
const bool bool_arr_3d[24] = {false, false, false, false, false, false,
                              false, false, false, false, false, false,
                              false, true,  false, true,  false, false,
                              false, false, false, true,  false, true};

const uword index_arr_1d[2] = {1, 3};
const uword index_arr_2d[4] = {1, 3, 9, 11};
const uword index_arr_3d[4] = {13, 15, 21, 23};

const uword index_dim1[3] = {0, 1, 3};
const uword index_dim2[2] = {0, 2};
const uword index_dim3[1] = {1};

template <typename T>
void test_print(const Matrix<T, 1>& m, std::string msg = "") {
  if (!msg.empty()) std::cout << msg << std::endl;
  for (uword r = 0; r != m.n_elem(); ++r) {
    std::cout << m(r);
    m.is_column_vector ? std::cout << '\n' : std::cout << ' ';
  }
  if (!m.is_column_vector) std::cout << std::endl;
}

template <typename T>
void test_print(const Matrix<T, 2>& m, std::string msg = "") {
  if (!msg.empty()) std::cout << msg << std::endl;
  for (uword r = 0; r != m.n_rows(); ++r) {
    for (uword c = 0; c != m.n_cols(); ++c) {
      std::cout << m(r, c) << ' ';
    }
    std::cout << '\n';
  }
}

template <typename T>
void test_print(const Matrix<T, 3>& m, std::string msg = "") {
  if (!msg.empty()) std::cout << msg << std::endl;
  for (uword s = 0; s != m.n_slices(); ++s) {
    std::cout << "slice " << s << ":" << std::endl;
    for (uword r = 0; r != m.n_rows(); ++r) {
      for (uword c = 0; c != m.n_cols(); ++c) {
        std::cout << m(r, c, s) << ' ';
      }
      std::cout << '\n';
    }
  }
}

void macro_test_init_array(bool print = false) {
  std::cout << "[TEST]: 01. Init Index Array (i.e., valarray<size_t>)"
            << std::endl;
  INIT_ARR1E(idx_arr1, 1);
  assert(idx_arr1.size() == 1);
  assert(idx_arr1[0] == 1);

  INIT_ARR2E(idx_arr2, 1, 2);
  assert(idx_arr2.size() == 2);
  assert(idx_arr2[0] == 1);
  assert(idx_arr2[1] == 2);

  INIT_ARR3E(idx_arr3, 1, 2, 3);
  assert(idx_arr3.size() == 3);
  assert(idx_arr3[0] == 1);
  assert(idx_arr3[1] == 2);
  assert(idx_arr3[2] == 3);
}

// ----- A1: Testing Matrix<T, N> Constructions -----

void matrix_test_constructor_a1_01(bool print = false) {
  std::cout << "[TEST]: A1_01. Default constructor\n"
            << "        vec: vec()\n"
            << "        mat: mat()\n"
            << "        cube: cube()" << std::endl;

  Matrix<double, 1> mat1d;
  if (print) test_print(mat1d, "mat1d =");
  assert(mat1d.n_elem() == 0);
  assert(mat1d.n_rows() == 0);
  assert(mat1d.n_cols() == 0);

  Matrix<double, 2> mat2d;
  if (print) test_print(mat2d, "mat2d =");
  assert(mat2d.n_elem() == 0);
  assert(mat2d.n_rows() == 0);
  assert(mat2d.n_cols() == 0);

  Matrix<double, 3> mat3d;
  if (print) test_print(mat3d, "mat3d =");
  assert(mat3d.n_elem() == 0);
  assert(mat3d.n_rows() == 0);
  assert(mat3d.n_cols() == 0);
  assert(mat3d.n_slices() == 0);
}

void matrix_test_constructor_a1_02(bool print = false) {
  std::cout << "[TEST]: A1_02. Dimensions only\n"
            << "        vec: explicit vec(n_rows)\n"
            << "        mat: mat(n_rows, n_cols)\n"
            << "        cube: cube(n_rows, n_cols, n_slices)" << std::endl;

  Matrix<double, 1> mat1d(4);
  if (print) test_print(mat1d, "mat1d =");
  assert(mat1d.n_elem() == 4);
  assert(mat1d.n_rows() == 4);
  assert(mat1d.n_cols() == 1);

  // It would not compile because Matrix<T, 1>(uword __n1) is explicit
  // Matrix<double, 1> mat1d_x = 4;
  // if (print) test_print(mat1d, "mat1d_x =");

  Matrix<double, 2> mat2d(4, 3);
  if (print) test_print(mat2d, "mat2d =");
  assert(mat2d.n_elem() == 12);
  assert(mat2d.n_rows() == 4);
  assert(mat2d.n_cols() == 3);

  Matrix<double, 3> mat3d(4, 3, 2);
  if (print) test_print(mat3d, "mat3d =");
  assert(mat3d.n_elem() == 24);
  assert(mat3d.n_rows() == 4);
  assert(mat3d.n_cols() == 3);
  assert(mat3d.n_slices() == 2);
}

void matrix_test_constructor_a1_03(bool print = false) {
  std::cout
      << "[TEST]: A1_03. Elements (single value) + Dimensions\n"
      << "        vec: vec(const elem_type& val, n_rows)\n"
      << "        mat: mat(const elem_type& val, n_rows, n_cols)\n"
      << "        cube: cube(const elem_type& val, n_rows, n_cols, n_slides)"
      << std::endl;

  Matrix<double, 1> mat1d(1.0, 4);
  if (print) test_print(mat1d, "mat1d =");
  assert(mat1d.n_elem() == 4);
  assert(mat1d.n_rows() == 4);
  assert(mat1d.n_cols() == 1);
  assert(mat1d(0) == 1.0);
  assert(mat1d(3) == 1.0);

  Matrix<double, 2> mat2d(2.0, 4, 3);
  if (print) test_print(mat2d, "mat2d =");
  assert(mat2d.n_elem() == 12);
  assert(mat2d.n_rows() == 4);
  assert(mat2d.n_cols() == 3);
  assert(mat2d(0, 0) == 2.0);
  assert(mat2d(3, 2) == 2.0);

  Matrix<double, 3> mat3d(3.0, 4, 3, 2);
  if (print) test_print(mat3d, "mat3d =");
  assert(mat3d.n_elem() == 24);
  assert(mat3d.n_rows() == 4);
  assert(mat3d.n_cols() == 3);
  assert(mat3d.n_slices() == 2);
  assert(mat3d(0, 0, 0) == 3.0);
  assert(mat3d(3, 2, 1) == 3.0);
}

void matrix_test_constructor_a1_04(bool print = false) {
  std::cout
      << "[TEST]: A1_04. Elements (raw array) + Dimension\n"
      << "        vec: vec(const elem_type* vals, n_rows)\n"
      << "        mat: mat(const elem_type* vals, n_rows, n_cols)\n"
      << "        cube: cube(const elem_type* vals, n_rows, n_cols, n_slides)"
      << std::endl;

  Matrix<double, 1> mat1d(arr_1d, 4);
  if (print) test_print(mat1d, "mat1d =");
  assert(mat1d.n_elem() == 4);
  assert(mat1d.n_rows() == 4);
  assert(mat1d.n_cols() == 1);
  assert(mat1d(0) == 1.0);
  assert(mat1d(1) == 2.0);
  assert(mat1d(2) == 3.0);
  assert(mat1d(3) == 4.0);

  Matrix<double, 2> mat2d(arr_2d, 4, 3);
  if (print) test_print(mat2d, "mat2d =");
  assert(mat2d.n_elem() == 12);
  assert(mat2d.n_rows() == 4);
  assert(mat2d.n_cols() == 3);
  assert(mat2d(0, 0) == 1.0);
  assert(mat2d(1, 0) == 2.0);
  assert(mat2d(2, 0) == 3.0);
  assert(mat2d(3, 0) == 4.0);
  assert(mat2d(0, 2) == 9.0);
  assert(mat2d(1, 2) == 10.0);
  assert(mat2d(2, 2) == 11.0);
  assert(mat2d(3, 2) == 12.0);

  Matrix<double, 3> mat3d(arr_3d, 4, 3, 2);
  if (print) test_print(mat3d, "mat3d =");
  assert(mat3d.n_elem() == 24);
  assert(mat3d.n_rows() == 4);
  assert(mat3d.n_cols() == 3);
  assert(mat3d.n_slices() == 2);
  assert(mat3d(0, 0, 0) == 1.0);
  assert(mat3d(1, 0, 0) == 2.0);
  assert(mat3d(2, 0, 0) == 3.0);
  assert(mat3d(3, 0, 0) == 4.0);
  assert(mat3d(0, 2, 1) == 21.0);
  assert(mat3d(1, 2, 1) == 22.0);
  assert(mat3d(2, 2, 1) == 23.0);
  assert(mat3d(3, 2, 1) == 24.0);
}

void matrix_test_constructor_a1_05(bool print = false) {
  std::cout << "[TEST]: A1_05. Copy constructor\n"
            << "        vec: vec(const vec& other)\n"
            << "        mat: mat(const mat& other)\n"
            << "        cube: cube(const cube& other)" << std::endl;

  Matrix<double, 1> mat1d_a(arr_1d, 4);
  Matrix<double, 1> mat1d_b(mat1d_a);
  if (print) test_print(mat1d_b, "mat1d_b =");
  assert(mat1d_b.n_elem() == 4);
  assert(mat1d_b.n_rows() == 4);
  assert(mat1d_b.n_cols() == 1);
  assert(mat1d_b(0) == 1.0);
  assert(mat1d_b(1) == 2.0);
  assert(mat1d_b(2) == 3.0);
  assert(mat1d_b(3) == 4.0);

  Matrix<double, 2> mat2d_a(arr_2d, 4, 3);
  Matrix<double, 2> mat2d_b = mat2d_a;
  if (print) test_print(mat2d_b, "mat2d_b =");
  assert(mat2d_b.n_elem() == 12);
  assert(mat2d_b.n_rows() == 4);
  assert(mat2d_b.n_cols() == 3);
  assert(mat2d_b(0, 0) == 1.0);
  assert(mat2d_b(1, 0) == 2.0);
  assert(mat2d_b(2, 0) == 3.0);
  assert(mat2d_b(3, 0) == 4.0);
  assert(mat2d_b(0, 2) == 9.0);
  assert(mat2d_b(1, 2) == 10.0);
  assert(mat2d_b(2, 2) == 11.0);
  assert(mat2d_b(3, 2) == 12.0);

  Matrix<double, 3> mat3d_a(arr_3d, 4, 3, 2), mat3d_b;
  mat3d_b = mat3d_a;
  if (print) test_print(mat3d_b, "mat3d_b =");
  assert(mat3d_b.n_elem() == 24);
  assert(mat3d_b.n_rows() == 4);
  assert(mat3d_b.n_cols() == 3);
  assert(mat3d_b.n_slices() == 2);
  assert(mat3d_b(0, 0, 0) == 1.0);
  assert(mat3d_b(1, 0, 0) == 2.0);
  assert(mat3d_b(2, 0, 0) == 3.0);
  assert(mat3d_b(3, 0, 0) == 4.0);
  assert(mat3d_b(0, 2, 1) == 21.0);
  assert(mat3d_b(1, 2, 1) == 22.0);
  assert(mat3d_b(2, 2, 1) == 23.0);
  assert(mat3d_b(3, 2, 1) == 24.0);
}

void matrix_test_constructor_a1_06(bool print = false) {
  std::cout << "[TEST]: A1_06. Move constructor\n"
            << "        vec: vec(vec&& other) noexcept\n"
            << "        mat: mat(mat&& other) noexcept\n"
            << "        cube: cube(cube&& other) noexcept" << std::endl;

#if defined MATRIX_LIB_USE_CPP11
  Matrix<double, 1> mat1d_a(arr_1d, 4);
  Matrix<double, 1> mat1d_b(std::move(mat1d_a));
  if (print) test_print(mat1d_b, "mat1d_b =");
  assert(mat1d_b.n_elem() == 4);
  assert(mat1d_b.n_rows() == 4);
  assert(mat1d_b.n_cols() == 1);
  assert(mat1d_b(0) == 1.0);
  assert(mat1d_b(1) == 2.0);
  assert(mat1d_b(2) == 3.0);
  assert(mat1d_b(3) == 4.0);

  Matrix<double, 2> mat2d_a(arr_2d, 4, 3);
  Matrix<double, 2> mat2d_b = std::move(mat2d_a);
  if (print) test_print(mat2d_b, "mat2d_b =");
  assert(mat2d_b.n_elem() == 12);
  assert(mat2d_b.n_rows() == 4);
  assert(mat2d_b.n_cols() == 3);
  assert(mat2d_b(0, 0) == 1.0);
  assert(mat2d_b(1, 0) == 2.0);
  assert(mat2d_b(2, 0) == 3.0);
  assert(mat2d_b(3, 0) == 4.0);
  assert(mat2d_b(0, 2) == 9.0);
  assert(mat2d_b(1, 2) == 10.0);
  assert(mat2d_b(2, 2) == 11.0);
  assert(mat2d_b(3, 2) == 12.0);

  Matrix<double, 3> mat3d_a(arr_3d, 4, 3, 2), mat3d_b;
  mat3d_b = std::move(mat3d_a);
  if (print) test_print(mat3d_b, "mat3d_b =");
  assert(mat3d_b.n_elem() == 24);
  assert(mat3d_b.n_rows() == 4);
  assert(mat3d_b.n_cols() == 3);
  assert(mat3d_b.n_slices() == 2);
  assert(mat3d_b(0, 0, 0) == 1.0);
  assert(mat3d_b(1, 0, 0) == 2.0);
  assert(mat3d_b(2, 0, 0) == 3.0);
  assert(mat3d_b(3, 0, 0) == 4.0);
  assert(mat3d_b(0, 2, 1) == 21.0);
  assert(mat3d_b(1, 2, 1) == 22.0);
  assert(mat3d_b(2, 2, 1) == 23.0);
  assert(mat3d_b(3, 2, 1) == 24.0);
#else
  std::cout << "[NOTE]: TEST A06 SKIPPED (C++11 NOT SUPPORTED)\n";
#endif
}

void matrix_test_constructor_a1_07(bool print = false) {
  std::cout << "[TEST]: A1_07. Constructor with Sub-Matrix (SliceMatrix)\n"
            << "        vec: vec(slice_vec)\n"
            << "        mat: mat(slice_vec)" << std::endl;

  Matrix<double, 2> mat2d(arr_2d, 4, 3);

  Matrix<double, 1> mat1d_row(mat2d.row(1));
  if (print) test_print(mat1d_row, "mat1d_row = ");
  assert(mat1d_row.n_elem() == 3);
  assert(mat1d_row.n_rows() == 1);
  assert(mat1d_row.n_cols() == 3);
  assert(mat1d_row(0) == 2);
  assert(mat1d_row(1) == 6);
  assert(mat1d_row(2) == 10);
  Matrix<double, 2> mat2d_row = mat2d.row(1);
  if (print) test_print(mat2d_row, "mat2d_row = ");
  assert(mat2d_row.n_elem() == 3);
  assert(mat2d_row.n_rows() == 1);
  assert(mat2d_row.n_cols() == 3);
  assert(mat2d_row(0, 0) == 2);
  assert(mat2d_row(0, 1) == 6);
  assert(mat2d_row(0, 2) == 10);

  Matrix<double, 1> mat1d_col(mat2d.col(2));
  if (print) test_print(mat1d_col, "mat1d_col = ");
  assert(mat1d_col.n_elem() == 4);
  assert(mat1d_col.n_rows() == 4);
  assert(mat1d_col.n_cols() == 1);
  assert(mat1d_col(0) == 9);
  assert(mat1d_col(1) == 10);
  assert(mat1d_col(2) == 11);
  assert(mat1d_col(3) == 12);
  Matrix<double, 2> mat2d_col = mat2d.col(2);
  if (print) test_print(mat2d_col, "mat2d_col = ");
  assert(mat2d_col.n_elem() == 4);
  assert(mat2d_col.n_rows() == 4);
  assert(mat2d_col.n_cols() == 1);
  assert(mat2d_col(0, 0) == 9);
  assert(mat2d_col(1, 0) == 10);
  assert(mat2d_col(2, 0) == 11);
  assert(mat2d_col(3, 0) == 12);
}

void matrix_test_constructor_a1_08(bool print = false) {
  std::cout << "[TEST]: A1_08. Constructor with Sub-Matrix (GsliceMatrix)\n"
            << "        vec: vec(gslice_vec)\n"
            << "        mat: mat(gslice_mat)\n"
            << "        cube: cube(gslice_cube)" << std::endl;
  Matrix<double, 1> mat1d_a(arr_1d, 4);
  Matrix<double, 1> mat1d_b(mat1d_a.subvec(1, 2));
  if (print) test_print(mat1d_b, "mat1d_b =");
  assert(mat1d_b.n_elem() == 2);
  assert(mat1d_b.n_rows() == 2);
  assert(mat1d_b.n_cols() == 1);
  assert(mat1d_b(0) == 2.0);
  assert(mat1d_b(1) == 3.0);

  Matrix<double, 2> mat2d_a(arr_2d, 4, 3);
  Matrix<double, 2> mat2d_b(mat2d_a.submat(2, 1, 3, 2));
  if (print) test_print(mat2d_b, "mat2d_b =");
  assert(mat2d_b.n_elem() == 4);
  assert(mat2d_b.n_rows() == 2);
  assert(mat2d_b.n_cols() == 2);
  assert(mat2d_b(0, 0) == 7.0);
  assert(mat2d_b(1, 0) == 8.0);
  assert(mat2d_b(0, 1) == 11.0);
  assert(mat2d_b(1, 1) == 12.0);

  Matrix<double, 3> mat3d_a(arr_3d, 4, 3, 2);
  Matrix<double, 3> mat3d_b(mat3d_a.subcube(2, 1, 1, 3, 2, 1));
  if (print) test_print(mat3d_b, "mat3d_b =");
  assert(mat3d_b.n_elem() == 4);
  assert(mat3d_b.n_rows() == 2);
  assert(mat3d_b.n_cols() == 2);
  assert(mat3d_b.n_slices() == 1);
  assert(mat3d_b(0, 0, 0) == 19.0);
  assert(mat3d_b(1, 0, 0) == 20.0);
  assert(mat3d_b(0, 1, 0) == 23.0);
  assert(mat3d_b(1, 1, 0) == 24.0);
}

void matrix_test_constructor_a1_09(bool print = false) {
  std::cout << "[TEST]: A1_09. Constructor with Sub-Matrix (MaskMatrix)\n"
            << "        vec: vec(mask_vec)\n"
            << "        mat: mat(mask_vec)" << std::endl;

  const bool idx[4] = {false, true, false, true};
  std::valarray<bool> bool_arr(idx, 4);

  Matrix<double, 1> mat1d_a(arr_1d, 4);
  Matrix<double, 1> mat1d_b = mat1d_a(bool_arr);
  if (print) test_print(mat1d_b, "mat1d_b = ");
  assert(mat1d_b.n_elem() == 2);
  assert(mat1d_b.n_rows() == 2);
  assert(mat1d_b.n_cols() == 1);
  assert(mat1d_b(0) == 2);
  assert(mat1d_b(1) == 4);

  Matrix<double, 2> mat2d = mat1d_a(bool_arr);
  if (print) test_print(mat2d, "mat2d = ");
  assert(mat2d.n_elem() == 2);
  assert(mat2d.n_rows() == 2);
  assert(mat2d.n_cols() == 1);
  assert(mat2d(0, 0) == 2);
  assert(mat2d(1, 0) == 4);
}

void matrix_test_constructor_a1_10(bool print = false) {
  std::cout << "[TEST]: A1_10. Constructor with Sub-Matrix (IndirectMatrix)\n"
            << "        vec: vec(indirect_vec)\n"
            << "        mat: mat(indirect_mat)\n"
            << "        cube: cube(indirect_cube)" << std::endl;

  const std::size_t idx1[2] = {1, 3};
  std::valarray<std::size_t> idx_arr1(idx1, 2);

  Matrix<double, 1> mat1d_a(arr_1d, 4);

  Matrix<double, 1> mat1d_b = mat1d_a(idx_arr1);
  if (print) test_print(mat1d_b, "mat1d_c = ");
  assert(mat1d_b.n_elem() == 2);
  assert(mat1d_b.n_rows() == 2);
  assert(mat1d_b.n_cols() == 1);
  assert(mat1d_b(0) == 2);
  assert(mat1d_b(1) == 4);

  const std::size_t idx2[2] = {0, 2};
  std::valarray<std::size_t> idx_arr2(idx2, 2);

  Matrix<double, 2> mat2d_a(arr_2d, 4, 3);

  Matrix<double, 2> mat2d_b = mat2d_a(idx_arr1, idx_arr2);
  if (print) test_print(mat2d_b, "mat2d_b = ");
  assert(mat2d_b.n_elem() == 4);
  assert(mat2d_b.n_rows() == 2);
  assert(mat2d_b.n_cols() == 2);
  assert(mat2d_b(0, 0) == 2);
  assert(mat2d_b(0, 1) == 10);
  assert(mat2d_b(1, 0) == 4);
  assert(mat2d_b(1, 1) == 12);

  const std::size_t idx3[3] = {1};
  std::valarray<std::size_t> idx_arr3(idx3, 1);

  Matrix<double, 3> mat3d_a(arr_3d, 4, 3, 2);

  Matrix<double, 3> mat3d_b = mat3d_a(idx_arr1, idx_arr2, idx_arr3);
  if (print) test_print(mat3d_b, "mat3d_b = ");
  assert(mat3d_b.n_elem() == 4);
  assert(mat3d_b.n_rows() == 2);
  assert(mat3d_b.n_cols() == 2);
  assert(mat3d_b.n_slices() == 1);
  assert(mat3d_b(0, 0, 0) == 14);
  assert(mat3d_b(0, 1, 0) == 22);
  assert(mat3d_b(1, 0, 0) == 16);
  assert(mat3d_b(1, 1, 0) == 24);
}

void matrix_test_constructor_a1_11(bool print = false) {
  std::cout << "[TEST]: A1_11. Nested initializer list\n"
            << "        vec: vec(initializer_list)\n"
            << "        mat: mat(nested initializer_list)\n"
            << "        cube: cube(nested initializer_list)" << std::endl;

#if defined MATRIX_LIB_USE_CPP11
  Matrix<double, 1> mat1d = {1, 2, 3, 4};
  if (print) test_print(mat1d, "mat1d =");
  assert(mat1d.n_elem() == 4);
  assert(mat1d.n_rows() == 4);
  assert(mat1d.n_cols() == 1);
  assert(mat1d(0) == 1.0);
  assert(mat1d(1) == 2.0);
  assert(mat1d(2) == 3.0);
  assert(mat1d(3) == 4.0);

  Matrix<double, 2> mat2d = {{1, 5, 9}, {2, 6, 10}, {3, 7, 11}, {4, 8, 12}};
  if (print) test_print(mat2d, "mat2d =");
  assert(mat2d.n_elem() == 12);
  assert(mat2d.n_rows() == 4);
  assert(mat2d.n_cols() == 3);
  assert(mat2d(0, 0) == 1.0);
  assert(mat2d(1, 0) == 2.0);
  assert(mat2d(2, 0) == 3.0);
  assert(mat2d(3, 0) == 4.0);
  assert(mat2d(0, 2) == 9.0);
  assert(mat2d(1, 2) == 10.0);
  assert(mat2d(2, 2) == 11.0);
  assert(mat2d(3, 2) == 12.0);

  Matrix<double, 3> mat3d = {
      {{1, 5, 9}, {2, 6, 10}, {3, 7, 11}, {4, 8, 12}},
      {{13, 17, 21}, {14, 18, 22}, {15, 19, 23}, {16, 20, 24}}};
  if (print) test_print(mat3d, "mat3d =");
  assert(mat3d.n_elem() == 24);
  assert(mat3d.n_rows() == 4);
  assert(mat3d.n_cols() == 3);
  assert(mat3d.n_slices() == 2);
  assert(mat3d(0, 0, 0) == 1.0);
  assert(mat3d(1, 0, 0) == 2.0);
  assert(mat3d(2, 0, 0) == 3.0);
  assert(mat3d(3, 0, 0) == 4.0);
  assert(mat3d(0, 2, 1) == 21.0);
  assert(mat3d(1, 2, 1) == 22.0);
  assert(mat3d(2, 2, 1) == 23.0);
  assert(mat3d(3, 2, 1) == 24.0);
#else
  std::cout << "[NOTE]: TEST A08 SKIPPED (C++11 NOT SUPPORTED)\n";
#endif
}

void matrix_test_constructor_a1_12(bool print = false) {
  std::cout << "[TEST]: A1_12. Elements (valarray) + Dimension\n"
            << "        vec: vec(valarray)\n"
            << "        mat: mat(valarray, n_rows, n_cols)\n"
            << "        cube: cube(valarray, n_rows, n_cols, n_slices)"
            << std::endl;

  const std::valarray<double> va1(arr_1d, 4);
  Matrix<double, 1> mat1d(va1);
  if (print) test_print(mat1d, "mat1d =");
  assert(mat1d.n_elem() == 4);
  assert(mat1d.n_rows() == 4);
  assert(mat1d.n_cols() == 1);
  assert(mat1d(0) == 1.0);
  assert(mat1d(1) == 2.0);
  assert(mat1d(2) == 3.0);
  assert(mat1d(3) == 4.0);

  const std::valarray<double> va2(arr_2d, 12);
  Matrix<double, 2> mat2d(va2, 4, 3);
  if (print) test_print(mat2d, "mat2d =");
  assert(mat2d.n_elem() == 12);
  assert(mat2d.n_rows() == 4);
  assert(mat2d.n_cols() == 3);
  assert(mat2d(0, 0) == 1.0);
  assert(mat2d(1, 0) == 2.0);
  assert(mat2d(2, 0) == 3.0);
  assert(mat2d(3, 0) == 4.0);
  assert(mat2d(0, 2) == 9.0);
  assert(mat2d(1, 2) == 10.0);
  assert(mat2d(2, 2) == 11.0);
  assert(mat2d(3, 2) == 12.0);

  const std::valarray<double> va3(arr_3d, 24);
  Matrix<double, 3> mat3d(va3, 4, 3, 2);
  if (print) test_print(mat3d, "mat3d =");
  assert(mat3d.n_elem() == 24);
  assert(mat3d.n_rows() == 4);
  assert(mat3d.n_cols() == 3);
  assert(mat3d.n_slices() == 2);
  assert(mat3d(0, 0, 0) == 1.0);
  assert(mat3d(1, 0, 0) == 2.0);
  assert(mat3d(2, 0, 0) == 3.0);
  assert(mat3d(3, 0, 0) == 4.0);
  assert(mat3d(0, 2, 1) == 21.0);
  assert(mat3d(1, 2, 1) == 22.0);
  assert(mat3d(2, 2, 1) == 23.0);
  assert(mat3d(3, 2, 1) == 24.0);
}

void matrix_test_constructor_a1_13(bool print = false) {
  std::cout << "[TEST]: A1_13. Constructs a Matrix <- from -> a Vector\n"
            << "        vec: vec(const mat&)\n"
            << "        mat: mat(const vec&)" << std::endl;

  Matrix<double, 1> mat1d_a1(arr_1d, 4);
  Matrix<double, 1> mat1d_a2(arr_1d, 4);

  mat1d_a2.is_column_vector = false;
  Matrix<double, 2> mat2d_a1(arr_1d, 4, 1);
  Matrix<double, 2> mat2d_a2(arr_1d, 1, 4);

  Matrix<double, 1> mat1d_b1(mat2d_a1);  // col
  Matrix<double, 1> mat1d_b2(mat2d_a2);  // row

  Matrix<double, 2> mat2d_b1(mat1d_a1);  // col
  Matrix<double, 2> mat2d_b2(mat1d_a2);  // row

  if (print) test_print(mat1d_b1, "mat1d_b1 =");
  if (print) test_print(mat1d_b2, "mat1d_b2 =");
  if (print) test_print(mat2d_b1, "mat2d_b1 =");
  if (print) test_print(mat2d_b2, "mat2d_b2 =");

  assert(mat1d_b1.n_elem() == 4);
  assert(mat1d_b1.n_rows() == 4);
  assert(mat1d_b1.n_cols() == 1);
  assert(mat1d_b1(0) == 1.0);
  assert(mat1d_b1(1) == 2.0);
  assert(mat1d_b1(2) == 3.0);
  assert(mat1d_b1(3) == 4.0);

  assert(mat1d_b2.n_elem() == 4);
  assert(mat1d_b2.n_rows() == 1);
  assert(mat1d_b2.n_cols() == 4);
  assert(mat1d_b2(0) == 1.0);
  assert(mat1d_b2(1) == 2.0);
  assert(mat1d_b2(2) == 3.0);
  assert(mat1d_b2(3) == 4.0);

  assert(mat2d_b1.n_elem() == 4);
  assert(mat2d_b1.n_rows() == 4);
  assert(mat2d_b1.n_cols() == 1);
  assert(mat2d_b1(0, 0) == 1.0);
  assert(mat2d_b1(1, 0) == 2.0);
  assert(mat2d_b1(2, 0) == 3.0);
  assert(mat2d_b1(3, 0) == 4.0);

  assert(mat2d_b2.n_elem() == 4);
  assert(mat2d_b2.n_rows() == 1);
  assert(mat2d_b2.n_cols() == 4);
  assert(mat2d_b2(0, 0) == 1.0);
  assert(mat2d_b2(0, 1) == 2.0);
  assert(mat2d_b2(0, 2) == 3.0);
  assert(mat2d_b2(0, 3) == 4.0);
}

void matrix_test_constructor_a1_14(bool print = false) {
  std::cout << "[TEST]: A1_14. Constructs a Matrix <- from -> a Vector (move)\n"
            << "        vec: vec(mat&&)\n"
            << "        mat: mat(vec&&)" << std::endl;

#if defined(MATRIX_LIB_USE_CPP11)
  Matrix<double, 1> mat1d_a1(arr_1d, 4);
  Matrix<double, 1> mat1d_a2(arr_1d, 4);

  mat1d_a2.is_column_vector = false;
  Matrix<double, 2> mat2d_a1(arr_1d, 4, 1);
  Matrix<double, 2> mat2d_a2(arr_1d, 1, 4);

  Matrix<double, 1> mat1d_b1(std::move(mat2d_a1));  // col
  Matrix<double, 1> mat1d_b2(std::move(mat2d_a2));  // row

  Matrix<double, 2> mat2d_b1(std::move(mat1d_a1));  // col
  Matrix<double, 2> mat2d_b2(std::move(mat1d_a2));  // row

  if (print) test_print(mat1d_b1, "mat1d_b1 =");
  if (print) test_print(mat1d_b2, "mat1d_b2 =");
  if (print) test_print(mat2d_b1, "mat2d_b1 =");
  if (print) test_print(mat2d_b2, "mat2d_b2 =");

  assert(mat1d_b1.n_elem() == 4);
  assert(mat1d_b1.n_rows() == 4);
  assert(mat1d_b1.n_cols() == 1);
  assert(mat1d_b1(0) == 1.0);
  assert(mat1d_b1(1) == 2.0);
  assert(mat1d_b1(2) == 3.0);
  assert(mat1d_b1(3) == 4.0);

  assert(mat1d_b2.n_elem() == 4);
  assert(mat1d_b2.n_rows() == 1);
  assert(mat1d_b2.n_cols() == 4);
  assert(mat1d_b2(0) == 1.0);
  assert(mat1d_b2(1) == 2.0);
  assert(mat1d_b2(2) == 3.0);
  assert(mat1d_b2(3) == 4.0);

  assert(mat2d_b1.n_elem() == 4);
  assert(mat2d_b1.n_rows() == 4);
  assert(mat2d_b1.n_cols() == 1);
  assert(mat2d_b1(0, 0) == 1.0);
  assert(mat2d_b1(1, 0) == 2.0);
  assert(mat2d_b1(2, 0) == 3.0);
  assert(mat2d_b1(3, 0) == 4.0);

  assert(mat2d_b2.n_elem() == 4);
  assert(mat2d_b2.n_rows() == 1);
  assert(mat2d_b2.n_cols() == 4);
  assert(mat2d_b2(0, 0) == 1.0);
  assert(mat2d_b2(0, 1) == 2.0);
  assert(mat2d_b2(0, 2) == 3.0);
  assert(mat2d_b2(0, 3) == 4.0);
#else
  std::cout << "[NOTE]: TEST A14 SKIPPED (C++11 NOT SUPPORTED)\n";
#endif
}

// ----- A2: Testing Matrix<T, N> Assignments -----

void matrix_test_assignment_a2_01(bool print = false) {
  std::cout << "[TEST]: A2_01. Copy assignment operator\n"
            << "        vec: vec& operator=(const vec&)\n"
            << "        mat: mat& operator=(const mat&)\n"
            << "        cube: cube& operator=(const cube&)" << std::endl;

  Matrix<double, 1> mat1d_a(arr_1d, 4);
  Matrix<double, 1> mat1d_b;
  mat1d_b = mat1d_a;
  if (print) test_print(mat1d_b, "mat1d_b =");
  assert(mat1d_b(0) == 1.0);
  assert(mat1d_b(1) == 2.0);
  assert(mat1d_b(2) == 3.0);
  assert(mat1d_b(3) == 4.0);

  Matrix<double, 2> mat2d_a(arr_2d, 4, 3);
  Matrix<double, 2> mat2d_b;
  mat2d_b = mat2d_a;
  if (print) test_print(mat2d_b, "mat2d_b =");
  assert(mat2d_b(0, 0) == 1.0);
  assert(mat2d_b(1, 0) == 2.0);
  assert(mat2d_b(2, 0) == 3.0);
  assert(mat2d_b(3, 0) == 4.0);
  assert(mat2d_b(0, 2) == 9.0);
  assert(mat2d_b(1, 2) == 10.0);
  assert(mat2d_b(2, 2) == 11.0);
  assert(mat2d_b(3, 2) == 12.0);

  Matrix<double, 3> mat3d_a(arr_3d, 4, 3, 2), mat3d_b;
  mat3d_b = mat3d_a;
  if (print) test_print(mat3d_b, "mat3d_b =");
  assert(mat3d_b(0, 0, 0) == 1.0);
  assert(mat3d_b(1, 0, 0) == 2.0);
  assert(mat3d_b(2, 0, 0) == 3.0);
  assert(mat3d_b(3, 0, 0) == 4.0);
  assert(mat3d_b(0, 2, 1) == 21.0);
  assert(mat3d_b(1, 2, 1) == 22.0);
  assert(mat3d_b(2, 2, 1) == 23.0);
  assert(mat3d_b(3, 2, 1) == 24.0);
}

void matrix_test_assignment_a2_02(bool print = false) {
  std::cout << "[TEST]: A2_02. Move assignment operator\n"
            << "        vec: vec& operator=(vec&&)\n"
            << "        mat: mat& operator=(mat&&)\n"
            << "        cube: cube& operator=(cube&&)" << std::endl;

#if defined MATRIX_LIB_USE_CPP11
  Matrix<double, 1> mat1d_a(arr_1d, 4);
  Matrix<double, 1> mat1d_b;
  mat1d_b = std::move(mat1d_a);
  if (print) test_print(mat1d_b, "mat1d_b =");
  assert(mat1d_b(0) == 1.0);
  assert(mat1d_b(1) == 2.0);
  assert(mat1d_b(2) == 3.0);
  assert(mat1d_b(3) == 4.0);

  Matrix<double, 2> mat2d_a(arr_2d, 4, 3);
  Matrix<double, 2> mat2d_b;
  mat2d_b = std::move(mat2d_a);
  if (print) test_print(mat2d_b, "mat2d_b =");
  assert(mat2d_b(0, 0) == 1.0);
  assert(mat2d_b(1, 0) == 2.0);
  assert(mat2d_b(2, 0) == 3.0);
  assert(mat2d_b(3, 0) == 4.0);
  assert(mat2d_b(0, 2) == 9.0);
  assert(mat2d_b(1, 2) == 10.0);
  assert(mat2d_b(2, 2) == 11.0);
  assert(mat2d_b(3, 2) == 12.0);

  Matrix<double, 3> mat3d_a(arr_3d, 4, 3, 2), mat3d_b;
  mat3d_b = std::move(mat3d_a);
  if (print) test_print(mat3d_b, "mat3d_b =");
  assert(mat3d_b(0, 0, 0) == 1.0);
  assert(mat3d_b(1, 0, 0) == 2.0);
  assert(mat3d_b(2, 0, 0) == 3.0);
  assert(mat3d_b(3, 0, 0) == 4.0);
  assert(mat3d_b(0, 2, 1) == 21.0);
  assert(mat3d_b(1, 2, 1) == 22.0);
  assert(mat3d_b(2, 2, 1) == 23.0);
  assert(mat3d_b(3, 2, 1) == 24.0);
#else
  std::cout << "[NOTE]: TEST B02 SKIPPED (C++11 NOT SUPPORTED)\n";
#endif
}

void matrix_test_assignment_a2_03(bool print = false) {
  std::cout << "[TEST]: A2_03. Assigns a Matrix with a value\n"
            << "        vec: vec& operator=(const elem_type& val)\n"
            << "        mat: mat& operator=(const elem_type& val)\n"
            << "        cube: cube& operator=(const elem_type& val)"
            << std::endl;

  Matrix<double, 1> mat1d(arr_1d, 4);
  if (print) test_print(mat1d, "mat1d =");
  mat1d = 1;
  if (print) test_print(mat1d, "mat1d =");
  assert(mat1d(0) == 1.0);
  assert(mat1d(1) == 1.0);
  assert(mat1d(2) == 1.0);
  assert(mat1d(3) == 1.0);

  Matrix<double, 2> mat2d(arr_2d, 4, 3);
  if (print) test_print(mat2d, "mat2d =");
  mat2d = 2;
  if (print) test_print(mat2d, "mat2d =");
  assert(mat2d(0, 0) == 2.0);
  assert(mat2d(1, 0) == 2.0);
  assert(mat2d(2, 0) == 2.0);
  assert(mat2d(3, 0) == 2.0);
  assert(mat2d(0, 2) == 2.0);
  assert(mat2d(1, 2) == 2.0);
  assert(mat2d(2, 2) == 2.0);
  assert(mat2d(3, 2) == 2.0);

  Matrix<double, 3> mat3d(arr_3d, 4, 3, 2);
  if (print) test_print(mat3d, "mat3d =");
  mat3d = 3;
  if (print) test_print(mat3d, "mat3d =");
  assert(mat3d(0, 0, 0) == 3.0);
  assert(mat3d(1, 0, 0) == 3.0);
  assert(mat3d(2, 0, 0) == 3.0);
  assert(mat3d(3, 0, 0) == 3.0);
  assert(mat3d(0, 2, 1) == 3.0);
  assert(mat3d(1, 2, 1) == 3.0);
  assert(mat3d(2, 2, 1) == 3.0);
  assert(mat3d(3, 2, 1) == 3.0);
}

void matrix_test_assignment_a2_04(bool print = false) {
  std::cout << "[TEST]: A2_04. Assignment with Sub-Matrix (SliceMatrix)\n"
            << "        vec: vec& operator=(slice_vec)\n"
            << "        mat: mat& operator=(slice_vec)" << std::endl;

  Matrix<double, 2> mat2d(arr_2d, 4, 3);

  Matrix<double, 1> mat1d_row;
  mat1d_row = mat2d.row(1);
  if (print) test_print(mat1d_row, "mat1d_row = ");
  assert(mat1d_row.n_elem() == 3);
  assert(mat1d_row.n_rows() == 1);
  assert(mat1d_row.n_cols() == 3);
  assert(mat1d_row(0) == 2);
  assert(mat1d_row(1) == 6);
  assert(mat1d_row(2) == 10);
  Matrix<double, 2> mat2d_row;
  mat2d_row = mat2d.row(1);
  if (print) test_print(mat2d_row, "mat2d_row = ");
  assert(mat2d_row.n_elem() == 3);
  assert(mat2d_row.n_rows() == 1);
  assert(mat2d_row.n_cols() == 3);
  assert(mat2d_row(0, 0) == 2);
  assert(mat2d_row(0, 1) == 6);
  assert(mat2d_row(0, 2) == 10);

  Matrix<double, 1> mat1d_col;
  mat1d_col = mat2d.col(2);
  if (print) test_print(mat1d_col, "mat1d_col = ");
  assert(mat1d_col.n_elem() == 4);
  assert(mat1d_col.n_rows() == 4);
  assert(mat1d_col.n_cols() == 1);
  assert(mat1d_col(0) == 9);
  assert(mat1d_col(1) == 10);
  assert(mat1d_col(2) == 11);
  assert(mat1d_col(3) == 12);
  Matrix<double, 2> mat2d_col;
  mat2d_col = mat2d.col(2);
  if (print) test_print(mat2d_col, "mat2d_col = ");
  assert(mat2d_col.n_elem() == 4);
  assert(mat2d_col.n_rows() == 4);
  assert(mat2d_col.n_cols() == 1);
  assert(mat2d_col(0, 0) == 9);
  assert(mat2d_col(1, 0) == 10);
  assert(mat2d_col(2, 0) == 11);
  assert(mat2d_col(3, 0) == 12);
}

void matrix_test_assignment_a2_05(bool print = false) {
  std::cout << "[TEST]: A2_05. Assignment with Sub-Matrix (GsliceMatrix)\n"
            << "        vec: vec& operator=(gslice_vec)\n"
            << "        mat: mat& operator=(gslice_mat)\n"
            << "        cube: cube& operator=(gslice_cube)" << std::endl;
  Matrix<double, 1> mat1d_a(arr_1d, 4);
  Matrix<double, 1> mat1d_b;
  mat1d_b = mat1d_a.subvec(1, 2);
  if (print) test_print(mat1d_b, "mat1d_b =");
  assert(mat1d_b.n_elem() == 2);
  assert(mat1d_b.n_rows() == 2);
  assert(mat1d_b.n_cols() == 1);
  assert(mat1d_b(0) == 2.0);
  assert(mat1d_b(1) == 3.0);

  Matrix<double, 2> mat2d_a(arr_2d, 4, 3);
  Matrix<double, 2> mat2d_b;
  mat2d_b = mat2d_a.submat(2, 1, 3, 2);
  if (print) test_print(mat2d_b, "mat2d_b =");
  assert(mat2d_b.n_elem() == 4);
  assert(mat2d_b.n_rows() == 2);
  assert(mat2d_b.n_cols() == 2);
  assert(mat2d_b(0, 0) == 7.0);
  assert(mat2d_b(1, 0) == 8.0);
  assert(mat2d_b(0, 1) == 11.0);
  assert(mat2d_b(1, 1) == 12.0);

  Matrix<double, 3> mat3d_a(arr_3d, 4, 3, 2);
  Matrix<double, 3> mat3d_b;
  mat3d_b = mat3d_a.subcube(2, 1, 1, 3, 2, 1);
  if (print) test_print(mat3d_b, "mat3d_b =");
  assert(mat3d_b.n_elem() == 4);
  assert(mat3d_b.n_rows() == 2);
  assert(mat3d_b.n_cols() == 2);
  assert(mat3d_b.n_slices() == 1);
  assert(mat3d_b(0, 0, 0) == 19.0);
  assert(mat3d_b(1, 0, 0) == 20.0);
  assert(mat3d_b(0, 1, 0) == 23.0);
  assert(mat3d_b(1, 1, 0) == 24.0);
}

void matrix_test_assignment_a2_06(bool print = false) {
  std::cout << "[TEST]: A2_06. Assignment with Sub-Matrix (MaskMatrix)\n"
            << "        vec: vec& operator=(mask_vec)\n"
            << "        mat: mat& operator=(mask_vec)" << std::endl;

  const bool idx[4] = {false, true, false, true};
  std::valarray<bool> bool_arr(idx, 4);

  Matrix<double, 1> mat1d_a(arr_1d, 4);
  Matrix<double, 1> mat1d_b;
  mat1d_b = mat1d_a(bool_arr);
  if (print) test_print(mat1d_b, "mat1d_b = ");
  assert(mat1d_b.n_elem() == 2);
  assert(mat1d_b.n_rows() == 2);
  assert(mat1d_b.n_cols() == 1);
  assert(mat1d_b(0) == 2);
  assert(mat1d_b(1) == 4);

  Matrix<double, 2> mat2d;
  mat2d = mat1d_a(bool_arr);
  if (print) test_print(mat2d, "mat2d = ");
  assert(mat2d.n_elem() == 2);
  assert(mat2d.n_rows() == 2);
  assert(mat2d.n_cols() == 1);
  assert(mat2d(0, 0) == 2);
  assert(mat2d(1, 0) == 4);
}

void matrix_test_assignment_a2_07(bool print = false) {
  std::cout << "[TEST]: A2_07. Assignment with Sub-Matrix (IndirectMatrix)\n"
            << "        vec: vec& operator=(indirect_vec)\n"
            << "        mat: mat& operator=(indirect_mat)\n"
            << "        cube: cube& operator=(indirect_cube)" << std::endl;

  const std::size_t idx1[2] = {1, 3};
  std::valarray<std::size_t> idx_arr1(idx1, 2);

  Matrix<double, 1> mat1d_a(arr_1d, 4);

  Matrix<double, 1> mat1d_b;
  mat1d_b = mat1d_a(idx_arr1);
  if (print) test_print(mat1d_b, "mat1d_c = ");
  assert(mat1d_b.n_elem() == 2);
  assert(mat1d_b.n_rows() == 2);
  assert(mat1d_b.n_cols() == 1);
  assert(mat1d_b(0) == 2);
  assert(mat1d_b(1) == 4);

  const std::size_t idx2[2] = {0, 2};
  std::valarray<std::size_t> idx_arr2(idx2, 2);

  Matrix<double, 2> mat2d_a(arr_2d, 4, 3);

  Matrix<double, 2> mat2d_b;
  mat2d_b = mat2d_a(idx_arr1, idx_arr2);
  if (print) test_print(mat2d_b, "mat2d_b = ");
  assert(mat2d_b.n_elem() == 4);
  assert(mat2d_b.n_rows() == 2);
  assert(mat2d_b.n_cols() == 2);
  assert(mat2d_b(0, 0) == 2);
  assert(mat2d_b(0, 1) == 10);
  assert(mat2d_b(1, 0) == 4);
  assert(mat2d_b(1, 1) == 12);

  const std::size_t idx3[3] = {1};
  std::valarray<std::size_t> idx_arr3(idx3, 1);

  Matrix<double, 3> mat3d_a(arr_3d, 4, 3, 2);

  Matrix<double, 3> mat3d_b;
  mat3d_b = mat3d_a(idx_arr1, idx_arr2, idx_arr3);
  if (print) test_print(mat3d_b, "mat3d_b = ");
  assert(mat3d_b.n_elem() == 4);
  assert(mat3d_b.n_rows() == 2);
  assert(mat3d_b.n_cols() == 2);
  assert(mat3d_b.n_slices() == 1);
  assert(mat3d_b(0, 0, 0) == 14);
  assert(mat3d_b(0, 1, 0) == 22);
  assert(mat3d_b(1, 0, 0) == 16);
  assert(mat3d_b(1, 1, 0) == 24);
}

void matrix_test_assignment_a2_08(bool print = false) {
  std::cout << "[TEST]: A2_08. Nested initializer list\n"
            << "        vec: vec& operator=(initializer_list)\n"
            << "        mat: mat& operator(nested initializer_list)\n"
            << "        cube: cube& operator(nested initializer_list)"
            << std::endl;

#if defined MATRIX_LIB_USE_CPP11
  Matrix<double, 1> mat1d;
  mat1d = {1, 2, 3, 4};
  if (print) test_print(mat1d, "mat1d =");
  assert(mat1d.n_elem() == 4);
  assert(mat1d.n_rows() == 4);
  assert(mat1d.n_cols() == 1);
  assert(mat1d(0) == 1.0);
  assert(mat1d(1) == 2.0);
  assert(mat1d(2) == 3.0);
  assert(mat1d(3) == 4.0);

  Matrix<double, 2> mat2d;
  mat2d = {{1, 5, 9}, {2, 6, 10}, {3, 7, 11}, {4, 8, 12}};
  if (print) test_print(mat2d, "mat2d =");
  assert(mat2d.n_elem() == 12);
  assert(mat2d.n_rows() == 4);
  assert(mat2d.n_cols() == 3);
  assert(mat2d(0, 0) == 1.0);
  assert(mat2d(1, 0) == 2.0);
  assert(mat2d(2, 0) == 3.0);
  assert(mat2d(3, 0) == 4.0);
  assert(mat2d(0, 2) == 9.0);
  assert(mat2d(1, 2) == 10.0);
  assert(mat2d(2, 2) == 11.0);
  assert(mat2d(3, 2) == 12.0);

  Matrix<double, 3> mat3d;
  mat3d = {{{1, 5, 9}, {2, 6, 10}, {3, 7, 11}, {4, 8, 12}},
           {{13, 17, 21}, {14, 18, 22}, {15, 19, 23}, {16, 20, 24}}};
  if (print) test_print(mat3d, "mat3d =");
  assert(mat3d.n_elem() == 24);
  assert(mat3d.n_rows() == 4);
  assert(mat3d.n_cols() == 3);
  assert(mat3d.n_slices() == 2);
  assert(mat3d(0, 0, 0) == 1.0);
  assert(mat3d(1, 0, 0) == 2.0);
  assert(mat3d(2, 0, 0) == 3.0);
  assert(mat3d(3, 0, 0) == 4.0);
  assert(mat3d(0, 2, 1) == 21.0);
  assert(mat3d(1, 2, 1) == 22.0);
  assert(mat3d(2, 2, 1) == 23.0);
  assert(mat3d(3, 2, 1) == 24.0);
#else
  std::cout << "[NOTE]: TEST A08 SKIPPED (C++11 NOT SUPPORTED)\n";
#endif
}

void matrix_test_assignment_a2_09(bool print = false) {
  std::cout << "[TEST]: A2_09. Assigns a Matrix <- with -> a Vector\n"
            << "        vec: vec& operator=(const mat&)\n"
            << "        mat: mat& operator=(const vec&)" << std::endl;

  Matrix<double, 1> mat1d_a1(arr_1d, 4);
  Matrix<double, 1> mat1d_a2(arr_1d, 4);

  mat1d_a2.is_column_vector = false;
  Matrix<double, 2> mat2d_a1(arr_1d, 4, 1);
  Matrix<double, 2> mat2d_a2(arr_1d, 1, 4);

  Matrix<double, 1> mat1d_b1;
  Matrix<double, 1> mat1d_b2;
  mat1d_b1 = mat2d_a1;  // col
  mat1d_b2 = mat2d_a2;  // row

  Matrix<double, 2> mat2d_b1;
  Matrix<double, 2> mat2d_b2;
  mat2d_b1 = mat1d_a1;  // col
  mat2d_b2 = mat1d_a2;  // row

  if (print) test_print(mat1d_b1, "mat1d_b1 =");
  if (print) test_print(mat1d_b2, "mat1d_b2 =");
  if (print) test_print(mat2d_b1, "mat2d_b1 =");
  if (print) test_print(mat2d_b2, "mat2d_b2 =");

  assert(mat1d_b1.n_elem() == 4);
  assert(mat1d_b1.n_rows() == 4);
  assert(mat1d_b1.n_cols() == 1);
  assert(mat1d_b1(0) == 1.0);
  assert(mat1d_b1(1) == 2.0);
  assert(mat1d_b1(2) == 3.0);
  assert(mat1d_b1(3) == 4.0);

  assert(mat1d_b2.n_elem() == 4);
  assert(mat1d_b2.n_rows() == 1);
  assert(mat1d_b2.n_cols() == 4);
  assert(mat1d_b2(0) == 1.0);
  assert(mat1d_b2(1) == 2.0);
  assert(mat1d_b2(2) == 3.0);
  assert(mat1d_b2(3) == 4.0);

  assert(mat2d_b1.n_elem() == 4);
  assert(mat2d_b1.n_rows() == 4);
  assert(mat2d_b1.n_cols() == 1);
  assert(mat2d_b1(0, 0) == 1.0);
  assert(mat2d_b1(1, 0) == 2.0);
  assert(mat2d_b1(2, 0) == 3.0);
  assert(mat2d_b1(3, 0) == 4.0);

  assert(mat2d_b2.n_elem() == 4);
  assert(mat2d_b2.n_rows() == 1);
  assert(mat2d_b2.n_cols() == 4);
  assert(mat2d_b2(0, 0) == 1.0);
  assert(mat2d_b2(0, 1) == 2.0);
  assert(mat2d_b2(0, 2) == 3.0);
  assert(mat2d_b2(0, 3) == 4.0);
}

void matrix_test_assignment_a2_10(bool print = false) {
  std::cout << "[TEST]: A2_10. Assigns a Matrix <- with -> a Vector (move)\n"
            << "        vec: vec& operator=(mat&&)\n"
               "        mat: mat& operator=(vec&&)"
            << std::endl;

#if defined(MATRIX_LIB_USE_CPP11)
  Matrix<double, 1> mat1d_a1(arr_1d, 4);
  Matrix<double, 1> mat1d_a2(arr_1d, 4);

  mat1d_a2.is_column_vector = false;
  Matrix<double, 2> mat2d_a1(arr_1d, 4, 1);
  Matrix<double, 2> mat2d_a2(arr_1d, 1, 4);

  Matrix<double, 1> mat1d_b1;
  Matrix<double, 1> mat1d_b2;
  mat1d_b1 = std::move(mat2d_a1);  // col
  mat1d_b2 = std::move(mat2d_a2);  // row

  Matrix<double, 2> mat2d_b1;
  Matrix<double, 2> mat2d_b2;
  mat2d_b1 = std::move(mat1d_a1);  // col
  mat2d_b2 = std::move(mat1d_a2);  // row

  if (print) test_print(mat1d_b1, "mat1d_b1 =");
  if (print) test_print(mat1d_b2, "mat1d_b2 =");
  if (print) test_print(mat2d_b1, "mat2d_b1 =");
  if (print) test_print(mat2d_b2, "mat2d_b2 =");

  assert(mat1d_b1.n_elem() == 4);
  assert(mat1d_b1.n_rows() == 4);
  assert(mat1d_b1.n_cols() == 1);
  assert(mat1d_b1(0) == 1.0);
  assert(mat1d_b1(1) == 2.0);
  assert(mat1d_b1(2) == 3.0);
  assert(mat1d_b1(3) == 4.0);

  assert(mat1d_b2.n_elem() == 4);
  assert(mat1d_b2.n_rows() == 1);
  assert(mat1d_b2.n_cols() == 4);
  assert(mat1d_b2(0) == 1.0);
  assert(mat1d_b2(1) == 2.0);
  assert(mat1d_b2(2) == 3.0);
  assert(mat1d_b2(3) == 4.0);

  assert(mat2d_b1.n_elem() == 4);
  assert(mat2d_b1.n_rows() == 4);
  assert(mat2d_b1.n_cols() == 1);
  assert(mat2d_b1(0, 0) == 1.0);
  assert(mat2d_b1(1, 0) == 2.0);
  assert(mat2d_b1(2, 0) == 3.0);
  assert(mat2d_b1(3, 0) == 4.0);

  assert(mat2d_b2.n_elem() == 4);
  assert(mat2d_b2.n_rows() == 1);
  assert(mat2d_b2.n_cols() == 4);
  assert(mat2d_b2(0, 0) == 1.0);
  assert(mat2d_b2(0, 1) == 2.0);
  assert(mat2d_b2(0, 2) == 3.0);
  assert(mat2d_b2(0, 3) == 4.0);
#else
  std::cout << "[NOTE]: TEST B10 SKIPPED (C++11 NOT SUPPORTED)\n";
#endif
}

void matrix_test_element_access_b1_02(bool print = false) {
  std::cout << "[TEST]: B1_02. Access elements with operator[]\n";

  const std::valarray<double> va1(arr_1d, 4);
  const Matrix<double, 1> mat1d(va1);
  if (print) test_print(mat1d, "mat1d =");
  assert(mat1d[0] == 1.0);
  assert(mat1d[1] == 2.0);
  assert(mat1d[2] == 3.0);
  assert(mat1d[3] == 4.0);

  const std::valarray<double> va2(arr_2d, 12);
  Matrix<double, 2> mat2d(va2, 4, 3);
  mat2d[0] = 0, mat2d[1] = 0, mat2d[2] = 0, mat2d[3] = 0;
  if (print) test_print(mat2d, "mat2d =");
  assert(mat2d(0, 0) == 0.0);
  assert(mat2d(1, 0) == 0.0);
  assert(mat2d(2, 0) == 0.0);
  assert(mat2d(3, 0) == 0.0);
  assert(mat2d(0, 2) == 9.0);
  assert(mat2d(1, 2) == 10.0);
  assert(mat2d(2, 2) == 11.0);
  assert(mat2d(3, 2) == 12.0);

  const std::valarray<double> va3(arr_3d, 24);
  const Matrix<double, 3> mat3d(va3, 4, 3, 2);
  if (print) test_print(mat3d, "mat3d =");
  assert(mat3d[0] == 1.0);
  assert(mat3d[1] == 2.0);
  assert(mat3d[2] == 3.0);
  assert(mat3d[3] == 4.0);
  assert(mat3d[20] == 21.0);
  assert(mat3d[21] == 22.0);
  assert(mat3d[22] == 23.0);
  assert(mat3d[23] == 24.0);
}

void matrix_2d_test_gslice_array(bool print = false) {
  std::cout
      << "[TEST]: 2D Matrix's std::gslice_array related member functions\n";

  const double arr[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
  Matrix<double, 2> mat2d_a(arr, 4, 3);

  if (print) test_print(mat2d_a, "mat2d_a =");
  if (print) std::cout << "Apply mat2d_a.submat(1,1,3,2) = 0\n";

  const uword sz[2] = {2, 3};
  const uword st[2] = {4, 1};

  const uword start = 4 * 1 + 1;
  const std::valarray<uword> size(sz, 2);
  const std::valarray<uword> stride(st, 2);
  std::gslice gs(start, size, stride);

  mat2d_a[gs] = 0;
  if (print) test_print(mat2d_a, "mat2d_a =");
  assert(mat2d_a(0, 0) == 1);
  assert(mat2d_a(0, 1) == 5);
  assert(mat2d_a(0, 2) == 9);
  assert(mat2d_a(3, 0) == 4);
  assert(mat2d_a(3, 1) == 0);
  assert(mat2d_a(3, 2) == 0);
}

void matrix_test_slicing_with_slicematrix_b2_01(bool print = false) {
  std::cout << "[TEST]: B2_01. Matrix<T, 1>::operator()(std::slice)"
            << " Non-const version\n";
  const double arr[] = {1, 2, 3, 4, 5};

  Matrix<double, 1> mat1d_a(arr, 5);
  mat1d_a(std::slice(1, 3, 1)) = 0;
  if (print) test_print(mat1d_a, "mat1d_a = ");
  assert(mat1d_a.n_elem() == 5);
  assert(mat1d_a.n_rows() == 5);
  assert(mat1d_a.n_cols() == 1);
  assert(mat1d_a(0) == 1);
  assert(mat1d_a(1) == 0);
  assert(mat1d_a(2) == 0);
  assert(mat1d_a(3) == 0);
  assert(mat1d_a(4) == 5);

  Matrix<double, 1> mat1d_b(arr, 5);
  Matrix<double, 1> mat1d_c = mat1d_b(std::slice(0, 3, 2));
  if (print) test_print(mat1d_c, "mat1d_c = ");
  assert(mat1d_c.n_elem() == 3);
  assert(mat1d_c.n_rows() == 3);
  assert(mat1d_c.n_cols() == 1);
  assert(mat1d_c(0) == 1);
  assert(mat1d_c(1) == 3);
  assert(mat1d_c(2) == 5);
}

void matrix_test_slicing_with_slicematrix_b2_02(bool print = false) {
  std::cout << "[TEST]: B2_02. Matrix<T, 1>::operator()(std::slice)"
            << " Const version\n";
  ;
  const double arr[] = {1, 2, 3, 4, 5};

  const Matrix<double, 1> mat1d_a(arr, 5);
  mat1d_a(std::slice(1, 3, 1)) = 0;  // Will Do Nothing!
  if (print) test_print(mat1d_a, "mat1d_a = ");
  assert(mat1d_a.n_elem() == 5);
  assert(mat1d_a.n_rows() == 5);
  assert(mat1d_a.n_cols() == 1);
  assert(mat1d_a(0) == 1);
  assert(mat1d_a(1) == 2);
  assert(mat1d_a(2) == 3);
  assert(mat1d_a(3) == 4);
  assert(mat1d_a(4) == 5);

  const Matrix<double, 1> mat1d_b(arr, 5);
  Matrix<double, 1> mat1d_c = mat1d_b(std::slice(0, 3, 2));
  if (print) test_print(mat1d_c, "mat1d_c = ");
  assert(mat1d_c.n_elem() == 3);
  assert(mat1d_c.n_rows() == 3);
  assert(mat1d_c.n_cols() == 1);
  assert(mat1d_c(0) == 1);
  assert(mat1d_c(1) == 3);
  assert(mat1d_c(2) == 5);
}

void matrix_test_slicing_with_slicematrix_b2_03(bool print = false) {
  std::cout << "[TEST]: B2_03. Matrix<T, 1>::subvec(i, j)"
            << " Non-const version\n";

  Matrix<double, 1> mat1d(arr_1d, 4);
  if (print) test_print(mat1d, "mat1d =");
  mat1d.subvec(2, 3) = 5.0;
  if (print) std::cout << "Apply mat1d.subvec(2, 3) = 5\n";
  if (print) test_print(mat1d, "mat1d =");
  assert(mat1d(0) == 1.0);
  assert(mat1d(1) == 2.0);
  assert(mat1d(2) == 5.0);
  assert(mat1d(3) == 5.0);

  // mat1d.subvec(0, 2) is a SliceMatrix<double>
  const Matrix<double, 1> mat1d_a(mat1d.subvec(0, 2));
  if (print) test_print(mat1d_a, "mat1d_a =");
  assert(mat1d_a.n_elem() == 3);
  assert(mat1d_a(0) == 1.0);
  assert(mat1d_a(1) == 2.0);
  assert(mat1d_a(2) == 5.0);

  // mat1d_a.subvec(0, 2) is a Matrix<T, 1>
  Matrix<double, 1> mat1d_b(mat1d_a.subvec(0, 1));
  if (print) test_print(mat1d_a, "mat1d_b =");
  assert(mat1d_b.n_elem() == 2);
  assert(mat1d_b(0) == 1.0);
  assert(mat1d_b(1) == 2.0);
}

void matrix_test_slicing_with_slicematrix_b2_04(bool print = false) {
  std::cout << "[TEST]: B2_04. Matrix<T, 1>::subvec(i, j)"
            << " Const version\n";

  Matrix<double, 1> mat1d(arr_1d, 4);
  if (print) test_print(mat1d, "mat1d =");
  mat1d.subvec(2, 3) = 5.0;
  if (print) std::cout << "Apply mat1d.subvec(2, 3) = 5\n";
  if (print) test_print(mat1d, "mat1d =");
  assert(mat1d(0) == 1.0);
  assert(mat1d(1) == 2.0);
  assert(mat1d(2) == 5.0);
  assert(mat1d(3) == 5.0);

  // mat1d.subvec(0, 2) is a SliceMatrix<double>
  const Matrix<double, 1> mat1d_a(mat1d.subvec(0, 2));
  if (print) test_print(mat1d_a, "mat1d_a =");
  assert(mat1d_a.n_elem() == 3);
  assert(mat1d_a(0) == 1.0);
  assert(mat1d_a(1) == 2.0);
  assert(mat1d_a(2) == 5.0);

  // mat1d_a.subvec(0, 2) is a Matrix<T, 1>
  Matrix<double, 1> mat1d_b(mat1d_a.subvec(0, 1));
  if (print) test_print(mat1d_a, "mat1d_b =");
  assert(mat1d_b.n_elem() == 2);
  assert(mat1d_b(0) == 1.0);
  assert(mat1d_b(1) == 2.0);
}

void matrix_test_slicing_with_slicematrix_b2_05(bool print = false) {
  std::cout << "[TEST]: B2_05. Matrix<T, 2>::row()/col()"
            << " Non-const version\n";

  Matrix<double, 2> mat2d_a(arr_2d, 4, 3);

  if (print) test_print(mat2d_a, "mat2d_a =");
  if (print) std::cout << "Apply mat2d_a.row(3) = 0\n";
  mat2d_a.row(3) = 0;
  if (print) test_print(mat2d_a, "mat2d_a =");
  assert(mat2d_a.n_elem() == 12);
  assert(mat2d_a.n_rows() == 4);
  assert(mat2d_a.n_cols() == 3);
  assert(mat2d_a(0, 0) == 1);
  assert(mat2d_a(1, 0) == 2);
  assert(mat2d_a(2, 0) == 3);
  assert(mat2d_a(3, 0) == 0);
  assert(mat2d_a(0, 1) == 5);
  assert(mat2d_a(1, 1) == 6);
  assert(mat2d_a(2, 1) == 7);
  assert(mat2d_a(3, 1) == 0);
  assert(mat2d_a(0, 2) == 9);
  assert(mat2d_a(1, 2) == 10);
  assert(mat2d_a(2, 2) == 11);
  assert(mat2d_a(3, 2) == 0);

  if (print) std::cout << "Apply mat2d_a.col(1) = 8\n";
  mat2d_a.col(1) = 8;
  if (print) test_print(mat2d_a, "mat2d_a =");
  assert(mat2d_a.n_elem() == 12);
  assert(mat2d_a.n_rows() == 4);
  assert(mat2d_a.n_cols() == 3);
  assert(mat2d_a(0, 0) == 1);
  assert(mat2d_a(1, 0) == 2);
  assert(mat2d_a(2, 0) == 3);
  assert(mat2d_a(3, 0) == 0);
  assert(mat2d_a(0, 1) == 8);
  assert(mat2d_a(1, 1) == 8);
  assert(mat2d_a(2, 1) == 8);
  assert(mat2d_a(3, 1) == 8);
  assert(mat2d_a(0, 2) == 9);
  assert(mat2d_a(1, 2) == 10);
  assert(mat2d_a(2, 2) == 11);
  assert(mat2d_a(3, 2) == 0);
}

void matrix_test_slicing_with_slicematrix_b2_06(bool print = false) {
  std::cout << "[TEST]: B2_06. Matrix<T, 2>::row()/col()"
            << " Const version\n";

  const Matrix<double, 2> mat2d_b(arr_2d, 4, 3);
  if (print) test_print(mat2d_b, "mat2d_b =");

  Matrix<double, 1> mat1d_a(mat2d_b.row(3));
  if (print) test_print(mat1d_a, "mat1d_a (mat2d_b.row(3)) =");
  assert(mat1d_a.n_elem() == 3);
  assert(mat1d_a.n_rows() == 1);
  assert(mat1d_a.n_cols() == 3);
  assert(mat1d_a(0) == 4.0);
  assert(mat1d_a(1) == 8.0);
  assert(mat1d_a(2) == 12.0);

  Matrix<double, 1> mat1d_b(mat2d_b.col(1));
  if (print) test_print(mat1d_b, "mat1d_b (mat2d_b.col(1)) =");
  assert(mat1d_b.n_elem() == 4);
  assert(mat1d_b.n_rows() == 4);
  assert(mat1d_b.n_cols() == 1);
  assert(mat1d_b(0) == 5.0);
  assert(mat1d_b(1) == 6.0);
  assert(mat1d_b(2) == 7.0);
  assert(mat1d_b(3) == 8.0);
}

void matrix_test_slicing_with_slicematrix_b2_07(bool print = false) {
  std::cout << "[TEST]: B2_07. Matrix<T, 2>::diag()"
            << " Non-const version\n";

  Matrix<double, 2> mat2d_a(arr_2d, 3, 4);
  mat2d_a.diag(0) = 0;
  mat2d_a.diag(1) = 1;
  mat2d_a.diag(-1) = -1;
  if (print) test_print(mat2d_a, "mat2d_a =");
  assert(mat2d_a.n_elem() == 12);
  assert(mat2d_a.n_rows() == 3);
  assert(mat2d_a.n_cols() == 4);
  assert(mat2d_a(0, 0) == 0);
  assert(mat2d_a(1, 0) == -1);
  assert(mat2d_a(2, 0) == 3);
  assert(mat2d_a(0, 1) == 1);
  assert(mat2d_a(1, 1) == 0);
  assert(mat2d_a(2, 1) == -1);
  assert(mat2d_a(0, 2) == 7);
  assert(mat2d_a(1, 2) == 1);
  assert(mat2d_a(2, 2) == 0);
  assert(mat2d_a(0, 3) == 10);
  assert(mat2d_a(1, 3) == 11);
  assert(mat2d_a(2, 3) == 1);

  Matrix<double, 2> mat2d_b(arr_2d, 4, 3);
  mat2d_b.diag(0) = 0;
  mat2d_b.diag(1) = 1;
  mat2d_b.diag(-1) = -1;
  if (print) test_print(mat2d_b, "mat2d_b =");
  assert(mat2d_b.n_elem() == 12);
  assert(mat2d_b.n_rows() == 4);
  assert(mat2d_b.n_cols() == 3);
  assert(mat2d_b(0, 0) == 0);
  assert(mat2d_b(1, 0) == -1);
  assert(mat2d_b(2, 0) == 3);
  assert(mat2d_b(3, 0) == 4);
  assert(mat2d_b(0, 1) == 1);
  assert(mat2d_b(1, 1) == 0);
  assert(mat2d_b(2, 1) == -1);
  assert(mat2d_b(3, 1) == 8);
  assert(mat2d_b(0, 2) == 9);
  assert(mat2d_b(1, 2) == 1);
  assert(mat2d_b(2, 2) == 0);
  assert(mat2d_b(3, 2) == -1);
}

void matrix_test_slicing_with_slicematrix_b2_08(bool print = false) {
  std::cout << "[TEST]: B2_08. Matrix<T, 2>::diag()"
            << " Const version\n";

  const double arr[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};

  const Matrix<double, 2> mat2d(arr, 3, 3);
  const Matrix<double, 1> mat1d_a = mat2d.diag(0);
  assert(mat1d_a.n_elem() == 3);
  assert(mat1d_a.n_rows() == 3);
  assert(mat1d_a.n_cols() == 1);
  assert(mat1d_a(0) == 1);
  assert(mat1d_a(1) == 5);
  assert(mat1d_a(2) == 9);

  const Matrix<double, 1> mat1d_b = mat2d.diag(1);
  assert(mat1d_b.n_elem() == 2);
  assert(mat1d_b.n_rows() == 2);
  assert(mat1d_b.n_cols() == 1);
  assert(mat1d_b(0) == 4);
  assert(mat1d_b(1) == 8);

  const Matrix<double, 1> mat1d_c = mat2d.diag(-1);
  assert(mat1d_c.n_elem() == 2);
  assert(mat1d_c.n_rows() == 2);
  assert(mat1d_c.n_cols() == 1);
  assert(mat1d_c(0) == 2);
  assert(mat1d_c(1) == 6);
}

void matrix_test_slicing_with_gslicematrix_b3_01(bool print = false) {
  std::cout << "[TEST]: B3_01. m(slice, slice) / m(slice, slice, slice)"
            << " Non-const version\n";
  Matrix<double, 2> mat2d(arr_2d, 4, 3);
  mat2d(std::slice(1, 2, 1), std::slice(0, 2, 2)) = 0;
  if (print) test_print(mat2d, "mat2d = ");
  assert(mat2d.n_elem() == 12);
  assert(mat2d.n_rows() == 4);
  assert(mat2d.n_cols() == 3);
  assert(mat2d(0, 0) == 1);
  assert(mat2d(1, 0) == 0);
  assert(mat2d(2, 0) == 0);
  assert(mat2d(3, 0) == 4);
  assert(mat2d(0, 1) == 5);
  assert(mat2d(1, 1) == 6);
  assert(mat2d(2, 1) == 7);
  assert(mat2d(3, 1) == 8);
  assert(mat2d(0, 2) == 9);
  assert(mat2d(1, 2) == 0);
  assert(mat2d(2, 2) == 0);
  assert(mat2d(3, 2) == 12);

  Matrix<double, 3> mat3d(arr_3d, 4, 3, 2);
  mat3d(std::slice(1, 2, 2), std::slice(0, 2, 2), std::slice(0, 2, 1)) = 0;
  if (print) test_print(mat3d, "mat3d = ");
  assert(mat3d.n_elem() == 24);
  assert(mat3d.n_rows() == 4);
  assert(mat3d.n_cols() == 3);
  assert(mat3d.n_slices() == 2);

  assert(mat3d(0, 0, 0) == 1);
  assert(mat3d(1, 0, 0) == 0);
  assert(mat3d(2, 0, 0) == 3);
  assert(mat3d(3, 0, 0) == 0);
  assert(mat3d(0, 1, 0) == 5);
  assert(mat3d(1, 1, 0) == 6);
  assert(mat3d(2, 1, 0) == 7);
  assert(mat3d(3, 1, 0) == 8);
  assert(mat3d(0, 2, 0) == 9);
  assert(mat3d(1, 2, 0) == 0);
  assert(mat3d(2, 2, 0) == 11);
  assert(mat3d(3, 2, 0) == 0);

  assert(mat3d(0, 0, 1) == 13);
  assert(mat3d(1, 0, 1) == 0);
  assert(mat3d(2, 0, 1) == 15);
  assert(mat3d(3, 0, 1) == 0);
  assert(mat3d(0, 1, 1) == 17);
  assert(mat3d(1, 1, 1) == 18);
  assert(mat3d(2, 1, 1) == 19);
  assert(mat3d(3, 1, 1) == 20);
  assert(mat3d(0, 2, 1) == 21);
  assert(mat3d(1, 2, 1) == 0);
  assert(mat3d(2, 2, 1) == 23);
  assert(mat3d(3, 2, 1) == 0);
}

void matrix_test_slicing_with_gslicematrix_b3_02(bool print = false) {
  std::cout << "[TEST]: B3_02. m(slice, slice) / m(slice, slice, slice)"
            << " Const version\n";
  const Matrix<double, 2> mat2d_a(arr_2d, 4, 3);
  Matrix<double, 2> mat2d_b = mat2d_a(std::slice(1, 2, 1), std::slice(0, 2, 2));
  if (print) test_print(mat2d_b, "mat2d_b = ");
  assert(mat2d_b.n_elem() == 4);
  assert(mat2d_b.n_rows() == 2);
  assert(mat2d_b.n_cols() == 2);
  assert(mat2d_b(0, 0) == 2);
  assert(mat2d_b(1, 0) == 3);
  assert(mat2d_b(0, 1) == 10);

  const Matrix<double, 3> mat3d_a(arr_3d, 4, 3, 2);
  Matrix<double, 3> mat3d_b =
      mat3d_a(std::slice(1, 2, 2), std::slice(0, 2, 2), std::slice(0, 2, 1));
  if (print) test_print(mat3d_b, "mat3d_b = ");
  assert(mat3d_b.n_elem() == 8);
  assert(mat3d_b.n_rows() == 2);
  assert(mat3d_b.n_cols() == 2);
  assert(mat3d_b.n_slices() == 2);
  assert(mat3d_b(0, 0, 0) == 2);
  assert(mat3d_b(1, 0, 0) == 4);
  assert(mat3d_b(0, 1, 0) == 10);
  assert(mat3d_b(1, 1, 0) == 12);
  assert(mat3d_b(0, 0, 1) == 14);
  assert(mat3d_b(1, 0, 1) == 16);
  assert(mat3d_b(0, 1, 1) == 22);
  assert(mat3d_b(1, 1, 1) == 24);
}

void matrix_test_slicing_with_gslicematrix_b3_03(bool print = false) {
  std::cout << "[TEST]: B3_03. m.submat() / m.subcube()"
            << " Non-const version\n";
  Matrix<double, 2> mat2d_a(arr_2d, 4, 3);
  if (print) test_print(mat2d_a, "mat2d_a =");
  if (print) std::cout << "Apply mat2d_a.submat(1,1,3,2) = 0\n";
  mat2d_a.submat(1, 1, 3, 2) = 0;
  if (print) test_print(mat2d_a, "mat2d_a =");
  assert(mat2d_a(0, 0) == 1);
  assert(mat2d_a(0, 1) == 5);
  assert(mat2d_a(0, 2) == 9);
  assert(mat2d_a(3, 0) == 4);
  assert(mat2d_a(3, 1) == 0);
  assert(mat2d_a(3, 2) == 0);

  Matrix<double, 3> mat3d_a(arr_3d_2, 4, 3, 3);

  if (print) test_print(mat3d_a, "mat3d_a =");
  if (print) std::cout << "Apply mat2d_a.submat(1,1,3,2) = 0\n";
  mat3d_a.subcube(1, 0, 1, 2, 2, 2) = 0;
  if (print) test_print(mat3d_a, "mat3d_a =");
  assert(mat3d_a(0, 0, 0) == 1);
  assert(mat3d_a(0, 1, 0) == 5);
  assert(mat3d_a(0, 2, 0) == 9);
  assert(mat3d_a(2, 0, 1) == 0);
  assert(mat3d_a(2, 1, 1) == 0);
  assert(mat3d_a(2, 2, 1) == 0);
}

void matrix_test_slicing_with_gslicematrix_b3_04(bool print = false) {
  std::cout << "[TEST]: B3_04. m.submat() / m.subcube()"
            << " Const version\n";
  const Matrix<double, 2> mat2d_b(arr_2d, 4, 3);
  Matrix<double, 2> mat2d_c(mat2d_b.submat(1, 1, 3, 2));
  if (print) test_print(mat2d_c, "mat2d_c = ");
  assert(mat2d_c.n_elem() == 6);
  assert(mat2d_c.n_rows() == 3);
  assert(mat2d_c.n_cols() == 2);
  assert(mat2d_c(0, 0) == 6);
  assert(mat2d_c(1, 0) == 7);
  assert(mat2d_c(2, 0) == 8);
  assert(mat2d_c(0, 1) == 10);
  assert(mat2d_c(1, 1) == 11);
  assert(mat2d_c(2, 1) == 12);

  const Matrix<double, 2> mat2d_d = mat2d_b.submat(0, 0, 2, 1);
  if (print) test_print(mat2d_d, "mat2d_d = ");
  assert(mat2d_d.n_elem() == 6);
  assert(mat2d_d.n_rows() == 3);
  assert(mat2d_d.n_cols() == 2);
  assert(mat2d_d(0, 0) == 1);
  assert(mat2d_d(1, 0) == 2);
  assert(mat2d_d(2, 0) == 3);
  assert(mat2d_d(0, 1) == 5);
  assert(mat2d_d(1, 1) == 6);
  assert(mat2d_d(2, 1) == 7);

  const Matrix<double, 3> mat3d_b(arr_3d_2, 4, 3, 3);
  Matrix<double, 3> mat3d_c(mat3d_b.subcube(0, 0, 0, 1, 2, 2));
  if (print) test_print(mat3d_c, "mat3d_c = ");
  assert(mat3d_c.n_elem() == 18);
  assert(mat3d_c.n_rows() == 2);
  assert(mat3d_c.n_cols() == 3);
  assert(mat3d_c.n_slices() == 3);
  assert(mat3d_c(0, 0, 0) == 1);
  assert(mat3d_c(0, 1, 0) == 5);
  assert(mat3d_c(0, 2, 0) == 9);
  assert(mat3d_c(1, 0, 2) == 26);
  assert(mat3d_c(1, 1, 2) == 30);
  assert(mat3d_c(1, 2, 2) == 34);

  const Matrix<double, 3> mat3d_d = mat3d_b.subcube(2, 0, 0, 3, 2, 2);
  if (print) test_print(mat3d_d, "mat3d_d = ");
  assert(mat3d_d.n_elem() == 18);
  assert(mat3d_d.n_rows() == 2);
  assert(mat3d_d.n_cols() == 3);
  assert(mat3d_c.n_slices() == 3);
  assert(mat3d_d(0, 0, 0) == 3);
  assert(mat3d_d(0, 1, 0) == 7);
  assert(mat3d_d(0, 2, 0) == 11);
  assert(mat3d_d(1, 0, 2) == 28);
  assert(mat3d_d(1, 1, 2) == 32);
  assert(mat3d_d(1, 2, 2) == 36);
}

void matrix_test_slicing_with_gslicematrix_b3_05(bool print = false) {
  std::cout << "[TEST]: B3_05. m.rows() / m.cols()"
            << " Non-Const version\n";

  Matrix<double, 2> mat2d_a(arr_2d, 4, 3);

  if (print) test_print(mat2d_a, "mat2d_a = ");
  if (print) std::cout << "Apply mat2d_a.rows(2, 3) = 0\n";
  mat2d_a.rows(2, 3) = 0;
  if (print) test_print(mat2d_a, "mat2d_a =");
  assert(mat2d_a.n_elem() == 12);
  assert(mat2d_a.n_rows() == 4);
  assert(mat2d_a.n_cols() == 3);
  assert(mat2d_a(0, 0) == 1);
  assert(mat2d_a(1, 0) == 2);
  assert(mat2d_a(2, 0) == 0);
  assert(mat2d_a(3, 0) == 0);
  assert(mat2d_a(0, 1) == 5);
  assert(mat2d_a(1, 1) == 6);
  assert(mat2d_a(2, 1) == 0);
  assert(mat2d_a(3, 1) == 0);
  assert(mat2d_a(0, 2) == 9);
  assert(mat2d_a(1, 2) == 10);
  assert(mat2d_a(2, 2) == 0);
  assert(mat2d_a(3, 2) == 0);

  Matrix<double, 2> mat2d_b(arr_2d, 4, 3);

  if (print) test_print(mat2d_b, "mat2d_b = ");
  if (print) std::cout << "Apply mat2d_b.cols(1, 2) = 0\n";
  mat2d_b.cols(1, 2) = 0;
  if (print) test_print(mat2d_b, "mat2d_b =");
  assert(mat2d_b.n_elem() == 12);
  assert(mat2d_b.n_rows() == 4);
  assert(mat2d_b.n_cols() == 3);
  assert(mat2d_b(0, 0) == 1);
  assert(mat2d_b(1, 0) == 2);
  assert(mat2d_b(2, 0) == 3);
  assert(mat2d_b(3, 0) == 4);
  assert(mat2d_b(0, 1) == 0);
  assert(mat2d_b(1, 1) == 0);
  assert(mat2d_b(2, 1) == 0);
  assert(mat2d_b(3, 1) == 0);
  assert(mat2d_b(0, 2) == 0);
  assert(mat2d_b(1, 2) == 0);
  assert(mat2d_b(2, 2) == 0);
  assert(mat2d_b(3, 2) == 0);
}

void matrix_test_slicing_with_gslicematrix_b3_06(bool print = false) {
  std::cout << "[TEST]: B3_06. m.rows() / m.cols()"
            << " Const version\n";

  const Matrix<double, 2> mat2d_a(arr_2d, 4, 3);

  const Matrix<double, 2> mat2d_b = mat2d_a.rows(0, 1);
  if (print) test_print(mat2d_b, "mat2d_b = ");
  assert(mat2d_b.n_elem() == 6);
  assert(mat2d_b.n_rows() == 2);
  assert(mat2d_b.n_cols() == 3);
  assert(mat2d_b(0, 0) == 1);
  assert(mat2d_b(0, 1) == 5);
  assert(mat2d_b(0, 2) == 9);
  assert(mat2d_b(1, 0) == 2);
  assert(mat2d_b(1, 1) == 6);
  assert(mat2d_b(1, 2) == 10);

  const Matrix<double, 2> mat2d_c = mat2d_a.cols(0, 1);
  if (print) test_print(mat2d_c, "mat2d_c = ");
  assert(mat2d_c.n_elem() == 8);
  assert(mat2d_c.n_rows() == 4);
  assert(mat2d_c.n_cols() == 2);
  assert(mat2d_c(0, 0) == 1);
  assert(mat2d_c(1, 0) == 2);
  assert(mat2d_c(2, 0) == 3);
  assert(mat2d_c(3, 0) == 4);
  assert(mat2d_c(0, 1) == 5);
  assert(mat2d_c(1, 1) == 6);
  assert(mat2d_c(2, 1) == 7);
  assert(mat2d_c(3, 1) == 8);
}

void matrix_test_slicing_with_gslicematrix_b3_07(bool print = false) {
  std::cout << "[TEST]: B3_07. m.row() / m.col() / m.slice()"
            << " Non-Const version\n";
  Matrix<double, 3> mat3d;
  mat3d = Matrix<double, 3>(arr_3d, 4, 3, 2);

  if (print) test_print(mat3d, "mat3d = ");

  if (print) std::cout << "Apply mat3d.row(1) = 1" << std::endl;
  mat3d.row(1) = 1;
  if (print) test_print(mat3d, "mat3d = ");

  if (print) std::cout << "Apply mat3d.col(1) = 2" << std::endl;
  mat3d.col(1) = 2;
  if (print) test_print(mat3d, "mat3d = ");

  if (print) std::cout << "Apply mat3d.slice(1) = 3" << std::endl;
  mat3d.slice(1) = 3;
  if (print) test_print(mat3d, "mat3d = ");

  assert(mat3d.n_elem() == 24);
  assert(mat3d.n_rows() == 4);
  assert(mat3d.n_cols() == 3);
  assert(mat3d.n_slices() == 2);

  assert(mat3d(0, 0, 0) == 1);
  assert(mat3d(1, 0, 0) == 1);
  assert(mat3d(2, 0, 0) == 3);
  assert(mat3d(3, 0, 0) == 4);
  assert(mat3d(0, 1, 0) == 2);
  assert(mat3d(1, 1, 0) == 2);
  assert(mat3d(2, 1, 0) == 2);
  assert(mat3d(3, 1, 0) == 2);
  assert(mat3d(0, 2, 0) == 9);
  assert(mat3d(1, 2, 0) == 1);
  assert(mat3d(2, 2, 0) == 11);
  assert(mat3d(3, 2, 0) == 12);

  assert(mat3d(0, 0, 1) == 3);
  assert(mat3d(1, 0, 1) == 3);
  assert(mat3d(2, 0, 1) == 3);
  assert(mat3d(3, 0, 1) == 3);
  assert(mat3d(0, 1, 1) == 3);
  assert(mat3d(1, 1, 1) == 3);
  assert(mat3d(2, 1, 1) == 3);
  assert(mat3d(3, 1, 1) == 3);
  assert(mat3d(0, 2, 1) == 3);
  assert(mat3d(1, 2, 1) == 3);
  assert(mat3d(2, 2, 1) == 3);
  assert(mat3d(3, 2, 1) == 3);
}

void matrix_test_slicing_with_gslicematrix_b3_08(bool print = false) {
  std::cout << "[TEST]: B3_08. m.row() / m.col() / m.slice()"
            << " Const version\n";
  const Matrix<double, 3> mat3d = Matrix<double, 3>(arr_3d, 4, 3, 2);
  Matrix<double, 2> mat2d_a = mat3d.row(1);
  Matrix<double, 2> mat2d_b = mat3d.col(1);
  Matrix<double, 2> mat2d_c = mat3d.slice(1);
  if (print) test_print(mat3d, "mat3d = ");
  if (print) test_print(mat2d_a, "mat2d_a = ");
  assert(mat2d_a.n_elem() == 6);
  assert(mat2d_a.n_rows() == 3);
  assert(mat2d_a.n_cols() == 2);
  assert(mat2d_a(0, 0) == 2);
  assert(mat2d_a(1, 0) == 6);
  assert(mat2d_a(2, 0) == 10);
  assert(mat2d_a(0, 1) == 14);
  assert(mat2d_a(1, 1) == 18);
  assert(mat2d_a(2, 1) == 22);

  if (print) test_print(mat2d_b, "mat2d_b = ");
  assert(mat2d_b.n_elem() == 8);
  assert(mat2d_b.n_rows() == 4);
  assert(mat2d_b.n_cols() == 2);
  assert(mat2d_b(0, 0) == 5);
  assert(mat2d_b(1, 0) == 6);
  assert(mat2d_b(2, 0) == 7);
  assert(mat2d_b(3, 0) == 8);
  assert(mat2d_b(0, 1) == 17);
  assert(mat2d_b(1, 1) == 18);
  assert(mat2d_b(2, 1) == 19);
  assert(mat2d_b(3, 1) == 20);

  if (print) test_print(mat2d_c, "mat2d_c = ");
  assert(mat2d_c.n_elem() == 12);
  assert(mat2d_c.n_rows() == 4);
  assert(mat2d_c.n_cols() == 3);
  assert(mat2d_c(0, 0) == 13);
  assert(mat2d_c(1, 0) == 14);
  assert(mat2d_c(2, 0) == 15);
  assert(mat2d_c(3, 0) == 16);
  assert(mat2d_c(0, 1) == 17);
  assert(mat2d_c(1, 1) == 18);
  assert(mat2d_c(2, 1) == 19);
  assert(mat2d_c(3, 1) == 20);
  assert(mat2d_c(0, 2) == 21);
  assert(mat2d_c(1, 2) == 22);
  assert(mat2d_c(2, 2) == 23);
  assert(mat2d_c(3, 2) == 24);
}

void matrix_test_slicing_with_gslicematrix_b3_09(bool print = false) {
  std::cout << "[TEST]: B3_09. m.rows() / m.cols() / m.slices()"
            << " Non-Const version\n";
  Matrix<double, 3> mat3d;
  mat3d = Matrix<double, 3>(arr_3d, 4, 3, 2);

  if (print) test_print(mat3d, "mat3d = ");

  if (print) std::cout << "Apply mat3d.rows(1,2) = 1" << std::endl;
  mat3d.rows(1, 2) = 1;
  if (print) test_print(mat3d, "mat3d = ");

  if (print) std::cout << "Apply mat3d.cols(1,2) = 2" << std::endl;
  mat3d.cols(1, 2) = 2;
  if (print) test_print(mat3d, "mat3d = ");

  if (print) std::cout << "Apply mat3d.slice(1,1) = 3" << std::endl;
  mat3d.slices(1, 1) = 3;
  if (print) test_print(mat3d, "mat3d = ");

  assert(mat3d.n_elem() == 24);
  assert(mat3d.n_rows() == 4);
  assert(mat3d.n_cols() == 3);
  assert(mat3d.n_slices() == 2);

  assert(mat3d(0, 0, 0) == 1);
  assert(mat3d(1, 0, 0) == 1);
  assert(mat3d(2, 0, 0) == 1);
  assert(mat3d(3, 0, 0) == 4);
  assert(mat3d(0, 1, 0) == 2);
  assert(mat3d(1, 1, 0) == 2);
  assert(mat3d(2, 1, 0) == 2);
  assert(mat3d(3, 1, 0) == 2);
  assert(mat3d(0, 2, 0) == 2);
  assert(mat3d(1, 2, 0) == 2);
  assert(mat3d(2, 2, 0) == 2);
  assert(mat3d(3, 2, 0) == 2);

  assert(mat3d(0, 0, 1) == 3);
  assert(mat3d(1, 0, 1) == 3);
  assert(mat3d(2, 0, 1) == 3);
  assert(mat3d(3, 0, 1) == 3);
  assert(mat3d(0, 1, 1) == 3);
  assert(mat3d(1, 1, 1) == 3);
  assert(mat3d(2, 1, 1) == 3);
  assert(mat3d(3, 1, 1) == 3);
  assert(mat3d(0, 2, 1) == 3);
  assert(mat3d(1, 2, 1) == 3);
  assert(mat3d(2, 2, 1) == 3);
  assert(mat3d(3, 2, 1) == 3);
}

void matrix_test_slicing_with_gslicematrix_b3_10(bool print = false) {
  std::cout << "[TEST]: B3_10. m.rows() / m.cols() / m.slices()"
            << " Const version\n";
  const Matrix<double, 3> mat3d = Matrix<double, 3>(arr_3d, 4, 3, 2);
  Matrix<double, 3> mat3d_a = mat3d.rows(1, 2);
  Matrix<double, 3> mat3d_b = mat3d.cols(1, 2);
  Matrix<double, 3> mat3d_c = mat3d.slices(1, 1);
  if (print) test_print(mat3d, "mat3d = ");
  if (print) test_print(mat3d_a, "mat3d_a = ");
  assert(mat3d_a.n_elem() == 12);
  assert(mat3d_a.n_rows() == 2);
  assert(mat3d_a.n_cols() == 3);
  assert(mat3d_a.n_slices() == 2);
  assert(mat3d_a(0, 0, 0) == 2);
  assert(mat3d_a(1, 0, 0) == 3);
  assert(mat3d_a(0, 1, 0) == 6);
  assert(mat3d_a(1, 1, 0) == 7);
  assert(mat3d_a(0, 2, 0) == 10);
  assert(mat3d_a(1, 2, 0) == 11);
  assert(mat3d_a(0, 0, 1) == 14);
  assert(mat3d_a(1, 0, 1) == 15);
  assert(mat3d_a(0, 1, 1) == 18);
  assert(mat3d_a(1, 1, 1) == 19);
  assert(mat3d_a(0, 2, 1) == 22);
  assert(mat3d_a(1, 2, 1) == 23);

  if (print) test_print(mat3d_b, "mat3d_b = ");
  assert(mat3d_b.n_elem() == 16);
  assert(mat3d_b.n_rows() == 4);
  assert(mat3d_b.n_cols() == 2);
  assert(mat3d_b.n_slices() == 2);
  assert(mat3d_b(0, 0, 0) == 5);
  assert(mat3d_b(1, 0, 0) == 6);
  assert(mat3d_b(2, 0, 0) == 7);
  assert(mat3d_b(3, 0, 0) == 8);
  assert(mat3d_b(0, 1, 0) == 9);
  assert(mat3d_b(1, 1, 0) == 10);
  assert(mat3d_b(2, 1, 0) == 11);
  assert(mat3d_b(3, 1, 0) == 12);
  assert(mat3d_b(0, 0, 1) == 17);
  assert(mat3d_b(1, 0, 1) == 18);
  assert(mat3d_b(2, 0, 1) == 19);
  assert(mat3d_b(3, 0, 1) == 20);
  assert(mat3d_b(0, 1, 1) == 21);
  assert(mat3d_b(1, 1, 1) == 22);
  assert(mat3d_b(2, 1, 1) == 23);
  assert(mat3d_b(3, 1, 1) == 24);

  if (print) test_print(mat3d_c, "mat3d_c = ");
  assert(mat3d_c.n_elem() == 12);
  assert(mat3d_c.n_rows() == 4);
  assert(mat3d_c.n_cols() == 3);
  assert(mat3d_c.n_slices() == 1);
  assert(mat3d_c(0, 0, 0) == 13);
  assert(mat3d_c(1, 0, 0) == 14);
  assert(mat3d_c(2, 0, 0) == 15);
  assert(mat3d_c(3, 0, 0) == 16);
  assert(mat3d_c(0, 1, 0) == 17);
  assert(mat3d_c(1, 1, 0) == 18);
  assert(mat3d_c(2, 1, 0) == 19);
  assert(mat3d_c(3, 1, 0) == 20);
  assert(mat3d_c(0, 2, 0) == 21);
  assert(mat3d_c(1, 2, 0) == 22);
  assert(mat3d_c(2, 2, 0) == 23);
  assert(mat3d_c(3, 2, 0) == 24);

  const Matrix<double, 3> mat3d_2(arr_3d_2, 4, 3, 3);
  const Matrix<double, 3> mat3d_d = mat3d_2.slices(1, 2);
  if (print) test_print(mat3d_d, "mat3d_d = ");
  assert(mat3d_d.n_elem() == 24);
  assert(mat3d_d.n_rows() == 4);
  assert(mat3d_d.n_cols() == 3);
  assert(mat3d_d(0, 0, 0) == 13);
  assert(mat3d_d(0, 1, 0) == 17);
  assert(mat3d_d(0, 2, 0) == 21);
  assert(mat3d_d(3, 0, 1) == 28);
  assert(mat3d_d(3, 1, 1) == 32);
  assert(mat3d_d(3, 2, 1) == 36);
}

void matrix_test_slicing_with_maskmatrix_b4_01(bool print = false) {
  std::cout << "[TEST]: B4_01. m(bool_array)"
            << " Non-Const version\n";

  std::valarray<bool> bool_1d(bool_arr_1d, 4);
  Matrix<double, 1> mat1d(arr_1d, 4);
  mat1d(bool_1d) = 0;
  if (print) test_print(mat1d, "mat1d = ");
  assert(mat1d.n_elem() == 4);
  assert(mat1d.n_rows() == 4);
  assert(mat1d.n_cols() == 1);
  assert(mat1d(0) == 1);
  assert(mat1d(1) == 0);
  assert(mat1d(2) == 3);
  assert(mat1d(3) == 0);

  std::valarray<bool> bool_2d(bool_arr_2d, 12);
  Matrix<double, 2> mat2d(arr_2d, 4, 3);
  mat2d(bool_2d) = 0;
  if (print) test_print(mat2d, "mat2d = ");
  assert(mat2d.n_elem() == 12);
  assert(mat2d.n_rows() == 4);
  assert(mat2d.n_cols() == 3);
  assert(mat2d(0, 0) == 1);
  assert(mat2d(1, 0) == 0);
  assert(mat2d(2, 0) == 3);
  assert(mat2d(3, 0) == 0);
  assert(mat2d(0, 1) == 5);
  assert(mat2d(1, 1) == 6);
  assert(mat2d(2, 1) == 7);
  assert(mat2d(3, 1) == 8);
  assert(mat2d(0, 2) == 9);
  assert(mat2d(1, 2) == 0);
  assert(mat2d(2, 2) == 11);
  assert(mat2d(3, 2) == 0);

  std::valarray<bool> bool_3d(bool_arr_3d, 24);
  Matrix<double, 3> mat3d(arr_3d, 4, 3, 2);
  mat3d(bool_3d) = 0;
  if (print) test_print(mat3d, "mat3d = ");
  assert(mat3d.n_elem() == 24);
  assert(mat3d.n_rows() == 4);
  assert(mat3d.n_cols() == 3);
  assert(mat3d.n_slices() == 2);
  assert(mat3d(0, 0, 0) == 1);
  assert(mat3d(1, 0, 0) == 2);
  assert(mat3d(2, 0, 0) == 3);
  assert(mat3d(3, 0, 0) == 4);
  assert(mat3d(0, 1, 0) == 5);
  assert(mat3d(1, 1, 0) == 6);
  assert(mat3d(2, 1, 0) == 7);
  assert(mat3d(3, 1, 0) == 8);
  assert(mat3d(0, 2, 0) == 9);
  assert(mat3d(1, 2, 0) == 10);
  assert(mat3d(2, 2, 0) == 11);
  assert(mat3d(3, 2, 0) == 12);
  assert(mat3d(0, 0, 1) == 13);
  assert(mat3d(1, 0, 1) == 0);
  assert(mat3d(2, 0, 1) == 15);
  assert(mat3d(3, 0, 1) == 0);
  assert(mat3d(0, 1, 1) == 17);
  assert(mat3d(1, 1, 1) == 18);
  assert(mat3d(2, 1, 1) == 19);
  assert(mat3d(3, 1, 1) == 20);
  assert(mat3d(0, 2, 1) == 21);
  assert(mat3d(1, 2, 1) == 0);
  assert(mat3d(2, 2, 1) == 23);
  assert(mat3d(3, 2, 1) == 0);
}

void matrix_test_slicing_with_maskmatrix_b4_02(bool print = false) {
  std::cout << "[TEST]: B4_02. m(bool_array)"
            << " Const version\n";
  Matrix<double, 1> mat1d_a;

  std::valarray<bool> bool_1d(bool_arr_1d, 4);
  const Matrix<double, 1> mat1d(arr_1d, 4);
  mat1d_a = mat1d(bool_1d);
  if (print) test_print(mat1d_a, "mat1d_a = ");
  assert(mat1d_a.n_elem() == 2);
  assert(mat1d_a.n_rows() == 2);
  assert(mat1d_a.n_cols() == 1);
  assert(mat1d_a(0) == 2);
  assert(mat1d_a(1) == 4);

  std::valarray<bool> bool_2d(bool_arr_2d, 12);
  const Matrix<double, 2> mat2d(arr_2d, 4, 3);
  mat1d_a = mat2d(bool_2d);
  if (print) test_print(mat1d_a, "mat1d_a = ");
  assert(mat1d_a.n_elem() == 4);
  assert(mat1d_a.n_rows() == 4);
  assert(mat1d_a.n_cols() == 1);
  assert(mat1d_a(0) == 2);
  assert(mat1d_a(1) == 4);
  assert(mat1d_a(2) == 10);
  assert(mat1d_a(3) == 12);

  std::valarray<bool> bool_3d(bool_arr_3d, 24);
  const Matrix<double, 3> mat3d(arr_3d, 4, 3, 2);
  mat1d_a = mat3d(bool_3d);
  if (print) test_print(mat1d_a, "mat1d_a = ");
  assert(mat1d_a.n_elem() == 4);
  assert(mat1d_a.n_rows() == 4);
  assert(mat1d_a.n_cols() == 1);
  assert(mat1d_a(0) == 14);
  assert(mat1d_a(1) == 16);
  assert(mat1d_a(2) == 22);
  assert(mat1d_a(3) == 24);
}

void matrix_test_slicing_with_indirectmatrix_b5_01(bool print = false) {
  std::cout << "[TEST]: B5_01. m(index_array)"
            << " Non-Const version\n";
  index_array index_1d(index_arr_1d, 2);
  Matrix<double, 1> mat1d(arr_1d, 4);
  mat1d(index_1d) = 0;
  if (print) test_print(mat1d, "mat1d = ");
  assert(mat1d.n_elem() == 4);
  assert(mat1d.n_rows() == 4);
  assert(mat1d.n_cols() == 1);
  assert(mat1d(0) == 1);
  assert(mat1d(1) == 0);
  assert(mat1d(2) == 3);
  assert(mat1d(3) == 0);

  index_array index_2d(index_arr_2d, 4);
  Matrix<double, 2> mat2d(arr_2d, 4, 3);
  mat2d(index_2d) = 0;
  assert(mat2d.n_elem() == 12);
  assert(mat2d.n_rows() == 4);
  assert(mat2d.n_cols() == 3);
  assert(mat2d(0, 0) == 1);
  assert(mat2d(1, 0) == 0);
  assert(mat2d(2, 0) == 3);
  assert(mat2d(3, 0) == 0);
  assert(mat2d(0, 1) == 5);
  assert(mat2d(1, 1) == 6);
  assert(mat2d(2, 1) == 7);
  assert(mat2d(3, 1) == 8);
  assert(mat2d(0, 2) == 9);
  assert(mat2d(1, 2) == 0);
  assert(mat2d(2, 2) == 11);
  assert(mat2d(3, 2) == 0);

  index_array index_3d(index_arr_3d, 4);
  Matrix<double, 3> mat3d(arr_3d, 4, 3, 2);
  mat3d(index_3d) = 0;
  assert(mat3d.n_elem() == 24);
  assert(mat3d.n_rows() == 4);
  assert(mat3d.n_cols() == 3);
  assert(mat3d.n_slices() == 2);
  assert(mat3d(0, 0, 0) == 1);
  assert(mat3d(1, 0, 0) == 2);
  assert(mat3d(2, 0, 0) == 3);
  assert(mat3d(3, 0, 0) == 4);
  assert(mat3d(0, 1, 0) == 5);
  assert(mat3d(1, 1, 0) == 6);
  assert(mat3d(2, 1, 0) == 7);
  assert(mat3d(3, 1, 0) == 8);
  assert(mat3d(0, 2, 0) == 9);
  assert(mat3d(1, 2, 0) == 10);
  assert(mat3d(2, 2, 0) == 11);
  assert(mat3d(3, 2, 0) == 12);
  assert(mat3d(0, 0, 1) == 13);
  assert(mat3d(1, 0, 1) == 0);
  assert(mat3d(2, 0, 1) == 15);
  assert(mat3d(3, 0, 1) == 0);
  assert(mat3d(0, 1, 1) == 17);
  assert(mat3d(1, 1, 1) == 18);
  assert(mat3d(2, 1, 1) == 19);
  assert(mat3d(3, 1, 1) == 20);
  assert(mat3d(0, 2, 1) == 21);
  assert(mat3d(1, 2, 1) == 0);
  assert(mat3d(2, 2, 1) == 23);
  assert(mat3d(3, 2, 1) == 0);
}

void matrix_test_slicing_with_indirectmatrix_b5_02(bool print = false) {
  std::cout << "[TEST]: B5_02. m(index_array)"
            << " Const version\n";
  Matrix<double, 1> mat1d_a;

  index_array index_1d(index_arr_1d, 2);
  const Matrix<double, 1> mat1d(arr_1d, 4);
  mat1d_a = mat1d(index_1d);
  if (print) test_print(mat1d_a, "mat1d_a = ");
  assert(mat1d_a.n_elem() == 2);
  assert(mat1d_a.n_rows() == 2);
  assert(mat1d_a.n_cols() == 1);
  assert(mat1d_a(0) == 2);
  assert(mat1d_a(1) == 4);

  index_array index_2d(index_arr_2d, 4);
  const Matrix<double, 2> mat2d(arr_2d, 4, 3);
  mat1d_a = mat2d(index_2d);
  if (print) test_print(mat1d_a, "mat1d_a = ");
  assert(mat1d_a.n_elem() == 4);
  assert(mat1d_a.n_rows() == 4);
  assert(mat1d_a.n_cols() == 1);
  assert(mat1d_a(0) == 2);
  assert(mat1d_a(1) == 4);
  assert(mat1d_a(2) == 10);
  assert(mat1d_a(3) == 12);

  index_array index_3d(index_arr_3d, 4);
  const Matrix<double, 3> mat3d(arr_3d, 4, 3, 2);
  mat1d_a = mat3d(index_3d);
  if (print) test_print(mat1d_a, "mat1d_a = ");
  assert(mat1d_a.n_elem() == 4);
  assert(mat1d_a.n_rows() == 4);
  assert(mat1d_a.n_cols() == 1);
  assert(mat1d_a(0) == 14);
  assert(mat1d_a(1) == 16);
  assert(mat1d_a(2) == 22);
  assert(mat1d_a(3) == 24);
}

void matrix_test_slicing_with_indirectmatrix_b5_03(bool print = false) {
  std::cout << "[TEST]: B5_03. m(ia1, ia2) / m(ia1, ia2, ia3)"
            << " Non-Const version\n";
  index_array idx_arr1(index_dim1, 3);
  index_array idx_arr2(index_dim2, 2);
  index_array idx_arr3(index_dim3, 1);

  Matrix<double, 2> mat2d(arr_2d, 4, 3);
  mat2d(idx_arr1, idx_arr2) = 0;
  if (print) test_print(mat2d, "mat2d = ");
  assert(mat2d.n_elem() == 12);
  assert(mat2d.n_rows() == 4);
  assert(mat2d.n_cols() == 3);
  assert(mat2d(0, 0) == 0);
  assert(mat2d(1, 0) == 0);
  assert(mat2d(2, 0) == 3);
  assert(mat2d(3, 0) == 0);
  assert(mat2d(0, 1) == 5);
  assert(mat2d(1, 1) == 6);
  assert(mat2d(2, 1) == 7);
  assert(mat2d(3, 1) == 8);
  assert(mat2d(0, 2) == 0);
  assert(mat2d(1, 2) == 0);
  assert(mat2d(2, 2) == 11);
  assert(mat2d(3, 2) == 0);

  Matrix<double, 3> mat3d(arr_3d, 4, 3, 2);
  mat3d(idx_arr1, idx_arr2, idx_arr3) = 0;
  if (print) test_print(mat3d, "mat3d = ");
  assert(mat3d.n_elem() == 24);
  assert(mat3d.n_rows() == 4);
  assert(mat3d.n_cols() == 3);
  assert(mat3d.n_slices() == 2);
  assert(mat3d(0, 0, 0) == 1);
  assert(mat3d(1, 0, 0) == 2);
  assert(mat3d(2, 0, 0) == 3);
  assert(mat3d(3, 0, 0) == 4);
  assert(mat3d(0, 1, 0) == 5);
  assert(mat3d(1, 1, 0) == 6);
  assert(mat3d(2, 1, 0) == 7);
  assert(mat3d(3, 1, 0) == 8);
  assert(mat3d(0, 2, 0) == 9);
  assert(mat3d(1, 2, 0) == 10);
  assert(mat3d(2, 2, 0) == 11);
  assert(mat3d(3, 2, 0) == 12);
  assert(mat3d(0, 0, 1) == 0);
  assert(mat3d(1, 0, 1) == 0);
  assert(mat3d(2, 0, 1) == 15);
  assert(mat3d(3, 0, 1) == 0);
  assert(mat3d(0, 1, 1) == 17);
  assert(mat3d(1, 1, 1) == 18);
  assert(mat3d(2, 1, 1) == 19);
  assert(mat3d(3, 1, 1) == 20);
  assert(mat3d(0, 2, 1) == 0);
  assert(mat3d(1, 2, 1) == 0);
  assert(mat3d(2, 2, 1) == 23);
  assert(mat3d(3, 2, 1) == 0);
}

void matrix_test_slicing_with_indirectmatrix_b5_04(bool print = false) {
  std::cout << "[TEST]: B5_04. m(ia1, ia2) / m(ia1, ia2, ia3)"
            << " Const version\n";
  index_array idx_arr1(index_dim1, 3);
  index_array idx_arr2(index_dim2, 2);
  index_array idx_arr3(index_dim3, 1);

  const Matrix<double, 2> mat2d(arr_2d, 4, 3);
  Matrix<double, 2> mat2d_a = mat2d(idx_arr1, idx_arr2);
  if (print) test_print(mat2d_a, "mat2d_a = ");
  assert(mat2d_a.n_elem() == 6);
  assert(mat2d_a.n_rows() == 3);
  assert(mat2d_a.n_cols() == 2);
  assert(mat2d_a(0, 0) == 1);
  assert(mat2d_a(1, 0) == 2);
  assert(mat2d_a(2, 0) == 4);
  assert(mat2d_a(0, 1) == 9);
  assert(mat2d_a(1, 1) == 10);
  assert(mat2d_a(2, 1) == 12);

  const Matrix<double, 3> mat3d(arr_3d, 4, 3, 2);
  Matrix<double, 3> mat3d_a = mat3d(idx_arr1, idx_arr2, idx_arr3);
  if (print) test_print(mat3d_a, "mat3d_a = ");
  assert(mat3d_a.n_elem() == 6);
  assert(mat3d_a.n_rows() == 3);
  assert(mat3d_a.n_cols() == 2);
  assert(mat3d_a.n_slices() == 1);
  assert(mat3d_a(0, 0, 0) == 13);
  assert(mat3d_a(1, 0, 0) == 14);
  assert(mat3d_a(2, 0, 0) == 16);
  assert(mat3d_a(0, 1, 0) == 21);
  assert(mat3d_a(1, 1, 0) == 22);
  assert(mat3d_a(2, 1, 0) == 24);
}

void matrix_test_scalar_operations_c2_01(bool print = false) {
  std::cout << "[TEST]: C2_01. scalar operations on matrix" << std::endl;
  Matrix<double, 1> mat1d_a(arr_1d, 4);
  Matrix<double, 2> mat2d_a(arr_2d, 4, 3);
  Matrix<double, 3> mat3d_a(arr_3d, 4, 3, 2);

  const Matrix<double, 1> mat1d_b(arr_1d, 4);
  const Matrix<double, 2> mat2d_b(arr_2d, 4, 3);
  const Matrix<double, 3> mat3d_b(arr_3d, 4, 3, 2);

  mat1d_a += 1;
  mat2d_a += 2;
  mat3d_a += 3;

  assert(all(mat1d_a == (mat1d_b + 1.0)));
  assert(all(mat1d_a == (1.0 + mat1d_b)));
  assert(all(mat2d_a == (mat2d_b + 2.0)));
  assert(all(mat2d_a == (2.0 + mat2d_b)));
  assert(all(mat3d_a == (mat3d_b + 3.0)));
  assert(all(mat3d_a == (3.0 + mat3d_b)));
}

void matrix_test_scalar_operations_c2_02(bool print = false) {
  std::cout << "[TEST]: C2_02. scalar operations on slicematrix" << std::endl;
  Matrix<double, 2> mat2d_a(arr_2d, 4, 3);
  if (print) test_print(mat2d_a, "mat2d = ");
  mat2d_a.row(0) += 1;
  mat2d_a.col(2) += 2;
  if (print) test_print(mat2d_a, "mat2d = ");
  assert(mat2d_a.n_elem() == 12);
  assert(mat2d_a.n_rows() == 4);
  assert(mat2d_a.n_cols() == 3);
  assert(mat2d_a(0, 0) == 2);
  assert(mat2d_a(1, 0) == 2);
  assert(mat2d_a(2, 0) == 3);
  assert(mat2d_a(3, 0) == 4);
  assert(mat2d_a(0, 1) == 6);
  assert(mat2d_a(1, 1) == 6);
  assert(mat2d_a(2, 1) == 7);
  assert(mat2d_a(3, 1) == 8);
  assert(mat2d_a(0, 2) == 12);
  assert(mat2d_a(1, 2) == 12);
  assert(mat2d_a(2, 2) == 13);
  assert(mat2d_a(3, 2) == 14);

  const Matrix<double, 2> mat2d_b(arr_2d, 4, 3);

  Matrix<double, 2> mat2d_c = mat2d_b.row(0) + 1.0;
  if (print) test_print(mat2d_c, "mat2d_c = ");
  assert(mat2d_c.n_elem() == 3);
  assert(mat2d_c.n_rows() == 1);
  assert(mat2d_c.n_cols() == 3);
  assert(mat2d_c(0, 0) == 2);
  assert(mat2d_c(0, 1) == 6);
  assert(mat2d_c(0, 2) == 10);

  Matrix<double, 2> mat2d_d = 2.0 + mat2d_b.col(2);
  if (print) test_print(mat2d_d, "mat2d_d = ");
  assert(mat2d_d.n_elem() == 4);
  assert(mat2d_d.n_rows() == 4);
  assert(mat2d_d.n_cols() == 1);
  assert(mat2d_d(0, 0) == 11);
  assert(mat2d_d(1, 0) == 12);
  assert(mat2d_d(2, 0) == 13);
  assert(mat2d_d(3, 0) == 14);
}

void matrix_test_scalar_operations_c2_03(bool print = false) {
  std::cout << "[TEST]: C2_03. scalar operations on gslicematrix" << std::endl;
  Matrix<double, 2> mat2d_a(arr_2d, 4, 3);
  if (print) test_print(mat2d_a, "mat2d = ");
  mat2d_a.rows(0, 1) += 1;
  mat2d_a.cols(1, 2) += 2;
  if (print) test_print(mat2d_a, "mat2d = ");
  assert(mat2d_a.n_elem() == 12);
  assert(mat2d_a.n_rows() == 4);
  assert(mat2d_a.n_cols() == 3);
  assert(mat2d_a(0, 0) == 2);
  assert(mat2d_a(1, 0) == 3);
  assert(mat2d_a(2, 0) == 3);
  assert(mat2d_a(3, 0) == 4);
  assert(mat2d_a(0, 1) == 8);
  assert(mat2d_a(1, 1) == 9);
  assert(mat2d_a(2, 1) == 9);
  assert(mat2d_a(3, 1) == 10);
  assert(mat2d_a(0, 2) == 12);
  assert(mat2d_a(1, 2) == 13);
  assert(mat2d_a(2, 2) == 13);
  assert(mat2d_a(3, 2) == 14);

  const Matrix<double, 2> mat2d_b(arr_2d, 4, 3);

  Matrix<double, 2> mat2d_c = mat2d_b.rows(0, 1) + 1.0;
  if (print) test_print(mat2d_c, "mat2d_c = ");
  assert(mat2d_c.n_elem() == 6);
  assert(mat2d_c.n_rows() == 2);
  assert(mat2d_c.n_cols() == 3);
  assert(mat2d_c(0, 0) == 2);
  assert(mat2d_c(1, 0) == 3);
  assert(mat2d_c(0, 1) == 6);
  assert(mat2d_c(1, 1) == 7);
  assert(mat2d_c(0, 2) == 10);
  assert(mat2d_c(1, 2) == 11);

  Matrix<double, 2> mat2d_d = 2.0 + mat2d_b.cols(1, 2);
  if (print) test_print(mat2d_d, "mat2d_d = ");
  assert(mat2d_d.n_elem() == 8);
  assert(mat2d_d.n_rows() == 4);
  assert(mat2d_d.n_cols() == 2);
  assert(mat2d_d(0, 0) == 7);
  assert(mat2d_d(1, 0) == 8);
  assert(mat2d_d(2, 0) == 9);
  assert(mat2d_d(3, 0) == 10);
  assert(mat2d_d(0, 1) == 11);
  assert(mat2d_d(1, 1) == 12);
  assert(mat2d_d(2, 1) == 13);
  assert(mat2d_d(3, 1) == 14);
}

void matrix_test_unary_add_minus_operator(bool print = false) {
  std::cout << "[TEST]: Applies unary add/minus operators to each element\n";

  const double a1[] = {1, 2, 3, 4};
  const std::valarray<double> va1(a1, 4);
  Matrix<double, 1> mat1d(va1);
  Matrix<double, 1> mat1d_a = +mat1d;
  Matrix<double, 1> mat1d_b = -mat1d;

  if (print) test_print(mat1d_a, "mat1d_a =");
  assert(mat1d_a(0) == 1.0);
  assert(mat1d_a(1) == 2.0);
  assert(mat1d_a(2) == 3.0);
  assert(mat1d_a(3) == 4.0);
  if (print) test_print(mat1d_b, "mat1d_b =");
  assert(mat1d_b(0) == -1.0);
  assert(mat1d_b(1) == -2.0);
  assert(mat1d_b(2) == -3.0);
  assert(mat1d_b(3) == -4.0);

  const double a2[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
  const std::valarray<double> va2(a2, 12);
  Matrix<double, 2> mat2d(va2, 4, 3);
  Matrix<double, 2> mat2d_a = +mat2d;
  Matrix<double, 2> mat2d_b = -mat2d;

  if (print) test_print(mat2d_a, "mat2d_a =");
  assert(mat2d_a(0, 0) == 1.0);
  assert(mat2d_a(1, 0) == 2.0);
  assert(mat2d_a(2, 0) == 3.0);
  assert(mat2d_a(3, 0) == 4.0);
  assert(mat2d_a(0, 2) == 9.0);
  assert(mat2d_a(1, 2) == 10.0);
  assert(mat2d_a(2, 2) == 11.0);
  assert(mat2d_a(3, 2) == 12.0);
  if (print) test_print(mat2d_b, "mat2d_b =");
  assert(mat2d_b(0, 0) == -1.0);
  assert(mat2d_b(1, 0) == -2.0);
  assert(mat2d_b(2, 0) == -3.0);
  assert(mat2d_b(3, 0) == -4.0);
  assert(mat2d_b(0, 2) == -9.0);
  assert(mat2d_b(1, 2) == -10.0);
  assert(mat2d_b(2, 2) == -11.0);
  assert(mat2d_b(3, 2) == -12.0);

  const double a3[] = {1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12,
                       13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24};
  const std::valarray<double> va3(a3, 24);
  Matrix<double, 3> mat3d(va3, 4, 3, 2);
  Matrix<double, 3> mat3d_a = +mat3d;
  Matrix<double, 3> mat3d_b = -mat3d;

  if (print) test_print(mat3d_a, "mat3d_a =");
  assert(mat3d_a(0, 0, 0) == 1.0);
  assert(mat3d_a(1, 0, 0) == 2.0);
  assert(mat3d_a(2, 0, 0) == 3.0);
  assert(mat3d_a(3, 0, 0) == 4.0);
  assert(mat3d_a(0, 2, 1) == 21.0);
  assert(mat3d_a(1, 2, 1) == 22.0);
  assert(mat3d_a(2, 2, 1) == 23.0);
  assert(mat3d_a(3, 2, 1) == 24.0);
  if (print) test_print(mat3d_b, "mat3d_b =");
  assert(mat3d_b(0, 0, 0) == -1.0);
  assert(mat3d_b(1, 0, 0) == -2.0);
  assert(mat3d_b(2, 0, 0) == -3.0);
  assert(mat3d_b(3, 0, 0) == -4.0);
  assert(mat3d_b(0, 2, 1) == -21.0);
  assert(mat3d_b(1, 2, 1) == -22.0);
  assert(mat3d_b(2, 2, 1) == -23.0);
  assert(mat3d_b(3, 2, 1) == -24.0);
}

void matrix_test_addition_assignment_operator(bool print = false) {
  std::cout
      << "[TEST]: Applies addition assignment operators to each element\n";

  const double a1[] = {1, 2, 3, 4};
  Matrix<double, 1> mat1d_a(a1, 4);
  Matrix<double, 1> mat1d_b(a1, 4);
  Matrix<double, 1> mat1d_c(a1, 4);

  if (print) test_print(mat1d_a, "mat1d_a = ");
  if (print) std::cout << "Apply mat1d_a += 1\n";
  mat1d_a += 1;
  if (print) test_print(mat1d_a, "mat1d_a = ");
  assert(mat1d_a(0) == 2.0);
  assert(mat1d_a(1) == 3.0);
  assert(mat1d_a(2) == 4.0);
  assert(mat1d_a(3) == 5.0);

  if (print) test_print(mat1d_b, "mat1d_b = ");
  if (print) std::cout << "Apply mat1d_b += mat1d_a\n";
  mat1d_b += mat1d_a;
  if (print) test_print(mat1d_b, "mat1d_b = ");
  assert(mat1d_b(0) == 3.0);
  assert(mat1d_b(1) == 5.0);
  assert(mat1d_b(2) == 7.0);
  assert(mat1d_b(3) == 9.0);

  if (print) test_print(mat1d_c, "mat1d_c = ");
  if (print) std::cout << "Apply mat1d_c = 1 + mat1d_c\n";
  mat1d_c = 1.0 + mat1d_c;
  if (print) test_print(mat1d_c, "mat1d_c = ");
  assert(mat1d_c(0) == 2.0);
  assert(mat1d_c(1) == 3.0);
  assert(mat1d_c(2) == 4.0);
  assert(mat1d_c(3) == 5.0);
}

void matrix_test_subtraction_assignment_operator(bool print = false) {
  std::cout
      << "[TEST]: Applies subtraction assignment operators to each element\n";

  const double a1[] = {1, 2, 3, 4};
  Matrix<double, 1> mat1d_a(a1, 4);
  Matrix<double, 1> mat1d_b(a1, 4);
  Matrix<double, 1> mat1d_c(a1, 4);

  if (print) test_print(mat1d_a, "mat1d_a = ");
  if (print) std::cout << "Apply mat1d_a -= 1\n";
  mat1d_a -= 1;
  if (print) test_print(mat1d_a, "mat1d_a = ");
  assert(mat1d_a(0) == 0.0);
  assert(mat1d_a(1) == 1.0);
  assert(mat1d_a(2) == 2.0);
  assert(mat1d_a(3) == 3.0);

  if (print) test_print(mat1d_b, "mat1d_b = ");
  if (print) std::cout << "Apply mat1d_b -= mat1d_a\n";
  mat1d_b -= mat1d_a;
  if (print) test_print(mat1d_b, "mat1d_b = ");
  assert(mat1d_b(0) == 1.0);
  assert(mat1d_b(1) == 1.0);
  assert(mat1d_b(2) == 1.0);
  assert(mat1d_b(3) == 1.0);

  if (print) test_print(mat1d_c, "mat1d_c = ");
  if (print) std::cout << "Apply mat1d_c = 1-mat1d_c\n";
  mat1d_c = 1.0 - mat1d_c;
  if (print) test_print(mat1d_c, "mat1d_c = ");
  assert(mat1d_c(0) == 0.0);
  assert(mat1d_c(1) == -1.0);
  assert(mat1d_c(2) == -2.0);
  assert(mat1d_c(3) == -3.0);
}

void matrix_test_multiplication_assignment_operator(bool print = false) {
  std::cout << "[TEST]: Applies multiplication assignment operators to each "
               "element\n";

  const double a1[] = {1, 2, 3, 4};
  Matrix<double, 1> mat1d_a(a1, 4);
  Matrix<double, 1> mat1d_b(a1, 4);
  Matrix<double, 1> mat1d_c(a1, 4);

  if (print) test_print(mat1d_a, "mat1d_a = ");
  if (print) std::cout << "Apply mat1d_a *= 2\n";
  mat1d_a *= 2;
  if (print) test_print(mat1d_a, "mat1d_a = ");
  assert(mat1d_a(0) == 2.0);
  assert(mat1d_a(1) == 4.0);
  assert(mat1d_a(2) == 6.0);
  assert(mat1d_a(3) == 8.0);

  if (print) test_print(mat1d_b, "mat1d_b = ");
  if (print) std::cout << "Apply mat1d_b *= mat1d_a\n";
  mat1d_b *= mat1d_a;
  if (print) test_print(mat1d_b, "mat1d_b = ");
  assert(mat1d_b(0) == 2.0);
  assert(mat1d_b(1) == 8.0);
  assert(mat1d_b(2) == 18.0);
  assert(mat1d_b(3) == 32.0);

  if (print) test_print(mat1d_c, "mat1d_c = ");
  if (print) std::cout << "Apply mat1d_c = 2.0 * mat1d_c\n";
  mat1d_c = 2.0 * mat1d_c;
  if (print) test_print(mat1d_c, "mat1d_c = ");
  assert(mat1d_c(0) == 2.0);
  assert(mat1d_c(1) == 4.0);
  assert(mat1d_c(2) == 6.0);
  assert(mat1d_c(3) == 8.0);
}

void matrix_test_division_assignment_operator(bool print = false) {
  std::cout
      << "[TEST]: Applies division assignment operators to each element\n";

  const double a1[] = {1, 2, 3, 4};
  Matrix<double, 1> mat1d_a(a1, 4);
  Matrix<double, 1> mat1d_b(a1, 4);
  Matrix<double, 1> mat1d_c(a1, 4);

  if (print) test_print(mat1d_a, "mat1d_a = ");
  if (print) std::cout << "Apply mat1d_a /= 2\n";
  mat1d_a /= 2;
  if (print) test_print(mat1d_a, "mat1d_a = ");
  assert(mat1d_a(0) == 0.5);
  assert(mat1d_a(1) == 1.0);
  assert(mat1d_a(2) == 1.5);
  assert(mat1d_a(3) == 2.0);

  if (print) test_print(mat1d_b, "mat1d_b = ");
  if (print) std::cout << "Apply mat1d_b /= mat1d_a\n";
  mat1d_b /= mat1d_a;
  if (print) test_print(mat1d_b, "mat1d_b = ");
  assert(mat1d_b(0) == 2.0);
  assert(mat1d_b(1) == 2.0);
  assert(mat1d_b(2) == 2.0);
  assert(mat1d_b(3) == 2.0);

  if (print) test_print(mat1d_c, "mat1d_c = ");
  if (print) std::cout << "Apply mat1d_c = 2 / mat1d_c\n";
  mat1d_c = 2.0 / mat1d_c;
  if (print) test_print(mat1d_c, "mat1d_c = ");
  assert(mat1d_c(0) == 2.0);
  assert(mat1d_c(1) == 1.0);
  assert(std::abs(mat1d_c(2) - 2.0 / 3) < 1e-5);
  assert(mat1d_c(3) == 0.5);
}

void matrix_test_member_function_sum(bool print = false) {
  std::cout << "[TEST]: Calculates the sum of all elements\n";

  const double a1[] = {1, 2, 3, 4};
  Matrix<double, 1> mat1d(a1, 4);
  if (print) {
    test_print(mat1d, "mat1d =");
    std::cout << "mat1d.sum() = " << mat1d.sum() << std::endl;
  }
  assert(mat1d.sum() == 10.0);
  assert(mat1d.subvec(1, 2).sum() == 5.0);

  const double a2[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
  Matrix<double, 2> mat2d(a2, 4, 3);
  if (print) {
    test_print(mat2d, "mat2d =");
    std::cout << "mat2d.sum() = " << mat2d.sum() << std::endl;
  }
  assert(mat2d.sum() == 78.0);
  assert(mat2d.col(0).sum() == 10.0);
  assert(mat2d.row(0).sum() == 15.0);

  const double a3[] = {1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12,
                       13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24};
  Matrix<double, 3> mat3d(a3, 4, 3, 2);
  if (print) {
    test_print(mat3d, "mat3d =");
    std::cout << "mat3d.sum() = " << mat3d.sum() << std::endl;
  }
  assert(mat3d.sum() == 300.0);
}

void matrix_test_member_function_min(bool print = false) {
  std::cout << "[TEST]: Gets the smallest element\n";

  const double a1[] = {1, 2, 3, 4};
  Matrix<double, 1> mat1d(a1, 4);
  if (print) {
    test_print(mat1d, "mat1d =");
    std::cout << "mat1d.min() = " << mat1d.min() << std::endl;
  }
  assert(mat1d.min() == 1.0);
  assert(mat1d.subvec(1, 3).min() == 2.0);

  const double a2[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
  Matrix<double, 2> mat2d(a2, 4, 3);
  if (print) {
    test_print(mat2d, "mat2d =");
    std::cout << "mat2d.min() = " << mat2d.min() << std::endl;
  }
  assert(mat2d.min() == 1.0);
  assert(mat2d.col(1).min() == 5.0);

  const double a3[] = {1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12,
                       13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24};
  Matrix<double, 3> mat3d(a3, 4, 3, 2);
  if (print) {
    test_print(mat3d, "mat3d =");
    std::cout << "mat3d.min() = " << mat3d.min() << std::endl;
  }
  assert(mat3d.min() == 1.0);
}

void matrix_test_member_function_max(bool print = false) {
  std::cout << "[TEST]: Gets the largest element\n";

  const double a1[] = {1, 2, 3, 4};
  Matrix<double, 1> mat1d(a1, 4);
  if (print) {
    test_print(mat1d, "mat1d =");
    std::cout << "mat1d.max() = " << mat1d.max() << std::endl;
  }
  assert(mat1d.max() == 4.0);
  assert(mat1d.subvec(0, 2).max() == 3.0);

  const double a2[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
  Matrix<double, 2> mat2d(a2, 4, 3);
  if (print) {
    test_print(mat2d, "mat2d =");
    std::cout << "mat2d.max() = " << mat2d.max() << std::endl;
  }
  assert(mat2d.max() == 12.0);
  assert(mat2d.col(1).max() == 8.0);

  const double a3[] = {1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12,
                       13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24};
  Matrix<double, 3> mat3d(a3, 4, 3, 2);
  if (print) {
    test_print(mat3d, "mat3d =");
    std::cout << "mat3d.max() = " << mat3d.max() << std::endl;
  }
  assert(mat3d.max() == 24.0);
}

#if defined MATRIX_LIB_USE_CPP11
auto fun = [](double& val) { val += 10.0; };
#else
void fun(double& val) { val += 10.0; }
#endif

void matrix_test_member_function_for_each(bool print = false) {
  std::cout << "[TEST]: Applys function to all element\n";

  const double a1[] = {1, 2, 3, 4};
  Matrix<double, 1> mat1d(a1, 4);
  if (print) test_print(mat1d, "mat1d =");
  mat1d.for_each(fun);
  if (print) test_print(mat1d, "mat1d =");
  assert(mat1d(0) == 11.0);
  assert(mat1d(1) == 12.0);
  assert(mat1d(2) == 13.0);
  assert(mat1d(3) == 14.0);

  const double a2[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
  Matrix<double, 2> mat2d(a2, 4, 3);
  if (print) test_print(mat2d, "mat2d =");
  mat2d.for_each(fun);
  if (print) test_print(mat2d, "mat2d =");
  assert(mat2d(0, 0) == 11.0);
  assert(mat2d(1, 0) == 12.0);
  assert(mat2d(2, 0) == 13.0);
  assert(mat2d(3, 0) == 14.0);
  assert(mat2d(0, 2) == 19.0);
  assert(mat2d(1, 2) == 20.0);
  assert(mat2d(2, 2) == 21.0);
  assert(mat2d(3, 2) == 22.0);

  const double a3[] = {1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12,
                       13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24};
  Matrix<double, 3> mat3d(a3, 4, 3, 2);
  if (print) test_print(mat3d, "mat3d =");
  mat3d.for_each(fun);
  if (print) test_print(mat3d, "mat3d =");
  assert(mat3d(0, 0, 0) == 11.0);
  assert(mat3d(1, 0, 0) == 12.0);
  assert(mat3d(2, 0, 0) == 13.0);
  assert(mat3d(3, 0, 0) == 14.0);
  assert(mat3d(0, 2, 1) == 31.0);
  assert(mat3d(1, 2, 1) == 32.0);
  assert(mat3d(2, 2, 1) == 33.0);
  assert(mat3d(3, 2, 1) == 34.0);
}

void matrix_1d_test_subvec(bool print = false) {
  std::cout << "[TEST]: 1D Matrix's member functions subvec\n";

  const double a1[] = {1, 2, 3, 4};
  Matrix<double, 1> mat1d(a1, 4);
  if (print) test_print(mat1d, "mat1d =");
  mat1d.subvec(2, 3) = 5.0;
  if (print) std::cout << "Apply mat1d.subvec(2, 3) = 5\n";
  if (print) test_print(mat1d, "mat1d =");
  assert(mat1d(0) == 1.0);
  assert(mat1d(1) == 2.0);
  assert(mat1d(2) == 5.0);
  assert(mat1d(3) == 5.0);

  // mat1d.subvec(0, 2) is a MatrixRef<T, 1>
  const Matrix<double, 1> mat1d_a(mat1d.subvec(0, 2));
  if (print) test_print(mat1d_a, "mat1d_a =");
  assert(mat1d_a.n_elem() == 3);
  assert(mat1d_a(0) == 1.0);
  assert(mat1d_a(1) == 2.0);
  assert(mat1d_a(2) == 5.0);

  // mat1d_a.subvec(0, 2) is a Matrix<T, 1>
  Matrix<double, 1> mat1d_b(mat1d_a.subvec(0, 1));
  if (print) test_print(mat1d_a, "mat1d_b =");
  assert(mat1d_b.n_elem() == 2);
  assert(mat1d_b(0) == 1.0);
  assert(mat1d_b(1) == 2.0);
}

void matrix_1d_test_transpose(bool print = false) {
  std::cout << "[TEST]: 1D Matrix's member functions transpose\n";

  const double arr[] = {1, 2, 3, 4};
  Matrix<double, 1> mat1d(arr, 4);
  Matrix<double, 2> mat2d = mat1d.t();

  if (print) test_print(mat1d, "mat1d = ");
  if (print) test_print(mat2d, "mat2d = ");

  assert(mat2d.n_elem() == 4);
  assert(mat2d.n_rows() == 1);
  assert(mat2d.n_cols() == 4);
  assert(mat2d(0, 0) == 1);
  assert(mat2d(0, 1) == 2);
  assert(mat2d(0, 2) == 3);
  assert(mat2d(0, 3) == 4);
}

void matrix_2d_test_transpose(bool print = false) {
  std::cout << "[TEST]: 2D Matrix's member functions transpose\n";

  const double arr[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
  Matrix<double, 2> mat2d_a(arr, 4, 3);
  Matrix<double, 2> mat2d_b = mat2d_a.t();

  if (print) test_print(mat2d_a, "mat2d_a = ");
  if (print) test_print(mat2d_b, "mat2d_b = ");

  assert(mat2d_b.n_elem() == 12);
  assert(mat2d_b.n_rows() == 3);
  assert(mat2d_b.n_cols() == 4);
  assert(mat2d_b(0, 0) == 1);
  assert(mat2d_b(1, 0) == 5);
  assert(mat2d_b(2, 0) == 9);
  assert(mat2d_b(0, 3) == 4);
  assert(mat2d_b(1, 3) == 8);
  assert(mat2d_b(2, 3) == 12);
}

void matrix_test_binary_addition_operator(bool print = false) {
  std::cout << "[TEST]: Applies binary addition operators to each element\n";

  const double a1[] = {1, 2, 3, 4};
  Matrix<double, 1> mat1d_a(a1, 4);
  Matrix<double, 1> mat1d_b(a1, 4);

  if (print) test_print(mat1d_a, "mat1d_a = ");
  if (print) std::cout << "Apply mat1d_a += 1\n";
  mat1d_a = mat1d_a + 1.0;
  if (print) test_print(mat1d_a, "mat1d_a = ");
  assert(mat1d_a(0) == 2.0);
  assert(mat1d_a(1) == 3.0);
  assert(mat1d_a(2) == 4.0);
  assert(mat1d_a(3) == 5.0);

  if (print) test_print(mat1d_b, "mat1d_b = ");
  if (print) std::cout << "Apply mat1d_b += mat1d_a\n";
  mat1d_b = mat1d_b + mat1d_a;
  if (print) test_print(mat1d_b, "mat1d_b = ");
  assert(mat1d_b(0) == 3.0);
  assert(mat1d_b(1) == 5.0);
  assert(mat1d_b(2) == 7.0);
  assert(mat1d_b(3) == 9.0);
}

void matrix_test_binary_subtraction_operator(bool print = false) {
  std::cout << "[TEST]: Applies binary subtraction operators to each element\n";

  const double a1[] = {1, 2, 3, 4};
  Matrix<double, 1> mat1d_a(a1, 4);
  Matrix<double, 1> mat1d_b(a1, 4);

  if (print) test_print(mat1d_a, "mat1d_a = ");
  if (print) std::cout << "Apply mat1d_a = mat1d_a - 1.0\n";
  mat1d_a = mat1d_a - 1.0;
  if (print) test_print(mat1d_a, "mat1d_a = ");
  assert(mat1d_a(0) == 0.0);
  assert(mat1d_a(1) == 1.0);
  assert(mat1d_a(2) == 2.0);
  assert(mat1d_a(3) == 3.0);

  if (print) test_print(mat1d_b, "mat1d_b = ");
  if (print) std::cout << "Apply mat1d_b = mat1d_b - mat1d_a\n";
  mat1d_b = mat1d_b - mat1d_a;
  if (print) test_print(mat1d_b, "mat1d_b = ");
  assert(mat1d_b(0) == 1.0);
  assert(mat1d_b(1) == 1.0);
  assert(mat1d_b(2) == 1.0);
  assert(mat1d_b(3) == 1.0);
}

void matrix_test_binary_equal_operator(bool print = false) {
  std::cout << "[TEST]: Matrix1 == Matrix2\n";
  const double a1[6] = {1, 2, 3, 4, 5, 6};
  const double a2[6] = {1, 2, 3, 7, 8, 9};
  Matrix<double, 2> m1(a1, 3, 2);
  Matrix<double, 2> m2(a1, 3, 2);
  Matrix<double, 2> m3(a2, 3, 2);
  bool_array ba1 = (m1 == m2);
  bool_array ba2 = (m1 == m3);
  assert(ba1[0] == true);
  assert(ba1[1] == true);
  assert(ba1[2] == true);
  assert(ba1[3] == true);
  assert(ba1[4] == true);
  assert(ba1[5] == true);

  assert(ba2[0] == true);
  assert(ba2[1] == true);
  assert(ba2[2] == true);
  assert(ba2[3] == false);
  assert(ba2[4] == false);
  assert(ba2[5] == false);
}

void matrix_test_matmul1(bool print = false) {
  std::cout << "[TEST]: Applies matrix multiplication 1\n";
  const double a1[] = {1, 2, 3};
  const double a2[] = {1, 2};
  Matrix<double, 1> m1(a1, 3);
  Matrix<double, 2> m2(a2, 1, 2);
  Matrix<double, 2> res = matmul(m1, m2);

  assert(res.n_elem() == 6);
  assert(res.n_rows() == 3);
  assert(res.n_cols() == 2);
  assert(res(0, 0) == 1);
  assert(res(0, 1) == 2);
  assert(res(2, 0) == 3);
  assert(res(2, 1) == 6);
}

void matrix_test_matmul2(bool print = false) {
  std::cout << "[TEST]: Applies matrix multiplication 2\n";
  const double a1[] = {1, 4, 2, 5, 3, 6};
  const double a2[] = {1, 2, 3};
  Matrix<double, 2> m1(a1, 2, 3);
  Matrix<double, 1> m2(a2, 3);
  Matrix<double, 1> res = matmul(m1, m2);

  assert(res(0) == 14.0);
  assert(res(1) == 32.0);
}

void matrix_test_matmul3(bool print = false) {
  std::cout << "[TEST]: Applies matrix multiplication 3\n";
  const double a1[] = {1, 4, 2, 5, 3, 6};
  const double a2[] = {7, 9, 11, 8, 10, 12};
  Matrix<double, 2> m1(a1, 2, 3);
  Matrix<double, 2> m2(a2, 3, 2);
  Matrix<double, 2> res = matmul(m1, m2);

  assert(res(0, 0) == 58.0);
  assert(res(0, 1) == 64.0);
  assert(res(1, 0) == 139.0);
  assert(res(1, 1) == 154.0);
}

void matrix_test_abs(bool print = false) {
  std::cout << "[TEST]: Applies the function abs to each element\n";

  const double a1[] = {1, -2, 3, -4};
  Matrix<double, 1> mat1d(a1, 4);
  if (print) test_print(mat1d, "mat1d =");
  Matrix<double, 1> mat1d_abs = abs(mat1d);
  if (print) test_print(mat1d_abs, "mat1d_abs =");
  assert(mat1d_abs(0) == 1.0);
  assert(mat1d_abs(1) == 2.0);
  assert(mat1d_abs(2) == 3.0);
  assert(mat1d_abs(3) == 4.0);

  const double a2[] = {1, 2, 3, 4, -5, -6, -7, -8, 9, 10, 11, 12};
  Matrix<double, 2> mat2d(a2, 4, 3);
  if (print) test_print(mat2d, "mat2d =");
  Matrix<double, 2> mat2d_abs = abs(mat2d);
  if (print) test_print(mat2d_abs, "mat2d_abs =");
  assert(mat2d_abs(0, 0) == 1.0);
  assert(mat2d_abs(1, 0) == 2.0);
  assert(mat2d_abs(2, 0) == 3.0);
  assert(mat2d_abs(3, 0) == 4.0);
  assert(mat2d_abs(0, 2) == 9.0);
  assert(mat2d_abs(1, 2) == 10.0);
  assert(mat2d_abs(2, 2) == 11.0);
  assert(mat2d_abs(3, 2) == 12.0);

  const double a3[] = {1,   2,   3,   4,   5,   6,   7,   8,
                       9,   10,  11,  12,  -13, -14, -15, -16,
                       -17, -18, -19, -20, -21, -22, -23, -24};
  Matrix<double, 3> mat3d(a3, 4, 3, 2);
  if (print) test_print(mat3d, "mat3d =");
  Matrix<double, 3> mat3d_abs = abs(mat3d);
  if (print) test_print(mat3d_abs, "mat3d_abs =");
  assert(mat3d_abs(0, 0, 0) == 1.0);
  assert(mat3d_abs(1, 0, 0) == 2.0);
  assert(mat3d_abs(2, 0, 0) == 3.0);
  assert(mat3d_abs(3, 0, 0) == 4.0);
  assert(mat3d_abs(0, 2, 1) == 21.0);
  assert(mat3d_abs(1, 2, 1) == 22.0);
  assert(mat3d_abs(2, 2, 1) == 23.0);
  assert(mat3d_abs(3, 2, 1) == 24.0);
}

void matrix_test_exp(bool print = false) {
  std::cout << "[TEST]: Applies the function exp to each element\n";

  const double a1[] = {1, 2, 3, 4};
  Matrix<double, 1> mat1d(a1, 4);
  if (print) test_print(mat1d, "mat1d =");
  Matrix<double, 1> mat1d_exp = exp(mat1d);
  if (print) test_print(mat1d_exp, "mat1d_exp =");
  assert(std::abs(mat1d_exp(0) - std::exp(1.0)) < 1e-5);
  assert(std::abs(mat1d_exp(1) - std::exp(2.0)) < 1e-5);
  assert(std::abs(mat1d_exp(2) - std::exp(3.0)) < 1e-5);
  assert(std::abs(mat1d_exp(3) - std::exp(4.0)) < 1e-5);

  const double a2[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
  Matrix<double, 2> mat2d(a2, 4, 3);
  if (print) test_print(mat2d, "mat2d =");
  Matrix<double, 2> mat2d_exp = exp(mat2d);
  if (print) test_print(mat2d_exp, "mat2d_exp =");
  assert(std::abs(mat2d_exp(0, 0) - std::exp(1.0)) < 1e-5);
  assert(std::abs(mat2d_exp(1, 0) - std::exp(2.0)) < 1e-5);
  assert(std::abs(mat2d_exp(2, 0) - std::exp(3.0)) < 1e-5);
  assert(std::abs(mat2d_exp(3, 0) - std::exp(4.0)) < 1e-5);
  assert(std::abs(mat2d_exp(0, 2) - std::exp(9.0)) < 1e-5);
  assert(std::abs(mat2d_exp(1, 2) - std::exp(10.0)) < 1e-5);
  assert(std::abs(mat2d_exp(2, 2) - std::exp(11.0)) < 1e-5);
  assert(std::abs(mat2d_exp(3, 2) - std::exp(12.0)) < 1e-5);

  const double a3[] = {1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12,
                       13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24};
  Matrix<double, 3> mat3d(a3, 4, 3, 2);
  if (print) test_print(mat3d, "mat3d =");
  Matrix<double, 3> mat3d_exp = exp(mat3d);
  if (print) test_print(mat3d_exp, "mat3d_exp =");
  assert(std::abs(mat3d_exp(0, 0, 0) - std::exp(1.0)) < 1e-5);
  assert(std::abs(mat3d_exp(1, 0, 0) - std::exp(2.0)) < 1e-5);
  assert(std::abs(mat3d_exp(2, 0, 0) - std::exp(3.0)) < 1e-5);
  assert(std::abs(mat3d_exp(3, 0, 0) - std::exp(4.0)) < 1e-5);
  assert(std::abs(mat3d_exp(0, 2, 1) - std::exp(21.0)) < 1e-5);
  assert(std::abs(mat3d_exp(1, 2, 1) - std::exp(22.0)) < 1e-5);
  assert(std::abs(mat3d_exp(2, 2, 1) - std::exp(23.0)) < 1e-5);
  assert(std::abs(mat3d_exp(3, 2, 1) - std::exp(24.0)) < 1e-5);
}

void matrix_test_dot(bool print = false) {
  std::cout << "[TEST]: Applies the function dot to two vectors\n";

  const double a1[] = {1, 2, 3, 4};
  Matrix<double, 1> mat1d_a(a1, 4);
  Matrix<double, 1> mat1d_b = mat1d_a;
  double res = dot(mat1d_a, mat1d_b);

  if (print) {
    test_print(mat1d_a, "mat1d_a =");
    test_print(mat1d_b, "mat1d_b =");
    std::cout << "dot(mat1d_a, mat1d_b) = " << res << std::endl;
  }
  assert(res == 30);
}

void matrix_et_add(bool print = false) {
  std::cout << "[TEST]: ET Addition\n";

  Matrix<double, 2> a(arr_2d, 4, 3);
  Matrix<double, 2> b(arr_2d, 4, 3);
  Matrix<double, 2> c(arr_2d, 4, 3);
  Matrix<double, 2> d(a.e() + b.e() + c.e(), 4, 3);
  a *= 3.0;
  bool_array ba = (a == d);

  if (print) {
    test_print(a, "a = ");
    test_print(d, "d = ");
    std::cout << "ba = " << std::endl;
    for (uword i = 0; i != ba.size(); ++i) {
      std::cout << i << " " << ba[i] << std::endl;
    }
    std::cout << "all(ba) = " << all(ba) << std::endl;
  }

  assert(a.n_elem() == d.n_elem());
  assert(a.n_rows() == d.n_rows());
  assert(a.n_cols() == d.n_cols());
  assert(all(ba) == true);
}

void print_msg(const std::string& msg) {
  std::cout << "\n----- " << msg << " -----\n" << std::endl;
}

int main() {
  vec a;
  mat b;
  cube c;

  bool print_flag = false;
  macro_test_init_array(print_flag);

  print_msg("A: Constructions and Assignments");

  print_msg("A1: Testing Matrix<T, N> Constructions");
  matrix_test_constructor_a1_01(print_flag);
  matrix_test_constructor_a1_02(print_flag);
  matrix_test_constructor_a1_03(print_flag);
  matrix_test_constructor_a1_04(print_flag);
  matrix_test_constructor_a1_05(print_flag);
  matrix_test_constructor_a1_06(print_flag);
  matrix_test_constructor_a1_07(print_flag);
  matrix_test_constructor_a1_08(print_flag);
  matrix_test_constructor_a1_09(print_flag);
  matrix_test_constructor_a1_10(print_flag);
  matrix_test_constructor_a1_11(print_flag);
  matrix_test_constructor_a1_12(print_flag);
  matrix_test_constructor_a1_13(print_flag);
  matrix_test_constructor_a1_14(print_flag);

  print_msg("A2: Testing Matrix<T, N> Assignments");
  matrix_test_assignment_a2_01(print_flag);
  matrix_test_assignment_a2_02(print_flag);
  matrix_test_assignment_a2_03(print_flag);
  matrix_test_assignment_a2_04(print_flag);
  matrix_test_assignment_a2_05(print_flag);
  matrix_test_assignment_a2_06(print_flag);
  matrix_test_assignment_a2_07(print_flag);
  matrix_test_assignment_a2_08(print_flag);
  matrix_test_assignment_a2_09(print_flag);
  matrix_test_assignment_a2_10(print_flag);

  print_msg("B: Subscripting and Slicing");
  print_msg("B1: Testing Matrix<T, N> Subscripting");
  // test for b1_01 (subscripting m(i, j, k)) is skipped
  matrix_test_element_access_b1_02(print_flag);
  matrix_2d_test_gslice_array(print_flag);

  print_msg("B2: Testing Matrix<T, N> Slicing with SliceMatrix<T>");
  matrix_test_slicing_with_slicematrix_b2_01(print_flag);
  matrix_test_slicing_with_slicematrix_b2_02(print_flag);
  matrix_test_slicing_with_slicematrix_b2_03(print_flag);
  matrix_test_slicing_with_slicematrix_b2_04(print_flag);
  matrix_test_slicing_with_slicematrix_b2_05(print_flag);
  matrix_test_slicing_with_slicematrix_b2_06(print_flag);
  matrix_test_slicing_with_slicematrix_b2_07(print_flag);
  matrix_test_slicing_with_slicematrix_b2_08(print_flag);

  print_msg("B3: Testing Matrix<T, N> Slicing with GSliceMatrix<T>");
  matrix_test_slicing_with_gslicematrix_b3_01(print_flag);
  matrix_test_slicing_with_gslicematrix_b3_02(print_flag);
  matrix_test_slicing_with_gslicematrix_b3_03(print_flag);
  matrix_test_slicing_with_gslicematrix_b3_04(print_flag);
  matrix_test_slicing_with_gslicematrix_b3_05(print_flag);
  matrix_test_slicing_with_gslicematrix_b3_06(print_flag);
  matrix_test_slicing_with_gslicematrix_b3_07(print_flag);
  matrix_test_slicing_with_gslicematrix_b3_08(print_flag);
  matrix_test_slicing_with_gslicematrix_b3_09(print_flag);
  matrix_test_slicing_with_gslicematrix_b3_10(print_flag);

  print_msg("B4: Testing Matrix<T, N> Slicing with MaskMatrix<T>");
  matrix_test_slicing_with_maskmatrix_b4_01(print_flag);
  matrix_test_slicing_with_maskmatrix_b4_02(print_flag);

  print_msg("B5: Testing Matrix<T, N> Slicing with IndirectMatrix<T>");
  matrix_test_slicing_with_indirectmatrix_b5_01(print_flag);
  matrix_test_slicing_with_indirectmatrix_b5_02(print_flag);
  matrix_test_slicing_with_indirectmatrix_b5_03(print_flag);
  matrix_test_slicing_with_indirectmatrix_b5_04(print_flag);

  print_msg("C: Testing Arithmetic Operations");

  print_msg("C2: Testing Scalar Operations");
  matrix_test_scalar_operations_c2_01(print_flag);
  matrix_test_scalar_operations_c2_02(print_flag);
  matrix_test_scalar_operations_c2_03(print_flag);

  matrix_test_unary_add_minus_operator(print_flag);
  matrix_test_addition_assignment_operator(print_flag);
  matrix_test_subtraction_assignment_operator(print_flag);
  matrix_test_multiplication_assignment_operator(print_flag);
  matrix_test_division_assignment_operator(print_flag);

  matrix_test_binary_addition_operator(print_flag);
  matrix_test_binary_subtraction_operator(print_flag);
  matrix_test_binary_equal_operator(print_flag);

  matrix_test_member_function_sum(print_flag);
  matrix_test_member_function_min(print_flag);
  matrix_test_member_function_max(print_flag);
  matrix_test_member_function_for_each(print_flag);

  matrix_test_abs(print_flag);
  matrix_test_exp(print_flag);
  matrix_test_dot(print_flag);

  std::cout << "\n----- Testing Matrix Transpose -----\n" << std::endl;

  matrix_1d_test_transpose(print_flag);
  matrix_2d_test_transpose(print_flag);

  std::cout << "\n----- Testing Matrix Multiplication -----\n" << std::endl;

  matrix_test_matmul1(print_flag);
  matrix_test_matmul2(print_flag);
  matrix_test_matmul3(print_flag);

  matrix_et_add(print_flag);

  std::cout << "\n----- All tests are done... -----\n" << std::endl;
  return 0;
}
