#include "matrix.h"

#include <cassert>
#include <iostream>
#include <string>
#include <valarray>

using namespace matrix_lib;

void test_print(const Matrix<double, 1>& m, std::string msg = "") {
  if (!msg.empty()) std::cout << msg << std::endl;
  for (uword r = 0; r != m.n_rows(); ++r) {
    std::cout << m(r);
    std::cout << '\n';
  }
}

void test_print(const Matrix<double, 2>& m, std::string msg = "") {
  if (!msg.empty()) std::cout << msg << std::endl;
  for (uword r = 0; r != m.n_rows(); ++r) {
    for (uword c = 0; c != m.n_cols(); ++c) {
      std::cout << m(r, c) << ' ';
    }
    std::cout << '\n';
  }
}

void test_print(const Matrix<double, 3>& m, std::string msg = "") {
  if (!msg.empty()) std::cout << msg << std::endl;
  for (uword s = 0; s != m.n_slices(); ++s) {
    std::cout << "slice " << s << ":" << std::endl;
    for (uword r = 0; r != m.n_rows(); ++r) {
      for (uword c = 0; c != m.n_cols(); ++c) {
        std::cout << m(r, c, s) << ' ';
      }
      std::cout << '\n';
    }
    std::cout << '\n';
  }
}

void matrix_test_constructor_01(bool print = false) {
  std::cout << "[TEST]: Constructs an empty Matrix\n";

  Matrix<double, 1> mat1d;
  if (print) test_print(mat1d, "mat1d =");
  assert(mat1d.n_elem() == 0);
  assert(mat1d.n_rows() == 0);

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

void matrix_test_constructor_02(bool print = false) {
  std::cout << "[TEST]: Constructs a Matrix with n_rows, n_cols and n_slices\n";

  Matrix<double, 1> mat1d(4);
  if (print) test_print(mat1d, "mat1d =");
  assert(mat1d.n_elem() == 4);
  assert(mat1d.n_rows() == 4);

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

void matrix_test_constructor_03(bool print = false) {
  std::cout
      << "[TEST]: Constructs a Matrix with val, n_rows, n_cols and n_slices\n";

  Matrix<double, 1> mat1d(1.0, 4);
  if (print) test_print(mat1d, "mat1d =");
  assert(mat1d(0) == 1.0);
  assert(mat1d(3) == 1.0);

  Matrix<double, 2> mat2d(2.0, 4, 3);
  if (print) test_print(mat2d, "mat2d =");
  assert(mat2d(0, 0) == 2.0);
  assert(mat2d(3, 2) == 2.0);

  Matrix<double, 3> mat3d(3.0, 4, 3, 2);
  if (print) test_print(mat3d, "mat3d =");
  assert(mat3d(0, 0, 0) == 3.0);
  assert(mat3d(3, 2, 1) == 3.0);
}

void matrix_test_constructor_04(bool print = false) {
  std::cout << "[TEST]: Constructs a Matrix with array, n_rows, n_cols and "
               "n_slices\n";

  const double a1[] = {1, 2, 3, 4};
  Matrix<double, 1> mat1d(a1, 4);
  if (print) test_print(mat1d, "mat1d =");
  assert(mat1d(0) == 1.0);
  assert(mat1d(1) == 2.0);
  assert(mat1d(2) == 3.0);
  assert(mat1d(3) == 4.0);

  const double a2[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
  Matrix<double, 2> mat2d(a2, 4, 3);
  if (print) test_print(mat2d, "mat2d =");
  assert(mat2d(0, 0) == 1.0);
  assert(mat2d(1, 0) == 2.0);
  assert(mat2d(2, 0) == 3.0);
  assert(mat2d(3, 0) == 4.0);
  assert(mat2d(0, 2) == 9.0);
  assert(mat2d(1, 2) == 10.0);
  assert(mat2d(2, 2) == 11.0);
  assert(mat2d(3, 2) == 12.0);

  const double a3[] = {1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12,
                       13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24};
  Matrix<double, 3> mat3d(a3, 4, 3, 2);
  if (print) test_print(mat3d, "mat3d =");
  assert(mat3d(0, 0, 0) == 1.0);
  assert(mat3d(1, 0, 0) == 2.0);
  assert(mat3d(2, 0, 0) == 3.0);
  assert(mat3d(3, 0, 0) == 4.0);
  assert(mat3d(0, 2, 1) == 21.0);
  assert(mat3d(1, 2, 1) == 22.0);
  assert(mat3d(2, 2, 1) == 23.0);
  assert(mat3d(3, 2, 1) == 24.0);
}

void matrix_test_constructor_05(bool print = false) {
  std::cout << "[TEST]: Constructs a Matrix with valarray, n_rows, n_cols and "
               "n_slices\n";

  const double a1[] = {1, 2, 3, 4};
  const std::valarray<double> va1(a1, 4);
  Matrix<double, 1> mat1d(a1, 4);
  if (print) test_print(mat1d, "mat1d =");
  assert(mat1d(0) == 1.0);
  assert(mat1d(1) == 2.0);
  assert(mat1d(2) == 3.0);
  assert(mat1d(3) == 4.0);

  const double a2[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
  const std::valarray<double> va2(a2, 12);
  Matrix<double, 2> mat2d(va2, 4, 3);
  if (print) test_print(mat2d, "mat2d =");
  assert(mat2d(0, 0) == 1.0);
  assert(mat2d(1, 0) == 2.0);
  assert(mat2d(2, 0) == 3.0);
  assert(mat2d(3, 0) == 4.0);
  assert(mat2d(0, 2) == 9.0);
  assert(mat2d(1, 2) == 10.0);
  assert(mat2d(2, 2) == 11.0);
  assert(mat2d(3, 2) == 12.0);

  const double a3[] = {1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12,
                       13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24};
  const std::valarray<double> va3(a3, 24);
  Matrix<double, 3> mat3d(va3, 4, 3, 2);
  if (print) test_print(mat3d, "mat3d =");
  assert(mat3d(0, 0, 0) == 1.0);
  assert(mat3d(1, 0, 0) == 2.0);
  assert(mat3d(2, 0, 0) == 3.0);
  assert(mat3d(3, 0, 0) == 4.0);
  assert(mat3d(0, 2, 1) == 21.0);
  assert(mat3d(1, 2, 1) == 22.0);
  assert(mat3d(2, 2, 1) == 23.0);
  assert(mat3d(3, 2, 1) == 24.0);
}

int main() {
  bool print_flag = false;
  matrix_test_constructor_01(print_flag);
  matrix_test_constructor_02(print_flag);
  matrix_test_constructor_03(print_flag);
  matrix_test_constructor_04(print_flag);
  matrix_test_constructor_05(print_flag);

  return 0;
}