#include "matrix.h"

#include <cassert>
#include <iostream>
#include <string>
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
  assert(mat1d(mat1d.n_elem() - 1) == 1.0);

  Matrix<double, 2> mat2d(2.0, 4, 3);
  if (print) test_print(mat2d, "mat2d =");
  assert(mat2d(0, 0) == 2.0);
  assert(mat2d(mat2d.n_rows() - 1, mat2d.n_cols() - 1) == 2.0);

  Matrix<double, 3> mat3d(3.0, 4, 3, 2);
  if (print) test_print(mat3d, "mat3d =");
  assert(mat3d(0, 0, 0) == 3.0);
  assert(mat3d(mat3d.n_rows() - 1, mat3d.n_cols() - 1, mat3d.n_slices() - 1) ==
         3.0);
}

int main() {
  bool print_flag = false;
  matrix_test_constructor_01(print_flag);
  matrix_test_constructor_02(print_flag);
  matrix_test_constructor_03(print_flag);

  return 0;
}