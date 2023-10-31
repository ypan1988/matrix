matrix.h: A small multidimensional matrix library
==

## Backgroud

**C++** does not come with its matrix library. This is partially true -- the Standard Template Library (STL) has introduced `std::valarray` for fast mathematical computations since **C++98**. The `std::valarray` by itself only acts like a 1D array, but it can be used to simulate a **N**-dimensional matrix quite easily with its helper classes (i.e., `std::slice_array`, `std::gslice_array`, `std::mask_array` and `std::indirect_array`) which have the reference semantics to a subset of the array.

This single header matrix library can be viewed as a wrapper of `std::valarray` and an extension of Bjarne Stroustrup's [matrix](https://www.stroustrup.com/Programming/Matrix/Matrix.h) implementation. It provides standard building blocks for performing basic vector, matrix and cube operations. Most of the APIs are coming from **Armadillo** (another **C++** library for linear algebra with syntax similar to **MATLAB**), but note that it will not be compatible for sure. While other APIs have their origins in a few different programming languages (e.g., **Fortran**/**Python**/**R**).

## Classes

`Matrix<T, N>` is an **N**-dimensional dense matrix of some value type **T** (**YP**: currently only 1D/2D/3D matrix is implemented). Internally it consists of a `std::valarray<T>` for storing its elements and a `std::size_t` array of size **N** to record its sizes on each dimension.

+ **Elements** are stored in [column-major ordering](https://en.wikipedia.org/wiki/Row-_and_column-major_order#Column-major_order) (like **Fortran**), and can be accessed by `get_elem()` (i.e., convert a `Matrix<T, N>` back to a `std::valarray<T>`):
  ```cpp
  // define an empty 2D matrix with T = double
  Matrix<double, 2> x;
  // ...
  // extract a read-only array of elements
  const std::valarray<double> &va1 = x.get_elem();
  // make a copy for the array of elements 
  std::valarray<double> va2 = x.get_elem();
  ```   

+ **Dimemsions** of `Matrix<T, N>` can be obtained through the following member functions (all the returned values have type `std::size_t`):
  | `Matrix<T, 1>` | `Matrix<T, 2>` | `Matrix<T, 3>` | Description |
  | :------------: | :------------: | :------------: | :---------- |
  | `x.n_elem()`   | `x.n_elem()`   | `x.n_elem()`   | total number of elements in x |
  | `x.n_rows()`   | `x.n_rows()`   | `x.n_rows()`   | number of rows in x           |
  |                | `x.n_cols()`   | `x.n_cols()`   | number of columns in x        |
  |                |                | `x.n_slices()` | number of slices in x         |
  ```cpp
  // define a 2x3 matrix
  Matrix<double, 2> x(2, 3);
  std::cout << x.n_elem() << std::endl;    // 6
  std::cout << x.n_rows() << std::endl;    // 2
  std::cout << x.n_cols() << std::endl;    // 3
  ```
+ For convenience the following typedefs have been defined for 1D/2D/3D matrix:

  | type  T  | `Matrix<T, 1>` | `Matrix<T, 2>` | `Matrix<T, 3>` |
  | :------: | :------------: | :------------: | :------------: |
  | `double` | `vec`          | `mat`          | `cube`         |
  | `double` | `dvec`         | `dmat`         | `dcube`        |
  | `float`  | `fvec`         | `fmat `        | `fcube`        |

  For example, we assume that **type T** is `double` in the first row, then obviously
  + `vec` is a 1D matrix (i.e., `vec = Matrix<double, 1>`)
  + `mat` is a 2D matrix (i.e., `mat = Matrix<double, 2>`)
  + `cube` is a 3D matrix (i.e., `cube = Matrix<double, 3>`)

  Note that in the following part of this document, we assume that `T = double` by default and `vec/mat/cube` are used for better readability, but it is possible to use other types.

  ## Constructors, Destructors and Assignments
  | `Matrix<T, 1>` | `Matrix<T, 2>` | `Matrix<T, 3>` | |
  | :------------  | :------------  | :------------  |-|
  | `vec()`                      | `mat()`                                 | `cube()`                                           |(1)|
  | `explicit vec(n_rows)`       | `mat(n_rows, n_cols)`                   | `cube(n_rows, n_cols, n_slices)`                   |(2)|
  | `vec(const T& val, n_rows)`  | `mat(const T& val, n_rows, n_cols)`     | `cube(const T& val, n_rows, n_cols, n_slides)`     |(3)|
  | `vec(const T* vals, n_rows)` | `mat(const T* vals, n_rows, n_cols)`    | `cube(const T* vals, n_rows, n_cols, n_slides)`    |(4)|
  | `vec(const vec& other)`      | `mat(const mat& other)`                 | `cube(const cube& other)`                          |(5)|
  | `vec(vec&& other) noexcept`  | `mat(mat&& other) noexcept`             | `cube(cube&& other) noexcept`                      |(6) C++11|
  | `vec(valarray)`              | `mat(valarray, n_rows, n_cols)`         | `cube(valarray, n_rows, n_cols, n_slides)`         |(7)|
  | `vec(initializer_list)`      | `mat(initializer_list, n_rows, n_cols)` | `cube(initializer_list, n_rows, n_cols, n_slices)` |(8) C++11|
  | `vec(mat)`                   | `mat(vec)`                              | NA                                                 |(9)|
  | `~vec()`                     | `~mat()`                                | `~cube()`                                          |(10)|

  The table above provides ways to construct new matrix from various sources:
  
  1) Default constructor. Constructs an empty `vec/mat/cube`
  2) Constructs a `vec/mat/cube` with the specified dimensions
  3) Constructs a `vec/mat/cube` with the specified dimensions, all elements set to `val`
  4) Constructs a `vec/mat/cube` with the contents of array pointed by `vals` and specified dimensions
  5) Copy constructor. Constructs an `vec/mat/cube` from the other one
  6) Move constructor. Constructs an `vec/mat/cube` from the other one using move semantics
  7) Constructs a `vec/mat/cube` with the contents of the `valarray` and specified dimensions
  8) Constructs a `vec/mat/cube` with the contents of the `initializer_list` and specified dimensions
  9) Constructs a `vec` from a **n x 1** `mat`; Constructs a **n x 1** `mat` from a `vec`
  10) Destructs the `vec/mat/cube`. The destructors of the elements (if **T** is a class) are called and the used storage is deallocated.


## Element Access and Matrix Slicing
