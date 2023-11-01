matrix.h: A small multidimensional matrix library
==

## Backgroud

**C++** does not come with its own matrix library. This is partially true -- the **Standard Template Library (STL)** has introduced `std::valarray` for fast mathematical computations since **C++98**. The `std::valarray` by itself only acts like a 1D array, but it can be used to simulate a **N**-dimensional matrix quite easily with its helper classes (i.e., `std::slice_array`, `std::gslice_array`, `std::mask_array` and `std::indirect_array`) which have the reference semantics to a subset of the array.

This single header matrix library can be viewed as a wrapper of `std::valarray` and an extension of Bjarne Stroustrup's [matrix](https://www.stroustrup.com/Programming/Matrix/Matrix.h) implementation. It provides standard building blocks for performing basic vector, matrix and cube operations. Most of the APIs are coming from **Armadillo** (another **C++** library for linear algebra with syntax similar to **MATLAB**), but note that it will not be compatible for sure. While other APIs have their origins in a few different programming languages (e.g., **Fortran**/**Python**/**R**).

## Classes

`Matrix<T, N>` is an **N**-dimensional dense matrix of some value type **T** (**YP**: currently only 1D/2D/3D matrix are implemented). Internally it consists of a `std::valarray<T>` for storing its elements in [column-major ordering](https://en.wikipedia.org/wiki/Row-_and_column-major_order#Column-major_order) (like **Fortran**) and a `std::size_t` array of size **N** to record its sizes on each dimension:

+ **Elements** (or the raw `std::valarray<T>`) can be accessed as a whole by `get_elem()` (i.e., convert a `Matrix<T, N>` back to a `std::valarray<T>`):
  ```cpp
  // define an empty 2D matrix with T = double
  Matrix<double, 2> x;
  // ...
  // make a copy for the array of elements 
  std::valarray<double> va2 = x.get_elem();
  // extract a read-only array of elements
  const std::valarray<double> &va1 = x.get_elem();
  ```
  See also the section **Element Access and Matrix Slicing**

+ **Dimemsions** of `Matrix<T, N>` can be accessed similarly by `get_dims()` while total size and size on each dimension can be obtained through the following member functions (all the returned values have type `std::size_t`):
  | `Matrix<T, 1>` | `Matrix<T, 2>` | `Matrix<T, 3>` | Description |
  | :------------: | :------------: | :------------: | :---------- |
  | `x.n_elem()`   | `x.n_elem()`   | `x.n_elem()`   | total number of elements in x |
  | `x.n_rows()`   | `x.n_rows()`   | `x.n_rows()`   | number of rows in x           |
  |                | `x.n_cols()`   | `x.n_cols()`   | number of columns in x        |
  |                |                | `x.n_slices()` | number of slices in x         |
  ```cpp
  // define a 2x3 matrix
  Matrix<double, 2> x(2, 3);
  // get the dimensions
  std::valarray<std::size_t> dims = x.get_dims();  // {2, 3}
  std::cout << x.n_elem() << std::endl;            // 6
  std::cout << x.n_rows() << std::endl;            // 2
  std::cout << x.n_cols() << std::endl;            // 3
  ```
For convenience the following typedefs have been provided for 1D/2D/3D matrix:

  | type  T  | `Matrix<T, 1>` | `Matrix<T, 2>` | `Matrix<T, 3>` |
  | :------: | :------------: | :------------: | :------------: |
  | `double` | `vec`          | `mat`          | `cube`         |
  | `double` | `dvec`         | `dmat`         | `dcube`        |
  | `float`  | `fvec`         | `fmat `        | `fcube`        |

  For example, we assume that type **T** is `double` in the first row, then obviously
  + `vec` is a 1D matrix (i.e., `vec = Matrix<double, 1>`)
  + `mat` is a 2D matrix (i.e., `mat = Matrix<double, 2>`)
  + `cube` is a 3D matrix (i.e., `cube = Matrix<double, 3>`)

  Note that in the following part of this document, we assume that
  + `T = double` by default, and `elem_type` should always be the same as type **T**
  + `vec/mat/cube` are used for better readability, but it is possible to use other types
  + Anything marked as **C++11** can be used only when the compiler supports **C++11** (**YP**: Since **std::valarray** is included in **STL** since **C++98**, this matrix library, as a wrapper of **std::valarray**, will support **pre-C++11** compilers) 

## Constructions and Assignments

### Constructors and Destructors

| `Matrix<double, 1>` | `Matrix<double, 2>` | `Matrix<double, 3>` | |
| :------------  | :------------  | :------------  |-|
| `vec()`                      | `mat()`                                 | `cube()`                                           |(1)|
| `explicit vec(n_rows)`       | `mat(n_rows, n_cols)`                   | `cube(n_rows, n_cols, n_slices)`                   |(2)|
| `vec(const elem_type& val, n_rows)`  | `mat(const elem_type& val, n_rows, n_cols)`  | `cube(const elem_type& val, n_rows, n_cols, n_slides)`  |(3)|
| `vec(const elem_type* vals, n_rows)` | `mat(const elem_type* vals, n_rows, n_cols)` | `cube(const elem_type* vals, n_rows, n_cols, n_slides)` |(4)|
| `vec(const vec& other)`      | `mat(const mat& other)`                 | `cube(const cube& other)`                          |(5)|
| `vec(vec&& other) noexcept`  | `mat(mat&& other) noexcept`             | `cube(cube&& other) noexcept`                      |(6) C++11|
| `vec(valarray)`              | `mat(valarray, n_rows, n_cols)`         | `cube(valarray, n_rows, n_cols, n_slides)`         |(7)|
| `vec(initializer_list)`      | `mat(initializer_list, n_rows, n_cols)` | `cube(initializer_list, n_rows, n_cols, n_slices)` |(8) C++11|
| `vec(mat)`                   | `mat(vec)`                              | NA                                                 |(9)|
| `~vec()`                     | `~mat()`                                | `~cube()`                                          |(10)|

The table above provides ways to construct new matrix from various sources:
  1) Default constructor. Constructs an empty `vec/mat/cube`.
  2) Constructs a `vec/mat/cube` with the specified dimensions.
  3) Constructs a `vec/mat/cube` with the specified dimensions, all elements set to `val`.
  4) Constructs a `vec/mat/cube` with the contents of array pointed by `vals` and specified dimensions.
  5) Copy constructor. Constructs an `vec/mat/cube` from the other one.
  6) Move constructor. Constructs an `vec/mat/cube` from the other one using move semantics.
  7) Constructs a `vec/mat/cube` with the contents of the `valarray` and specified dimensions.
  8) Constructs a `vec/mat/cube` with the contents of the `initializer_list` and specified dimensions.
  9) Constructs a `vec` from a **n x 1** `mat`; Constructs a **n x 1** `mat` from a `vec`.
  10) Destructs the `vec/mat/cube`. The destructors of the elements (if **T** is a class) are called and the used storage is deallocated.

### Assignments
| `Matrix<T, 1>` | `Matrix<T, 2>` | `Matrix<T, 3>` | |
| :------------  | :------------  | :------------  |-|
| `vec& operator=(const vec& other)` | `mat& operator=(const mat& other)` | `cube& operator=(const cube& other)` |(1)|
| `vec& operator=(vec&& other)`      | `mat& operator=(mat&& other)`      | `cube& operator=(cube&& other)`      |(2) C++11|
| `vec& operator=(const mat& other)` | `mat& operator=(const vec& other)` | NA                                   |(3)|
| `vec& operator=(mat&& other)`      | `mat& operator=(vec&& other)`      | NA                                   |(4) C++11|
| `vec& operator=(const elem_type& val)` | `mat& operator=(const elem_type& val)` | `cube& operator=(const elem_type& val)` |(5)|

The table above provides ways to replace the contents of the matrix:
  1) Copy assignment operator.
  2) Move assignment operator.
  3) Copy assignment operator for the assignment between a `vec` and a **n x 1** `mat`.
  4) Move assignment operator for the assignment between a `vec` and a **n x 1** `mat`.
  5) Replaces each value in `*this` with a copy of `val`.

## Element Access and Matrix Slicing
