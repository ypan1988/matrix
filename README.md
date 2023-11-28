matrix.h: A small multidimensional matrix library
==

## Introduction

**C++** does not come with its own matrix library. This is partially true -- the **Standard Template Library (STL)** has introduced [`std::valarray`](https://en.cppreference.com/w/cpp/numeric/valarray) for fast mathematical computations since **C++98**. The `std::valarray` by itself only acts like a 1D array, but it can be used to simulate a **N**-dimensional matrix quite easily with its helper classes (i.e., `std::slice_array`, `std::gslice_array`, `std::mask_array` and `std::indirect_array`) which have the reference semantics to a subset of the array.

This single header matrix library can be viewed as a wrapper of `std::valarray` and its helper classes, and an extension of Bjarne Stroustrup's [matrix](https://www.stroustrup.com/Programming/Matrix/Matrix.h) implementation. It provides standard building blocks for performing basic vector, matrix and cube operations. Most of the APIs are coming from **Armadillo** (a popular **C++** library for linear algebra with syntax similar to **MATLAB**), but note that it will not be compatible for sure. While other APIs have their origins in a few different programming languages (e.g., **Fortran**).

## A Matrix Template

### Matrix<T, N>

`Matrix<T, N>` is an **N**-dimensional dense matrix of some value type **T** (**YP**: currently only 1D/2D/3D matrix are implemented). Internally it consists of a `std::valarray<T>` for storing its elements in [column-major ordering](https://en.wikipedia.org/wiki/Row-_and_column-major_order#Column-major_order) (like **Fortran**) and a `std::size_t` array of size **N** to record its sizes on each dimension:

+ **Elements** (or the raw `std::valarray<T>`) can be accessed as a whole by `get_elem()` (i.e., convert a `Matrix<T, N>` back to a `std::valarray<T>`):
  ```cpp
  // define an empty 2D matrix with T = double
  Matrix<double, 2> x;
  // ...
  // make a copy for the array of elements 
  std::valarray<double> va1 = x.get_elem();
  // extract a read-only array of elements
  const std::valarray<double> &va2 = x.get_elem();
  ```

+ **Dimensions** of `Matrix<T, N>` can be accessed similarly by `get_dims()` while total size and size on each dimension can be obtained through the following member functions (all the returned values have type `std::size_t`):

+ | `Matrix<T, 1>` | `Matrix<T, 2>` | `Matrix<T, 3>` | Description |
  | :------------: | :------------: | :------------: | :---------- |
  | `x.n_elem()`   | `x.n_elem()`   | `x.n_elem()`   | total number of elements in x |
  | `x.n_rows()`   | `x.n_rows()`   | `x.n_rows()`   | number of rows in x           |
  | `x.n_cols()`   | `x.n_cols()`   | `x.n_cols()`   | number of columns in x        |
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

For convenience the following typedefs of `Matrix<T, N>` have been provided for 1D/2D/3D matrix:

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

### sub-Matrix<T, N>
Similar to `Matrix<T, N>`, four wrappers for the corresponding helper classes of `std::valarray<T>` are included in this library, namely `SliceMatrix<T>`, `GsliceMatrix<T,N>`, `IndirectMatrix<T,N>` and `MaskMatrix<T>` (see table below for a short description).

  | STL Class                | Wrapper Class              | Description of the Wrapper Class                                              |
  | :----------------------  | :------------------------  | :--------------------------------------------------------------------------   |
  | `std::slice_array<T>`    | `SliceMatrix<T>`           | a sub-Matrix described by `std::slice`                                        |
  | `std::gslice_array<T>`   | `GsliceMatrix<T,N>`        | a sub-Matrix described by `std::gslice`                                       |
  | `std::indirect_array<T>` | `IndirectMatrix<T,N>`      | a sub-Matrix described by an `index_array` (i.e., `std::valarray<std::size>`) |
  | `std::mask_array<T>`     | `MaskMatrix<T>`            | a sub-Matrix described by a `bool_array` (i.e., `std::valarray<bool>`)        |

These wrapper classes behave just like the `Matrix` except they refer to a `Matrix`, rather than owning their own elements. They can be viewed as a reference to a sub-`Matrix` which will be discussed in more detail in section **Subscripting and Slicing**. Most of the time, it is safe to ignore the existence of these three helper classes when dealing with sub-Matrix. But for completeness, we provide following typedefs for 1D/2D/3D sub-Matrix:

  | type  T  | `SliceMatrix<T>` <br> `GsliceMatrix<T, 1>` <br> `IndirectMatrix<T, 1>` <br> `MaskMatrix<T>` | `GsliceMatrix<T, 2>` <br> `IndirectMatrix<T, 2>`  | `Matrix<T, 3>` <br> `IndirectMatrix<T, 3>` |
  | :------: | :-------------------------------------------------------------------: | :--------------------------------: | :---------------------------------: |
  | `double` | `slice_vec`  <br> `gslice_vec`  <br> `indirect_vec`  <br> `mask_vec`  | `gslice_mat`  <br> `indirect_mat`  | `gslice_cube` <br> `indirect_cube`  |
  | `double` | `slice_dvec` <br> `gslice_dvec` <br> `indirect_dvec` <br> `mask_dvec` | `gslice_dmat` <br> `indirect_dmat` | `gslice_dcube`<br> `indirect_dcube` |
  | `float`  | `slice_fvec` <br> `gslice_fvec` <br> `indirect_vec`  <br> `mask_fvec` | `gslice_fmat` <br> `indirect_fmat` | `gslice_fcube`<br> `indirect_fcube` |


## Constructions and Assignments

### Constructors and Destructors

| Matrix<T, N>::Matrix /  Matrix<T, N>::~Matrix <br> (with T = double, N = 1/2/3)                                                                     |   |
|:----------------------------------------------------------------------------------------------------------------------------------------------------|---|
| `vec()`                              <br> `mat()`                                      <br> `cube()`                                                |(1)|
| `explicit vec(n_rows)`               <br> `mat(n_rows, n_cols)`                        <br> `cube(n_rows, n_cols, n_slices)`                        |(2)|
| `vec(const elem_type& val, n_rows)`  <br> `mat(const elem_type& val, n_rows, n_cols)`  <br> `cube(const elem_type& val, n_rows, n_cols, n_slides)`  |(3)|
| `vec(const elem_type* vals, n_rows)` <br> `mat(const elem_type* vals, n_rows, n_cols)` <br> `cube(const elem_type* vals, n_rows, n_cols, n_slides)` |(4)|
| `vec(const vec& other)`              <br> `mat(const mat& other)`                      <br> `cube(const cube& other)`                               |(5)|
| `vec(vec&& other) noexcept`          <br> `mat(mat&& other) noexcept`                  <br> `cube(cube&& other) noexcept`                           |(6) C++11|
| `vec(slice_vec)`                     <br> `mat(slice_vec)`                             <br>                                                         |(7)|
| `vec(gslice_vec)`                    <br> `mat(gslice_mat)`                            <br> `cube(gslice_cube)`                                     |(8)|
| `vec(indirect_vec)`                  <br> `mat(indirect_mat)`                          <br> `cube(indirect_cube)`                                   |(9)|
| `vec(mask_vec)`                      <br> `mat(mask_vec)`                              <br>                                                         |(10)|
| `vec(initializer_list)`              <br> `mat(initializer_list, n_rows, n_cols)`      <br> `cube(initializer_list, n_rows, n_cols, n_slices)`      |(11) C++11|
| `vec(valarray)`                      <br> `mat(valarray, n_rows, n_cols)`              <br> `cube(valarray, n_rows, n_cols, n_slides)`              |(12)|
| `vec(mat)`                           <br> `mat(vec)`                                   <br>                                                         |(9)|
| `~vec()`                             <br> `~mat()`                                     <br> `~cube()`                                               |(10)|

The table above provides ways to construct new matrix from various sources:
  1) Default constructor. Constructs an empty `vec/mat/cube`.
  2) Constructs a `vec/mat/cube` with the specified number of elements in each dimension.
  3) Constructs a `vec/mat/cube` with all elements set to `val` and the specified number of elements in each dimension.
  4) Constructs a `vec/mat/cube` with elements set to the contents of array pointed by `vals` and the specified number of elements in each dimension. If this array contains less than total number of elements (i.e., products of the specified number of elements in each dimension), the behavior is undefined.
  5) Copy constructor. Constructs a `vec/mat/cube` from another one using copy semantics.
  6) Move constructor. Constructs a `vec/mat/cube` from another one using move semantics.
  7) Constructs a `vec/mat/cube` from a `slice_vec` (i.e., a `sub-Matrix` described by `std::slice`).
  8) Constructs a `vec/mat/cube` from a `gslice_vec/gslice_mat/gslice_cube` (i.e., a `sub-Matrix` described by `std::gslice`).
  9) Constructs a `vec/mat/cube` from a `indirect_vec/indirect_mat/indirect_cube` (i.e., a `sub-Matrix` described by `index_array`).
  10) Constructs a `vec/mat/cube` from a `mask_vec` (i.e., a `sub-Matrix` described by `bool_array`).
  11) Constructs a `vec/mat/cube` with elements set to the contents of the `initializer_list` and the specified number of elements in each dimension.
  12) Constructs a `vec/mat/cube` with elements set to the contents of the `valarray` and the specified number of elements in each dimension.
  9) Constructs a `vec` from a **n x 1** `mat`; Constructs a **n x 1** `mat` from a `vec`.
  10) Destructs the `vec/mat/cube`. The destructors of the elements (if **T** is a class) are called and the used storage is deallocated.

### Assignments
| Matrix<T, N>::operator= <br> (with T = double, N = 1/2/3)                                                                       |   |
| :-----------------------------------------------------------------------------------------------------------------------------  |---|
| `vec& operator=(const vec& other)`     <br> `mat& operator=(const mat& other)`     <br> `cube& operator=(const cube& other)`    |(1)|
| `vec& operator=(vec&& other)`          <br> `mat& operator=(mat&& other)`          <br> `cube& operator=(cube&& other)`         |(2) C++11 |
| `vec& operator=(const mat& other)`     <br> `mat& operator=(const vec& other)`     <br>                                         |(3)|
| `vec& operator=(mat&& other)`          <br> `mat& operator=(vec&& other)`          <br>                                         |(4) C++11 |
| `vec& operator=(const elem_type& val)` <br> `mat& operator=(const elem_type& val)` <br> `cube& operator=(const elem_type& val)` |(5)|

The table above provides ways to replace the contents of the matrix:
  1) Copy assignment operator.
  2) Move assignment operator.
  3) Copy assignment operator for the assignment between a `vec` and a **n x 1** `mat`.
  4) Move assignment operator for the assignment between a `vec` and a **n x 1** `mat`.
  5) Replaces each value in `*this` with a copy of `val`.

## Subscripting and Slicing
### Matrix Subscripting
| Matrix<T, N>::operator()  (with T = double, N = 1/2/3)                                                                                                    |   |
| :-------------------------------------------------------------------------------------------------------------------------------------------------------  |---|
| vec : `const elem_type& operator()(i) const` <br> mat : `const elem_type& operator()(i, j) const` <br> cube: `const elem_type& operator()(i, j, k) const` |(1)|
| vec : `elem_type& operator()(i)`             <br> mat : `elem_type& operator()(i, j)`             <br> cube: `elem_type& operator()(i, j, k)`             |(2)|

### Matrix Slicing with SliceMatrix
| Matrix<T, N>'s member functions  (with T = double, N = 1/2/3)                                                                                             |   |
|:----------------------------------------------------------------------------------------------------------------------------------------------------------|---|
| mat : `vec row(i) const / vec col(i) const`                                                                                                               |(1)|
| mat : `slice_vec row(i) / slice_vec col(i)`                                                                                                               |(2)|

### Matrix Slicing with GsliceMatrix
| Matrix<T, N>'s member functions  (with T = double, N = 1/2/3)                                                                                                                                                              |   |
|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---|
| cube: `mat slice(i) const`                                                                                                                                                                                                 |(1)|
| cube: `gslice_mat slice(i)`                                                                                                                                                                                                |(2)|
| mat : `mat rows(first_row, last_row) const / mat cols(first_col, last_col) const`                                      <br> cube: `cube slices(first_slice, last_slice) const`                                             |(3)|
| mat : `gslice_mat rows(first_row, last_row) / gslice_mat cols(first_col, last_col)`                                    <br> cube: `gslice_cube slices(first_slice, last_slice)`                                            |(4)|
| vec : `vec subvec(first_row, last_row) const` <br> mat : `mat submat(first_row, first_col, last_row, last_col) const`  <br> cube: `cube subcube(first_row, first_col, first_slice, last_row, last_col, last_slice) const`  |(5)|
| vec : `gslice_vec subvec(first_row, last_row)`<br> mat : `gslice_mat submat(first_row, first_col, last_row, last_col)` <br> cube: `gslice_cube subcube(first_row, first_col, first_slice, last_row, last_col, last_slice)` |(6)|
| vec : `vec operator()(std::slice s1) const`   <br> mat : `mat operator()(std::slice s1, std::slice s2) const`          <br> cube: `cube operator()(std::slice s1, std::slice s2, std::slice s3) const`                     |(7)|
| vec : `gslice_vec operator()(std::slice s1)`  <br> mat : `gslice_mat operator()(std::slice s1, std::slice s2)`         <br> cube: `gslice_cube operator()(std::slice s1, std::slice s2, std::slice s3)`                    |(8)|

### Matrix Slicing with IndirectMatrix
| Matrix<T, N>'s member functions  (with T = double, N = 1/2/3)                                                                                                                                                                                         |   |
| :---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  |---|
| vec : `vec elem(const index_array& ia) const`         <br> mat: `vec elem(const index_array& ia) const`                                  <br> cube: `vec elem(const index_array& ia) const`                                                           |(1)|
| vec : `indirect_vec elem(const index_array& ia)`      <br> mat: `indirect_vec elem(const index_array& ia)`                               <br> cube: `indirect_vec elem(const index_array& ia)`                                                        |(2)|
| vec : `vec operator()(const index_array& ia) const`   <br> mat: `mat operator()(const index_array& ia1, const index_array& ia2) const`   <br> cube: `cube operator()(const index_array& ia1, const index_array& ia2, const index_array& ia3) const`   |(3)|
| vec : `indirect_vec operator()(const index_array& ia)`<br> mat: `indirect_mat operator()(const index_array& ia1, const index_array& ia2)`<br> cube: `indirect_cube operator()(const index_array& ia1, const index_array& ia2, const index_array& ia3)`|(4)|

### Matrix Slicing with MaskMatrix
| Matrix<T, N>'s member functions  (with T = double, N = 1/2/3)                                                                                                                                                                                         |   |
| :-----------------------------------------------------------------------------------------------------------------------------------------------------------------|---|
| vec : `vec operator()(const bool_array& ba) const` <br> mat: `vec operator()(const bool_array& ba) const` <br> cube: `vec operator()(const bool_array& ba) const` |(1)|
| vec : `mask_vec operator()(const bool_array& ba)`  <br> mat: `mask_vec operator()(const bool_array& ba)`  <br> cube: `mask_vec operator()(const bool_array& ba)`  |(2)|