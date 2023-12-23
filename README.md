`matrix.h`: A small multidimensional matrix library
==

## Introduction

**C++** does not come with its own matrix library. This is partially true -- the **Standard Template Library (STL)** has introduced [`std::valarray`](https://en.cppreference.com/w/cpp/numeric/valarray) for fast mathematical computations since **C++98**. The `std::valarray` by itself only acts like a 1D array, but it can be used to simulate a **N**-dimensional matrix quite easily with its helper classes (i.e., [`std::slice_array`](https://en.cppreference.com/w/cpp/numeric/valarray/slice_array), [`std::gslice_array`](https://en.cppreference.com/w/cpp/numeric/valarray/gslice_array), [`std::mask_array`](https://en.cppreference.com/w/cpp/numeric/valarray/mask_array) and [`std::indirect_array`](https://en.cppreference.com/w/cpp/numeric/valarray/mask_array)) which have the reference semantics to a subset of the array.

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

  | `Matrix<T, 1>` | `Matrix<T, 2>` | `Matrix<T, 3>` | Description                   |
  | :------------: | :------------: | :------------: |:------------------------------|
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

### Sub-Matrix<T>
Similar to `Matrix<T, N>`, four wrappers for the corresponding helper classes of `std::valarray<T>` are included in this library, namely `SliceMatrix<T>`, `GsliceMatrix<T>`, `MaskMatrix<T>` and `IndirectMatrix<T>` (see table below for a short description).

  | Wrapper Class       | Description of the Wrapper Class                                              |
  |:--------------------|:------------------------------------------------------------------------------|
  | `SliceMatrix<T>`    | Reference to a subset of `Matrix<T,N>` specified by the `std::slice` object.  |
  | `GsliceMatrix<T>`   | Reference to a subset of `Matrix<T,N>` specified by the `std::gslice` object. |
  | `MaskMatrix<T>`     | Reference to a subset of `Matrix<T,N>` specified by the `bool_array` object.  |
  | `IndirectMatrix<T>` | Reference to a subset of `Matrix<T,N>` specified by the `index_array` object. |

These wrapper classes behave just like the `Matrix` except they refer to a `Matrix`, rather than owning their own elements. They can be viewed as a reference to a sub-`Matrix` which will be discussed in more detail in section **Subscripting and Slicing**. Most of the time, it is safe to ignore the existence of these three helper classes when dealing with sub-Matrix. But for better readability of this doc, we provide following typedefs for 1D/2D/3D sub-Matrix:

  |       type  T       | <div style="width:290px">double</div> |
  |:-------------------:|:-------------------------------------:|
  |  `SliceMatrix<T>`   |            `submat_slice`             |
  |  `GsliceMatrix<T>`  |            `submat_gslice`            |
  |   `MaskMatrix<T>`   |             `submat_mask`             |
  | `IndirectMatrix<T>` |           `submat_indirect`           |



## Constructions and Assignments

### Constructors and Destructors

| Matrix<T, N>::Matrix /  Matrix<T, N>::~Matrix <br> (with T = double, N = 1/2/3)                                                                                     |              |
|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------|
| vec: `vec()`                              <br> mat: `mat()`                                      <br> cube: `cube()`                                                | (1)          |
| vec: `explicit vec(n_rows)`               <br> mat: `mat(n_rows, n_cols)`                        <br> cube: `cube(n_rows, n_cols, n_slices)`                        | (2)          |
| vec: `vec(const elem_type& val, n_rows)`  <br> mat: `mat(const elem_type& val, n_rows, n_cols)`  <br> cube: `cube(const elem_type& val, n_rows, n_cols, n_slides)`  | (3)          |
| vec: `vec(const elem_type* vals, n_rows)` <br> mat: `mat(const elem_type* vals, n_rows, n_cols)` <br> cube: `cube(const elem_type* vals, n_rows, n_cols, n_slides)` | (4)          |
| vec: `vec(const vec&)`                    <br> mat: `mat(const mat&)`                            <br> cube: `cube(const cube&)`                                     | (5)          |
| vec: `vec(vec&&) noexcept`                <br> mat: `mat(mat&&) noexcept`                        <br> cube: `cube(cube&&) noexcept`                                 | (6) C++11    |
| vec: `vec(submat_slice)`                  <br> mat: `mat(submat_slice)`                          <br>                                                               | (7)          |
| vec: `vec(submat_gslice)`                 <br> mat: `mat(submat_gslice)`                         <br> cube: `cube(submat_gslice)`                                   | (8)          |
| vec: `vec(submat_mask)`                   <br> mat: `mat(submat_mask)`                           <br>                                                               | (9)          |
| vec: `vec(submat_indirect)`               <br> mat: `mat(submat_indirect)`                       <br> cube: `cube(submat_indirect)`                                 | (10)         |
| vec: `vec(initializer_list)`              <br> mat: `mat(initializer_list, n_rows, n_cols)`      <br> cube: `cube(initializer_list, n_rows, n_cols, n_slices)`      | (11v1) C++11 |
| vec: `vec(initializer_list)`              <br> mat: `mat(nested initializer_list)`               <br> cube: `cube(nested initializer_list)`                         | (11v2) C++11 |
| vec: `vec(valarray)`                      <br> mat: `mat(valarray, n_rows, n_cols)`              <br> cube: `cube(valarray, n_rows, n_cols, n_slides)`              | (12)         |
| vec: `vec(const mat&)`                    <br> mat: `mat(const vec&)`                            <br>                                                               | (13)         |
| vec: `vec(mat&&)`                         <br> mat: `mat(vec&&)`                                 <br>                                                               | (14) C++11   |
| vec: `~vec()`                             <br> mat: `~mat()`                                     <br> cube: `~cube()`                                               | (15)         |

The table above provides ways to construct new matrix from various sources:
  1) Default constructor. Constructs an empty `vec/mat/cube`.
  2) Constructs a `vec/mat/cube` with the specified number of elements in each dimension.
  3) Constructs a `vec/mat/cube` with all elements set to `val` and the specified number of elements in each dimension.
  4) Constructs a `vec/mat/cube` with elements set to the contents of array pointed by `vals` and the specified number of elements in each dimension. If this array contains less than total number of elements (i.e., products of the specified number of elements in each dimension), the behavior is undefined.
  5) Copy constructor. Constructs a `vec/mat/cube` from another one using copy semantics.
  6) Move constructor. Constructs a `vec/mat/cube` from another one using move semantics.
  7) Constructs a `vec/mat/cube` from a `submat_slice` (i.e., a `sub-Matrix` described by `std::slice`).
  8) Constructs a `vec/mat/cube` from a `submat_gslice` (i.e., a `sub-Matrix` described by `std::gslice`).
  9) Constructs a `vec/mat/cube` from a `submat_mask` (i.e., a `sub-Matrix` described by `bool_array`).
  10) Constructs a `vec/mat/cube` from a `submat_indirect` (i.e., a `sub-Matrix` described by `index_array`).
  11) Constructs a `vec/mat/cube` with elements set to the contents of the `initializer_list` and the specified number of elements in each dimension.
  12) Constructs a `vec/mat/cube` with elements set to the contents of the `valarray` and the specified number of elements in each dimension.
  13) Constructs a `vec` from a **n x 1** `mat` using copy semantics; Constructs a n x 1 or 1 x n `mat` from a `vec` using copy semantics.
  14) Constructs a `vec` from a **n x 1** `mat` using move semantics; Constructs a n x 1 or 1 x n `mat` from a `vec` using move semantics.
  15) Destructs the `vec/mat/cube`. The destructors of the elements (if **T** is a class) are called and the used storage is deallocated.

### Assignments
| Matrix<T, N>::operator= <br> (with T = double, N = 1/2/3)                                                                                             |            |
|:------------------------------------------------------------------------------------------------------------------------------------------------------|------------|
| vec: `vec& operator=(const vec&)`           <br> mat: `mat& operator=(const mat&)`              <br> cube: `cube& operator=(const cube&)`             | (1)        |
| vec: `vec& operator=(vec&&)`                <br> mat: `mat& operator=(mat&&)`                   <br> cube: `cube& operator=(cube&&)`                  | (2) C++11  |
| vec: `vec& operator=(const elem_type& val)` <br> mat: `mat& operator=(const elem_type& val)`    <br> cube: `cube& operator=(const elem_type& val)`    | (3)        |
| vec: `vec& operator=(submat_slice)`         <br> mat: `mat& operator=(submat_slice)`            <br>                                                  | (4)        |
| vec: `vec& operator=(submat_gslice)`        <br> mat: `mat& operator=(submat_gslice)`           <br> cube: `cube& operator=(submat_gslice)`           | (5)        |
| vec: `vec& operator=(submat_mask)`          <br> mat: `mat& operator=(submat_mask)`             <br>                                                  | (6)        |
| vec: `vec& operator=(submat_indirect)`      <br> mat: `mat& operator=(submat_indirect)`         <br> cube: `cube& operator=(submat_indirect)`         | (7)        |
| vec: `vec& operator=(initializer_list)`     <br> mat: `mat& operator=(nested_initializer_list)` <br> cube: `cube& operator=(nested_initializer_list)` | (8) C++11  |
| vec: `vec& operator=(const mat&)`           <br> mat: `mat& operator=(const vec&)`              <br>                                                  | (9)        |
| vec: `vec& operator=(mat&&)`                <br> mat: `mat& operator=(vec&&)`                   <br>                                                  | (10) C++11 |

The table above provides ways to replace the contents of the matrix:
  1) Copy assignment operator.
  2) Move assignment operator.
  3) Replaces each value in `*this` with a copy of `val`.
  4) Assignment with Sub-Matrix (SliceMatrix).
  5) Assignment with Sub-Matrix (GSliceMatrix).
  6) Assignment with Sub-Matrix (MaskMatrix).
  7) Assignment with Sub-Matrix (IndirectMatrix).
  8) Assignment with nested initializer list.
  9) Copy assignment operator for the assignment between a `vec` and a n x 1 or 1 x n `mat`.
  10) Move assignment operator for the assignment between a `vec` and a n x 1 or 1 x n `mat`.

## Subscripting and Slicing
### Matrix Subscripting
| Matrix<T, N>::operator()  (with T = double, N = 1/2/3)                                                                                                                       |      |
|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------|
| vec : `const elem_type& operator()(i) const` <br> mat : `const elem_type& operator()(i, j) const` <br> cube: `const elem_type& operator()(i, j, k) const`                    | (1)  |
| vec : `elem_type& operator()(i)`             <br> mat : `elem_type& operator()(i, j)`             <br> cube: `elem_type& operator()(i, j, k)`                                | (2)  |
| vec : `const elem_type& operator[](i) const` <br> mat : `const elem_type& operator[](i) const`    <br> cube: `const elem_type& operator[](i) const`                          | (3)  |
| vec : `elem_type& operator[](i)`             <br> mat : `elem_type& operator[](i)`                <br> cube: `elem_type& operator[](i)`                                      | (4)  |


### Matrix Slicing with SliceMatrix
One subset of a `std::valarray` is a `std::slice`, which selects every nth element of a `std::valarray` for some integer n. As we shall see, this in turn makes it possible to select elements from a row/col/diag of 2D matrix.
A declaration of a `std::slice` has the form `std::slice s(start, size, stride);` which specifies the indices `start, start + stride, start + 2*stride, ...` in a `std::valarray`.

| `SliceMatrix<T>-related member function                                                                                                                                    ` |      |
|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------|
| vec : `vec operator()(std::slice s1) const`                                                                                                                                  | (1)  |
| vec : `submat_slice operator()(std::slice s1)`                                                                                                                               | (2)  |
| vec : `vec subvec(first_row, last_row) const`                                                                                                                                | (3)  |
| vec : `submat_slice subvec(first_row, last_row)`                                                                                                                             | (4)  |
| mat : `vec row(i) const`     <br> mat: `vec col(i) const`                                                                                                                    | (5)  |
| mat : `submat_slice row(i)`  <br> mat: `submat_slice col(i)`                                                                                                                 | (6)  |
| mat : `vec diag(int k) const`                                                                                                                                                | (7)  |
| mat : `submat_slice diag(int k)`                                                                                                                                             | (8)  |

### Matrix Slicing with GsliceMatrix
| <div style="width:750px">member function</div>                                                                                                                               |      |
|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------|
| mat : `mat operator()(std::slice s1, std::slice s2) const`             <br> cube: `cube operator()(std::slice s1, std::slice s2, std::slice s3) const`                       | (1a) |
| mat : `submat_gslice operator()(std::slice s1, std::slice s2)`         <br> cube: `submat_gslice operator()(std::slice s1, std::slice s2, std::slice s3)`                    | (1b) |
| mat : `mat submat(first_row, first_col, last_row, last_col) const`     <br> cube: `cube subcube(first_row, first_col, first_slice, last_row, last_col, last_slice) const`    | (2a) |
| mat : `submat_gslice submat(first_row, first_col, last_row, last_col)` <br> cube: `submat_gslice subcube(first_row, first_col, first_slice, last_row, last_col, last_slice)` | (2b) |
| mat : `mat rows(first_row, last_row) const / mat cols(first_col, last_col) const`                                                                                            | (3a) |
| mat : `submat_gslice rows(first_row, last_row) / submat_gslice cols(first_col, last_col)`                                                                                    | (3b) |
| cube: `mat row(i) const` <br> cube: `mat col(i) const` <br> cube: `mat slice(i) const`                                                                                       | (4a) |
| cube: `submat_gslice row(i)` <br> cube: `submat_gslice col(i)` <br> cube: `submat_gslice slice(i)`                                                                           | (4b) |
| cube: `cube rows(first_slice, last_slice) const` <br> cube: `cube cols(first_slice, last_slice) const` <br> cube: `cube slices(first_slice, last_slice) const`               | (5a) |
| cube: `submat_cube rows(first_slice, last_slice)` <br> cube: `submat_cube cols(first_slice, last_slice)` <br> cube: `submat_cube slices(first_slice, last_slice)`            | (5a) |

### Matrix Slicing with MaskMatrix
| <div style="width:750px">member function</div>                                                                                                                               |      |
|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------|
| vec : `vec operator()(const bool_array& ba) const`   <br> mat: `vec operator()(const bool_array& ba) const` <br> cube: `vec operator()(const bool_array& ba) const`          | (1)  |
| vec : `submat_mask operator()(const bool_array& ba)` <br> mat: `submat_mask operator()(const bool_array& ba)`  <br> cube: `submat_mask operator()(const bool_array& ba)`     | (2)  |

### Matrix Slicing with IndirectMatrix
| <div style="width:750px">member function</div>                                                                                                                                                                                                                 |   |
|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---|
| vec : `vec elem(const index_array& ia) const`            <br> mat: `vec elem(const index_array& ia) const`                                     <br> cube: `vec elem(const index_array& ia) const`                                                              |(1)|
| vec : `submat_indirect elem(const index_array& ia)`      <br> mat: `submat_indirect elem(const index_array& ia)`                               <br> cube: `submat_indirect elem(const index_array& ia)`                                                        |(2)|
| vec : `vec operator()(const index_array& ia) const`      <br> mat: `mat operator()(const index_array& ia1, const index_array& ia2) const`      <br> cube: `cube operator()(const index_array& ia1, const index_array& ia2, const index_array& ia3) const`      |(3)|
| vec : `submat_indirect operator()(const index_array& ia)`<br> mat: `submat_indirect operator()(const index_array& ia1, const index_array& ia2)`<br> cube: `submat_indirect operator()(const index_array& ia1, const index_array& ia2, const index_array& ia3)` |(4)|
