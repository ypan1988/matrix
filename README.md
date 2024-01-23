`matrix.h`: A multidimensional matrix library
==

## 1. Introduction

The **C++** language does not come with its own matrix library. This
is partially true -- **C++98** introduced a special container named
[`std::valarray`](https://en.cppreference.com/w/cpp/numeric/valarray)
for the processing of arrays of numeric values as an attempt to copy
**Fortran**'s model of multidimensional array.

The `std::valarray` by itself only acts like a 1D array, but it can be
used to simulate a N-dimensional matrix with its helper classes (i.e.,
[`std::slice_array`](https://en.cppreference.com/w/cpp/numeric/valarray/slice_array),
[`std::gslice_array`](https://en.cppreference.com/w/cpp/numeric/valarray/gslice_array),
[`std::mask_array`](https://en.cppreference.com/w/cpp/numeric/valarray/mask_array)
and
[`std::indirect_array`](https://en.cppreference.com/w/cpp/numeric/valarray/mask_array))
which have the reference semantics to a subset of the array.

However, to make more than trivial use of `std:valarray`, you have to
learn a lot about these ancillary classes, some of which are pretty
complex and none of which seems very well documented. Hence, we write
this single header matrix library as a wrapper and/or an extension of
`std::valarray`:

+ The code is based on Bjarne Stroustrup's
[matrix](https://www.stroustrup.com/Programming/Matrix/Matrix.h)
implementation and provides standard building blocks for performing
basic vector, matrix and cube operations;

+ Most of the APIs are coming from `std::valarray` itself and
[**Armadillo**](https://arma.sourceforge.net/docs.html) (a popular
**C++** library for linear algebra with syntax similar to **MATLAB**),
but note that it will not be compatible for sure. While other APIs
have their origins in a few different programming languages (e.g.,
**Fortran**);

+ This is not a C++11 only library -- Since **std::valarray** is
included in **STL** since **C++98**(and it is not likely to be removed
in the near future), this matrix library, which can be based on
`std::valarray`, can be used with **pre-C++11** compilers.


## 2. A Matrix Template

### 2.1 Matrix<T, N>

`Matrix<T, N>` is an **N**-dimensional dense matrix (**N = 1, 2 and
3**) of some value type **T**. Internally it consists of a
`std::valarray<T>` for storing its elements in [column-major
ordering](https://en.wikipedia.org/wiki/Row-_and_column-major_order#Column-major_order)
(like **Fortran**) and a `std::size_t` array of size **N** to record
its sizes on each dimension:

+ **Elements** (or the raw `std::valarray<T>`) can be accessed as a
  whole by `elem()/e()` (i.e., convert a `Matrix<T, N>` back to a
  `std::valarray<T>`):

  ```cpp
  // define an empty 2D matrix with T = double
  Matrix<double, 2> x;
  // ...
  // extract a read-only array of elements
  const std::valarray<double> &va2 = x.elem(); // Or x.e()
  // make a copy for the array of elements
  std::valarray<double> va1 = x.elem(); // Or x.e() again
  ```

+ **Dimensions** of `Matrix<T, N>` can be accessed similarly by
  `get_dims()` which returns a `std::valarray<std::size_t>`. While
  total size and size on each dimension can be obtained through the
  member functions in the following table with the return type
  `std::size_t`:

  | `Matrix<T, 1>` | `Matrix<T, 2>` | `Matrix<T, 3>` | Description              |
  | :------------: | :------------: | :------------: |:-------------------------|
  | `m1.n_elem()`  | `m2.n_elem()`  | `m3.n_elem()`  | total number of elements |
  | `m1.n_rows()`  | `m2.n_rows()`  | `m3.n_rows()`  | number of rows           |
  | `m1.n_cols()`  | `m2.n_cols()`  | `m3.n_cols()`  | number of columns        |
  |                |                | `m3.n_slices()`| number of slices         |

  ```cpp
  // define a 2x3 matrix
  Matrix<double, 2> m2(2, 3);
  // get the dimensions
  std::valarray<std::size_t> dims = m2.get_dims();  // {2, 3}
  std::cout << m2.n_elem() << std::endl;            // 6
  std::cout << m2.n_rows() << std::endl;            // 2
  std::cout << m2.n_cols() << std::endl;            // 3
  ```

### 2.2 Sub-Matrices

Similar to `Matrix<T, N>`, four types of Sub-Matrix classes for the
corresponding helper classes of `std::valarray<T>` are included in
this library, namely `SliceMatrix<T>`, `GsliceMatrix<T>`,
`MaskMatrix<T>` and `IndirectMatrix<T>` (see table below for a short
description).

  | Helper Classes           | Wrapper Classes     | Description of the Wrapper Class                                                             |
  |:-------------------------|:--------------------|:---------------------------------------------------------------------------------------------|
  | `std::slice_array<T>`    | `SliceMatrix<T>`    | Reference to a subset of `Matrix<T,N>` specified by the `std::slice` object.                 |
  | `std::gslice_array<T>`   | `GsliceMatrix<T>`   | Reference to a subset of `Matrix<T,N>` specified by the `std::gslice` object.                |
  | `std::mask_array<T>`     | `MaskMatrix<T>`     | Reference to a subset of `Matrix<T,N>` specified by the `std::valarray<bool>` object.        |
  | `std::indirect_array<T>` | `IndirectMatrix<T>` | Reference to a subset of `Matrix<T,N>` specified by the `std::valarray<std::size_t>` object. |

These helper classes provides writing access to sub-part of
`Matrix<T, N>` and they will be discussed in more detail in section
**Subscripting and Slicing**. Most of the time, it is safe to ignore
the existence of these three helper classes when dealing with
sub-Matrices.

For convenience, in the following part of this document we assume that:
  + `uword = std::size_t` / `bool_array = std::valarray<bool>` /
    `index_array = std::valarray<std::size_t>`
  + `T = double` by default, `elem_type` refers to `Matrix<T,
    N>::elem_type`(hence it should always be the same as type **T**)
  + We use `vec = Matrix<double, 1>` / `mat = Matrix<double, 2>` /
    `cube = Matrix<double, 3>` for better readability, but it is
    possible to use other types
  + Similarly we denote `slice_view = SliceMatrix<double>` /
    `gslice_view = GsliceMatrix<double>` / `mask_view =
    MaskMatrix<double>` / `indirect_view = MaskMatrix<double>`
  + Any API marked as **C++11** can be used only when the compiler
    supports **C++11**.

## 3.Constructions and Assignments

### 3.1 Constructors and Destructors

| `std::valarray<T>::valarray`<br> `std::valarray<T>::~valarray` | `Matrix<T, N>::Matrix` <br> `Matrix<T, N>::~Matrix` <br> `(T = double, N = 1/2/3)                                                   `     | `  ID  `   |
|:---------------------------------------------------------------|:------------------------------------------------------------------------------------------------------------------------------------------|:-----------|
| `valarray()`                                                   | vec: `vec()`                          <br> mat: `mat()`                              <br> cube: `cube()`                                  | 1          |
| `explicit valarray(count)`                                     | vec: `explicit vec(nr)`               <br> mat: `mat(nr, nc)`                        <br> cube: `cube(nr, nc, ns)`                        | 2          |
| `valarray(const T& val, count)`                                | vec: `vec(const elem_type& val, nr)`  <br> mat: `mat(const elem_type& val, nr, nc)`  <br> cube: `cube(const elem_type& val, nr, nc, ns)`  | 3          |
| `valarray(const T* vals, count)`                               | vec: `vec(const elem_type* vals, nr)` <br> mat: `mat(const elem_type* vals, nr, nc)` <br> cube: `cube(const elem_type* vals, nr, nc, ns)` | 4          |
| `valarray(const valarray& other)`                              | vec: `vec(const vec& other)`          <br> mat: `mat(const mat& other)`              <br> cube: `cube(const cube& other)`                 | 5          |
| `valarray(valarray&& other) noexcept`                          | vec: `vec(vec&& other) noexcept`      <br> mat: `mat(mat&& other) noexcept`          <br> cube: `cube(cube&& other) noexcept`             | 6 (C++11)  |
| `valarray(const slice_array<T>& sa)`                           | vec: `vec(const slice_view& other)`   <br> mat: `mat(const slice_view& other)`       <br>                                                 | 7          |
| `valarray(const gslice_array<T>& gsa)`                         | vec: `vec(const gslice_view& other)`  <br> mat: `mat(const gslice_view& other)`      <br> cube: `cube(const gslice_view& other)`          | 8          |
| `valarray(const mask_array<T>& ma)`                            | vec: `vec(const mask_view& other)`    <br> mat: `mat(const mask_view& other)`        <br>                                                 | 9          |
| `valarray(const indirect_array<T>& ia)`                        | vec: `vec(const indirect_view&other)` <br> mat: `mat(const indirect_view& other)`    <br> cube: `cube(const indirect_view& other)`        | 10         |
| `valarray(initializer_list<T> il)`                             | vec: `vec(initializer_list il)`       <br> mat: `mat(nested_initializer_list)`       <br> cube: `cube(nested_initializer_list)`           | 11 (C++11) |
|                                                                | vec: `vec(const valarray& other)`     <br> mat: `mat(const valarray& other, nr, nc)` <br> cube: `cube(const valarray& other, nr, nc, ns)` | 12         |
|                                                                | vec: `vec(const mat& other)`          <br> mat: `mat(const vec& other)`              <br>                                                 | 13         |
|                                                                | vec: `vec(mat&& other)`               <br> mat: `mat(vec&& other)`                   <br>                                                 | 14 (C++11) |
| `~valarray()`                                                  | vec: `~vec()`                         <br> mat: `~mat()`                             <br> cube: `~cube()`                                 | 15         |

A `Matrix<T, N>` (**N = 1,2,3**) can be created from various sources:

  1. Default constructor. Constructs an empty `vec/mat/cube`.
  2. Constructs a `vec/mat/cube` with specified number of elements in
     each dimension.
  3. Constructs a `vec/mat/cube` with all elements set to `val` and
     the specified number of elements in each dimension.
  4. Constructs a `vec/mat/cube` with elements set to the contents of
     array pointed by `vals` and the specified number of elements in
     each dimension. If this array contains less than total number of
     elements (i.e., products of the specified number of elements in
     each dimension), the behavior is undefined.
  5. Copy constructor. Constructs a `vec/mat/cube` from another one
     using copy semantics.
  6. Move constructor. Constructs a `vec/mat/cube` from another one
     using move semantics.
  7. Constructs a `vec/mat/cube` from a `slice_view` (i.e., a
     `sub-Matrix` described by `std::slice`).
  8. Constructs a `vec/mat/cube` from a `gslice_view` (i.e., a
     `sub-Matrix` described by `std::gslice`).
  9. Constructs a `vec/mat/cube` from a `mask_view` (i.e., a
     `sub-Matrix` described by `bool_array`).
  10. Constructs a `vec/mat/cube` from a `indirect_view` (i.e., a
      `sub-Matrix` described by `index_array`).
  11. Constructs a `vec/mat/cube` with elements set to the contents of
      the `initializer_list` and the specified number of elements in
      each dimension.
  12. Constructs a `vec/mat/cube` with elements set to the contents of
      the `valarray` and the specified number of elements in each
      dimension.
  13. Constructs a `vec` from a **n x 1** or **1 x n** `mat` using
      copy semantics; Constructs a **n x 1** or **1 x n** `mat` from a
      `vec` using copy semantics.
  14. Constructs a `vec` from a **n x 1** or ***1 x n* `mat` using
      move semantics; Constructs a **n x 1** or **1 x n** `mat` from a
      `vec` using move semantics.
  15. Destructs the `vec/mat/cube`. The destructors of the elements
      (if **T** is a class) are called and the used storage is
      deallocated.

### 3.2 Assignments
| `std::valarray<T>::operator=`                            | `Matrix<T, N>::operator= (T = double, N = 1/2/3)                                          `                                                           | `   ID   ` |
|:---------------------------------------------------------|:------------------------------------------------------------------------------------------------------------------------------------------------------|------------|
| `valarray<T>& operator=(const valarray<T>& other)`       | vec: `vec& operator=(const vec&)`           <br> mat: `mat& operator=(const mat&)`              <br> cube: `cube& operator=(const cube&)`             | 1          |
| `valarray<T>& operator=(valarray<T>&& other) noexcept`   | vec: `vec& operator=(vec&&)`                <br> mat: `mat& operator=(mat&&)`                   <br> cube: `cube& operator=(cube&&)`                  | 2 (C++11)  |
| `valarray<T>& operator=(const T& val)`                   | vec: `vec& operator=(const elem_type& val)` <br> mat: `mat& operator=(const elem_type& val)`    <br> cube: `cube& operator=(const elem_type& val)`    | 3          |
| `valarray<T>& operator=(const slice_array<T>& other)`    | vec: `vec& operator=(slice_view)`           <br> mat: `mat& operator=(slice_view)`              <br>                                                  | 4          |
| `valarray<T>& operator=(const gslice_array<T>& other)`   | vec: `vec& operator=(gslice_view)`          <br> mat: `mat& operator=(gslice_view)`             <br> cube: `cube& operator=(gslice_view)`             | 5          |
| `valarray<T>& operator=(const mask_array<T>& other)`     | vec: `vec& operator=(mask_view)`            <br> mat: `mat& operator=(mask_view)`               <br>                                                  | 6          |
| `valarray<T>& operator=(const indirect_array<T>& other)` | vec: `vec& operator=(indirect_view)`        <br> mat: `mat& operator=(indirect_view)`           <br> cube: `cube& operator=(indirect_view)`           | 7          |
| `valarray<T>& operator=(initializer_list<T> il)`         | vec: `vec& operator=(initializer_list)`     <br> mat: `mat& operator=(nested_initializer_list)` <br> cube: `cube& operator=(nested_initializer_list)` | 8 (C++11)  |
|                                                          | vec: `vec& operator=(const mat&)`           <br> mat: `mat& operator=(const vec&)`              <br>                                                  | 9          |
|                                                          | vec: `vec& operator=(mat&&)`                <br> mat: `mat& operator=(vec&&)`                   <br>                                                  | 10 (C++11) |

The table above provides ways to replace the contents of the `Matrix<T, N>`:
  1. Copy assignment operator.
  2. Move assignment operator.
  3. Replaces each value in `*this` with a copy of `val`.
  4. Assignment with a `slice_view` (i.e., a `sub-Matrix` described
     by `std::slice`).
  5. Assignment with a `gslice_view` (i.e., a `sub-Matrix` described
     by `std::gslice`).
  6. Assignment with a `mask_view` (i.e., a `sub-Matrix` described
     by `bool_array`).
  7. Assignment with a `indirect_view` (i.e., a `sub-Matrix`
      described by `index_array`).
  8. Assignment with nested initializer list.
  9. Copy assignment operator for the assignment between a `vec` and a
     **n x 1** or **1 x n** `mat`.
  10. Move assignment operator for the assignment between a `vec` and
      a **n x 1** or **1 x n** `mat`.

## 4.Subscripting and Slicing
### 4.1 Matrix Subscripting
| `Matrix<T, N> (T = double, N = 1/2/3) Subscription                                       `                                                                | ` ID `|
|:----------------------------------------------------------------------------------------------------------------------------------------------------------|-------|
| vec : `const elem_type& operator()(i) const` <br> mat : `const elem_type& operator()(i, j) const` <br> cube: `const elem_type& operator()(i, j, k) const` | 1     |
| vec : `elem_type& operator()(i)`             <br> mat : `elem_type& operator()(i, j)`             <br> cube: `elem_type& operator()(i, j, k)`             | 2     |
| vec : `const elem_type& operator[](i) const` <br> mat : `const elem_type& operator[](i) const`    <br> cube: `const elem_type& operator[](i) const`       | 3     |
| vec : `elem_type& operator[](i)`             <br> mat : `elem_type& operator[](i)`                <br> cube: `elem_type& operator[](i)`                   | 4     |


### 4.2 Matrix Slicing with SliceMatrix
One subset of a `std::valarray` is a `std::slice`, which selects every
nth element of a `std::valarray` for some integer n. As we shall see,
this in turn makes it possible to select elements from a row/col/diag
of 2D matrix.  A declaration of a `std::slice` has the form
`std::slice s(start, size, stride);` which specifies the indices
`start, start + stride, start + 2*stride, ...` in a `std::valarray`.

| `SliceMatrix<T>-related member function                                                  ` |` ID `|
|:-------------------------------------------------------------------------------------------|------|
| vec : `vec operator()(std::slice s1) const`                                                | 1  |
| vec : `slice_view operator()(std::slice s1)`                                               | 2  |
| vec : `vec subvec(first_row, last_row) const`                                              | 3  |
| vec : `slice_view subvec(first_row, last_row)`                                             | 4  |
| mat : `vec row(i) const / vec col(i) const`                                                | 5  |
| mat : `slice_view row(i) / slice_view col(i)`                                              | 6  |
| mat : `vec diag(int k) const`                                                              | 7  |
| mat : `slice_view diag(int k)`                                                             | 8  |

### 4.3 Matrix Slicing with GsliceMatrix
| `GsliceMatrix<T>-related member function                                                 `                                                        |` ID `|
|:--------------------------------------------------------------------------------------------------------------------------------------------------|------|
| mat : `mat operator()(std::slice s1, std::slice s2) const`     <br> cube: `cube operator()(std::slice s1, std::slice s2, std::slice s3) const`    | 1  |
| mat : `gslice_view operator()(std::slice s1, std::slice s2)`   <br> cube: `gslice_view operator()(std::slice s1, std::slice s2, std::slice s3)`   | 2  |
| mat : `mat submat(fr, fc, lr, lc) const`                       <br> cube: `cube subcube(fr, fc, fs, lr, lc, ls) const`                            | 3  |
| mat : `gslice_view submat(fr, fc, lr, lc)`                     <br> cube: `gslice_view subcube(fr, fc, fs, lr, lc, ls)`                           | 4  |
| mat : `mat rows(fr, lr) const`                                 <br> mat : `mat cols(fc, lc) const`                                                | 5  |
| mat : `gslice_view rows(fr, lr)`                               <br> mat : `gslice_view cols(fc, lc)`                                              | 6  |
| cube: `mat row(i) const`                                       <br> cube: `mat col(i) const`         <br> cube: `mat slice(i) const`              | 7  |
| cube: `gslice_view row(i)`                                     <br> cube: `gslice_view col(i)`       <br> cube: `gslice_view slice(i)`            | 8  |
| cube: `cube rows(fr, lr) const`                                <br> cube: `cube cols(fc, lc) const`  <br> cube: `cube slices(fs, ls) const`       | 9  |
| cube: `gslice_view rows(fr, lr)`                               <br> cube: `gslice_view cols(fc, lc)` <br> cube: `gslice_view slices(fs, ls)`      | 10 |

### 4.4 Matrix Slicing with MaskMatrix
| `MaskMatrix<T>-related member function                                                   `                                                                               |` ID `|
|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------|
| vec : `vec operator()(const bool_array& ba) const`   <br> mat: `vec operator()(const bool_array& ba) const`    <br> cube: `vec operator()(const bool_array& ba) const`   | 1  |
| vec : `mask_view operator()(const bool_array& ba)`   <br> mat: `mask_view operator()(const bool_array& ba)`    <br> cube: `mask_view operator()(const bool_array& ba)`   | 2  |

### 4.5 Matrix Slicing with IndirectMatrix
| `IndirectMatrix<T>-related member function                                               `                                                                                                        |` ID `|
|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------|
| vec : `vec operator()(const index_array& ia) const`            <br> mat : `vec operator()(const index_array& ia) const`            <br> cube: `vec operator()(const index_array& ia) const`       | 1  |
| vec : `indirect_view operator()(const index_array& ia)`        <br> mat : `indirect_view operator()(const index_array& ia)`        <br> cube: `indirect_view operator()(const index_array& ia)`   | 2  |
| mat : `mat operator()(const index_array& ia1, const index_array& ia2) const`      <br> cube: `cube operator()(const index_array& ia1, const index_array& ia2, const index_array& ia3) const`      | 3  |
| mat : `indirect_view operator()(const index_array& ia1, const index_array& ia2)`  <br> cube: `indirect_view operator()(const index_array& ia1, const index_array& ia2, const index_array& ia3)`   | 4  |
