/**
 * The MIT License (MIT)
 *
 * Copyright (c) 2015 Frank Astier
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "../include/quaternion.h"
#include "../include/quaternion_io.h"

#include <set>
#include <unordered_set>
#include <chrono> // definitely needed on some platforms, but maybe not with clang
#include <random>
#include <iomanip>

#include <boost/math/quaternion.hpp>
#include <boost/rational.hpp>
#include <map>

using namespace std;
using namespace boost;
using namespace quaternion;

//----------------------------------------------------------------------------------------------------------------------
// Test utilities
//----------------------------------------------------------------------------------------------------------------------
/**
 * Prints out arrays to a stream - sometimes useful when debugging unit tests
 */
template <typename T, size_t n>
inline ostream& operator<<(ostream& out, const array<T,n>& x) {
  out << "{";
  for (size_t i = 0; i < n; ++i) {
    out << x[i];
    if (i < n-1)
      out << ",";
  }
  return out << "}";
}

/**
 * Prints out a vector of elements of type T to a stream
 */
template <typename T>
inline ostream& operator<<(ostream& out, const std::vector<T>& x) {
  out << "(";
  for (size_t i = 0; i < x.size(); ++i) {
    out << x[i];
    if (i < x.size()-1)
      out << ",";
  }
  return out << ")";
}

/**
 * Compare a boost quaternion to quaternion, within epsilon.
 * Boost is required only for testing, to compare against a known implementation,
 * so this method is not provided alongside the other nearly_equal in quaternion.h.
 */
template <typename T>
inline bool operator==(const quaternion::Quaternion<T>& us, const boost::math::quaternion<T>& them) {
  return us.a() == them.R_component_1()
         && us.b() == them.R_component_2()
         && us.c() == them.R_component_3()
         && us.d() == them.R_component_4();
}

template <typename T, typename T1>
inline bool nearly_equal(const quaternion::Quaternion<T>& us, const boost::math::quaternion<T>& them, T1 eps) {
  return is_nearly_equal(us.a(), them.R_component_1(), eps)
         && is_nearly_equal(us.b(), them.R_component_2(), eps)
         && is_nearly_equal(us.c(), them.R_component_3(), eps)
         && is_nearly_equal(us.d(), them.R_component_4(), eps);
}

/**
 * This used to test e.g. the polar representation of the quaternions.
 */
template <typename T, size_t n>
inline bool nearly_equal(const array<T,n>& x, const array<T,n>& y, T eps) {
  for (size_t i = 0; i < n; ++i)
    if (!is_nearly_equal(x[i], y[i], eps))
      return false;
  return true;
};

/**
 * This used to test e.g. the rotation matrices.
 */
template <typename T, size_t n1, size_t n2>
inline bool
nearly_equal(const array<array<T,n1>,n2>& x, const array<array<T,n1>,n2>& y, T eps) {
  for (size_t i = 0; i < n1; ++i)
    for (size_t j = 0; j < n2; ++j)
      if (!is_nearly_equal(x[i][j], y[i][j], eps))
        return false;
  return true;
};

/**
 * Transpose of a square matrix.
 */
template <typename T, size_t n>
inline array<array<T,n>,n> transpose(const array<array<T,n>,n>& a) {
  array<array<T,n>,n> r;
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < n; ++j)
      r[i][j] = a[j][i];
  return r;
}

template <typename T>
inline T det(const array<array<T,2>,2>& m) {
  return m[0][0] * m[1][1] - m[0][1] * m[1][0];
}

template <typename T, size_t n>
inline T det(const array<array<T,n>,n>& m) {
  T d = 0;
  for (size_t j = 0; j < n; ++j) {
    array<array<T,n-1>,n-1> sub;
    for (size_t k = 1; k < n; ++k)
      for (size_t l = 0; l < n; ++l)
        if (l != j)
          sub[k-1][l - (l > j)] = m[k][l];
    d += ((j % 2) ? -1 : 1) * m[0][j] * det(sub);
  }
  return d;
}

/**
 * Conjugate of a square matrix.
 * Using "conj" because it is the std:: name for the conjugate of a std::complex.
 * TODO: need 1 function for both real and complex matrices.
 */
template <typename T, size_t n>
inline array<array<T,n>,n> conj(const array<array<T,n>,n>& a) {
  array<array<T,n>,n> r;
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < n; ++j)
      r[i][j] = conj(a[i][j]); // std::conj will transform reals into complex here
  return r;
}

/**
 * Initialize a m x n matrix from a list of values of size m * n.
 */
template <size_t m, size_t n = m, typename T =double>
inline array<array<T,n>,m> make_mat(std::initializer_list<T> l) {
  array<array<T,n>,m> mr;
  const double* it = l.begin();
  for (size_t i = 0; i < m; ++i)
    for (size_t j = 0; j < n; ++j)
      mr[i][j] = *it++;
  return mr;
};

/**
 * Square matrix addition.
 */
template <typename T, size_t n>
inline array<array<T,n>,n>
operator+(const array<array<T,n>,n>& a, const array<array<T,n>,n>& b) {
  array<array<T,n>,n> r;
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < n; ++j)
      r[i][j] = a[i][j] + b[i][j];
  return r;
}

/**
 * Square matrix multiplication.
 * TODO: move to mat utils lib
 */
template <typename T, size_t n>
inline array<array<T,n>,n>
operator*(const array<array<T,n>,n>& a, const array<array<T,n>,n>& b) {
  array<array<T,n>,n> r;
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < n; ++j) {
      T val = 0;
      for (size_t k = 0; k < n; ++k)
        val += a[i][k] * b[k][j];
      r[i][j] = val;
  }
  return r;
}

/**
 * A random number generator, for the speed tests.
 */
std::default_random_engine generator;
std::uniform_real_distribution<float> distribution(0.0,1.0);
auto rng = bind(distribution, generator);

/**
 * A function to generate random quaternions.
 */
template <typename T, typename G>
inline Quaternion<T> random_quaternion(G& g) {
  return {g(), g(), g(), g()};
}

/**
 * For convenience and ease of reading the testing code.
 */
typedef std::complex<float> Cf;
typedef std::complex<double> Cd;
typedef std::complex<long double> Cld;

typedef boost::math::quaternion<float> qf;
typedef boost::math::quaternion<double> qd;

void test_nearly_equal() {
  cout << "Testing nearly equal for array" << endl;

  {
    array<array<float,3>,3> A;
    A[0] = {{1, 0, 0}};
    A[1] = {{0, 1, 0}};
    A[2] = {{0, 0, 1}};
    assert(nearly_equal(A, A, 0.0f));
  }

  {
    array<array<float,3>,3> A, B;
    A[0] = {{1, 0, 0}};
    A[1] = {{0, 1, 0}};
    A[2] = {{0, 0, 1}};
    B[0] = {{1, 0, 0}};
    B[1] = {{0, 1, 0}};
    B[2] = {{0, 0, 2}};
    assert(!nearly_equal(A, B, 1e-6f));
  }
}

void test_det() {
  cout << "Testing det" << endl;
  {
    auto mat22 = make_mat<2,2,double>;
    assert(det(mat22({0,0,0,0})) == 0);
    assert(det(mat22({1,0,0,0})) == 0);
    assert(det(mat22({1,1,0,0})) == 0);
    assert(det(mat22({1,0,0,1})) == 1);
    assert(det(mat22({1,0,1,0})) == 0);
    assert(det(mat22({1,0,1,1})) == 1);
    assert(det(mat22({0,1,1,0})) == -1);
    assert(det(mat22({1,2,3,4})) == -2);
  }
  {
    auto mat33 = make_mat<3,3,double>;
    assert(det(mat33({0,0,0, 0,0,0, 0,0,0})) == 0);
    assert(det(mat33({1,0,0, 0,1,0, 0,0,1})) == 1);
    assert(det(mat33({1,0,0, 0,1,2, 0,3,4})) == -2);
    assert(det(mat33({1,2,3, 4,5,6, 7,8,9})) == 0);
    assert(det(mat33({1,8,3, 2,5,7, 4,6,9})) == 59);
  }
}

//----------------------------------------------------------------------------------------------------------------------
// Unit tests
//----------------------------------------------------------------------------------------------------------------------
void test_constructors() {
  cout << "Testing constructors" << endl;
  {
    Qf x;
    assert(x.a() == 0 && x.b() == 0 && x.c() == 0 && x.d() == 0);
    assert(x.real() == 0);
    assert(x.unreal() == Qf_0);
  }

  {
    Qf x(0,0,0,0);
    assert(x.a() == 0 && x.b() == 0 && x.c() == 0 && x.d() == 0);
  }

  {
    Qf x(1);
    assert(x.a() == 1 && x.b() == 0 && x.c() == 0 && x.d() == 0);
  }

  {
    Qf x(1, 2);
    assert(x.a() == 1 && x.b() == 2 && x.c() == 0 && x.d() == 0);
  }

  {
    Qf x(1, 2, 3);
    assert(x.a() == 1 && x.b() == 2 && x.c() == 3 && x.d() == 0);
  }

  {
    Qf x(1, 2, 3, 4);
    assert(x.a() == 1 && x.b() == 2 && x.c() == 3 && x.d() == 4);

  }

  {
    Qf x(1,2);
    assert(x.a() == 1 && x.b() == 2 && x.c() == 0 && x.d() == 0);
  }

  {
    Qf x(1,2,3,4);
    assert(x.a() == 1 && x.b() == 2 && x.c() == 3 && x.d() == 4);
  }

  {
    Qf x(Cf(1,2));
    assert(x.a() == 1 && x.b() == 2 && x.c() == 0 && x.d() == 0);
  }

  {
    Qf x(Cf(1,2), Cf(3,4));
    assert(x.a() == 1 && x.b() == 2 && x.c() == 3 && x.d() == 4);
    assert(x.c1() == Cf(1,2));
    assert(x.c2() == Cf(3,4));
  }

  {
    Qd x(Cf(1,2), Cf(3,4));
    assert(x.a() == 1 && x.b() == 2 && x.c() == 3 && x.d() == 4);
  }

  {
    Qd y(Cf(1,2), Cf(3,4));
    Qd x(y);
    assert(x.a() == 1 && x.b() == 2 && x.c() == 3 && x.d() == 4);
  }

  {
    Qf y(1,2,3,4);
    Qd x(y);
    assert(x.a() == 1 && x.b() == 2 && x.c() == 3 && x.d() == 4);
  }

  {
    Qd y(Cf(1,2), Cf(3,4));
    Qd x = y;
    assert(x.a() == 1 && x.b() == 2 && x.c() == 3 && x.d() == 4);
  }

  {
    Qf y(1,2,3,4);
    Qd x = y;
    assert(x.a() == 1 && x.b() == 2 && x.c() == 3 && x.d() == 4);
    assert(x.c1() == complex<double>(1,2));
    assert(x.c2() == complex<double>(3,4));
  }

  {
    int d[4] = {1,2,3,4};
    Qf x(d);
    assert(x.a() == 1 && x.b() == 2 && x.c() == 3 && x.d() == 4);
  }

  {
    int d[] = {1,2,3,4, 5,6,7,8};
    int* it = d;
    Qf x(it), y(it + 4);
    assert(x.a() == 1 && x.b() == 2 && x.c() == 3 && x.d() == 4);
    assert(y.a() == 5 && y.b() == 6 && y.c() == 7 && y.d() == 8);
  }

  {
    vector<int> d{1,2,3,4};
    Qf x(d.begin());
    assert(x.a() == 1 && x.b() == 2 && x.c() == 3 && x.d() == 4);
  }

  {
    set<int> d{1,2,3,4};
    Qf x(d.begin());
    assert(x.a() == 1 && x.b() == 2 && x.c() == 3 && x.d() == 4);
  }

  {
    Qf x{1,0,2,0};
    assert(x.a() == 1 && x.b() == 0 && x.c() == 2 && x.d() == 0);
  }

  {
    Qd x({1,2,3,4});
    assert(x.a() == 1 && x.b() == 2 && x.c() == 3 && x.d() == 4);
  }

  {
    Qf x(1,2,3,4);
    Qf y;
    y += 3;
    y = x;
    assert(x == y);
  }
}

void test_trigonometric_constructors() {
  cout << "Testing trigonometric constructors" << endl;
  {
    Qf x = spherical(0.0f,0.0f,0.0f,0.0f);
    assert(x.a() == 0 && x.b() == 0 && x.c() == 0 && x.d() == 0);
  }

  {
    Qd x = spherical(10.0, 3.1415, 3.1415/2.0, 3.1415/4.0);
    assert(nearly_equal(x, boost::math::spherical(10.0, 3.1415, 3.1415 / 2, 3.1415 / 4), 1e-6));
  }

  {
    Qd x = semipolar(10.0, 3.1415, 3.1415/2, 3.1415/4);
    assert(nearly_equal(x, boost::math::semipolar(10.0, 3.1415, 3.1415 / 2, 3.1415 / 4), 1e-6));
  }

  {
    Qd x = multipolar(10.0, 3.1415, 3.1415/2, 3.1415/4);
    assert(nearly_equal(x, boost::math::multipolar(10.0, 3.1415, 3.1415 / 2, 3.1415 / 4), 1e-6));
  }

  {
    Qd x = cylindrospherical(-2.0, -3.0, 3.1415/2, 3.1415/4);
    assert(nearly_equal(x, boost::math::cylindrospherical(-2.0, -3.0, 3.1415 / 2, 3.1415 / 4), 1e-6));
  }

  {
    Qd x = cylindrical(-2.0, 3.1415/2, 3.0, 4.0);
    assert(nearly_equal(x, boost::math::cylindrical(-2.0, 3.1415 / 2, 3.0, 4.0), 1e-6));
  }
}

void test_IJK() {
  cout << "Testing IJK" << endl;
  assert(Qf_0 == 0);
  assert(Qf_1 == 1);
  assert(Qf_i * Qf_i == -1);
  assert(Qf_j * Qf_j == -1);
  assert(Qf_k * Qf_k == -1);
  assert(Qf_i * Qf_j * Qf_k == -1);

  assert(Qd_0 == 0);
  assert(Qd_1 == 1);
  assert(Qd_i * Qd_i == -1);
  assert(Qd_j * Qd_j == -1);
  assert(Qd_k * Qd_k == -1);
  assert(Qd_i * Qd_j * Qd_k == -1);

  assert(Qld_0 == 0);
  assert(Qld_1 == 1);
  assert(Qld_i * Qld_i == -1);
  assert(Qld_j * Qld_j == -1);
  assert(Qld_k * Qld_k == -1);
  assert(Qld_i * Qld_j * Qld_k == -1);

  assert(Qf_i/Qf_j == -Qf_k);
  assert(Qf_j/Qf_k == -Qf_i);
  assert(Qf_k/Qf_i == -Qf_j);

  assert(1.0f/Qf_i == -Qf_i);
  assert(1.0f/Qf_j == -Qf_j);
  assert(1.0f/Qf_k == -Qf_k);

  assert(Qf_i*(Qf_k/Qf_i) == -Qf_k);
  assert((Qf_k/Qf_i)*Qf_i == Qf_k);
}

void test_accessors() {
  cout << "Testing accessors" << endl;
  {
    Qd x;
    assert(x.real() == 0);
    assert(x.unreal() == 0);
    assert(x.unreal() == Qd());
    assert(x.is_zero());
    assert(x.is_zero(1e-6));
    assert(!x.is_non_zero());
    assert(!x.is_non_zero(1e-6));
    assert(x.is_finite());
    assert(!x.is_inf());
    assert(!x.is_nan());
    assert(x.is_real());
    assert(x.is_complex());
    assert(!x.is_unreal());
    assert(!x.is_unit());
    assert(-x == x);
  }
  {
    Qd x(1);
    assert(x.real() == 1);
    assert(x.unreal() == 0);
    assert(x.unreal() == Qd());
    assert(!x.is_zero());
    assert(!x.is_zero(1e-6));
    assert(x.is_non_zero());
    assert(x.is_non_zero(1e-6));
    assert(x.is_finite());
    assert(!x.is_inf());
    assert(!x.is_nan());
    assert(x.is_real());
    assert(x.is_complex());
    assert(!x.is_unreal());
    assert(x.is_unit());
    assert(-x == Qd(-1));
  }
  {
    Qd x(0,1);
    assert(x.real() == 0);
    assert(x.unreal() == Qd(0,1));
    assert(!x.is_zero());
    assert(!x.is_zero(1e-6));
    assert(x.is_non_zero());
    assert(x.is_non_zero(1e-6));
    assert(x.is_finite());
    assert(!x.is_inf());
    assert(!x.is_nan());
    assert(!x.is_real());
    assert(x.is_complex());
    assert(x.is_unreal());
    assert(x.is_unit());
    assert(-x == Qd(0,-1));
  }
  {
    Qd x(1e-12);
    assert(x.real() == 1e-12);
    assert(x.unreal() == 0);
    assert(!x.is_zero());
    assert(x.is_zero(1e-6));
    assert(x.is_non_zero());
    assert(!x.is_non_zero(1e-6));
    assert(x.is_finite());
    assert(!x.is_inf());
    assert(!x.is_nan());
    assert(x.is_real());
    assert(x.is_complex());
    assert(!x.is_unreal());
    assert(!x.is_unit());
  }
  {
    Qf x(1, 2, 3, 4);
    assert(x.real() == 1);
    assert(x.unreal() != 0);
    assert(x.unreal() == Qf(0, 2, 3, 4));
    assert(!x.is_zero());
    assert(!x.is_zero(1e-6));
    assert(x.is_non_zero());
    assert(x.is_non_zero(1e-6));
    assert(x.is_finite());
    assert(!x.is_inf());
    assert(!x.is_nan());
    assert(!x.is_real());
    assert(!x.is_complex());
    assert(!x.is_unreal());
    assert(!x.is_unit());
    assert(-x == Qf(-1,-2,-3,-4));
  }
  {
    assert(is_inf(Qd_i/0));
    assert(is_nan(Qd_i/Qd_0)); // TODO: this should be inf, but it slows down operator/
    assert(is_nan(Qd_0/Qd_0));
  }

  {
    Qf x(3.14);
    assert(x.real() == 3.14f);
  }

  {
    Qf x(3.14, 2.71);
    assert(x.c1().real() == 3.14f && x.c1().imag() == 2.71f);
  }

  {
    Qf x(1,2,3,4);
    array<float,4> a = x.to_array();
    array<float,4> r{{1,2,3,4}};
    assert(a == r); // exact equality expected (operator== is in std)
  }
}

void test_conjugate() {
  cout << "Testing conjugate" << endl;

  {
    assert(conj(Qf(0)) == 0);
    assert(conj(Qf(1)) == 1);
    assert(conj(Qf_i) == -Qf_i);
    assert(conj(Qf_j) == -Qf_j);
    assert(conj(Qf_k) == -Qf_k);
  }
  {
    Qf x(1, 2, 3, 4);
    assert(conj(x) == Qf(1, -2, -3, -4));
    assert(conj(conj(x)) == x);
    assert(conj(x) + x == Qf(2, 0, 0, 0));
    assert(x - conj(x) == Qf(0, 4, 6, 8));
  }
}

void test_complex_matrix_2d_representation() {
  cout << "Testing complex matrix 2d representation" << endl;
  using MR = complex_matrix_2d<double>;

  {
    assert(to_complex_matrix_2d(Qd()) == MR());
    assert(from_complex_matrix_2d(MR()) == 0);
  }
  {
    MR r;
    r[0] = {{Cd(1,0), Cd()}};
    r[1] = {{Cd(), Cd(1,0)}};
    assert(to_complex_matrix_2d(Qd_1) == r);
    assert(from_complex_matrix_2d(r) == Qd_1);
  }
  {
    MR r;
    r[0] = {{Cd(0,1), Cd()}};
    r[1] = {{Cd(), Cd(0,-1)}};
    assert(to_complex_matrix_2d(Qd_i) == r);
    assert(from_complex_matrix_2d(r) == Qd_i);
  }
  {
    MR r;
    r[0] = {{Cd(), Cd(1,0)}};
    r[1] = {{Cd(-1,0), Cd()}};
    assert(to_complex_matrix_2d(Qd_j) == r);
    assert(from_complex_matrix_2d(r) == Qd_j);
  }
  {
    MR r;
    r[0] = {{Cd(), Cd(0,1)}};
    r[1] = {{Cd(0,1), Cd()}};
    assert(to_complex_matrix_2d(Qd_k) == r);
    assert(from_complex_matrix_2d(r) == Qd_k);
  }
  {
    MR r;
    r[0] = {{Cd(1, 2), Cd(3, 4)}};
    r[1] = {{Cd(-3, 4), Cd(1, -2)}};
    assert(to_complex_matrix_2d(Qd(1,2,3,4)) == r);
    assert(from_complex_matrix_2d(r) == Qd(1,2,3,4));
  }
  {
    Qd x(1,2,3,4); MR m = to_complex_matrix_2d(x);
    assert(from_complex_matrix_2d(to_complex_matrix_2d(x)) == x);
    assert(to_complex_matrix_2d(from_complex_matrix_2d(m)) == m);
    assert(norm_squared(x) == det(m));
    assert(to_complex_matrix_2d(conj(x)) == transpose(conj(m)));
  }
  {
    Qd x(1,2,3,4), y(5,6,7,8);
    MR mx = to_complex_matrix_2d(x);
    MR my = to_complex_matrix_2d(y);
    assert(to_complex_matrix_2d(x + y) == mx + my);
    assert(to_complex_matrix_2d(x * y) == mx * my);
  }
}


void test_real_matrix_4d_representation() {
  cout << "Testing real matrix 4d representation" << endl;
  using MR = real_matrix_4d<double>;
  auto mat4d = make_mat<4,4,double>;
  {
    assert(to_real_matrix_4d(Qd_0) == mat4d({0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0}));
    assert(to_real_matrix_4d(Qd_1) == mat4d({1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1}));
    assert(to_real_matrix_4d(Qd_i) == mat4d({0,1,0,0, -1,0,0,0, 0,0,0,-1, 0,0,1,0}));
    assert(to_real_matrix_4d(Qd_j) == mat4d({0,0,1,0, 0,0,0,1, -1,0,0,0, 0,-1,0,0}));
    assert(to_real_matrix_4d(Qd_k) == mat4d({0,0,0,1, 0,0,-1,0, 0,1,0,0, -1,0,0,0}));
  }
  {
    Qd x(1,2,3,4); MR m = to_real_matrix_4d(x);
    assert(from_real_matrix_4d(to_real_matrix_4d(x)) == x);
    assert(to_real_matrix_4d(from_real_matrix_4d(m)) == m);
    assert(norm_squared(x) == std::sqrt(det(m)));
    assert(to_real_matrix_4d(conj(x)) == transpose(m)); // std::conj puts us in complex, but we are reals
  }
  {
    Qd x(1,2,3,4), y(5,6,7,8);
    MR mx = to_real_matrix_4d(x);
    MR my = to_real_matrix_4d(y);
    assert(to_real_matrix_4d(x + y) == mx + my);
    assert(to_real_matrix_4d(x * y) == mx * my);
  }
}

void test_to_polar_representation() {
  cout << "Testing polar representation" << endl;
  typedef polar_representation<double> PR;
  {
    PR expected{{0,0,0,0,0}};
    assert(to_polar_representation(Qd()) == expected);

    expected = {{1,1,1,1,1}};
    assert(!nearly_equal(to_polar_representation(Qd(1)), expected, 1e-6));

    expected = {{3,3.141592,0,0,0}};
    assert(nearly_equal(to_polar_representation(Qd(-3)), expected, 1e-6));

    expected = {{1,3.1415926f/2,1,0,0}};
    assert(nearly_equal(to_polar_representation(Qd(0, 1)), expected, 1e-6));

    expected = {{3,3.1415926f/2,-1,0,0}};
    assert(nearly_equal(to_polar_representation(Qd(0, -3)), expected, 1e-6));

    expected = {{1,3.1415926f/2,0,1,0}};
    assert(nearly_equal(to_polar_representation(Qd(0, 0, 1)), expected, 1e-6));

    expected = {{2.5,3.1415926f/2,0,-1,0}};
    assert(nearly_equal(to_polar_representation(Qd(0, 0, -2.5)), expected, 1e-6));

    expected = {{1,3.1415926f/2,0,0,1}};
    assert(nearly_equal(to_polar_representation(Qd(0, 0, 0, 1)), expected, 1e-6));

    expected = {{3.5,3.1415926f/2,0,0,-1}};
    assert(nearly_equal(to_polar_representation(Qd(0, 0, 0, -3.5)), expected, 1e-6));
  }
}

void test_to_rotation_matrix() {
  cout << "Testing rotation matrix" << endl;
  typedef rotation_matrix<double> RM;
  auto rm = make_mat<3,3,double>;

  {
    Qd x(1);
    RM expected = rm({1,0,0, 0,1,0, 0,0,1});
    assert(nearly_equal(to_rotation_matrix(x), expected, 1e-16));
    assert(from_rotation_matrix(expected) == x);
    assert(from_rotation_matrix(to_rotation_matrix(x)) == x);
    assert(to_rotation_matrix(from_rotation_matrix(expected)) == expected);
  }
  {
    Qd x(0,1);
    RM expected = rm({1,0,0, 0,-1,0, 0,0,-1});
    assert(nearly_equal(to_rotation_matrix(x), expected, 1e-16));
    assert(from_rotation_matrix(expected) == x);
    assert(from_rotation_matrix(to_rotation_matrix(x)) == x);
    assert(to_rotation_matrix(from_rotation_matrix(expected)) == expected);
  }
  {
    Qd x(0,0,1);
    RM expected = rm({-1,0,0, 0,1,0, 0,0,-1});
    assert(nearly_equal(to_rotation_matrix(x), expected, 1e-16));
    assert(from_rotation_matrix(expected) == x);
    assert(from_rotation_matrix(to_rotation_matrix(x)) == x);
    assert(to_rotation_matrix(from_rotation_matrix(expected)) == expected);
  }
  {
    Qd x(0,0,0,1);
    RM expected = rm({-1,0,0, 0,-1,0, 0,0,1});
    assert(nearly_equal(to_rotation_matrix(x), expected, 1e-16));
    assert(from_rotation_matrix(expected) == x);
    assert(from_rotation_matrix(to_rotation_matrix(x)) == x);
    assert(to_rotation_matrix(from_rotation_matrix(expected)) == expected);
  }
  {
    // TODO: more tests here
  }
}

void test_euler_angles() {
  cout << "Testing Euler angles" << endl;
  const double pi = 3.14159265358979323846;
  using A = array<double,3>;
  {
    //assert((to_euler(Qd_0) == array<double,3>{{0,0,0}})); only unit quaternions
    assert((to_euler(Qd_i) == A{{0,0,pi}}));
    assert((to_euler(Qd_j) == A{{pi,0,0}}));
    assert((to_euler(Qd_k) == A{{pi,0,pi}}));
  }
  {
    assert(nearly_equal(from_euler(A{{pi,0,0}}), Qd_i, 1e-16));
    assert(nearly_equal(from_euler(A{{0,pi,0}}), Qd_j, 1e-16));
    assert(nearly_equal(from_euler(A{{0,0,pi}}), -Qd_i + Qd_k, 1e-16));
    assert(nearly_equal(from_euler(A{{pi,0,pi}}), Qd_j, 1e-1));
  }
}

void test_norms() {
  cout << "Testing norms" << endl;
  {
    assert(abs(Qld_0) == 0);
    assert(abs(Qld_1) == 1);
    assert(abs(Qld_i) == 1);
    assert(abs(Qld_j) == 1);
    assert(abs(Qld_k) == 1);
  }

  {
    Qf x;
    assert(norm_squared(x) == 0);
    assert(std::abs(abs(x)) < 1e-6);
    assert(norm_l0(x) == 0);
    assert(norm_l1(x) == 0);
    assert(norm_sup(x) == 0);
    assert(norm_lk(x, 0.5) == 0);
    assert(unreal_norm_squared(x) == 0);
    assert(is_real(x));
    assert(is_complex(x));
    assert(!is_unreal(x));
  }

  {
    Qd x{1,2,3,4};
    assert(norm_squared(x) == 1 + 4 + 9 + 16);
    assert(std::abs(abs(x) - std::sqrt(1 + 4 + 9 + 16)) < 1e-6);
    assert(normalize(x).is_unit(1e-6));
    assert(nearly_equal(normalize(x), x / std::sqrt(1 + 4 + 9 + 16), 1e-6));
    assert(unreal_norm_squared(x) == 29);
    assert(norm_l0(x) == 4);
    assert(norm_l1(x) == 10);
    assert(norm_sup(x) == 4);
    assert(is_nearly_equal(norm_lk(x, .5), pow(1 + sqrt(2) + sqrt(3) + sqrt(4), 2), 1e-6));
  }

  {
    Qd x{1,-4,3,2};
    assert(norm_l0(x) == 4);
    assert(norm_l1(x) == 1+2+3+4);
    assert(norm_sup(x) == 4);
    assert(unreal_norm_squared(x) == 29);
    assert(!is_unit(x));
    assert(is_unit(normalize(x), 1e-6));
    assert(!is_unreal(x));
    assert(is_nearly_equal(norm_lk(x, .5), pow(1 + sqrt(2) + sqrt(3) + sqrt(4), 2), 1e-6));
    x -= 1;
    assert(is_unreal(x, 1e-6));
  }
}

void test_equality() {
  cout << "Testing equality" << endl;
  assert(Qf(1) == 1);
  assert(1 == Qf(1));
  assert(Qf(1) != 2);
  assert(2 != Qf(1));
  assert(Qf(1,2) != 1);
  assert(Qf(1,2) == Cf(1,2));
  assert(Cf(1,2) == Qf(1,2));
  assert(Qf(1,2) != Cf(3,4));
  assert(Cf(3,4) != Qf(1,2));
  assert(Qf(1,2,3) == Qf(1,2,3));
  assert(Qf(1,2,3) != Qf(4,5,6));
  assert(Qf(1,2,3,4) == Qf(1, 2, 3, 4));
  assert(Qf(1,2,3,4) != Qf(4,3,2,1));
  assert(nearly_equal(Qf(1, 2, 3, 4), Qf(1, 2, 3, 4), 0));
  assert(nearly_equal(Qf(1, 2, 3, 4), Qf(1, 2, 3, 4), 1e-6));
  assert(nearly_equal(Qf(1, 2, 3, 4), Qf(1, 2, 3, 3.9999999f), 1e-6));
  assert(!nearly_equal(Qf(1, 2, 3, 4), Qf(0, 2, 3, 3.9999999f), 1e-6));
  assert(!nearly_equal(Qf(1, 2, 3, 4), Qf(1), 1e-6));
  assert(!nearly_equal(Qf(1, 2, 3, 4), Qf(1, 2), 1e-6));
  assert(!nearly_equal(Qf(1, 2, 3, 4), Qf(1, 2, 3), 1e-6));
  assert(!nearly_equal(Qf(1, 2, 3, 4), Cf(1), 1e-6));
  assert(!nearly_equal(Qf(1, 2, 3, 4), Cf(1, 2), 1e-6));
  assert(nearly_equal(Qf(1), Cf(1), 1e-6));
  assert(!nearly_equal(Qf(1), Cf(2), 1e-6));
  assert(nearly_equal(Qf(1, 2), Cf(1, 2), 1e-6));
  assert(!nearly_equal(Qf(1), Cf(1, 2), 1e-6));
  assert(!nearly_equal(Qf(3, 4), Cf(1, 2), 1e-6));
}

void test_plus_minus() {
  cout << "Testing +/-" << endl;
  assert(-Qf(1,2,3,4) == Qf(-1,-2,-3,-4));
  assert(- -Qf(1,2,3,4) == Qf(1,2,3,4));
  assert(- - -Qf(1,2,3,4) == Qf(-1,-2,-3,-4));
  assert(+Qf(1,2,3,4) == Qf(1,2,3,4));
  assert(+ +Qf(1,2,3,4) == Qf(1,2,3,4));
  assert(+ + +Qf(1,2,3,4) == Qf(1,2,3,4));
}

void test_unary_w_scalar() {
  cout << "Testing unary operators with scalar" << endl;
  {
    Qd x(1,2,3,4);
    x += (long double) 3;
    assert(x == Qd(4,2,3,4));
  }

  {
    Qd x(1,2,3,4);
    x -= 3;
    assert(x == Qd(-2,2,3,4));
  }

  {
    Qd x(1,2,3,4);
    x *= -3.5;
    assert(x == Qd(-3.5,-3.5*2,-3.5*3,-3.5*4));
  }

  {
    Qd x(1,2,3,4);
    x /= 2;
    assert(x == Qd(1.0/2,2.0/2,3.0/2,4.0/2));
  }

  {
    Qd x(1,2);
    x /= 3;
    assert(x == complex<double>(1.0/3, 2.0/3));
  }
}

void test_unary_w_complex() {
  cout << "Testing unary operators with complex" << endl;
  {
    Qf x(1); Cf y(3.14, 2.718);
    x += y;
    assert(nearly_equal(x, Cf(4.14, 2.718), 1e-6));
  }

  {
    Qf x(1,2), y(3.14, 2.718);
    x += y;
    assert(nearly_equal(x, Cf(4.14, 4.718), 1e-6));
  }

  {
    Qf x(1,2), y(3.14, 2.718);
    x -= y;
    assert(nearly_equal(x, Cf(-2.14f, -0.718f), 1e-6));
  }

  {
    Qf x(1,2,3,4); Cf y(5,6);
    x *= y;
    assert(x == Qf(1,2,3,4) * Qf(5,6));
  }

  {
    Qf x(1,2,3,4); Cf y(5,6);
    x /= y;
    assert(nearly_equal(x, Qf(1, 2, 3, 4) / Qf(5, 6), 1e-6));
  }
}

void test_unary_w_quaternion() {
  cout << "Testing unary operators with quaternion" << endl;
  {
    Qf x(1), y(3.14);
    x += y;
    assert(std::abs(x.real() - 4.14f) < 1e-6);
  }

  {
    Qf x(1,2), y(3.14, 2.718);
    x += y;
    assert(nearly_equal(x, Cf(4.14, 4.718), 1e-6));
  }

  {
    Qf x(1,2,3,4), y(5,6,7,8);
    x += y;
    assert(x == Qf(6,8,10,12));
    Qd z(1,1,1,1);
    x += z;
    assert(x == Qf(7,9,11,13));
  }

  {
    Qf x(1,2,3,4), y(5,6,7,8);
    qf q1(1,2,3,4), q2(5,6,7,8);
    x *= y;
    assert(x == Qf(-60,12,30,24));
    assert(x == q1 * q2);
    Qd z(1,1,1,1);
    x *= z;
    assert(x == Qf(-126,-42,-18,-54));
  }

  {
    Qf x(1,2,3,4), y(5,6,7,8);
    qf q1(1,2,3,4), q2(5,6,7,8);
    x /= y;
    assert(nearly_equal(x, q1 / q2, 1e-6));
    Qd z(1,1,1,1);
    x /= z;
    assert(nearly_equal(x, (q1 / q2) / qf(1, 1, 1, 1), 1e-6));
  }
}

void test_operators() {
  cout << "Testing operators" << endl;
  {
    assert(Qd(1,2,3,4) + 4 == Qd(5,2,3,4));
    assert(Qd(1,2,3,4) + 4 == 4 + Qd(1,2,3,4));
  }
  {
    assert(Qd(1,2,3,4) - 4 == Qd(-3,2,3,4));
    assert(Qd(1,2,3,4) - 4 == -4 + Qd(1,2,3,4));
  }
  {
    assert(Qd(1,2,3,4) * 4 == Qd(4,8,12,16));
    assert(Qd(1,2,3,4) * 4 == 4 * Qd(1,2,3,4));
  }
  {
    assert(Qd(1,2,3,4) / 4 == Qd(1.0/4,2.0/4,3.0/4,1.0));
    assert(4 / Qd(1,2,3,4) == 4 * inverse(Qd(1,2,3,4)));
  }

  {
    assert(Qd(1,2,3,4) + Cd(1,2) == Qd(2,4,3,4));
    assert(Qd(1,2,3,4) + Cd(1,2) == Cd(1,2) + Qd(1,2,3,4));
  }
  {
    assert(Qd(1,2,3,4) - Cd(1,2) == Qd(0,0,3,4));
    assert(Qd(1,2,3,4) - Cd(1,2) == -Cd(1,2) + Qd(1,2,3,4));
  }
  {

    assert(Qd(1,2,3,4) * Cd(1,2) == qd(1,2,3,4) * Cd(1,2));
    assert(Qd(1,2,3,4) * Cd(1,2) == Qd(-3,4,11,-2));
    assert(Qd(1,2,3,4) * Cd(1,2) == Cd(1,2) * Qd(1,2,3,4));
  }
  {
    assert(Qd(1,2,3,4) / Cd(1,2) == qd(1,2,3,4) / Cd(1,2));
    assert(Qd(1,2,3,4) / Cd(1,2) == Qd(1,0,-1,2));
    assert(Cd(1,2) / Qd(1,2,3,4) == Cd(1,2) * inverse(Qd(1,2,3,4)));
  }
  {
    assert(Qf(1,2,3,4) + Qf(4,5,6,7) == Qf(5,7,9,11));
    assert(Qf(1,2,3,4) - Qf(4,5,6,7) == Qf(-3,-3,-3,-3));
  }
  {
    assert(Qf(1,2,3,4) * Qf(4,-5,6,-7) == qf(1,2,3,4) * qf(4,-5,6,-7));
    assert(nearly_equal(Qf(1,2,3,4) / Qf(4,-5,6,-7), qf(1,2,3,4) / qf(4,-5,6,-7), 1e-6));
  }
}

void test_pow2() {
  cout << "Testing pow2" << endl;

  assert(pow2(Qf_0) == 0);
  assert(pow2(Qf_1) == 1);
  assert(pow2(Qf_1) == Qf_1 * Qf_1);
  assert(pow2(Qf_i) == -1);
  assert(pow2(Qf_j) == -1);
  assert(pow2(Qf_k) == -1);
  assert(pow2(Qf(1,2,3,4)) == Qf(1,2,3,4) * Qf(1,2,3,4));
  assert(pow2(Qf(3.14)) == 3.14f * 3.14f);
  assert(pow2(Qf(1,2)) == Cf(1,2) * Cf(1,2));
}

void test_pow3() {
  cout << "Testing pow3" << endl;

  assert(pow3(Qd_0) == 0);
  assert(pow3(Qd_1) == 1);
  assert(pow3(Qd_i) == -Qd_i);
  assert(pow3(Qd_j) == -Qd_j);
  assert(pow3(Qd_k) == -Qd_k);

  assert(pow3(Qld_0) == 0);
  assert(pow3(Qld_1) == 1);
  assert(pow3(Qld_i) == -Qld_i);
  assert(pow3(Qld_j) == -Qld_j);
  assert(pow3(Qld_k) == -Qld_k);
  assert(pow3(Qf(1,2,3,4)) == Qf(1,2,3,4) * Qf(1,2,3,4) * Qf(1,2,3,4));
}

void test_pow() {
  cout << "Testing pow" << endl;

  { // pow should work as is for various combinations of floats/ints.
    // The library should figure out which actual method to use.
    // We don't want to have to write pow(x, (int) 3) or pow(x, (float)3.5),
    // as it is not convenient.
    Quaternion<int> x(1,2,3,4);
    Quaternion<long> y(1,2,3,4);
    Quaternion<float> z(1,2,3,4);
    Quaternion<double> t(1,2,3,4);
    assert(pow(x,3) == pow(y,3));
    assert(pow(x,3) == pow(z,3));
    assert(pow(x,3) == pow(t,3));
    assert(pow(x,3.5) == pow(y,3.5));
    assert(pow(x,3.5) == pow(z,3.5));
    assert(pow(x,3.5) == pow(t,3.5));
  }

  assert(pow(Qf_0,2) == 0);
  assert(pow(Qf_1,2) == 1);
  assert(pow(Qf_i,2) == -1);
  assert(pow(Qf_j,2) == -1);
  assert(pow(Qf_k,2) == -1);

  assert(pow(Qf_1, 0) == 1);
  assert(pow(Qf_i, 0) == 1);
  assert(pow(Qf_j, 0) == 1);
  assert(pow(Qf_k, 0) == 1);

  assert(nearly_equal(pow(Qf(1), -1), Qf(1), 1e-6));
  assert(nearly_equal(pow(Qf(2), -3), Qf(1.0f / 8), 1e-6));
  assert(nearly_equal(pow(Qf(-2), -3), Qf(-1.0f / 8), 1e-6));
  assert(nearly_equal(pow(Qd_i, -2), pow(Cd(0, 1), -2.0f), 1e-15));
  assert(pow(Qd_i, -2) == -1);
  assert(pow(Qd_j, -2) == -1);
  assert(pow(Qd_k, -2) == -1);

  assert(pow(Qf_1, 0.5f) == 1);
  assert(nearly_equal(pow(-Qf_1, 0.5f), Qf_i, 1e-6));
  assert(nearly_equal(pow(Qf_i, 0.5f), sqrt(Cf(0, 1)), 1e-6));
  assert(nearly_equal(pow(-Qf_i, 0.5f), sqrt(Cf(0, -1)), 1e-6));
  assert(nearly_equal(pow(Qf_j, 0.5f), Qf(1.0f / sqrt(2.0f), 0, 1.0f / sqrt(2.0f)), 1e-6));
  assert(nearly_equal(pow(-Qf_j, 0.5f), Qf(1.0f / sqrt(2.0f), 0, -1.0f / sqrt(2.0f)), 1e-6));
  assert(nearly_equal(pow(Qf_k, 0.5f), Qf(1.0f / sqrt(2.0f), 0, 0, 1.0f / sqrt(2.0f)), 1e-6));
  assert(nearly_equal(pow(-Qf_k, 0.5f), Qf(1.0f / sqrt(2.0f), 0, 0, -1.0f / sqrt(2.0f)), 1e-6));
  assert(pow(Qf_1, -0.33f) == 1);
  assert(nearly_equal(pow(-Qf_1, -0.33f), pow(Cf(-1, 0), -0.33f), 1e-6));
  assert(nearly_equal(pow(Qf_i, -0.33f), pow(Cf(0, 1), -0.33f), 1e-6));
  assert(nearly_equal(pow(-Qf_i, -0.33f), pow(Cf(0, -1), -0.33f), 1e-6));

  for (size_t i = 0; i < 1000; ++i) {
    int n = (int) random() % 20;
    Qld x(1+rand()%5,rand()%5,rand()%5,rand()%5);
    Qld y = Qld_1;
    for (int j = 0; j < n; ++j)
      y *= x;
    assert(norm_squared(y - pow(x, n)) < 1e-10);
  }

  // Need Qld for precision, and even with Qld, not better than 1e-5...
  for (size_t i = 0; i < 1000; ++i) {
    long double n = ((long double)(random() % 20))/(1 + (long double)(random() % 15));
    Qld x(1+rand()%5,rand()%5,rand()%5,rand()%5);
    assert(nearly_equal(pow(x, n), exp(n * log(x)), 1e-5));
  }
}

void test_q_pow() {
  cout << "Testing q pow" << endl;

  assert(pow(Qd_0,Qd(2)) == 0);
  assert(pow(Qd_1,Qd(2)) == 1);
  assert(pow(Qd_i,Qd(2)) == -1);
  assert(pow(Qd_j,Qd(2)) == -1);
  assert(pow(Qd_k,Qd(2)) == -1);

  for (size_t i = 0; i < 10; ++i) {
    Qld x(rand()%5,rand()%5,rand()%5,rand()%5);
    Qld y(rand()%5,rand()%5,rand()%5,rand()%5);
    assert(nearly_equal(pow(x, y), exp(y * log(x)), 1e-6));
  }
}

void test_log() {
  cout << "Testing log" << endl;

  { // Make sure it works for reals
    Qld x(2,0,0,0);
    assert(std::abs(log(x).real() - std::log(2)) < 1e-10);
  }

  size_t N = 10;
  for (size_t i = 0; i < N; ++i) {
    Qld x(rand() % 5, rand() % 5, rand() % 5, rand() % 5);
    assert(nearly_equal(x, exp(log(x)), 1e-6)); // but not the other way around!
  }
}

void test_exp() {
  cout << "Testing exp" << endl;

  { // Make sure it works for reals
    Qld x(2,0,0,0);
    assert(std::abs(exp(x).real() - std::exp(2)) < 1e-10);
  }

  {
    size_t N = 1000;
    for (size_t i = 0; i < N; ++i) {
      Qd x(rand() % 5, rand() % 5, rand() % 5, rand() % 5);
      qd bx(x.a(), x.b(), x.c(), x.d());
      Qd y = exp(x);
      qd by = exp(bx);
      assert(nearly_equal(y, by, 1e-6));
    }
  }
}

void test_dot() {
  cout << "Testing dot product" << endl;
  assert(dot(Qf(1, 2, 3, 4), Qf(2, 2, 2, 2)) == 20);
  assert(dot(Qf(-1, 2, -3, 4), Qf(-1, 2, -3, 4)) == norm_squared(Qf(1, 2, 3, 4)));
  // TODO: verify this works in general
  Qf a(1,2,3,4), b(5,6,7,8);
  Qf d1 = dot(a,b);
  Qf d2 = .5*(conj(b) * a + conj(a) * b);
  assert(d1 == d2);
}

void test_cross() {
  cout << "Testing cross product" << endl;
  // TODO: verify
  assert(Qf(0,-2,4,-2) == cross(Qf(1, 2, 3, 4), Qf(2, 2, 2, 2)));
  // TODO: verify this works in general
  Qf a(0,2,3,4), b(0,6,7,8);
  Qf p1 = cross(a,b);
  Qf p2 = .5*(a * b - conj(b) * conj(a));
  assert(p1 == p2);
}

void test_commutator() {
  cout << "Testing commutator" << endl;
  Qf c = commutator(Qf(1,2,3,4), Qf(2,2,2,2));
  assert(c == Qf(0,-4,8,-4));
  assert(c == Qf(1,2,3,4) * Qf(2,2,2,2) - Qf(2,2,2,2) * Qf(1,2,3,4));
  // TODO: verify this works in general
  Qf a(1,2,3,4), b(5,6,7,8);
  Qf c1 = commutator(a,b);
  Qf c2 = 2 * cross(a,b);
  assert(c1 == c2);
}

void test_trigo() {
  cout << "Testing trigonometric functions" << endl;
  // TODO: verify ordinary real and complex trigo against std
  {
    Qd x{1,2,3,4}; qd qx{1,2,3,4};
    assert(nearly_equal(sin(x), sin(qx), 1e-6));
  }
  {
    Qd x{1,2,3,4}; qd qx{1,2,3,4};
    assert(nearly_equal(cos(x), cos(qx), 1e-6));
  }
  {
    Qd x{1,2,3,4};
    assert(nearly_equal(tan(x), sin(x) / cos(x), 1e-6));
  }
  {
    for (double t = 0; t < 6.28; t += .01) {
      assert(nearly_equal(sin(Qd(t)), sin(t), 1e-6));
      assert(nearly_equal(cos(Qd(t)), cos(t), 1e-6));
      assert(nearly_equal(tan(Qd(t)), tan(t), 1e-6));
    }
  }
  {
    for (double t = 0; t < 6.28; t += .01) {
      assert(nearly_equal(sin(Qd(0, t)), sin(Cd(0, t)), 1e-6));
      assert(nearly_equal(cos(Qd(0, t)), cos(Cd(0, t)), 1e-6));
      assert(nearly_equal(tan(Qd(0, t)), tan(Cd(0, t)), 1e-6));
    }
  }
  {
    for (size_t i = 0; i < 1000; ++i) {
      Qd x = random_quaternion<float>(rng);
      qd qx(x.a(), x.b(), x.c(), x.d());
      assert(nearly_equal(sin(x), sin(qx), 1e-6));
    }
  }
  {
    for (size_t i = 0; i < 1000; ++i) {
      Qd x = random_quaternion<float>(rng);
      qd qx(x.a(), x.b(), x.c(), x.d());
      assert(nearly_equal(cos(x), cos(qx), 1e-6));
    }
  }
  {
    for (size_t i = 0; i < 1000; ++i) {
      Qd x = random_quaternion<float>(rng);
      qd qx(x.a(), x.b(), x.c(), x.d());
      assert(nearly_equal(tan(x), tan(qx), 1e-6));
    }
  }
}

void test_hyper_trigo() {
  cout << "Testing hyperbolic trigonometric functions" << endl;
  {
    Qd x{1,2,3,4};
    assert(nearly_equal(sinh(x), (exp(x) - exp(-x)) / 2, 1e-6));
  }
  {
    Qd x{1,2,3,4};
    assert(nearly_equal(cosh(x), (exp(x) + exp(-x)) / 2, 1e-6));
  }
  {
    Qd x{1,2,3,4};
    assert(nearly_equal(tanh(x), sinh(x) / cosh(x), 1e-6));
  }
  {
    for (size_t i = 0; i < 1000; ++i) {
      Qd x = random_quaternion<float>(rng);
      qd qx(x.a(), x.b(), x.c(), x.d());
      assert(nearly_equal(sinh(x), sinh(qx), 1e-6));
    }
  }
  {
    for (size_t i = 0; i < 1000; ++i) {
      Qd x = random_quaternion<float>(rng);
      qd qx(x.a(), x.b(), x.c(), x.d());
      assert(nearly_equal(cosh(x), cosh(qx), 1e-6));
    }
  }
  {
    for (size_t i = 0; i < 1000; ++i) {
      Qd x = random_quaternion<float>(rng);
      qd qx(x.a(), x.b(), x.c(), x.d());
      assert(nearly_equal(tanh(x), tanh(qx), 1e-6));
    }
  }
  { // TODO: relative precision around 0 is not too good
    for (double t = -4; t < 4; t += .1) {
      assert(nearly_equal(sinh(Qd(t)), sinh(t), 1e-1));
      assert(nearly_equal(cosh(Qd(t)), cosh(t), 1e-1));
      assert(nearly_equal(tanh(Qd(t)), tanh(t), 1e-1));
    }
  }
  {
    for (double t = -4; t < 4; t += .1) {
      assert(nearly_equal(sinh(Qd(0, t)), sinh(Cd(0, t)), 1e-6));
      assert(nearly_equal(cosh(Qd(0, t)), cosh(Cd(0, t)), 1e-6));
      assert(nearly_equal(tanh(Qd(0, t)), tanh(Cd(0, t)), 1e-6));
    }
  }
}

void test_axby() {
  cout << "Testing axby" << endl;
  {
    Qf x(1,2,3,4), y(5,6,7,8);
    assert(axby(0,x,0,y) == 0);
  }
  {
    Qf x(1,2,3,4), y(5,6,7,8);
    assert(axby(1,x,0,y) == x);
  }
  {
    Qf x(1,2,3,4), y(5,6,7,8);
    assert(axby(0,x,1,y) == y);
  }
  {
    Qf x(1,2,3,4), y(5,6,7,8);
    assert(axby(-1,x,3,y) == -x + 3 * y);
  }
}

void test_swap() {
  cout << "Testing swap" << endl;
  {
    Qd x(1,2,3,4), y(5,6,7,8);
    std::swap(x,y);
    assert(x == Qd(5,6,7,8));
    assert(y == Qd(1,2,3,4));
  }
}

void test_lexicographic_order() {
  cout << "Testing lexicographic order" << endl;
  quaternion::lexicographic_order<double> lt;
  {
    assert(lt(Qd_0, Qd_1));
    assert(lt(Qd_1, Qd(2)));
    assert(!lt(Qd_1, Qd_i));
    assert(lt(Qd_i, Qd_1));
    assert(lt(Qd_j, Qd_i));
    assert(lt(Qd_j, Qd_1));
    assert(lt(Qd_k, Qd_j));
    assert(lt(Qd_k, Qd_i));
    assert(lt(Qd_k, Qd_1));
    assert(!lt(Qd(1,2,3,4), Qd(1,2,3,4)));
    assert(lt(Qd(1,2,3,4), Qd(5,6,7,8)));
    assert(lt(Qd(1,2,3,4), Qd(1,5,7,8)));
    assert(lt(Qd(1,2,3,4), Qd(1,2,7,8)));
    assert(lt(Qd(1,2,3,4), Qd(1,2,3,8)));
  }
}

void test_hash() {
  cout << "Testing hash" << endl;
  {
    Qf x(1,2,3,4);
    quaternion::hash<Qf::value_type> h;
    assert(h(x) == 4040899601354814717);
  }
}

void test_stl() {
  cout << "Testing STL" << endl;
  {
    vector<Qf> qs{{1, 2, 3, 4},
                  {5, 6, 7, 8},
                  {1, 3, 5, 7},
                  {2, 4, 6, 8}};
    auto v = accumulate(qs.begin(), qs.end(), Qf_1, multiplies<Qf>());
    assert(v == Qf_1 * Qf(1, 2, 3, 4) * Qf(5, 6, 7, 8) * Qf(1, 3, 5, 7) * Qf(2, 4, 6, 8));
  }

  {
    unordered_set<Qd, quaternion::hash<double>> q_set = {{1, 2, 3, 4}, {5, 6, 7, 8}, {1, 2, 3, 4}};
    assert(q_set.size() == 2);
    auto v = accumulate(q_set.begin(), q_set.end(), Qd_0, plus<Qd>());
    assert(v == Qd(1,2,3,4) + Qd(5,6,7,8));
  }
  {
    using Qi = Quaternion<int>;
    map<Qi, int, quaternion::lexicographic_order<int>> m = {{Qi(2), 3}, {Qi(3,4,5), 1}, {Qi(3,4,6), 2}};
    assert(m.size() == 3);
    auto t = accumulate(m.begin(), m.end(), 0, [](int tt, const pair<Qi,int>& q){ return tt + q.second; });
    assert(t == 6);
  }
}

void test_binary_quaternions() {
  cout << "Testing binary quaternions" << endl;
  {
    using Qb = Quaternion<bool>;
    Qb x(0,1,0,1), y(1,0,1,0);
    assert(x * y == Qb(0,0,0,1));
    //cout << pow(x, 5) << endl; this doesn't compile, requires the "negative" of a binary number
    //cout << exp(x) << endl; // this doesn't compile for now
  }
}

void test_integer_quaternions() {
  cout << "Testing integer quaternions" << endl;
  {
    using Qi = Quaternion<int>;
    Qi x(1, 2, 3, 4), y(5, 6, 7, 8);
    assert(x + y == Qi(6, 8, 10, 12));
    //cout << pow(x, 5) << endl;
    //cout << pow(x, 3.5) << endl;
    //cout << exp(x) << endl; // this doesn't compile for now
  }
  {
    Cd x(1,1);
  }
  // TODO: more tests here
}

// TODO: try quaternion with complex anyway

void test_quaternion_matrix() {
  cout << "Testing quaternion matrix: ";
  {
    constexpr size_t n = 8; // increasing n quickly becomes a perf test
    array<array<Qd,n>,n> M;
    for (size_t i = 0; i < n; ++i)
      for (size_t j = 0; j < n; ++j)
        M[i][j] = random_quaternion<double>(rng);
    cout << det(M);
  }
  cout << " OK" << endl;
}

void test_io() {
  stringstream s;
  s << Qf() << Qf(1) << Qf(-1) << Qf(0,1) << Qf(0,-1) << Qf(0,0,1) << Qf(0,0,-1);
  s << Qf(0,0,0,1) << Qf(0,0,0,-1) << Qf(1,2,3,4) << Qf(-1,2,-3,4);
  assert(s.str() == "01-1i-ij-jk-k1+2i+3j+4k-1+2i-3j+4k");
}

void test_io_eps() {
  cout << "Testing io/eps" << endl;
  Qd x(1e-1, 1e-2, 1e-3, 1e-4);
  stringstream s;
  s << set_display_eps(1e-6) << x;
  s << set_display_eps(1e-4) << x;
  s << set_display_eps(1e-3) << x;
  s << set_display_eps(1e-2) << x;
  s << set_display_eps(1e-1) << x;
  s << set_display_eps(1) << x;
  s << set_display_eps(1e-6);
  assert(s.str() == "0.1+0.01i+0.001j+0.0001k0.1+0.01i+0.001j0.1+0.01i0.1");
}

void test_io_style() {
  cout << "Testing io/style" << endl;
  Qd x(1,2,3,4);
  stringstream s;
  s << set_display_eps(1e-6);
  s << set_display_style(q_nice) << x << " ";
  s << set_display_style(q_compact) << x ;
  assert(s.str() == "1+2i+3j+4k {1,2,3,4}");
}

/**
 * Some operations don't compile with boost rational, in particular
 * those that involve e.g. std::fabs (like operator<< does!).
 * This works fine, but I'm disabling it for now so that Quaternion.h does
 * not depend on boost at all.
 */
//void test_boost_rational() {
//  cout << "Testing with boost::rational" << endl;
//  typedef rational<long> Rl;
//  typedef Quaternion<Rl> Qrl;
//  {
//    Qrl x({1,2},{3,4},{5,6},{7,8}), y({2,3},{4,5},{6,7},{8,9});
//    Qrl z = x * y;
//    assert(z.a() == Rl(-554,315));
//    assert(z.b() == Rl(481,540));
//    assert(z.c() == Rl(641,630));
//    assert(z.d() == Rl(253,252));
//  }
//}

void test_precision() {
  cout << "Testing precision" << endl;
  {
    Qd x(0,1,0,0);
    x = Qd_i * (2.0 + exp(log(x)) / x - 2.0);
    assert(nearly_equal(x, Qd_i, 1e-16));
  }

  {
    Qd x(0.0001,0.002,3.5,3.1415926535);
    assert(nearly_equal(exp(log(x)), x, 1e-8));
  }

  {
    Qf x(0.0001,0.002,3.5,3.1415926535);
    assert(abs(exp(log(x)) - x) < 1e-6);
  }
}

/**
 * Between standard C++, boost and vectorclass which uses intrinsics, no difference.
 * The compiler is probably optimizing well enough, or the intrinsics are not used properly.
 * Whoever is in first position (boot or quaternion) in this micro-benchmark, is twice as
 * fast as the second.... But only on the Mac/clang or gcc. That doesn't happen on Linux/gcc.
 */
void test_multiplication_speed() {
  cout << "Testing multiplication speed" << endl;
  size_t N = 100000;

  Qf q1 = random_quaternion<float>(rng), q2 = random_quaternion<float>(rng);

  {
    float certificate = 0.0;
    auto start = std::chrono::system_clock::now();
    for (size_t i = 0; i < N; ++i) {
      Qf r = q1 * (q2 + (float)i);
      certificate += r.a() + r.b() + r.c() + r.d();
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::nanoseconds diff = end - start;
    cout << "Quaternion: " << (diff.count() / N) << "ns" << endl;
    cout << "Certificate=" << certificate << endl;
  }

  { // With Boost
    qf a(q1.a(),q1.b(),q1.c(),q1.d()), b(q2.a(),q2.b(),q2.c(),q2.d());
    float certificate = 0.0;
    auto start = std::chrono::system_clock::now();
    for (size_t i = 0; i < N; ++i) {
      qf r = a * (b + (float)i);
      certificate += r.R_component_1() + r.R_component_2() + r.R_component_3() + r.R_component_4();
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::nanoseconds diff = end - start;
    cout << "Boost: " << (diff.count() / N) << "ns" << endl;
    cout << "Certificate=" << certificate << endl;
  }
}

// TODO: warm up tests correctly, and get some variance
void test_pow_speed() {
  cout << "Testing pow speed" << endl;
  size_t N = 1000000;

  std::default_random_engine generator;
  std::lognormal_distribution<float> distribution(0.0,1.0);
  auto rng = bind(distribution, generator);
  Qf q1 = random_quaternion<float>(rng);


  { // With Boost
    qf a(q1.a(),q1.b(),q1.c(),q1.d());
    float certificate = 0.0;
    auto start = std::chrono::system_clock::now();
    for (size_t i = 0; i < N; ++i) {
      qf r = pow(a, 15);
      certificate += r.R_component_1() + r.R_component_2() + r.R_component_3() + r.R_component_4();
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::nanoseconds diff = end - start;
    cout << "Boost: " << (diff.count() / N) << "ns" << endl;
    cout << "Certificate=" << certificate << endl;
  }

  {
    float certificate = 0.0;
    auto start = std::chrono::system_clock::now();
    for (size_t i = 0; i < N; ++i) {
      Qf r = pow(q1, 15);
      certificate += r.a() + r.b() + r.c() + r.d();
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::nanoseconds diff = end - start;
    cout << "quaternion: " << (diff.count() / N) << "ns" << endl;
    cout << "Certificate=" << certificate << endl;
  }

  // TODO: match results
}

void test_exp_speed() {
  cout << "Testing exp speed" << endl;

  size_t N = 100000;

  {
    float certificate = 0.0;
    auto start = std::chrono::system_clock::now();
    for (size_t i = 0; i < N; ++i) {
      Qld x(rand() % 5, rand() % 5, rand() % 5, rand() % 5);
      Qld r = exp(x);
      certificate += r.a() + r.b() + r.c() + r.d();
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::nanoseconds diff = end - start;
    cout << "quaternion: " << (diff.count() / N) << "ns" << endl;
    cout << "Certificate=" << certificate << endl;
  }

  {
    float certificate = 0.0;
    auto start = std::chrono::system_clock::now();
    for (size_t i = 0; i < N; ++i) {
      qd x(rand() % 5, rand() % 5, rand() % 5, rand() % 5);
      qd r = exp(x);
      certificate += r.R_component_1() + r.R_component_2() + r.R_component_3() + r.R_component_4();
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::nanoseconds diff = end - start;
    cout << "Boost: " << (diff.count() / N) << "ns" << endl;
    cout << "Certificate=" << certificate << endl;
  }
}

void test_axby_speed() {
  cout << "Testing axby speed" << endl;
  size_t N = 100000;

  Qf q1 = random_quaternion<float>(rng), q2 = random_quaternion<float>(rng);

  { // With Boost
    qf a(q1.a(),q1.b(),q1.c(),q1.d()), b(q2.a(),q2.b(),q2.c(),q2.d());
    float certificate = 0.0;
    auto start = std::chrono::system_clock::now();
    for (size_t i = 0; i < N; ++i) {
      qf r = ((float) i) * a + ((float) i+1) * b;
      certificate += r.R_component_1() + r.R_component_2() + r.R_component_3() + r.R_component_4();
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::nanoseconds diff = end - start;
    cout << "Boost: " << (diff.count() / N) << "ns" << endl;
    cout << "Certificate=" << certificate << endl;
  }

  {
    float certificate = 0.0;
    auto start = std::chrono::system_clock::now();
    for (size_t i = 0; i < N; ++i) {
      Qf r = axby(i, q1, i + 1, q2);
      certificate += r.a() + r.b() + r.c() + r.d();
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::nanoseconds diff = end - start;
    cout << "Quaternion: " << (diff.count() / N) << "ns" << endl;
    cout << "Certificate=" << certificate << endl;
  }
}

void test_tan_speed() {
  cout << "Testing tan speed" << endl;
  size_t N = 100000;

  Qf q1 = random_quaternion<float>(rng);

  { // With Boost
    qf a(q1.a(),q1.b(),q1.c(),q1.d());
    float certificate = 0.0;
    auto start = std::chrono::system_clock::now();
    for (size_t i = 0; i < N; ++i) {
      qf r = tan(a);
      certificate += r.R_component_1() + r.R_component_2() + r.R_component_3() + r.R_component_4();
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::nanoseconds diff = end - start;
    cout << "Boost: " << (diff.count() / N) << "ns" << endl;
    cout << "Certificate=" << certificate << endl;
  }

  {
    float certificate = 0.0;
    auto start = std::chrono::system_clock::now();
    for (size_t i = 0; i < N; ++i) {
      Qf r = tan(q1);
      certificate += r.a() + r.b() + r.c() + r.d();
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::nanoseconds diff = end - start;
    cout << "Quaternion: " << (diff.count() / N) << "ns" << endl;
    cout << "Certificate=" << certificate << endl;
  }
}

int main(int argc, char** argv) {

  test_nearly_equal();
  test_det();

  test_constructors();
  test_trigonometric_constructors();
  test_IJK();
  test_accessors();
  test_conjugate();
  test_complex_matrix_2d_representation();
  test_real_matrix_4d_representation();
  test_to_polar_representation();
  test_to_rotation_matrix();
  test_euler_angles();
  test_norms();
  test_equality();
  test_plus_minus();
  test_unary_w_scalar();
  test_unary_w_complex();
  test_unary_w_quaternion();
  test_operators();
  //test_boost_rational(); disabled for now, but keep, works fine
  test_pow2();
  test_pow3();
  test_pow();
  test_q_pow();
  test_log();
  test_exp();
  test_dot();
  test_cross();
  test_commutator();
  test_trigo();
  test_hyper_trigo();
  test_axby();
  test_swap();
  test_lexicographic_order();
  test_hash();
  test_stl();
  test_binary_quaternions(); // really requires specialized implementation
  test_integer_quaternions(); // some operations don't work (transcendentals)
  test_quaternion_matrix();
  test_io();
  test_io_eps();
  test_io_style();

  cout << endl;
  test_precision();

  cout << endl;
  test_exp_speed();
  test_multiplication_speed();
  test_pow_speed();
  test_axby_speed();
  test_tan_speed();

  return 0;
}
