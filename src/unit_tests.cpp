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

#include "quaternion.h"
#include "quaternion_io.h"

#include <set>
#include <random>
#include <chrono>
#include <iomanip>

#include <boost/math/quaternion.hpp>
#include <boost/rational.hpp>

using namespace std;
using namespace boost;
using namespace boost::math;

//----------------------------------------------------------------------------------------------------------------------
// Test utilities
//----------------------------------------------------------------------------------------------------------------------
///**
// * Prints out arrays to a stream
// */
//template <typename T, size_t n>
//inline ostream& operator<<(ostream& out, const array<T,n>& x) {
//  out << "{";
//  for (size_t i = 0; i < n; ++i) {
//    out << x[i];
//    if (i < n-1)
//      out << ",";
//  }
//  return out << "}";
//}
//
///**
// * Prints out a vector of elements of type T to a stream
// */
//template <typename T>
//inline ostream& operator<<(ostream& out, const vector<T>& x) {
//  out << "(";
//  for (size_t i = 0; i < x.size(); ++i) {
//    out << x[i];
//    if (i < x.size()-1)
//      out << ",";
//  }
//  return out << ")";
//}

/**
 * Compare a boost quaternion to quaternion, within epsilon.
 * Boost is required only for testing, to compare against a known implementation,
 * so this method is not provided alongside the other nearly_equal in quaternion.h.
 */
template <typename T, typename T1>
inline bool nearly_equal(const Quaternion<T>& us, const boost::math::quaternion<T>& them, T1 eps) {
  return is_near_equal_relative(us.a(), them.R_component_1(), eps)
         && is_near_equal_relative(us.b(), them.R_component_2(), eps)
         && is_near_equal_relative(us.c(), them.R_component_3(), eps)
         && is_near_equal_relative(us.d(), them.R_component_4(), eps);
}

/**
 * This method useful in unit tests, to compare against boost, which has
 * a (supposedly) well tested quaternion library.
 * NOTE: hard to compare large values, requires floating point relative comparison.
 * TODO: maybe remove this and keep only a nearly_equal for uniformity
 */
template <typename T>
inline bool operator==(const Quaternion<T>& x, const boost::math::quaternion<T>& boost_y) {
  return nearly_equal(x, boost_y, 1e-6);
}

template <typename T, typename T2>
inline bool nearly_equal(T x, T y, T2 eps) {
  return is_near_equal_relative(x, y, eps);
}

/**
 * This used to test e.g. the polar representation of the quaternions.
 * TODO: this could go into utils because of generality
 */
template <typename T, size_t n>
inline bool nearly_equal(const std::array<T,n>& x, const std::array<T,n>& y, T eps) {
  for (size_t i = 0; i < n; ++i)
    if (!is_near_equal_relative(x[i], y[i], eps))
      return false;
  return true;
};

/**
 * This used to test e.g. the rotation matrices.
 */
template <typename T, size_t n1, size_t n2>
inline bool nearly_equal(const std::array<std::array<T,n1>,n2>& x, const std::array<std::array<T,n1>,n2>& y, T eps) {
  for (size_t i = 0; i < n1; ++i)
    for (size_t j = 0; j < n2; ++j)
      if (!is_near_equal_relative(x[i][j], y[i][j], eps))
        return false;
  return true;
};

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
  cout << "Testing nearly equal" << endl;

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

//----------------------------------------------------------------------------------------------------------------------
// Unit tests
//----------------------------------------------------------------------------------------------------------------------
void test_constructors() {
  cout << "Testing constructors" << endl;
  {
    Qf x;
    assert(x.a() == 0 && x.b() == 0 && x.c() == 0 && x.d() == 0);
    assert(x.real() == 0);
    assert(x.unreal() == 0);
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
    Qf x; x.spherical(0,0,0,0);
    assert(x.a() == 0 && x.b() == 0 && x.c() == 0 && x.d() == 0);
  }

  {
    Qd x = Qd::spherical(10, 3.1415, 3.1415/2, 3.1415/4);
    assert(x == boost::math::spherical(10.0, 3.1415, 3.1415/2, 3.1415/4));
  }

  {
    Qd x = Qd::semipolar(10, 3.1415, 3.1415/2, 3.1415/4);
    assert(x == boost::math::semipolar(10.0, 3.1415, 3.1415/2, 3.1415/4));
  }

  {
    Qd x = Qd::multipolar(10, 3.1415, 3.1415/2, 3.1415/4);
    assert(x == boost::math::multipolar(10.0, 3.1415, 3.1415/2, 3.1415/4));
  }

  {
    Qd x = Qd::cylindrospherical(-2.0, -3.0, 3.1415/2, 3.1415/4);
    assert(x == boost::math::cylindrospherical(-2.0, -3.0, 3.1415/2, 3.1415/4));
  }

  {
    Qd x = Qd::cylindrical(-2.0, 3.1415/2, 3.0, 4.0);
    assert(x == boost::math::cylindrical(-2.0, 3.1415/2, 3.0, 4.0));
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
    Qf x(1, 2, 3, 4);
    assert(x.real() == 1);
    assert(x.unreal() != 0); // TODO: check what is happening here exactly
    assert(x.unreal() == Qf(0, 2, 3, 4));
    assert(!x.is_real());
    assert(!x.is_complex());
    assert(!x.is_unreal());
    assert(!x.is_unit());
    assert(-x == Qf(-1,-2,-3,-4));
  }

  {
    Qf x(3.14);
    float a = x;
    assert(a == 3.14f);
  }

  {
    Qf x(3.14, 2.71);
    Cf c = x;
    assert(c.real() == 3.14f && c.imag() == 2.71f);
  }

  {
    Qf x(1,2,3,4);
    array<float,4> a = x; // can cast to array directly
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

void test_to_matrix_representation() {
  cout << "Testing matrix representation" << endl;
  typedef typename Qd::matrix_representation MR;

  {
    assert(Qd().to_matrix_representation() == MR());
  }
  {
    MR r;
    r[0] = {{Cd(1,0), Cd()}};
    r[1] = {{Cd(), Cd(1,0)}};
    assert(Qd(1).to_matrix_representation() == r);
  }
  {
    MR r;
    r[0] = {{Cd(1, 2), Cd(3, 4)}};
    r[1] = {{Cd(-3, 4), Cd(1, -2)}};
    assert(Qd(1,2,3,4).to_matrix_representation() == r);
  }
}

void test_to_polar_representation() {
  cout << "Testing polar representation" << endl;
  typedef typename Qd::polar_representation PR;
  {
    PR expected{{0,0,0,0,0}};
    assert(Qd().to_polar_representation() == expected);

    expected = {{1,1,1,1,1}};
    assert(!nearly_equal(Qd(1).to_polar_representation(), expected, 1e-6));

    expected = {{3,3.141592,0,0,0}};
    assert(nearly_equal(Qd(-3).to_polar_representation(), expected, 1e-6));

    expected = {{1,3.1415926f/2,1,0,0}};
    assert(nearly_equal(Qd(0,1).to_polar_representation(), expected, 1e-6));

    expected = {{3,3.1415926f/2,-1,0,0}};
    assert(nearly_equal(Qd(0,-3).to_polar_representation(), expected, 1e-6));

    expected = {{1,3.1415926f/2,0,1,0}};
    assert(nearly_equal(Qd(0,0,1).to_polar_representation(), expected, 1e-6));

    expected = {{2.5,3.1415926f/2,0,-1,0}};
    assert(nearly_equal(Qd(0,0,-2.5).to_polar_representation(), expected, 1e-6));

    expected = {{1,3.1415926f/2,0,0,1}};
    assert(nearly_equal(Qd(0,0,0,1).to_polar_representation(), expected, 1e-6));

    expected = {{3.5,3.1415926f/2,0,0,-1}};
    assert(nearly_equal(Qd(0,0,0,-3.5).to_polar_representation(), expected, 1e-6));
  }
}

void test_to_rotation_matrix() {
  cout << "Testing rotation matrix" << endl;
  typedef typename Qd::rotation_matrix RM;
  {
    Qd x(1);
    RM expected;
    expected[0] = {{1,0,0}};
    expected[1] = {{0,1,0}};
    expected[2] = {{0,0,1}};
    assert(nearly_equal(x.to_rotation_matrix(), expected, 1e-10));
    RM r = x.to_rotation_matrix();
    Qd y; y.from_rotation_matrix(r);
    assert(y == Qd_1);
  }
  {
    Qd x(0,1);
    RM expected;
    expected[0] = {{1,0,0}};
    expected[1] = {{0,-1,0}};
    expected[2] = {{0,0,-1}};
    assert(nearly_equal(x.to_rotation_matrix(), expected, 1e-10));
    RM r = x.to_rotation_matrix();
    Qd y; y.from_rotation_matrix(r);
    assert(y == Qd_i);
  }
  {
    Qd x(0,0,1);
    RM expected;
    expected[0] = {{-1,0,0}};
    expected[1] = {{0,1,0}};
    expected[2] = {{0,0,-1}};
    assert(nearly_equal(x.to_rotation_matrix(), expected, 1e-10));
    RM r = x.to_rotation_matrix();
    Qd y; y.from_rotation_matrix(r);
    assert(y == Qd_j);
  }
  {
    Qd x(0,0,0,1);
    RM expected;
    expected[0] = {{-1,0,0}};
    expected[1] = {{0,-1,0}};
    expected[2] = {{0,0,1}};
    assert(nearly_equal(x.to_rotation_matrix(), expected, 1e-10));
    RM r = x.to_rotation_matrix();
    Qd y; y.from_rotation_matrix(r);
    assert(y == Qd_k);
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
    assert(nearly_equal(normalize(x), x/std::sqrt(1+4+9+16), 1e-6));
    assert(unreal_norm_squared(x) == 29);
    assert(norm_l0(x) == 4);
    assert(norm_l1(x) == 10);
    assert(norm_sup(x) == 4);
    assert(nearly_equal(norm_lk(x, .5), pow(1 + sqrt(2) + sqrt(3) + sqrt(4), 2), 1e-6));
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
    assert(nearly_equal(norm_lk(x, .5), pow(1 + sqrt(2) + sqrt(3) + sqrt(4),2), 1e-6));
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
  assert(nearly_equal(Qf(1,2,3,4), Qf(1,2,3,4), 0));
  assert(nearly_equal(Qf(1,2,3,4), Qf(1,2,3,4), 1e-6));
  assert(nearly_equal(Qf(1,2,3,4), Qf(1,2,3,3.9999999f), 1e-6));
  assert(!nearly_equal(Qf(1,2,3,4), Qf(0,2,3,3.9999999f), 1e-6));
  assert(!nearly_equal(Qf(1,2,3,4), Qf(1), 1e-6));
  assert(!nearly_equal(Qf(1,2,3,4), Qf(1,2), 1e-6));
  assert(!nearly_equal(Qf(1,2,3,4), Qf(1,2,3), 1e-6));
  assert(!nearly_equal(Qf(1,2,3,4), Cf(1), 1e-6));
  assert(!nearly_equal(Qf(1,2,3,4), Cf(1,2), 1e-6));
  assert(nearly_equal(Qf(1), Cf(1), 1e-6));
  assert(!nearly_equal(Qf(1), Cf(2), 1e-6));
  assert(nearly_equal(Qf(1,2), Cf(1,2), 1e-6));
  assert(!nearly_equal(Qf(1), Cf(1,2), 1e-6));
  assert(!nearly_equal(Qf(3,4), Cf(1,2), 1e-6));
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
    assert(nearly_equal(x, Cf(4.14,4.718), 1e-6));
  }

  {
    Qf x(1,2), y(3.14, 2.718);
    x -= y;
    assert(nearly_equal(x, Cf(-2.14f,-0.718f), 1e-6));
  }

  {
    Qf x(1,2,3,4); Cf y(5,6);
    x *= y;
    assert(x == Qf(1,2,3,4) * Qf(5,6));
  }

  {
    Qf x(1,2,3,4); Cf y(5,6);
    x /= y;
    assert(nearly_equal(x, Qf(1,2,3,4) / Qf(5,6), 1e-6));
  }
}

void test_unary_w_quaternion() {
  cout << "Testing unary operators with quaternion" << endl;
  {
    Qf x(1), y(3.14);
    x += y;
    assert(std::abs(x - 4.14f) < 1e-6);
  }

  {
    Qf x(1,2), y(3.14, 2.718);
    x += y;
    assert(nearly_equal(x, Cf(4.14,4.718), 1e-6));
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
    assert(q1 * q2 == x);
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
    assert(nearly_equal(x, (q1/q2)/qf(1,1,1,1), 1e-6));
  }
}

void test_operators() {
  cout << "Testing operators" << endl;

  {
    assert(Qf(1,2,3,4) * 3 == Qf(3,6,9,12));
  }

  {
    assert(3 / Qf(1,2) == 3.0f / Cf(1,2));
    assert(3 / Qf(1,2,3,4) == 3.0f / qf(1,2,3,4));
  }

  {
    assert(Qf(1,2,3,4) + Qf(4,5,6,7) == Qf(5,7,9,11));
    assert(Qf(1,2,3,4) - Qf(4,5,6,7) == Qf(-3,-3,-3,-3));
  }
  {
    assert(Qf(1,2,3,4) * Qf(4,-5,6,-7) == qf(1,2,3,4) * qf(4,-5,6,-7));
    assert(Qf(1,2,3,4) / Qf(4,-5,6,-7) == qf(1,2,3,4) / qf(4,-5,6,-7));
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

  assert(pow(Qf_0,2) == 0);
  assert(pow(Qf_1,2) == 1);
  assert(pow(Qf_i,2) == -1);
  assert(pow(Qf_j,2) == -1);
  assert(pow(Qf_k,2) == -1);

  assert(pow(Qf_1, 0) == 1);
  assert(pow(Qf_i, 0) == 1);
  assert(pow(Qf_j, 0) == 1);
  assert(pow(Qf_k, 0) == 1);

  assert(nearly_equal(pow(Qf(1),-1),Qf(1),1e-6));
  assert(nearly_equal(pow(Qf(2),-3),Qf(1.0f/8),1e-6));
  assert(nearly_equal(pow(Qf(-2),-3),Qf(-1.0f/8),1e-6));
  assert(nearly_equal(pow(Qd_i,-2),pow(Cd(0,1),-2.0f),1e-15));
  // TODO: verify
  assert(nearly_equal(pow(Qd_j,-2),Qd(-1),1e-10));
  assert(nearly_equal(pow(Qd_k,-2),Qd(-1),1e-10));

  assert(pow(Qf_1, 0.5f) == 1);
  assert(nearly_equal(pow(-Qf_1, 0.5f), Qf_i, 1e-6));
  assert(nearly_equal(pow(Qf_i, 0.5f), sqrt(Cf(0,1)), 1e-6));
  assert(nearly_equal(pow(-Qf_i, 0.5f), sqrt(Cf(0,-1)), 1e-6));
  assert(nearly_equal(pow(Qf_j, 0.5f), Qf(1.0f/sqrt(2.0f), 0, 1.0f/sqrt(2.0f)), 1e-6));
  assert(nearly_equal(pow(-Qf_j, 0.5f), Qf(1.0f/sqrt(2.0f), 0, -1.0f/sqrt(2.0f)), 1e-6));
  assert(nearly_equal(pow(Qf_k, 0.5f), Qf(1.0f/sqrt(2.0f), 0, 0, 1.0f/sqrt(2.0f)), 1e-6));
  assert(nearly_equal(pow(-Qf_k, 0.5f), Qf(1.0f/sqrt(2.0f), 0, 0, -1.0f/sqrt(2.0f)), 1e-6));
  assert(pow(Qf_1, -0.33f) == 1);
  assert(nearly_equal(pow(-Qf_1, -0.33f), pow(Cf(-1,0), -0.33f), 1e-6));
  assert(nearly_equal(pow(Qf_i, -0.33f), pow(Cf(0,1), -0.33f), 1e-6));
  assert(nearly_equal(pow(-Qf_i, -0.33f), pow(Cf(0,-1), -0.33f), 1e-6));

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
    assert(nearly_equal(pow(x,n), exp(n * log(x)), 1e-5));
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
    assert(nearly_equal(pow(x,y), exp(y * log(x)), 1e-6));
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
      Qld x(rand() % 5, rand() % 5, rand() % 5, rand() % 5);
      quaternion<long double> bx(x.a(), x.b(), x.c(), x.d());
      Qld y = exp(x);
      quaternion<long double> by = exp(bx);
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
  cout << "Testing trigonometic functions" << endl;
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
    assert(nearly_equal(tan(x), sin(x)/cos(x), 1e-6));
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
    assert(nearly_equal(sinh(x), (exp(x) - exp(-x))/2, 1e-6));
  }
  {
    Qd x{1,2,3,4};
    assert(nearly_equal(cosh(x), (exp(x) + exp(-x))/2, 1e-6));
  }
  {
    Qd x{1,2,3,4};
    assert(nearly_equal(tanh(x), sinh(x)/cosh(x), 1e-6));
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

void test_dot_sum_product() {
  cout << "Testing dot sum product" << endl;
  {
    Qf q1(1,2,3,4), q2(5,6,7,8), q3(2,4,6,8);
    Qf q[] = {q1,q2,q3};
    float s[] = {-1,2,-3};
    // TODO: just accumulate? Qf(). is clumsy
    assert(Qf().dot_sum_product(s, s + 3, q) == -1 * q1 + 2 * q2 - 3 * q3);
  }
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

void test_stl() {
  cout << "Testing STL" << endl;
  vector<Qf> qs{{1,2,3,4},{5,6,7,8},{1,3,5,7},{2,4,6,8}};
  auto v = accumulate(qs.begin(), qs.end(), Qf_1, multiplies<Qf>());
  assert(v == Qf_1 * Qf(1,2,3,4) * Qf(5,6,7,8) * Qf(1,3,5,7) * Qf(2,4,6,8));
}

/**
 * Some operations don't compile with boost rational, in particular
 * those that involve e.g. std::fabs (like operator<< does!).
 */
void test_boost_rational() {
  cout << "Testing with boost::rational" << endl;
  typedef rational<long> Rl;
  typedef Quaternion<Rl> Qrl;
  {
    Qrl x({1,2},{3,4},{5,6},{7,8}), y({2,3},{4,5},{6,7},{8,9});
    Qrl z = x * y;
    assert(z.a() == Rl(-554,315));
    assert(z.b() == Rl(481,540));
    assert(z.c() == Rl(641,630));
    assert(z.d() == Rl(253,252));
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

void test_pow_speed() {
  cout << "Testing pow speed" << endl;
  size_t N = 100000;

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
      quaternion<long double> x(rand() % 5, rand() % 5, rand() % 5, rand() % 5);
      quaternion<long double> r = exp(x);
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
      Qld x(rand() % 5, rand() % 5, rand() % 5, rand() % 5);
      Qld r = exp(x);
      certificate += r.a() + r.b() + r.c() + r.d();
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::nanoseconds diff = end - start;
    cout << "quaternion: " << (diff.count() / N) << "ns" << endl;
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

int main() {

  test_nearly_equal();

  test_constructors();
  test_trigonometric_constructors();
  test_IJK();
  test_accessors();
  test_conjugate();
  test_to_matrix_representation();
  test_to_polar_representation();
  test_to_rotation_matrix();
  test_norms();
  test_equality();
  test_plus_minus();
  test_unary_w_scalar();
  test_unary_w_complex();
  test_unary_w_quaternion();
  test_operators();
  test_stl();
  test_boost_rational();
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
  test_dot_sum_product();
  test_io();
  test_io_eps();
  test_io_style();

  cout << endl;
  test_exp_speed();
  test_multiplication_speed();
  test_pow_speed();
  test_axby_speed();
  test_tan_speed();

  return 0;
}
