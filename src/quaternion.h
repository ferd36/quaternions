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

/**
 * A Quaternion class.
 */

#ifndef QuaternionS_Quaternion_H
#define QuaternionS_Quaternion_H

#include <limits>
#include <array>
#include <complex>
#include <iterator>
#include <assert.h>

#include "Quaternion_utils.h"

// TODO: add namespace

/**
 * A Quaternion class.
 * TODO: provide same operations as boost
 * TODO: provide round to zero, floor, ceil
 * TODO: remove copies/constructions in expressions
 * TODO: fast sinus? (or combined operationrs)
 * TODO: IEEE NaN?
 * TODO: study matrix representation and isomorphism
 * TODO: check references to make sure functionality covered
 * TODO: preconditions
 * TODO: sort out when to provide member functions v external functions
 * TODO: integer Quaternions, binary Quaternions
 * TODO: http://www.gamedev.net/page/resources/_/technical/math-and-physics/Quaternion-powers-r1095
 */
template<typename T =double> // assert operations for numeric is_specialized??
// T has to be real or integer for exp, log, can't accept e.g. complex
// if custom type, check requirements - can't have complex
class Quaternion {
public:
  /**
   * The value of each component of a Quaternion.
   */
  typedef T value_type;

  /**
   * 2x2 matrices are related to Quaternions and Pauli matrices.
   */
  typedef std::array<std::array<std::complex<T>, 2>, 2> matrix_2_2;

  /**
   * The polar representation of a Quaternion.
   */
  typedef std::array<T, 5> polar_representation;

  /**
   * The type used for matrix representations of Quaternions.
   * TODO: array or valarray? - maybe make template parameters, and polar_representation too
   */
  typedef std::array<std::array<std::complex<T>, 2>, 2> matrix_representation;

  /**
   * A 3D rotation matrix.
   */
  typedef std::array<std::array<T, 3>, 3> rotation_matrix;

  /**
  * Construct a Quaternion from at most 4 components of type T.
  * Specifying only a != 0 makes the Quaternion a real.
  * Specifying only a != and b != 0 makes the Quaternion an ordinary complex number.
  */
  Quaternion(T a = 0, T b = 0, T c = 0, T d = 0)
      : _a(a), _b(b), _c(c), _d(d) { }

  /**
   * Construct a Quaternion from at most 4 components of type T.
   * Specifying only a != 0 makes the Quaternion a real.
   * Specifying only a != and b != 0 makes the Quaternion an ordinary complex number.
   */
  template<typename T1 = T, IS_NOT_ITERATOR(T1)>
  Quaternion(T1 a = 0, T1 b = 0, T1 c = 0, T1 d = 0)
      : _a(a), _b(b), _c(c), _d(d) { }

  /**
   * Construct a Quaternion from 2 complex<T>.
   * This sets all 4 components of the Quaternion.
   */
  template<typename T1>
  Quaternion(const std::complex<T1>& x, const std::complex<T1>& y = std::complex<T1>(0, 0))
      : _a(x.real()), _b(x.imag()), _c(y.real()), _d(y.imag()) { }

  /**
   * Construct from a pointer to a range of 4 elements ("float[4]").
   * TODO: make sure we can't use the next constructor here
   */
  template<typename T1 = T>
  Quaternion(T1* it)
      : _a(*it), _b(*++it), _c(*++it), _d(*++it) { }

  /**
   * Construct from an iterator to a range of 4 elements.
   */
  template<typename It, IS_ITERATOR(It)>
  Quaternion(It it)
      : _a(*it), _b(*++it), _c(*++it), _d(*++it) { }

  /**
   * Copy constructor.
   */
  template<typename T1>
  Quaternion(const Quaternion<T1>& y)
      : _a(y.a()), _b(y.b()), _c(y.c()), _d(y.d()) { }

  /**
   * Assignment operator.
   */
  template<typename T1>
  Quaternion &operator=(const Quaternion<T1>& other) {
    _a = other.a();
    _b = other.b();
    _c = other.c();
    _d = other.d();
    return *this;
  }

  // TODO: move to factory?
  static Quaternion spherical(T rho, T theta, T phi1, T phi2) {

    T d = std::sin(phi2);
    T cr = std::cos(phi2);
    T c = cr * std::sin(phi1);
    cr *= std::cos(phi1);
    T b = cr * std::sin(theta);
    T a = cr * std::cos(theta);

    return {rho * a, rho * b, rho * c, rho * d};
  }

  static Quaternion semipolar(T rho, T alpha, T theta1, T theta2) {

    T ca = std::cos(alpha);
    T sa = std::sin(alpha);
    T a = ca * std::cos(theta1);
    T b = ca * std::sin(theta1);
    T c = sa * std::cos(theta2);
    T d = sa * std::sin(theta2);

    return {rho * a, rho * b, rho * c, rho * d};
  }

  static Quaternion multipolar(T rho1, T theta1, T rho2, T theta2) {

    T a = rho1 * std::cos(theta1);
    T b = rho1 * std::sin(theta1);
    T c = rho2 * std::cos(theta2);
    T d = rho2 * std::sin(theta2);

    return {a, b, c, d};
  }

  static Quaternion cylindrospherical(T t, T radius, T longitude, T latitude) {

    T cl = std::cos(latitude);
    T b = radius * cl * std::cos(longitude);
    T c = radius * cl * std::sin(longitude);
    T d = radius * std::sin(latitude);

    return {t, b, c, d};
  }

  static Quaternion cylindrical(T r, T angle, T h1, T h2) {

    T a = r * std::cos(angle);
    T b = r * std::sin(angle);

    return {a, b, h1, h2};
  }

  // TODO: copy to valarray, array, vector...

  /**
   * Accessors for all 4 components of the Quaternion.
   */
  T a() const { return _a; }
  T b() const { return _b; }
  T c() const { return _c; }
  T d() const { return _d; }

  /**
   * The complex components of this Quaternion.
   */
  std::complex<T> c1() const { return {_a,_b}; }
  std::complex<T> c2() const { return {_c,_d}; }

  /**
   * Only for real Quaternions, will assert if not real.
   */
  operator T() const {
    assert(_b == 0 && _c == 0 && _d == 0);
    return _a;
  }

  /**
   * Only for complex, will assert if not complex.
   */
  operator std::complex<T>() const {
    assert(_c == 0 && _d == 0);
    return std::complex<T>(_a,_b);
  }

  /**
   * Ordered list form.
   */
  operator std::array<T,4>() const {
    return {{_a, _b, _c, _d}};
  }

  /**
   * The polar representation of a Quaternion.
   * Returns 5 numbers:
   * - the Euclidean norm of the Quaternion,
   * - the polar angle theta,
   * - and each of the components of the "unreal unit direction".
   */
  polar_representation to_polar_representation() const {
    T nu = unreal_norm_squared();
    T n = std::sqrt(nu + _a*_a);
    assert(nu >= 0);
    if (nu > 0) {
      T theta = std::acos(_a / n);
      T ns = sqrt(nu);
      return {{n, theta, _b / ns, _c / ns, _d / ns}};
    }
    const T pi = std::atan2(+0., -0.);
    // theta = 0 or pi, because n = +/- a().
    return {{n, n == _a ? 0 : pi, 0, 0, 0}};
  }

  /**
   * Returns a matrix representation of a Quaternion.
   */
  matrix_representation to_matrix_representation() const {
    std::array<std::complex<T>, 2> r0{{std::complex<T>(a(), b()), std::complex<T>(c(), d())}};
    std::array<std::complex<T>, 2> r1{{std::complex<T>(-c(), d()), std::complex<T>(a(), -b())}};
    return {{r0, r1}};
  }

  /**
   * Returns a 3D rotation matrix.
   */
  rotation_matrix to_rotation_matrix() const {
    // 21 operations?
    T a2 = _a*_a, b2 = _b*_b, c2 = _c*_c, d2 = _d*_d;
    T ab = _a*_b, ac = _a*_c, ad = _a*_d;
    T bc = _b*_c, bd = _b*_d;
    T cd = _c*_d;
    std::array<T, 3> r0{{a2+b2-c2-d2,2*(bc-ad),2*(bd+ac)}};
    std::array<T, 3> r1{{2*(bc+ad),a2-b2+c2-d2,2*(cd-ab)}};
    std::array<T, 3> r2{{2*(bd-ac),2*(cd+ab),a2-b2-c2+d2}};
    return {{r0, r1, r2}};
  }

  // TODO: clumsy: maybe standalone factory?
  void from_rotation_matrix(const rotation_matrix& rm) {
    T t = rm[0][0] + rm[1][1] + rm[2][2];
    if (t > 0) {
      T s = 0.5 / std::sqrt(t + 1);
      _a = 0.25 / s;
      _b = (rm[2][1] - rm[1][2]) * s;
      _c = (rm[0][2] - rm[2][0]) * s;
      _d = (rm[1][0] - rm[0][1]) * s;
    } else {
      if (rm[0][0] > rm[1][1] && rm[0][0] > rm[2][2]) {
        T s = 2.0 * std::sqrt(1.0 + rm[0][0] - rm[1][1] - rm[2][2]);
        _a = (rm[2][1] - rm[1][2]) / s;
        _b = 0.25 * s;
        _c = (rm[0][1] + rm[1][0]) / s;
        _d = (rm[0][2] + rm[2][0]) / s;
      } else if (rm[1][1] > rm[2][2]) {
        T s = 2.0 * std::sqrt(1.0 + rm[1][1] - rm[0][0] - rm[2][2]);
        _a = (rm[0][2] - rm[2][0]) / s;
        _b = (rm[0][1] + rm[1][0]) / s;
        _c = 0.25 * s;
        _d = (rm[1][2] + rm[2][1]) / s;
      } else {
        T s = 2.0 * std::sqrt(1.0 + rm[2][2] - rm[0][0] - rm[1][1]);
        _a = (rm[1][0] - rm[0][1]) / s;
        _b = (rm[0][2] + rm[2][0]) / s;
        _c = (rm[1][2] + rm[2][1]) / s;
        _d = 0.25 * s;
      }
    }
  }

  /**
   * The real and "unreal" parts of the Quaternion.
   */
  T real() const { return _a; }
  Quaternion unreal() const { return {0, _b, _c, _d}; }

  /**
   * The squared of the norm of the Quaternion.
   * (The square is sometimes useful, and it avoids paying for a sqrt).
   */
  T norm_squared() const {
    return _a * _a + _b * _b + _c * _c + _d * _d;
  }

  /**
   * The norm of the Quaternion (the l2 norm).
   */
  T abs() const {
    return std::sqrt(norm_squared());
  }

  /**
   * The L2 norm of the "unreal" components of the Quaternion,
   * comes back often in computations.
   */
  T unreal_norm_squared() const {
    return _b * _b + _c * _c + _d * _d;
  }

  /**
   * Return true if this Quaternion has norm 1, false otherwise.
   */
  template <typename T1 =T>
  bool is_unit(T1 eps =0) const {
    return is_scalar_zero(norm_squared() - T(1), eps);
  }

  /**
   * Return true if this Quaternion is real, false otherwise.
   */
  template <typename T1 =T>
  bool is_real(T1 eps =0) const {
    return is_scalar_zero(_b, eps)
           && is_scalar_zero(_c, eps)
           && is_scalar_zero(_d, eps);
  }

  /**
   * Return true if this Quaternion is complex, false otherwise.
   */
  template <typename T1 =T>
  bool is_complex(T1 eps =0) const {
    return is_scalar_zero(_c, eps) && is_scalar_zero(_d, eps);
  }

  /**
   * Return true if this Quaternion is real, false otherwise.
   */
  template <typename T1 =T>
  bool is_unreal(T1 eps =0) const {
    return is_scalar_zero(_a, eps)
           && !(is_scalar_zero(_b, eps)
                && is_scalar_zero(_c, eps)
                && is_scalar_zero(_d, eps));
  }

  /**
   * Unary minus.
   */
  Quaternion operator-() const {
    return {-_a, -_b, -_c, -_d};
  }

  /**
  * Unary plus.
  */
  Quaternion operator+() const {
    return *this;
  }

  /**
   * Unary +=.
   */
  Quaternion operator+=(T y) {
    _a += y;
    return *this;
  }

  /**
  * Unary +=.
  */
  Quaternion operator-=(T y) {
    _a -= y;
    return *this;
  }

  /**
  * Scaling by a constant.
  */
  Quaternion operator*=(T k) {
    _a = k * _a;
    _b = k * _b;
    _c = k * _c;
    _d = k * _d;
    return *this;
  }

  /**
   * Dividing by a constant.
   */
  Quaternion operator/=(T k) {
    _a /= k;
    _b /= k;
    _c /= k;
    _d /= k;
    return *this;
  }

  /**
   * Unary +=.
   */
  template<typename T1>
  Quaternion operator+=(const std::complex<T1>& y) {
    _a += y.real();
    _b += y.imag();
    return *this;
  }

  /**
  * Unary -=.
  */
  template<typename T1>
  Quaternion operator-=(const std::complex<T1>& y) {
    _a -= y.real();
    _b -= y.imag();
    return *this;
  }

  /**
   * Unary *=.
   */
  template<typename T1>
  Quaternion operator*=(const std::complex<T1>& y) {

    T at = _a * y.real() - _b * y.imag();
    T bt = _a * y.imag() + _b * y.real();
    T ct = _c * y.real() + _d * y.imag();
    T dt = - _c * y.imag() + _d * y.real();

    _a = at;
    _b = bt;
    _c = ct;
    _d = dt;

    return *this;
  }

  /**
   * Unary /=.
   */
  template<typename T1>
  Quaternion operator/=(const std::complex<T1>& y) {

    T n2 = y.real() * y.real() + y.imag() * y.imag();
    T at = _a * y.real() + _b * y.imag();
    T bt = - _a * y.imag() + _b * y.real();
    T ct = _c * y.real() - _d * y.imag();
    T dt = _c * y.imag() + _d * y.real();

    _a = at / n2;
    _b = bt / n2;
    _c = ct / n2;
    _d = dt / n2;

    return *this;
  }

  /**
   * Unary +=.
   */
  template<typename T1>
  Quaternion operator+=(const Quaternion<T1>& y) {
    _a += y.a();
    _b += y.b();
    _c += y.c();
    _d += y.d();
    return *this;
  }

  /**
   * Unary -=.
   */
  template <typename T1>
  Quaternion operator-=(const Quaternion<T1>& y) {
    _a -= y._a;
    _b -= y._b;
    _c -= y._c;
    _d -= y._d;
    return *this;
  }

  /**
   * Unary multiplication.
   * 28 operations
   */
  template <typename T1>
  Quaternion operator*=(const Quaternion<T1>& y) {

    T at = _a * y.a() - _b * y.b() - _c * y.c() - _d * y.d();
    T bt = _a * y.b() + _b * y.a() + _c * y.d() - _d * y.c();
    T ct = _a * y.c() - _b * y.d() + _c * y.a() + _d * y.b();
    T dt = _a * y.d() + _b * y.c() - _c * y.b() + _d * y.a();

    _a = at;
    _b = bt;
    _c = ct;
    _d = dt;

    return *this;
  }

  /**
   * Unary division with other Quaternion.
   */
  template <typename T1>
  Quaternion operator/=(const Quaternion<T1>& y) {

    T n2 = y.norm_squared();
    T at = _a * y.a() + _b * y.b() + _c * y.c() + _d * y.d();
    T bt = - _a * y.b() + _b * y.a() - _c * y.d() + _d * y.c();
    T ct = - _a * y.c() + _b * y.d() + _c * y.a() - _d * y.b();
    T dt = - _a * y.d() - _b * y.c() + _c * y.b() + _d * y.a();

    _a = at / n2;
    _b = bt / n2;
    _c = ct / n2;
    _d = dt / n2;

    return *this;
  }

  /**
   * k1 * this Quaternion + k2 * y
   * Improves performance by reducing number of constructions/copies.
   * TODO: move to Quaternion algorithms? but faster here?
   */
  template<typename K>
  Quaternion axby(K k1, K k2, const Quaternion &y) {
    _a = k1 * _a + k2 * y._a;
    _b = k1 * _b + k2 * y._b;
    _c = k1 * _c + k2 * y._c;
    _d = k1 * _d + k2 * y._d;
    return *this;
  }

  /**
   * TODO: does this have merit versus std::accumulate? faster?
   */
  template <typename S_it, typename Q_it>
  Quaternion dot_sum_product(S_it s_begin, S_it s_end, Q_it q_begin) {
    for (; s_begin != s_end; ++s_begin, ++q_begin) {
      _a += *s_begin * q_begin->a();
      _b += *s_begin * q_begin->b();
      _c += *s_begin * q_begin->c();
      _d += *s_begin * q_begin->d();
    }
    return *this;
  }

private:
  T _a, _b, _c, _d; // the full state for a Quaternion
};

/**
 * Predefined Quaternions on floats.
 */
typedef Quaternion<float> Qf;
const Qf Qf_0;
const Qf Qf_1(1);
const Qf Qf_i(0, 1);
const Qf Qf_j(0, 0, 1);
const Qf Qf_k(0, 0, 0, 1);

/**
 * Predefined Quaternions on doubles.
 */
typedef Quaternion<double> Qd;
const Qd Qd_0;
const Qd Qd_1(1);
const Qd Qd_i(0, 1);
const Qd Qd_j(0, 0, 1);
const Qd Qd_k(0, 0, 0, 1);

/**
 * Predefined Quaternions on long doubles.
 */
typedef Quaternion<long double> Qld;
const Qld Qld_0;
const Qld Qld_1(1);
const Qld Qld_i(0, 1);
const Qld Qld_j(0, 0, 1);
const Qld Qld_k(0, 0, 0, 1);

/**
 * Multiplication by a constant on the left.
 */
template<typename T, typename T1>
inline Quaternion<T> operator*(T1 k, const Quaternion<T>& x) {
  return Quaternion<T>(x) *= k;
}

// Same, swapping the lhs and rhs.
template<typename T, typename T1>
inline Quaternion<T> operator*(const Quaternion<T>& x, T1 k) {
  return k * x;
}

// TODO: multiplication by a complex number

/**
 * Division by a number.
 */
template<typename T, typename T1>
inline Quaternion<T> operator/(const Quaternion<T>& x, T1 k) {
  return Quaternion<T>(x) /= k;
}

/** +
 * Returns the conjugate of x, as a new Quaternion (x is unchanged).
 */
template<typename T>
inline Quaternion<T> conj(const Quaternion<T> &x) {
  return {x.a(), -x.b(), -x.c(), -x.d()};
}

/** +
 * Norms on a Quaternion.
 */
template<typename T>
inline T norm_squared(const Quaternion<T> &x) {
  return x.norm_squared();
}

// abs = l2 norm = euclidean norm
template<typename T>
inline T abs(const Quaternion<T> &x) {
  return x.abs();
}

template<typename T>
inline T unreal_norm_squared(const Quaternion<T> &x) {
  return x.unreal_norm_squared();
}

// Hamming
template <typename T>
inline T norm_l0(const Quaternion<T>& x) {
  return (x.a() != 0) + (x.b() != 0) + (x.c() != 0) + (x.d() != 0);
}

// l1 norm = taxicab = manhattan
template<typename T>
inline T norm_l1(const Quaternion<T>& x) {
  return std::abs(x.a()) + std::abs(x.b()) + std::abs(x.c()) + std::abs(x.d());
}

template<typename T, typename T1>
inline T norm_lk(const Quaternion<T> &x, T1 k) {
  return std::pow(std::pow(std::abs(x.a()),k)
                  + std::pow(std::abs(x.b()),k)
                  + std::pow(std::abs(x.c()),k)
                  + std::pow(std::abs(x.d()),k),1.0/k);
}

// norm sup = max norm = norm inf
template<typename T>
inline T norm_sup(const Quaternion<T>& x) {
  return std::max(std::max(std::abs(x.a()), std::abs(x.b())),
                  std::max(std::abs(x.c()), std::abs(x.d())));
}

/**
 * Quaternion tests.
 */
template<typename T, typename T1 =T>
inline bool is_unit(const Quaternion<T>& x, T1 eps =0) {
  return x.is_unit(eps);
}

template<typename T, typename T1 =T>
inline bool is_real(const Quaternion<T>& x, T1 eps =0) {
  return x.is_real(eps);
}

template<typename T, typename T1 =T>
inline bool is_complex(const Quaternion<T>& x, T1 eps =0) {
  return x.is_complex(eps);
}

template<typename T, typename T1 =T>
inline bool is_unreal(const Quaternion<T>& x, T1 eps =0) {
  return x.is_unreal(eps);
}

/**
 * Equality.
 */
template<typename T>
inline bool operator==(const Quaternion<T>& x, const Quaternion<T>& y) {
  return x.a() == y.a() && x.b() == y.b() && x.c() == y.c() && x.d() == y.d();
}

template<typename T>
inline bool operator!=(const Quaternion<T>& x, const Quaternion<T>& y) {
  return !(x == y);
}

/**
 * This is more costly, but very useful in practice: Quaternions transcendentals
 * require a lot of floating point operations, so accuracy degrades quickly.
 */
template <typename T, typename T1> // T1 allows eps to be the default type double
inline bool nearly_equal(const Quaternion<T>& x, const Quaternion<T>& y, T1 eps) {
  return is_near_equal_relative(x.a(), y.a(), eps)
         && is_near_equal_relative(x.b(), y.b(), eps)
         && is_near_equal_relative(x.c(), y.c(), eps)
         && is_near_equal_relative(x.d(), y.d(), eps);
}

template<typename T, typename T2>
inline bool operator==(const Quaternion<T>& x, T2 y) {
  return is_real(x) && x.a() == y;
}

template<typename T, typename T2>
inline bool operator!=(const Quaternion<T>& x, T2 y) {
  return !(x == y);
}

// Same, swapping the lhs and rhs.
template<typename T, typename T2>
inline bool operator==(T2 y, const Quaternion<T>& x) {
  return x == y;
}

template<typename T, typename T2>
inline bool operator!=(T2 y, const Quaternion<T>& x) {
  return x != y;
}

template <typename T, typename  T1>
inline bool nearly_equal(const Quaternion<T>& x, T1 y, T1 eps) {
  return is_real(x, eps) && is_near_equal_relative(x.a(), y, eps);
}

template<typename T, typename T2>
inline bool operator==(const Quaternion<T>& x, const std::complex<T2>& y) {
  return is_complex(x) && x.a() == y.real() && x.b() == y.imag();
}

template<typename T, typename T2>
inline bool operator!=(const Quaternion<T>& x, const std::complex<T2>& y) {
  return !(x == y);
}

// TODO: make name uniform with is_near_equal_relative
template <typename T, typename T1>
inline bool nearly_equal(const Quaternion<T>& x, const std::complex<T>& y, T1 eps) {
  return is_complex(x, eps)
         && is_near_equal_relative(x.a(), y.real(), eps)
         && is_near_equal_relative(x.b(), y.imag(), eps);
}

// Same, swapping the lhs and rhs.
template<typename T, typename T2>
inline bool operator==(const std::complex<T2>& y, const Quaternion<T>& x) {
  return x == y;
}

template<typename T, typename T2>
inline bool operator!=(const std::complex<T2>& y, const Quaternion<T>& x) {
  return x != y;
}

// TODO: equality of Quaternion and complex, of Quaternion and array/container

template<typename T>
inline Quaternion<T> operator+(const Quaternion<T>& x, T y) {
  return {x.a() + y, x.b(), x.c(), x.d()};
}

template<typename T>
inline Quaternion<T> operator+(const Quaternion<T>& x, const Quaternion<T>& y) {
  return Quaternion<T>(x) += y;
}

template<typename T>
inline Quaternion<T> operator-(const Quaternion<T>& x, const Quaternion<T>& y) {
  return Quaternion<T>(x) -= y;
}

/**
 * SSE operations: tried 2 implementations (SO and vectorclass): not faster.
 * Boost: as fast as boost implementation.
 */
template<typename T>
inline Quaternion<T> operator*(const Quaternion<T>& x, const Quaternion<T>& y) {
  return Quaternion<T>(x) *= y;
}

template<typename T>
inline Quaternion<T> inverse(const Quaternion<T>& x) {
  return conj(x) / norm_squared(x);
}

template<typename T, typename T1>
inline Quaternion<T> operator/(T1 k, const Quaternion<T>& x) {
  return k * inverse(x);
}

template<typename T>
inline Quaternion<T> operator/(const Quaternion<T>& x, const Quaternion<T>& y) {
  return x * inverse(y);
}

template<typename T>
inline T dot(const Quaternion<T>& x, const Quaternion<T>& y) {
  return x.a() * y.a() + x.b() * y.b() + x.c() * y.c() + x.d() * y.d();
}

/**
 * 9 operations
 */
template<typename T>
inline Quaternion<T> cross(const Quaternion<T>& x, const Quaternion<T>& y) {
  return {0,
          x.c() * y.d() - x.d() * y.c(),
          x.d() * y.b() - x.b() * y.d(),
          x.b() * y.c() - x.c() * y.b()};
}

template<typename T>
inline Quaternion<T> commutator(const Quaternion<T>& x, const Quaternion<T>& y) {
  return x * y - y * x;
}

template<typename T>
inline Quaternion<T> normalize(const Quaternion<T>& x) {
  assert(abs(x) > 0); // or this is not normalizable
  return x / abs(x);
}

/**
 * Exponential of a Quaternion.
 * This code seems to be quite a bit faster than boost, while giving
 * the same results. Boost uses a Taylor approximation for sinc,
 * which *might* (not sure) be why they are slightly slower here.
 *
 * exp(log(x)) == x always, but log(exp(x)) != x is already not true
 * for complex number, because the log is multi-valued.
 *
 * NOTE: the precision is not great with so many floating point operations
 */
template<typename T>
inline Quaternion<T> exp(const Quaternion<T>& x) {
  T un = x.unreal_norm_squared();
  if (un == 0)
    return {std::exp(x.a())};
  assert(un > 0);
  T n1 = std::sqrt(un); // un > 0, no problem
  T ea = std::exp(x.a());
  T n2 = ea * std::sin(n1) / n1;
  return {ea * std::cos(n1), n2 * x.b(), n2 * x.c(), n2 * x.d()};
}

/**
 * Log of a Quaternion.
 * exp(log(x)) == x always, but log(exp(x)) != x is already not true
 * for complex number, because the log is multi-valued.
 *
 * NOTE: the precision is not great with so many floating point operations
 */
template<typename T>
inline Quaternion<T> log(const Quaternion<T>& x) {
  T nu2 = x.unreal_norm_squared();
  if (nu2 == 0) {
    if (x.a() > 0)
      return {std::log(x.a())};
    else { // TODO: reduce number of instructions
      std::complex<T> l = log(std::complex<T>(x.a(),0)); // TODO: is that correct?
      return {l.real(), l.imag()};
    }
  }
  T a = x.a();
  assert(nu2 > 0);
  T n = std::sqrt(a * a + nu2); // n > 0
  T th = std::acos(a / n) / std::sqrt(nu2); // -a <= a/n <= 1
  return {std::log(n), th * x.b(), th * x.c(), th * x.d()}; // n > 0
}

/**
 * 10 operations:
 * a^2 - b^2 - c^2 - d^2
 * 2 a b
 * 2 a c
 * 2 a d
 */
template<typename T>
inline Quaternion<T> pow2(const Quaternion<T>& x) {
  T aa = 2 * x.a();
  return {x.a() * x.a() - x.unreal_norm_squared(),
          aa * x.b(),
          aa * x.c(),
          aa * x.d()};
}

/**
 * 14 operations:
 * a (a^2 - 3 (b^2 + c^2 + d^2))
 * -b (-3 a^2 + b^2 + c^2 + d^2)
 * -c (-3 a^2 + b^2 + c^2 + d^2)
 * -d (-3 a^2 + b^2 + c^2 + d^2)
 */
template<typename T>
inline Quaternion<T> pow3(const Quaternion<T>& x) {
  T a2 = x.a() * x.a();
  T n1 = x.unreal_norm_squared();
  T n2 = 3 * a2 - n1;
  return {x.a() * (a2 - 3 * n1),
          x.b() * n2,
          x.c() * n2,
          x.d() * n2};
}

/**
 * 18 operations:
 * a^4 - 6 a^2 (b^2 + c^2 + d^2) + (b^2 + c^2 + d^2)^2
 * -4 a b (-a^2 + b^2 + c^2 + d^2)
 * -4 a c (-a^2 + b^2 + c^2 + d^2)
 * -4 a d (-a^2 + b^2 + c^2 + d^2)
 */
template<typename T>
inline Quaternion<T> pow4(const Quaternion<T>& x) {
  T a2 = x.a() * x.a();
  T n1 = x.unreal_norm_squared();
  T n2 = 4 * x.a() * (a2 - n1);
  return {a2 * a2 - 6 * a2 * n1 + n1 * n1,
          x.b() * n2,
          x.c() * n2,
          x.d() * n2};
}

/**
 * I benchmarked that method written via the polar representation,
 * and it turned out to be much slower, and less numerically stable,
 * than this implementation. This implementation is also much faster
 * than the boost implementation. However, via the polar representation
 * I could compute pow for any real exponent, whereas this method is
 * limited to integer exponents.
 */
template<typename T>
inline Quaternion<T> pow(const Quaternion<T>& x, int expt) {

  if (expt < 0)
    return inverse(pow(x, -expt));
  if (expt == 0)
    return {1};
  if (expt == 1)
    return x;
  if (expt == 2)
    return pow2(x);
  if (expt == 3)
    return pow3(x);
  if (expt == 4)
    return pow4(x);

  Quaternion<T> x4 = pow4(x), y = x4;
  for (size_t i = 1; i < expt / 4; ++i)
    y *= x4;
  if (expt % 4 == 3)
    y *= pow3(x);
  if (expt % 4 == 2)
    y *= pow2(x);
  if (expt % 4 == 1)
    y *= x;
  return y;
}

/**
 * Real power of a Quaternion.
 */
template<typename T>
inline Quaternion<T> pow(const Quaternion<T>& x, T a) {
  if (std::floor(a) == a) // TODO: worth it for speed? helps with stability
    return pow(x, (int)a);
  return exp(a * log(x));
}

/**
 * Quaternion power of a Quaternion.
 * (that should cover all the other cases for the exponent...)
 * TODO: test against pow just above
 */
template<typename T>
inline Quaternion<T> pow(const Quaternion<T>& x, const Quaternion<T>& a) {
  if (a.is_real())
    return pow(x, a.a());
  return exp(a * log(x));
}

// TODO: sqrt
// TOOD: sin, cos, tan ...
template<typename T>
inline Quaternion<T> cos(const Quaternion<T>& x)
{
  T z = abs(x.unreal());
  T w = -std::sin(x.real()) * std::sinh(z) / z;

  return {std::cos(x.real()) * std::cosh(z), w * x.b(), w * x.c(), w * x.d()};
}

template<typename T>
inline Quaternion<T> sin(const Quaternion<T>& x)
{
  T z = abs(x.unreal());
  T w = std::cos(x.real()) * std::sinh(z) / z;

  return {std::sin(x.real()) * std::cosh(z), w * x.b(), w * x.c(), w * x.d()};
}

/**
 * TODO: optimize instruction count
 * TODO: reciprocals
 * TODO: checkout /0
 */
template<typename T>
inline Quaternion<T> tan(const Quaternion<T>& x)
{
  T z = abs(x.unreal());
  T n = std::sinh(2*z);
  T d = std::cos(2*x.a()) + std::cosh(2*z);
  T r = n/(z*d);
  return {std::sin(2*x.a())/d, r*x.b(), r*x.c(), r*x.d()};
}

template<typename T>
inline Quaternion<T> cosh(const Quaternion<T>& x)
{
  return (exp(x) + exp(-x))/2;
}

template<typename T>
inline Quaternion<T> sinh(const Quaternion<T>& x)
{
  return (exp(x) - exp(-x))/2;
}

template<typename T>
inline Quaternion<T> tanh(const Quaternion<T>& x)
{
  return sinh(x)/cosh(x);
}

/**
 * result = a*x + b*y
 */
template<typename T, typename K>
inline Quaternion<T> axby(K k1, const Quaternion<T>& x, K k2, const Quaternion<T>& y) {
  return Quaternion<T>(x).axby(k1, k2, y);
}

#endif //QuaternionS_Quaternion_H
