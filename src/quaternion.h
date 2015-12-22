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
 * A quaternion class.
 */

#ifndef QUATERNIONS_QUATERNION_H
#define QUATERNIONS_QUATERNION_H

#include <limits>
#include <array>
#include <complex>
#include <iterator>
#include <assert.h>

#include "quaternion_utils.h"

// TODO: add namespace

/**
 * A quaternion class.
 * TODO: combine operations: axby...
 * TODO: provide same operations as boost
 * TODO: provide round to zero, floor, ceil
 * TODO: check if 0 detection is worth doing in *
 * TODO: check with std algos
 * TODO: remove copies/constructions in expressions
 * TODO: fast sinus? Sinus in CPU instruction?
 * TODO: study matrix representation and isomorphism
 * TODO: do we need the static_cast?
 * TODO: check references to make sure functionality covered
 * TODO: preconditions
 */
template<typename T =double> // assert operations for numeric is_specialized??
// T has to be real or integer for exp, log, can't accept e.g. complex
// if custom type, check requirements
class Quaternion {
public:
  /**
   * The value of each component of a quaternion.
   */
  typedef T value_type;

  /**
   * The polar representation of a quaternion.
   */
  typedef std::array<T, 5> polar_representation;

  /**
   * The type used for matrix representations of quaternions.
   * TODO: array or valarray? - maybe make template parameters, and polar_representation too
   */
  typedef std::array<std::array<std::complex<T>, 2>, 2> matrix_representation;

  /**
  * Construct a quaternion from at most 4 components of type T.
  * Specifying only a != 0 makes the quaternion a real.
  * Specifying only a != and b != 0 makes the quaternion an ordinary complex number.
  */
  Quaternion(T a = 0, T b = 0, T c = 0, T d = 0)
      : _a(a), _b(b), _c(c), _d(d) { }

  /**
   * Construct a quaternion from at most 4 components of type T.
   * Specifying only a != 0 makes the quaternion a real.
   * Specifying only a != and b != 0 makes the quaternion an ordinary complex number.
   */
  template<typename T1 = T, IS_NOT_ITERATOR(T1)>
  Quaternion(T1 a = 0, T1 b = 0, T1 c = 0, T1 d = 0)
      : _a(a), _b(b), _c(c), _d(d) { }

  /**
   * Construct a quaternion from 2 complex<T>.
   * This sets all 4 components of the quaternion.
   */
  template<typename T1>
  Quaternion(const std::complex<T1>& x, const std::complex<T1> &y = std::complex<T1>(0, 0))
      : _a(x.real()), _b(x.imag()), _c(y.real()), _d(y.imag()) { }

  /**
   * Construct from a pointer to a range of 4 elements ("float[4]").
   * TODO: make sure we can't use the next constructor here
   */
  template<typename T1 = T>
  Quaternion(T1 *it)
      : _a(*it), _b(*++it), _c(*++it), _d(*++it) { }

  /**
   * Construct from an iterator to a range of 4 elements.
   */
  template<typename It, IS_ITERATOR(It)>
  Quaternion(It it)
      : _a(static_cast<T>(*it)),
        _b(static_cast<T>(*++it)),
        _c(static_cast<T>(*++it)),
        _d(static_cast<T>(*++it)) { }

  /**
   * Copy constructor.
   */
  template<typename T1>
  Quaternion(const Quaternion<T1> &y)
      : _a(static_cast<T>(y.a())),
        _b(static_cast<T>(y.b())),
        _c(static_cast<T>(y.c())),
        _d(static_cast<T>(y.d())) { }

  /**
   * Assignment operator.
   */
  template<typename T1>
  Quaternion &operator=(const Quaternion<T1> &other) {
    _a = static_cast<T>(other.a());
    _b = static_cast<T>(other.b());
    _c = static_cast<T>(other.c());
    _d = static_cast<T>(other.d());
    return *this;
  }

  template<typename T1>
  static Quaternion spherical(T1 rho, T1 theta, T1 phi1, T1 phi2) {

    T d = std::sin(phi2);
    T cr = std::cos(phi2);
    T c = cr * std::sin(phi1);
    cr *= std::cos(phi1);
    T b = cr * std::sin(theta);
    T a = cr * std::cos(theta);

    return {rho * a, rho * b, rho * c, rho * d};
  }

  template<typename T1>
  static Quaternion semipolar(T1 rho, T1 alpha, T1 theta1, T1 theta2) {

    T ca = std::cos(alpha);
    T sa = std::sin(alpha);
    T a = ca * std::cos(theta1);
    T b = ca * std::sin(theta1);
    T c = sa * std::cos(theta2);
    T d = sa * std::sin(theta2);

    return {rho * a, rho * b, rho * c, rho * d};
  }

  template<typename T1>
  static Quaternion multipolar(T1 rho1, T1 theta1, T1 rho2, T1 theta2) {

    T a = rho1 * std::cos(theta1);
    T b = rho1 * std::sin(theta1);
    T c = rho2 * std::cos(theta2);
    T d = rho2 * std::sin(theta2);

    return {a, b, c, d};
  }

  template<typename T1>
  static Quaternion cylindrospherical(T1 t, T1 radius, T1 longitude, T1 latitude) {

    T cl = std::cos(latitude);
    T b = radius * cl * std::cos(longitude);
    T c = radius * cl * std::sin(longitude);
    T d = radius * std::sin(latitude);

    return {t, b, c, d};
  }

  template<typename T1>
  static Quaternion cylindrical(T1 r, T1 angle, T1 h1, T1 h2) {

    T a = r * std::cos(angle);
    T b = r * std::sin(angle);

    return {a, b, h1, h2};
  }

  // TODO: copy to valarray, array, vector...

  /**
   * Accessors for all 4 components of the quaternion.
   */
  T a() const { return _a; }
  T b() const { return _b; }
  T c() const { return _c; }
  T d() const { return _d; }

  /**
   * The complex components of this quaternion.
   */
  std::complex<T> c1() const { return {_a,_b}; }
  std::complex<T> c2() const { return {_c,_d}; }

  operator T() const {
    assert(_b == 0 && _c == 0 && _d == 0);
    return _a;
  }

  operator std::complex<T>() const {
    assert(_c == 0 && _d == 0);
    return std::complex<T>(_a,_b);
  }

  operator std::array<T,4>() const {
    return {{_a, _b, _c, _d}};
  }

  /**
   * The polar representation of a quaternion.
   * Returns 5 numbers:
   * - the Euclidean norm of the quaternion,
   * - the polar angle theta,
   * - and each of the components of the "unreal unit direction".
   */
  polar_representation to_polar_representation() const {
    T n = norm();
    T theta = std::acos(_a / n);
    T nu = _b * _b + _c * _c + _d * _d;
    if (nu != 0) {
      T ns = 1.0 / std::sqrt(nu); //n*std::sin(theta);
      return {{n, theta, _b / ns, _c / ns, _d / ns}};
    }
    const T pi = std::atan2(+0., -0.);
    // theta = 0 or pi, because n = +/- a().
    return {{n, n == _a ? 0 : pi, 0, 0, 0}};
  }

  /**
   * Returns a matrix representation of a quaternion.
   */
  matrix_representation to_matrix_representation() const {
    std::array<std::complex<T>, 2> r0{{std::complex<T>(a(), b()), std::complex<T>(c(), d())}};
    std::array<std::complex<T>, 2> r1{{std::complex<T>(-c(), d()), std::complex<T>(a(), -b())}};
    return {{r0, r1}};
  }

  /**
   * The real and "unreal" parts of the quaternion.
   */
  T real() const { return _a; }
  Quaternion unreal() const { return {0, _b, _c, _d}; }

  /**
   * The conjugate of this quaternion.
   */
  Quaternion conjugate() const {
    return {_a, -_b, -_c, -_d};
  }

  /**
   * The squared of the norm of the quaternion.
   * (The square is sometimes useful, and it avoids paying for a sqrt).
   */
  T norm2() const {
    return _a * _a + _b * _b + _c * _c + _d * _d;
  }

  /**
   * The norm of the quaternion (the l2 norm).
   */
  T norm() const {
    return std::sqrt(norm2());
  }

  /**
   * The l1 norm of the quaternion.
   */
  T norm_l1() const {
    return std::abs(_a) + std::abs(_b) + std::abs(_c) + std::abs(_d);
  }

  /**
   * The value of the largest components of the quaternion.
   */
  T norm_sup() const {
    return std::max(std::max(std::abs(_a), std::abs(_b)),
                    std::max(std::abs(_c), std::abs(_d)));
  }

  /**
   * The L2 norm of the "unreal" components of the quaternion,
   * comes back often in computations.
   */
  T unreal_norm2() const {
    return _b * _b + _c * _c + _d * _d;
  }

  /**
   * Return true if this quaternion has norm 1, false otherwise.
   */
  template <typename T1 =T>
  bool is_unit(T1 eps =0) const {
    return is_scalar_zero(norm2() - T(1), eps);
  }

  /**
   * Return true if this quaternion is real, false otherwise.
   */
  template <typename T1 =T>
  bool is_real(T1 eps =0) const {
    return is_scalar_zero(_b, eps)
           && is_scalar_zero(_c, eps)
           && is_scalar_zero(_d, eps);
  }

  /**
   * Return true if this quaternion is complex, false otherwise.
   */
  template <typename T1 =T>
  bool is_complex(T1 eps =0) const {
    return is_scalar_zero(_c, eps) && is_scalar_zero(_d, eps);
  }

  /**
   * Return true if this quaternion is real, false otherwise.
   */
  template <typename T1 =T>
  bool is_unreal(T1 eps =0) const {
    return !is_scalar_zero(_b, eps)
           && is_scalar_zero(_c, eps)
           && is_scalar_zero(_d, eps);
  }

  /**
   * Unary minus.
   */
  Quaternion operator-() const {
    return Quaternion(-_a, -_b, -_c, -_d);
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
  Quaternion operator+=(const std::complex<T1> &y) {
    _a += y.real();
    _b += y.imag();
    return *this;
  }

  /**
  * Unary +=.
  */
  template<typename T1>
  Quaternion operator-=(const std::complex<T1> &y) {
    _a -= y.real();
    _b -= y.imag();
    return *this;
  }

  // TODO *=, /= with complex

  /**
   * Unary +=.
   */
  template<typename T1>
  Quaternion operator+=(const Quaternion<T1> &y) {
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
  Quaternion operator-=(const Quaternion<T1> &y) {
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
   * TODO: inline and reduce number of operations
   */
  template <typename T1>
  Quaternion operator/=(const Quaternion<T1>& y) {
    return operator*=(y.conjugate()/y.norm2());
  }

  /**
   * k1 * this quaternion + k2 * y
   * Improves performance by reducing number of constructions/copies.
   */
  template<typename K>
  Quaternion axby(K k1, K k2, const Quaternion &y) {
    _a = k1 * _a + k2 * y._a;
    _b = k1 * _b + k2 * y._b;
    _c = k1 * _c + k2 * y._c;
    _d = k1 * _d + k2 * y._d;
    return *this;
  }


private:
  T _a, _b, _c, _d; // the full state for a Quaternion
};

/**
 * Predefined quaternions on floats.
 */
typedef Quaternion<float> Qf;
const Qf Qf_0;
const Qf Qf_1(1);
const Qf Qf_i(0, 1);
const Qf Qf_j(0, 0, 1);
const Qf Qf_k(0, 0, 0, 1);

/**
 * Predefined quaternions on doubles.
 */
typedef Quaternion<double> Qd;
const Qd Qd_0;
const Qd Qd_1(1);
const Qd Qd_i(0, 1);
const Qd Qd_j(0, 0, 1);
const Qd Qd_k(0, 0, 0, 1);

/**
 * Predefined quaternions on long doubles.
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

/**
 * Returns the conjugate of x, as a new quaternion (x is unchanged).
 */
template<typename T>
inline Quaternion<T> conjugate(const Quaternion<T>& x) {
  return Quaternion<T>(x).conjugate();
}

/**
 * The norms on a quaternion.
 */
template<typename T>
inline T norm2(const Quaternion<T>& x) {
  return x.norm2();
}

template<typename T>
inline T norm(const Quaternion<T>& x) {
  return x.norm();
}

template<typename T>
inline T nom_l1(const Quaternion<T>& x) {
  return x.norm_l1();
}

template<typename T>
inline T norm_sup(const Quaternion<T>& x) {
  return x.norm_sup();
}

/**
 * quaternion tests.
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
inline bool operator==(const Quaternion<T>& x, const Quaternion<T> &y) {
  return x.a() == y.a() && x.b() == y.b() && x.c() == y.c() && x.d() == y.d();
}

template<typename T>
inline bool operator!=(const Quaternion<T>& x, const Quaternion<T> &y) {
  return !(x == y);
}

/**
 * This is more costly, but very useful in practice: quaternions transcendentals
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

// TODO: equality of quaternion and complex, of quaternion and array/container

template<typename T>
inline Quaternion<T> operator+(const Quaternion<T>& x, T y) {
  return {x.a() + y, x.b(), x.c(), x.d()};
}

template<typename T>
inline Quaternion<T> operator+(const Quaternion<T>& x, const Quaternion<T> &y) {
  return Quaternion<T>(x) += y;
}

template<typename T>
inline Quaternion<T> operator-(const Quaternion<T>& x, const Quaternion<T> &y) {
  return Quaternion<T>(x) -= y;
}

/**
 * SSE operations: tried 2 implementations (SO and vectorclass): not faster.
 * Boost: as fast as boost implementation.
 * TODO: micro-threads?
 * TODO: on Mac/clang and only there, this is twice as slow as boost?
 */
template<typename T>
inline Quaternion<T> operator*(const Quaternion<T>& x, const Quaternion<T> &y) {
  return Quaternion<T>(x) *= y;
}

template<typename T>
inline Quaternion<T> inverse(const Quaternion<T>& x) {
  return conjugate(x) / norm2(x);
}

template<typename T, typename T1>
inline Quaternion<T> operator/(T1 k, const Quaternion<T>& x) {
  return k * inverse(x);
}

template<typename T>
inline Quaternion<T> operator/(const Quaternion<T>& x, const Quaternion<T> &y) {
  return x * inverse(y);
}

template<typename T>
inline T dot(const Quaternion<T>& x, const Quaternion<T> &y) {
  return x.a() * y.a() + x.b() * y.b() + x.c() * y.c() + x.d() * y.d();
}

/**
 * 9 operations
 */
template<typename T>
inline Quaternion<T> cross(const Quaternion<T>& x, const Quaternion<T> &y) {
  return {0,
          x.c() * y.d() - x.d() * y.c(),
          x.d() * y.b() - x.b() * y.d(),
          x.b() * y.c() - x.c() * y.b()};
}

template<typename T>
inline Quaternion<T> commutator(const Quaternion<T>& x, const Quaternion<T> &y) {
  return x * y - y * x;
}

template<typename T>
inline Quaternion<T> normalize(const Quaternion<T>& x) {
  return x / norm(x);
}

/**
 * Exponential of a quaternion.
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
  T un = x.unreal_norm2();
  if (un == 0)
    return {std::exp(x.a())};
  assert(un > 0);
  T n1 = std::sqrt(un); // un > 0, no problem
  T ea = std::exp(x.a());
  T n2 = ea * std::sin(n1) / n1;
  return {ea * std::cos(n1), n2 * x.b(), n2 * x.c(), n2 * x.d()};
}

/**
 * Log of a quaternion.
 * exp(log(x)) == x always, but log(exp(x)) != x is already not true
 * for complex number, because the log is multi-valued.
 *
 * NOTE: the precision is not great with so many floating point operations
 */
template<typename T>
inline Quaternion<T> log(const Quaternion<T>& x) {
  T nu2 = x.unreal_norm2();
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
  return {x.a() * x.a() - x.unreal_norm2(),
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
  T n1 = x.unreal_norm2();
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
  T n1 = x.unreal_norm2();
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
 * Real power of a quaternion.
 */
template<typename T>
inline Quaternion<T> pow(const Quaternion<T>& x, T a) {
  if (std::floor(a) == a) // TODO: worth it for speed? helps with stability
    return pow(x, (int)a);
  return exp(a * log(x));
}

/**
 * quaternion power of a quaternion.
 * (that should cover all the other cases for the exponent...)
 * TODO: test against pow just above
 */
template<typename T>
inline Quaternion<T> pow(const Quaternion<T>& x, const Quaternion<T>& a) {
  if (a.is_real())
    return pow(x, a.a());
  // if (a.is_complex()) // TODO: worth it?
  return exp(a * log(x));
}

// TODO: sqrt
// TOOD: sin, cos, tan ...

/**
 * result = a*x + b*y
 */
template<typename T, typename K>
inline Quaternion<T> axby(K k1, const Quaternion<T>& x, K k2, const Quaternion<T> &y) {
  return Quaternion<T>(x).axby(k1, k2, y);
}

#endif //QUATERNIONS_QUATERNION_H
