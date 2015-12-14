//
// Created by Frank Astier on 2015/12/13.
//

#ifndef QUATERNIONS_QUATERNION_H
#define QUATERNIONS_QUATERNION_H

#include <iostream>
#include <cmath>
#include <limits>
#include <array>
#include <complex>

using namespace std;

/**
 * Two utility functions to work with numbers approximately equal to zero.
 * If eps == 0, does a "hard" comparison to 0.
 * Otherwise, uses a ball of radius eps around 0. If the scalar is inside
 * that ball, it is equivalent to 0.
 */
template <typename T>
inline bool is_scalar_zero(T x, T eps = 0) {
  if (eps > 0)
    return std::abs(x) < eps;
  else
    return x == 0;
}

template <typename T>
inline T round_to_zero(T x) {
  return is_scalar_zero(x) ? 0 : x;
}

/**
 * A quaternion class.
 * TODO: add precise requirements on numeric types
 * TODO: policy to compress
 * TODO: byte type, as policy, to optimize mult
 * TODO: threads, as policy
 * TODO: combine operations
 * TODO: provide same operations as boost
 * TODO: provide round to zero as optional
 */
template <typename T =double> // assert operations for numeric
class Quaternion {
public:
  /**
   * Construct a quaternion from at most 4 components of type T.
   * Specifying only a != 0 makes the quaternion a real.
   * Specifying only a != and b != 0 makes the quaternion an ordinary complex number.
   */
  Quaternion(T a =0, T b =0, T c =0, T d =0)
      : _a(a), _b(b), _c(c), _d(d)
  {}

  /**
   * Construct a quaternion from a single complex<T>.
   * The j and k components are 0.
   */
  Quaternion(const std::complex<T>& x)
      : _a(const_cast<std::complex<T>&>(x).real()),
        _b(const_cast<std::complex<T>&>(x).imag()),
        _c(0), _d(0)
  {}

  /**
   * Construct a quaternion from 2 complex<T>.
   * This will set all 4 components of the quaternion.
   */
  Quaternion(const std::complex<T>& x, const std::complex<T>& y)
      : _a(const_cast<std::complex<T>&>(x).real()),
        _b(const_cast<std::complex<T>&>(x).imag()),
        _c(const_cast<std::complex<T>&>(x).real()),
        _d(const_cast<std::complex<T>&>(x).imag())
  {}

  /**
   * Copy constructor.
   */
  Quaternion(const Quaternion& other)
      : _a(other.a()),
        _b(other.b()),
        _c(other.c()),
        _d(other.d())
  {}

  /**
   * Accessors for all 4 components of the quaternion.
   */
  T a() const { return _a; }
  T b() const { return _b; }
  T c() const { return _c; }
  T d() const { return _d; }

  /**
   * The real and "unreal" parts of the quaternion.
   */
  T real() const { return _a; }
  Quaternion unreal() const { return Quaternion(0,b(),c(),d()); }

  /**
   * Transforms this quaternion into its conjugate.
   */
  void conjugate() {
    _b = -b();
    _c = -c();
    _d = -d();
  }

  /**
   * The squared of the norm of the quaternion.
   * (The square is sometimes useful, and it avoids paying for a sqrt).
   */
  T norm2() const {
    return a()*a() + b()*b() + c()*c() + d()*d();
  }

  /**
   * The norm of the quaternion (the l2 norm).
   */
  T norm() const {
    return sqrt(norm2());
  }

  /**
   * The l1 norm of the quaternion.
   */
  T l1() const {
    return std::abs(a()) + std::abs(b()) + std::abs(c()) + std::abs(d());
  }

  /**
   * The l2 norm of the quaternion.
   */
  T l2() const {
    return norm();
  }

  /**
   * The larget of the components of the quaternion.
   */
  T sup() const {
    return std::max(std::max(a(),b()), std::max(c(),d()));
  }

  /**
   * Return true if this quaternion is zero, false otherwise.
   * TODO: introduce epsilon
   */
  bool is_zero() const {
    return a() == 0 && b() == 0 && c() == 0 && d() == 0;
  }

  /**
   * Return true if this quaternion is one, false otherwise.
   * TODO: introduce epsilon
   */
  bool is_one() const {
    return a() == 1 && b() == 0 && c() == 0 && d() == 0;
  }

  /**
   * Return true if this quaternion has norm 1, false otherwise.
   * TODO: introduce epsilon
   */
  bool is_unit() const {
    return norm() == 1;
  }

  /**
   * Return true if this quaternion is real, false otherwise.
   * TODO: introduce epsilon
   */
  bool is_real() const {
    return b() == 0 && c() == 0 && d() == 0;
  }

  /**
   * Return true if this quaternion is complex, false otherwise.
   * TODO: introduce epsilon
   */
  bool is_complex() const {
    return c() == 0 && d() == 0;
  }

  /**
   * Return true if this quaternion is real, false otherwise.
   * TODO: introduce epsilon
   */
  bool is_unreal() const {
    return b() != 0 || c() == 0 || d() == 0;
  }

  /**
   * Unary minus.
   */
  Quaternion operator -() const {
    return Quaternion(-a(), -b(), -c(), -d());
  }

  /**
   * Unary +=.
   */
  Quaternion operator +=(const Quaternion& y)
  {
    _a = a()+y.a();
    _b = b()+y.b();
    _c = c()+y.c();
    _d = d()+y.d();
    return *this;
  }

  /**
  * Unary -=.
  */
  Quaternion operator -=(const Quaternion& y)
  {
    _a = a()-y.a();
    _b = b()-y.b();
    _c = c()-y.c();
    _d = d()-y.d();
    return *this;
  }

  /**
   * Scaling by a constant.
   */
  Quaternion operator *=(T k)
  {
    if (k == 0) {
      _a = _b = _c = _d = 0;
    }
    if (k != 1) {
      _a = k * _a;
      _b = k * _b;
      _c = k * _c;
      _d = k * _d;
    }
    return *this;
  }

  /**
  * Scaling by a constant.
  */
  Quaternion operator /=(T k)
  {
    if (k != 1) {
      _a = _a / k;
      _b = _b / k;
      _c = _c / k;
      _d = _d / k;
    }
    return *this;
  }

  /**
   * Unary multiplication.
   * 28 operations
   */
  Quaternion operator *=(const Quaternion& y)
  {
    // TODO: if 0 == y, if 1  == y, optional
    T ta = a() * y.a() - b() * y.b() - c() * y.c() - d() * y.d();
    T tb = a() * y.b() + b() * y.a() + c() * y.d() - d() * y.c();
    T tc = a() * y.c() - b() * y.d() + c() * y.a() + d() * y.b();
    T td = a() * y.d() + b() * y.c() - c() * y.b() + d() * y.a();
    _a = ta;
    _b = tb;
    _c = tc;
    _d = td;
    return *this;
  }

  /**
   * The type used for matrix representations of quaternions.
   */
  typedef array<array<std::complex<T>,2>,2> MatrixRepresentation;

  /**
   * Returns a matrix representation of a quaternion.
   */
  MatrixRepresentation to_matrix_representation() const {
    array<std::complex<T>,2> r0{{std::complex<T>(a(),b()), std::complex<T>(c(),d())}};
    array<std::complex<T>,2> r1{{std::complex<T>(-c(),d()), std::complex<T>(a(),-b())}};
    return MatrixRepresentation{{r0,r1}};
  }

  /**
   * Print format control flags.
   */
  static T scalar_zero_threshold;

  /**
   * Print a quaternion to a stream in various formats.
   */
  std::ostream& print(std::ostream& out) const {
      if (is_zero())
        return out << 0;
      if (is_one())
        return out << 1;
      if (*this == Quaternion<T>(-1))
        return out << -1;
      if (*this == Quaternion<T>(0,1))
        return out << "i";
      if (*this == Quaternion<T>(0,-1))
        return out << "-i";
      if (*this == Quaternion<T>(0,0,1))
        return out << "j";
      if (*this == Quaternion<T>(0,0,-1))
        return out << "-j";
      if (*this == Quaternion<T>(0,0,0,1))
        return out << "k";
      if (*this == Quaternion<T>(0,0,0,-1))
        return out << "-k";
      auto s = [](T x) { return x < 0 ? "" : "+"; };
      if (!is_scalar_zero(a(), scalar_zero_threshold))
        out << a();
      if (!is_scalar_zero(b(), scalar_zero_threshold))
        out << s(b()) << b() << "i";
      if (!is_scalar_zero(c(), scalar_zero_threshold))
        out << s(c()) << c() << "j";
      if (!is_scalar_zero(d(), scalar_zero_threshold))
        out << s(d()) << d() << "k";
      return out;
  }

private:
  T _a,_b,_c,_d;
};

template <typename T>
T Quaternion<T>::scalar_zero_threshold = std::numeric_limits<T>::epsilon();

/**
 * Predefined quaternions on floats.
 */
typedef Quaternion<float> Qf;
const Qf Qf_0 = Qf();
const Qf Qf_1 = Qf(1);
const Qf Qf_i = Qf(0,1);
const Qf Qf_j = Qf(0,0,1);
const Qf Qf_k = Qf(0,0,0,1);

/**
 * Predefined quaternions on doubles.
 */
typedef Quaternion<double> Qd;
const Qd Qd_0 = Qd();
const Qd Qd_1 = Qd(1);
const Qd Qd_i = Qd(0,1);
const Qd Qd_j = Qd(0,0,1);
const Qd Qd_k = Qd(0,0,0,1);

/**
 * Predefined quaternions on long doubles.
 */
typedef Quaternion<long double> Qld;
const Qld Qld_0 = Qld();
const Qld Qld_1 = Qld(1);
const Qld Qld_i = Qld(0,1);
const Qld Qld_j = Qld(0,0,1);
const Qld Qld_k = Qld(0,0,0,1);

/**
 * This streaming operator made me wonder if I should sneak "smart" code
 * in the quaternion arithmetic, in order to optimize it for space, but that
 * turned out not worthwhile (see CQuaternion).
 * TODO: control format for file or human readable. Also write operator>>
 */
template <typename T>
inline ostream& operator<<(ostream& out, const Quaternion<T>& q) {
  return q.print(out);
}

/**
 * Multiplication by a constant on the left.
 */
template <typename T, typename T1>
inline Quaternion<T> operator*(T1 k, const Quaternion<T>& x) {
  return Quaternion<T>(x) *= k;
}

// Same, swapping the lhs and rhs.
template <typename T, typename T1>
inline Quaternion<T> operator*(const Quaternion<T>& x, T1 k) {
  return k * x;
}

// TODO: multiplication by a complex number

/**
 * Division by a number.
 */
template <typename T, typename T1>
inline Quaternion<T> operator/(const Quaternion<T>& x, T1 k) {
  return x * T(1)/k;
}

/**
 * Returns the conjugate of x, as a new quaternion (x is unchanged).
 */
template <typename T>
inline Quaternion<T> conjugate(const Quaternion<T>& x) {
  Quaternion<T> r(x);
  r.conjugate();
  return r;
}

/**
 * The norms on a quaternion.
 */
template <typename T>
inline T norm2(const Quaternion<T>& x) {
  return x.norm2();
}

template <typename T>
inline T norm(const Quaternion<T>& x) {
  return x.norm();
}

template <typename T>
inline T l1(const Quaternion<T>& x) {
  return x.l1();
}

template <typename T>
inline T l2(const Quaternion<T>& x) {
  return x.l2();
}

template <typename T>
inline T sup(const Quaternion<T>& x) {
  return x.sup();
}

/**
 * Quaternion tests.
 */
template <typename T>
inline bool is_zero(const Quaternion<T>& x) {
  return x.is_zero();
}

template <typename T>
inline bool is_one(const Quaternion<T>& x) {
  return x.is_one();
}

template <typename T>
inline bool is_unit(const Quaternion<T>& x) {
  return x.is_unit();
}

template <typename T>
inline bool is_real(const Quaternion<T>& x) {
  return x.is_real();
}

template <typename T>
inline bool is_complex(const Quaternion<T>& x) {
  return x.is_complex();
}

template <typename T>
inline bool is_unreal(const Quaternion<T>& x) {
  return x.is_unreal();
}

/**
 * Equality of two quaternions.
 * TODO: for very large values, equality is a classic floating point problem
 */
template <typename T>
inline bool operator==(const Quaternion<T>& x, const Quaternion<T>& y) {
  return x.a() == y.a() && x.b() == y.b() && x.c() == y.c() && x.d() == y.d();
}

template <typename T>
inline bool operator!=(const Quaternion<T>& x, const Quaternion<T>& y) {
  return !(x == y);
}

/**
 * Equality of a quaternion and a real, or at least a number that can converted to a real.
 */
template <typename T, typename T2>
inline bool operator==(const Quaternion<T>& x, T2 y) {
  return is_real(x) && x.a() == static_cast<T>(y);
}

template <typename T, typename T2>
inline bool operator!=(const Quaternion<T>& x, T2 y) {
  return !(x == y);
}

// Same, swapping the lhs and rhs.
template <typename T, typename T2>
inline bool operator==(T2 y, const Quaternion<T>& x) {
  return x == y;
}

template <typename T, typename T2>
inline bool operator!=(T2 y, const Quaternion<T>& x) {
  return x != y;
}

// TODO: equality of quaternion and complex, of quaternion and array/container

template <typename T>
inline Quaternion<T> operator+(const Quaternion<T>& x, const Quaternion<T>& y) {
  return Quaternion<T>(x.a()+y.a(),x.b()+y.b(),x.c()+y.c(),x.d()+y.d());
}

template <typename T>
inline Quaternion<T> operator-(const Quaternion<T>& x, const Quaternion<T>& y) {
  return Quaternion<T>(x.a()-y.a(),x.b()-y.b(),x.c()-y.c(),x.d()-y.d());
}

/**
 * SSE operations: tried 2 implementations (SO and vectorclass): not faster.
 * Boost: as fast as boost implementation.
 * TODO: threads?
 */
template <typename T>
inline Quaternion<T> operator*(const Quaternion<T>& x, const Quaternion<T>& y) {
  return Quaternion<T>(x) *= y;
}

template <typename T>
inline Quaternion<T> inverse(const Quaternion<T>& x) {
  return conjugate(x) / norm2(x);
}

template <typename T, typename T1>
inline Quaternion<T> operator/(T1 k, const Quaternion<T>& x) {
  return k * inverse(x);
}

template <typename T>
inline Quaternion<T> operator/(const Quaternion<T>& x, const Quaternion<T>& y) {
  return x * inverse(y);
}

template <typename T>
inline T dot(const Quaternion<T>& x, const Quaternion<T>& y) {
  return x.b()*y.b() + x.c()*y.c() + x.d() * y.d();
}

/**
 * 9 operations
 */
template <typename T>
inline Quaternion<T> cross(const Quaternion<T>& x, const Quaternion<T>& y) {
  return Quaternion<T>(0,
                       x.c()*y.d() - x.d()*y.c(),
                       x.d()*y.b() - x.b()*y.d(),
                       x.b()*y.c() - x.c()*y.b());
}

template <typename T>
inline Quaternion<T> commutator(const Quaternion<T>& x, const Quaternion<T>& y) {
  return x * y - y * x;
}

template <typename T>
inline Quaternion<T> normalize(const Quaternion<T>& x) {
  return Quaternion<T>(x) / norm(x);
}

/**
 * 10 operations:
 * a^2 - b^2 - c^2 - d^2
 * 2 a b
 * 2 a c
 * 2 a d
 */
template <typename T>
inline Quaternion<T> pow2(const Quaternion<T>& x) {
  T aa = 2*x.a();
  return Quaternion<T>(x.a()*x.a() - x.b()*x.b() - x.c()*x.c() - x.d()*x.d(),
                       aa*x.b(),
                       aa*x.c(),
                       aa*x.d());
}

/**
 * 14 operations:
 * a (a^2 - 3 (b^2 + c^2 + d^2))
 * -b (-3 a^2 + b^2 + c^2 + d^2)
 * -c (-3 a^2 + b^2 + c^2 + d^2)
 * -d (-3 a^2 + b^2 + c^2 + d^2)
 */
template <typename T>
inline Quaternion<T> pow3(const Quaternion<T>& x) {
  T a2 = x.a()*x.a();
  T n1 = x.b()*x.b() + x.c()*x.c() + x.d()*x.d();
  T n2 = 3 * a2 - n1;
  return Quaternion<T>(x.a() * (a2 - 3 * n1),
                       x.b()*n2,
                       x.c()*n2,
                       x.d()*n2);
}

/**
 * 18 operations:
 * a^4 - 6 a^2 (b^2 + c^2 + d^2) + (b^2 + c^2 + d^2)^2
 * -4 a b (-a^2 + b^2 + c^2 + d^2)
 * -4 a c (-a^2 + b^2 + c^2 + d^2)
 * -4 a d (-a^2 + b^2 + c^2 + d^2)
 */
template <typename T>
inline Quaternion<T> pow4(const Quaternion<T>& x) {
  T a2 = x.a()*x.a();
  T n1 = x.b()*x.b() + x.c()*x.c() + x.d()*x.d();
  T n2 = 4 * x.a() * (a2 - n1);
  return Quaternion<T>(a2*a2 - 6 * a2 * n1 + n1 * n1,
                       x.b()*n2,
                       x.c()*n2,
                       x.d()*n2);
}

template <typename T, typename I>
inline Quaternion<T> pow(const Quaternion<T>& x, I exponent) {

//  T n1 = norm(x);
//
  if (exponent < 0)
    return inverse(pow(x, -exponent));
  if (exponent == 0)
    return Quaternion<T>(1);
  if (exponent == 1)
    return x;
  if (exponent == 2)
    return pow2(x);
  if (exponent == 3)
    return pow3(x);
  if (exponent == 4)
    return pow4(x);

  Quaternion<T> x4 = pow4(x), y = x4;
  for (size_t i = 1; i < exponent / 4; ++i)
    y *= x4;
  if (exponent % 4 == 3)
    y *= pow3(x);
  if (exponent % 4 == 2)
    y *= pow2(x);
  if (exponent % 4 == 1)
    y *= x;

  return y;
}

template <typename T>
inline Quaternion<T> exp(const Quaternion<T>& x) {
  T n1 = sqrt(x.b()*x.b() + x.c()*x.c() + x.d()*x.d());
  T n2 = exp(x.a())*sin(n1)/n1;
  return Quaternion<T>(exp(x.a())*cos(n1),
                       n2*x.b(),
                       n2*x.c(),
                       n2*x.d());
}


#endif //QUATERNIONS_QUATERNION_H
