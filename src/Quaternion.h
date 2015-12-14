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

template <typename T>
inline bool is_zero(T x) {
  return std::abs(x) < std::numeric_limits<T>::epsilon();
}

template <typename T>
inline T round_to_zero(T x) {
  return is_zero(x) ? 0 : x;
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
  Quaternion(T a =0, T b =0, T c =0, T d =0)
      : _a(a), _b(b), _c(c), _d(d)
  {}

  Quaternion(const std::complex<T>& x)
      : _a(x.real()), _b(x.imag()), _c(0), _d(0)
  {}

  Quaternion(const std::complex<T>& x, const std::complex<T>& y)
      : _a(x.real()), _b(x.imag()), _c(y.real()), _d(y.imag())
  {}

  Quaternion(const Quaternion& other)
      : _a(other.a()),
        _b(other.b()),
        _c(other.c()),
        _d(other.d())
  {}

  T a() const { return _a; }
  T b() const { return _b; }
  T c() const { return _c; }
  T d() const { return _d; }

  Quaternion operator-() const {
    return Quaternion(-a(), -b(), -c(), -d());
  }

  Quaternion operator +=(const Quaternion& y)
  {
    // TODO: if 0 == y, if 1  == y, optional
    _a = a()*y.a();
    _b = b()*y.b();
    _c = c()*y.c();
    _d = d()*y.d();
    return *this;
  }

  /**
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

  typedef array<array<std::complex<T>,2>,2> MatrixRepresentation;

  MatrixRepresentation to_matrix() const {
    array<std::complex<T>,2> r0{{std::complex<T>(a(),b()), std::complex<T>(c(),d())}};
    array<std::complex<T>,2> r1{{std::complex<T>(-c(),d()), std::complex<T>(a(),-b())}};
    return MatrixRepresentation{{r0,r1}};
  }

private:
  T _a,_b,_c,_d;
};

typedef Quaternion<float> Qf;
const Qf Qf_0 = Qf();
const Qf Qf_1 = Qf(1);
const Qf Qf_i = Qf(0,1);
const Qf Qf_j = Qf(0,0,1);
const Qf Qf_k = Qf(0,0,0,1);

typedef Quaternion<double> Qd;
const Qd Qd_0 = Qd();
const Qd Qd_1 = Qd(1);
const Qd Qd_i = Qd(0,1);
const Qd Qd_j = Qd(0,0,1);
const Qd Qd_k = Qd(0,0,0,1);

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
  if (q == Quaternion<T>())
    return out << 0;
  if (q == Quaternion<T>(1))
    return out << 1;
  if (q == Quaternion<T>(-1))
    return out << -1;
  if (q == Quaternion<T>(0,1))
    return out << "i";
  if (q == Quaternion<T>(0,-1))
    return out << "-i";
  if (q == Quaternion<T>(0,0,1))
    return out << "j";
  if (q == Quaternion<T>(0,0,-1))
    return out << "-j";
  if (q == Quaternion<T>(0,0,0,1))
    return out << "k";
  if (q == Quaternion<T>(0,0,0,-1))
    return out << "-k";
  auto s = [](T x) { return x < 0 ? "" : "+"; };
  if (!is_zero(q.a()))
    out << q.a();
  if (!is_zero(q.b()))
    out << s(q.b()) << q.b() << "i";
  if (!is_zero(q.c()))
    out << s(q.c()) << q.c() << "j";
  if (!is_zero(q.d()))
    out << s(q.d()) << q.d() << "k";
  return out;
}

template <typename T, typename T1>
inline Quaternion<T> operator*(T1 k, const Quaternion<T>& x) {
//  if (is_zero(k)) // TODO: make optional
//    return Q_0<T>;
//  if (is_zero(k-1))
//    return x;
  return Quaternion<T>(k*x.a(), k*x.b(), k*x.c(), k*x.d());
}

template <typename T, typename T1>
inline Quaternion<T> operator*(const Quaternion<T>& x, T1 k) {
  return k * x;
}

template <typename T, typename T1>
inline Quaternion<T> operator/(const Quaternion<T>& x, T1 k) {
//  if (is_zero(k-1))
//    return x;
  return Quaternion<T>(x.a()/k, x.b()/k, x.c()/k, x.d()/k);
}

template <typename T>
inline Quaternion<T> conjugate(const Quaternion<T>& x) {
  return Quaternion<T>(x.a(), -x.b(), -x.c(), -x.d());
}

template <typename T>
inline T norm2(const Quaternion<T>& x) {
  return x.a()*x.a() + x.b()*x.b() + x.c()*x.c() + x.d()*x.d();
}

template <typename T>
inline T norm(const Quaternion<T>& x) {
  return std::sqrt(norm2(x));
}

template <typename T>
inline bool is_zero(const Quaternion<T>& x) {
  return is_zero(norm2(x));
}

template <typename T>
inline bool is_unit(const Quaternion<T>& x) {
  return is_zero(norm2(x) - 1);
}

template <typename T>
inline bool is_real(const Quaternion<T>& x) {
  return is_zero(x.b()) && is_zero(x.c()) && is_zero(x.d());
}

/**
 * TODO: for very large values, equality is a classic floating point problem
 */
template <typename T>
inline bool operator==(const Quaternion<T>& x, const Quaternion<T>& y) {
  return is_zero(norm2(x - y));
}

template <typename T, typename T2>
inline bool operator==(const Quaternion<T>& x, T2 y) {
  T yy = static_cast<T>(y);
  return is_real(x) && is_zero(x.a() - yy);
}

template <typename T, typename T2>
inline bool operator==(T2 y, const Quaternion<T>& x) {
  return x == y;
}

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
 * TODO: threads
 */
template <typename T>
inline Quaternion<T> operator*(const Quaternion<T>& x, const Quaternion<T>& y) {
  Quaternion<T> r(x);
  return r *= y;
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