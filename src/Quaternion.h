/**
 * A quaternion class.
 */

#ifndef QUATERNIONS_QUATERNION_H
#define QUATERNIONS_QUATERNION_H

#include <iostream>
#include <cmath>
#include <limits>
#include <array>
#include <complex>
#include <assert.h>

/**
 * Utility function to work with numbers approximately equal to zero.
 * If eps == 0, does a "hard" comparison to 0.
 * Otherwise, uses a ball of radius eps around 0. If the scalar is inside
 * that ball, it is equivalent to 0.
 * TODO: refine for large floating point numbers
 */
template <typename T>
inline bool is_scalar_zero(T x, T eps =0) {
  if (eps > 0)
    return std::fabs(x) < eps;
  else
    return x == 0;
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
   * The value of each component of a quaternion.
   */
  typedef T value_type;

  /**
   * The polar representation of a quaternion.
   */
  typedef std::array<T,5> polar_representation;

  /**
   * The type used for matrix representations of quaternions.
   */
  typedef std::array<std::array<std::complex<T>,2>,2> matrix_representation;

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
      : _a(other._a),
        _b(other._b),
        _c(other._c),
        _d(other._d)
  {}

  /**
   * Casting constructor.
   */
  template <typename T1>
  explicit Quaternion(Quaternion<T1>& y)
      : _a(static_cast<T>(y._a)),
        _b(static_cast<T>(y._b)),
        _c(static_cast<T>(y._c)),
        _d(static_cast<T>(y._d))
  {}

  /**
   * Assingment operator.
   */
  Quaternion& operator=(const Quaternion& other)
  {
    if (&other != this) {
      _a = other._a;
      _b = other._b;
      _c = other._c;
      _d = other._d;
    }
    return *this;
  }

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
  Quaternion unreal() const { return Quaternion(0,_b,_c,_d); }

  /**
   * Transforms this quaternion into its conjugate.
   */
  void conjugate() {
    _b = -_b;
    _c = -_c;
    _d = -_d;
  }

  /**
   * The squared of the norm of the quaternion.
   * (The square is sometimes useful, and it avoids paying for a sqrt).
   */
  T norm2() const {
    return _a*_a + _b*_b + _c*_c + _d*_d;
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
  T l1() const {
    return std::abs(_a) + std::abs(_b) + std::abs(_c) + std::abs(_d);
  }

  /**
   * The l2 norm of the quaternion.
   */
  T l2() const {
    return norm();
  }

  /**
   * The value of the largest components of the quaternion.
   */
  T sup() const {
    return std::max(std::max(std::abs(_a),std::abs(_b)),
                    std::max(std::abs(_c),std::abs(_d)));
  }

  /**
   * Return true if this quaternion has norm 1, false otherwise.
   */
  bool is_unit() const {
    return is_scalar_zero(norm() - T(1), scalar_zero_threshold);
  }

  /**
   * Return true if this quaternion is real, false otherwise.
   */
  bool is_real() const {
    return is_scalar_zero(_b, scalar_zero_threshold)
           && is_scalar_zero(_c, scalar_zero_threshold)
           && is_scalar_zero(_d, scalar_zero_threshold);
  }

  /**
   * Return true if this quaternion is complex, false otherwise.
   */
  bool is_complex() const {
    return is_scalar_zero(_c, scalar_zero_threshold)
           && is_scalar_zero(_d, scalar_zero_threshold);
  }

  /**
   * Return true if this quaternion is real, false otherwise.
   */
  bool is_unreal() const {
    return !is_scalar_zero(_b, scalar_zero_threshold)
           && is_scalar_zero(_c, scalar_zero_threshold)
           && is_scalar_zero(_d, scalar_zero_threshold);
  }

  /**
   * Unary minus.
   */
  Quaternion operator -() const {
    return Quaternion(-_a, -_b, -_c, -_d);
  }

  /**
   * Unary +=.
   */
  Quaternion operator +=(const Quaternion& y)
  {
    _a += y._a;
    _b += y._b;
    _c += y._c;
    _d += y._d;
    return *this;
  }

  /**
   * Unary -=.
   */
  Quaternion operator -=(const Quaternion& y)
  {
    _a -= y._a;
    _b -= y._b;
    _c -= y._c;
    _d -= y._d;
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
      _a /= k;
      _b /= k;
      _c /= k;
      _d /= k;
    }
    return *this;
  }

  /**
   * Unary multiplication.
   * 28 operations
   */
  Quaternion operator *=(const Quaternion& y)
  {
    T ar = static_cast<T>(y._a);
    T br = static_cast<T>(y._b);
    T cr = static_cast<T>(y._c);
    T dr = static_cast<T>(y._d);

    T at = _a*ar-_b*br-_c*cr-_d*dr;
    T bt = _a*br+_b*ar+_c*dr-_d*cr;
    T ct = _a*cr-_b*dr+_c*ar+_d*br;
    T dt = _a*dr+_b*cr-_c*br+_d*ar;

    _a = at;
    _b = bt;
    _c = ct;
    _d = dt;

    return *this;
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
    T nu = _b*_b + _c*_c + _d*_d;
    if (nu != 0) {
      T ns = 1.0/std::sqrt(nu); //n*std::sin(theta);
      return {{n, theta, _b/ns, _c/ns, _d/ns}};
    }
    const T pi = std::atan2(+0., -0.);
    // theta = 0 or pi, because n = +/- a().
    return {{n, n == _a ? 0 : pi, 0, 0, 0}};
  }

  /**
   * Returns a matrix representation of a quaternion.
   */
  matrix_representation to_matrix_representation() const {
    std::array<std::complex<T>,2> r0{{std::complex<T>(a(),b()), std::complex<T>(c(),d())}};
    std::array<std::complex<T>,2> r1{{std::complex<T>(-c(),d()), std::complex<T>(a(),-b())}};
    return matrix_representation{{r0, r1}};
  }

  /**
   * Print format control flags.
   */
  static T scalar_zero_threshold; // if 0, does "hard" equality tests for zero
  static int print_style;

  /**
   * Print a quaternion to a stream in various formats.
   * TODO: introduce eps and make faster with constants?
   */
  std::ostream& print(std::ostream& out) const {
    if (print_style == 0) {
      if (*this == 0)
        return out << 0;
      if (*this == Quaternion<T>(1))
        return out << 1;
      if (*this == Quaternion<T>(-1))
        return out << -1;
      if (*this == Quaternion<T>(0, 1))
        return out << "i";
      if (*this == Quaternion<T>(0, -1))
        return out << "-i";
      if (*this == Quaternion<T>(0, 0, 1))
        return out << "j";
      if (*this == Quaternion<T>(0, 0, -1))
        return out << "-j";
      if (*this == Quaternion<T>(0, 0, 0, 1))
        return out << "k";
      if (*this == Quaternion<T>(0, 0, 0, -1))
        return out << "-k";
      auto s = [](T x) { return x < 0 ? "" : "+"; }; // print out the sign correctly
      if (!is_scalar_zero(a(), scalar_zero_threshold))
        out << a();
      if (!is_scalar_zero(b(), scalar_zero_threshold))
        out << s(b()) << b() << "i";
      if (!is_scalar_zero(c(), scalar_zero_threshold))
        out << s(c()) << c() << "j";
      if (!is_scalar_zero(d(), scalar_zero_threshold))
        out << s(d()) << d() << "k";
    } else if (print_style == 1) {
      out << "{" << a() << "," << b() << "," << c() << "," << d() << "}";
    }
    return out;
  }

  // TODO: read from stream

private:
  T _a,_b,_c,_d;
};

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
 * IO manipulators to control the format when printing quaternions out to a stream.
 */
template <typename T>
T Quaternion<T>::scalar_zero_threshold;

template <typename T>
int Quaternion<T>::print_style;

template <typename T>
struct SetScalarZeroThreshold { T eps = 0; };

template <typename T>
inline SetScalarZeroThreshold<T> set_eps(T eps) {
  SetScalarZeroThreshold<T> sszt;
  sszt.eps = eps;
  return sszt;
}

template <typename T>
inline std::ostream& operator<<(std::ostream& out, SetScalarZeroThreshold<T> sszt) {
  Quaternion<T>::scalar_zero_threshold = sszt.eps;
  return out;
}

template <typename T>
struct SetPrintStyle { int style; };

template <typename T>
inline SetPrintStyle<T> set_style_nice() {
  SetPrintStyle<T> sps;
  sps.style = 0;
  return sps;
}

template <typename T>
inline SetPrintStyle<T> set_style_compact() {
  SetPrintStyle<T> sps;
  sps.style = 1;
  return sps;
}

template <typename T>
inline std::ostream& operator<<(std::ostream& out, SetPrintStyle<T> sps) {
  Quaternion<T>::print_style = sps.style;
  return out;
}

/**
 * This streaming operator made me wonder if I should sneak "smart" code
 * in the quaternion arithmetic, in order to optimize it for space, but that
 * turned out not worthwhile (see CQuaternion).
 * TODO: control format for file or human readable. Also write operator>>
 */
template <typename T>
inline std::ostream& operator<<(std::ostream& out, const Quaternion<T>& q) {
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
  const T eps = Quaternion<T>::scalar_zero_threshold;
  if (eps)
    return x.a() == y.a() && x.b() == y.b() && x.c() == y.c() && x.d() == y.d();
  else
    return is_scalar_zero(x.a()-y.a(), eps)
           && is_scalar_zero(x.b()-y.b(), eps)
           && is_scalar_zero(x.c()-y.c(), eps)
           && is_scalar_zero(x.d()-y.d(), eps);
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
  return Quaternion<T>(x) += y;
}

template <typename T>
inline Quaternion<T> operator-(const Quaternion<T>& x, const Quaternion<T>& y) {
  return Quaternion<T>(x) -= y;
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
  return {a2*a2 - 6 * a2 * n1 + n1 * n1,
          x.b()*n2,
          x.c()*n2,
          x.d()*n2};
}

/**
 * I benchmarked that method written via the polar representation,
 * and it turned out to be much slower, and less numerically stable,
 * than this implementation. This implementation is also much faster
 * than the boost implementation. However, via the polar representation
 * I could compute pow for any real exponent, whereas this method is
 * limited to integer exponents.
 */
template <typename T>
inline Quaternion<T> pow(const Quaternion<T>& x, int expt) {

  if (expt < 0)
    return inverse(pow(x, -expt));
  if (expt == 0)
    return Quaternion<T>(1);
  if (expt == 1)
    return x;
  if (expt == 2)
    return pow2(x);
  if (expt == 3)
    return pow3(x);
  if (expt == 4)
    return pow4(x);

  assert(0 < expt);
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
