# Quaternions  [![Build status](https://travis-ci.org/FrankAstier/quaternions.svg?branch=master)](https://travis-ci.org/FrankAstier/quaternions) [![Coverage Status](https://coveralls.io/repos/FrankAstier/quaternions/badge.svg?branch=master&service=github&bust=1)](https://coveralls.io/github/FrankAstier/quaternions?branch=master)

A C++11 library to work with quaternions, as a single header file.

## Introduction
- In mathematics, the quaternions are a number system that extends the complex numbers.
  They were first described by Irish mathematician William Rowan Hamilton in 1843 and applied to mechanics
  in three-dimensional space. A feature of quaternions is that multiplication of two quaternions is noncommutative.
  Hamilton defined a quaternion as the quotient of two directed lines in a three-dimensional space or equivalently
  as the quotient of two vectors. Quaternions find uses in both theoretical and applied mathematics, in particular
  for calculations involving three-dimensional rotations such as in three-dimensional computer graphics, computer
  vision and crystallographic texture analysis. In practical applications, they can be used alongside other methods,
  such as Euler angles and rotation matrices, or as an alternative to them, depending on the application.

## Design objectives
- This library was designed to be simple, fast and convenient.
- For simplicity, there is a single(parametric) class, quaternion::Quaternion, that lives in a single header file.
- Computing the power of a quaternion in particular has been optimized to be significantly faster than boost.
- The methods provided have been designed to make it as natural as possible to use quaternion::Quaternion.

## Usage
- include \<quaternions/src/quaternion.h\>
- using namespace quaternion;
- To build the unit tests: cmake . ; make OR make -f Makefile.mk
- To make the doc: doxygen Doxyfile.
- The file unit_tests.cpp shows usage examples for all the function in the quaternion namespace.
- The file unit_tests.cpp contains various benchmarks.

## Notes
- Boost provides quaternions at: http://www.boost.org/doc/libs/1_59_0/libs/math/doc/html/quaternions.html
  This implementation is as fast as or faster than Boost. In particular pow is much faster than Boost. The speedups
  in this implementation were obtained by factorizing the a few low powers, and re-using those factorizations for
  higher powers.
- I have tried to use intrinsics(SSE), but didn't find it to be faster than "naive" code.
- Expression templates: unless the expression templates do some serious work, gcc can optimize the code to get very
  good performance, essentially equal to expression templates. So expression templates seem to be an older technique
  that's no longer required by modern compilers(although clang seems to lag compared to gcc in terms of optimization).

## IEEE-754/IEC-559 compliance
- If the compiler supports IEC-559, which can be verified by calling std::numeric_limits::is_iec559,
  then quaternion::Quaternion will leverage this.
- However the Quaternion algorithms do not currently take extra steps to comply with IEC-559's prescribed way of handling
  NaNs, underflows, overflows...

## Tested with:
- clang
- gcc 4.8

## Requirements
- Boost is required for the unit tests, but not for the Quaternion class itself.

## Others/References
- https://en.wikipedia.org/wiki/Quaternion
- http://www.boost.org/doc/libs/1_59_0/libs/math/doc/html/quaternions.html
- *vectorclass* from Agner Fog, can be found at: http://www.agner.org/optimize/, and provides quaternion classes.
- http://www.geometrictools.com/index.html

## Quaternion synopsis

### Quaternion<\T\> members
#### Constructors
    Quaternion(T a=0, T b=0, T c=0, T d=0)

    template<typename T1 = T, IS_NOT_ITERATOR(T1)>
    Quaternion(T1 a=0, T1 b=0, T1 c=0, T1 d=0)

    template<typename T1>
    Quaternion(const std::complex<T1> &x, const std::complex<T1> &y=std::complex<T1>(0, 0))

    template<typename T1 = T>
    Quaternion(T1 *it)

    template<typename It, IS_ITERATOR(It)>
    Quaternion(It it)

    template<typename T1>
    Quaternion(const Quaternion<T1> &y)

    template<typename T1>
    Quaternion & operator=(const Quaternion<T1> &other)

#### Accessors
    T a() const
    T b() const
    T c() const
    T d() const

    std::complex<T> c1() const
    std::complex<T> c2() const

    T to_real() const
    std::complex<T> to_complex() const
    std::array<T, 4> to_array() const

    T real() const
    Quaternion unreal() const

    T norm_squared() const
    T abs() const
    T unreal_norm_squared() const

#### Tests
    template<typename T1 = T>
    bool is_unit(T1 eps=0) const

    template<typename T1 = T>
    bool is_real(T1 eps=0) const

    template<typename T1 = T>
    bool is_complex(T1 eps=0) const

    template<typename T1 = T>
    bool is_unreal(T1 eps=0) const

#### Arithmetic
    Quaternion operator+() const
    Quaternion operator-() const
    Quaternion operator+=(T y)
    Quaternion operator-=(T y)
    Quaternion operator*=(T k)
    Quaternion operator/=(T k)

    template<typename T1>
    Quaternion operator+=(const std::complex<T1> &y)

    template<typename T1>
    Quaternion operator-=(const std::complex<T1> &y)

    template<typename T1>
    Quaternion operator*=(const std::complex<T1> &y)

    template<typename T1>
    Quaternion operator/=(const std::complex<T1> &y)

    template<typename T1>
    Quaternion operator+=(const Quaternion<T1> &y)

    template<typename T1>
    Quaternion operator-=(const Quaternion<T1> &y)

    template<typename T1>
    Quaternion operator*=(const Quaternion<T1> &y)

    template<typename T1>
    Quaternion operator/=(const Quaternion<T1> &y)


### Typedef
    typedef Quaternion<float> Qf
    typedef Quaternion<double> Qd
    typedef Quaternion<long double> Qld

    template<typename T>
    using polar_representation = std::array<T, 5>

    template<typename T>
    using matrix_representation = std::array<std::array<std::complex<T>, 2>, 2>

    template<typename T>
    using rotation_matrix = std::array<std::array<T, 3>, 3>

### Enumerations
    enum  DisplayStyle { q_nice, q_compact }

### Predefined constants
    const Qf Qf_1(1)
    const Qf Qf_i(0, 1)
    const Qf Qf_j(0, 0, 1)
    const Qf Qf_k(0, 0, 0, 1)

    const Qd Qd_1(1)
    const Qd Qd_i(0, 1)
    const Qd Qd_j(0, 0, 1)
    const Qd Qd_k(0, 0, 0, 1)

    const Qld Qld_1(1)
    const Qld Qld_i(0, 1)
    const Qld Qld_j(0, 0, 1)
    const Qld Qld_k(0, 0, 0, 1)

### Functions
#### Constructors
    template<typename T>
    Quaternion<T> spherical(T rho, T theta, T phi1, T phi2)

    template<typename T>
    Quaternion<T> semipolar(T rho, T alpha, T theta1, T theta2)

    template<typename T>
    Quaternion<T> multipolar(T rho1, T theta1, T rho2, T theta2)

    template<typename T>
    Quaternion<T> cylindrospherical(T t, T radius, T longitude, T latitude)

    template<typename T>
    Quaternion<T> cylindrical(T r, T angle, T h1, T h2)

    template<typename T>
    polar_representation<T> to_polar_representation(const Quaternion<T> &x)

    template<typename T>
    matrix_representation<T> to_matrix_representation(const Quaternion<T> &x)

    template<typename T>
    rotation_matrix<T> to_rotation_matrix(const Quaternion<T> &x)

    template<typename T>
    Quaternion<T> from_rotation_matrix(const rotation_matrix<T> &rm)

#### Various
    template<typename T>
    Quaternion<T> conj(const Quaternion<T> &x)

    template<typename T>
    T norm_squared(const Quaternion<T> &x)

    template<typename T>
    T abs(const Quaternion<T> &x)

    template<typename T>
    T unreal_norm_squared(const Quaternion<T> &x)

    template<typename T>
    T norm_l0(const Quaternion<T> &x)

    template<typename T>
    T norm_l1(const Quaternion<T> &x)

    template<typename T , typename T1>
    T norm_lk(const Quaternion<T> &x, T1 k)

    template<typename T>
    T norm_sup(const Quaternion<T> &x)

    template<typename T>
    Quaternion<T> normalize(const Quaternion<T> &x)

#### Tests
    template<typename T , typename T1 = T>
    bool is_unit(const Quaternion<T> &x, T1 eps=0)

    template<typename T , typename T1 = T>
    bool is_real(const Quaternion<T> &x, T1 eps=0)

    template<typename T , typename T1 = T>
    bool is_complex(const Quaternion<T> &x, T1 eps=0)

    template<typename T , typename T1 = T>
    bool is_unreal(const Quaternion<T> &x, T1 eps=0)

#### Equality
    template<typename T , typename T2 , IS_CONVERTIBLE(T2, T)>
    bool operator==(const Quaternion<T> &x, T2 y)

    template<typename T , typename T2 , IS_CONVERTIBLE(T2, T)>
    bool operator==(T2 y, const Quaternion<T> &x)

    template<typename T , typename T2 , IS_CONVERTIBLE(T2, T)>
    bool operator!=(const Quaternion<T> &x, T2 y)

    template<typename T , typename T2 , IS_CONVERTIBLE(T2, T)>
    bool operator!=(T2 y, const Quaternion<T> &x)

    template<typename T , typename T2 , typename T3 , IS_CONVERTIBLE(T2, T) , IS_CONVERTIBLE(T3, T)>
    bool nearly_equal(const Quaternion<T> &x, T2 y, T3 eps)

    template<typename T , typename T2 , typename T3 , IS_CONVERTIBLE(T2, T) , IS_CONVERTIBLE(T3, T)>
    bool nearly_equal(T2 y, const Quaternion<T> &x, T3 eps)

    template<typename T , typename T2 , IS_CONVERTIBLE(T2, T)>
    bool operator==(const Quaternion<T> &x, const std::complex<T2> &y)

    template<typename T , typename T2 , IS_CONVERTIBLE(T2, T)>
    bool operator!=(const Quaternion<T> &x, const std::complex<T2> &y)

    template<typename T , typename T2 , IS_CONVERTIBLE(T2, T)>
    bool operator==(const std::complex<T2> &y, const Quaternion<T> &x)

    template<typename T , typename T2 , IS_CONVERTIBLE(T2, T)>
    bool operator!=(const std::complex<T2> &y, const Quaternion<T> &x)

    template<typename T , typename T2 , typename T3 , IS_CONVERTIBLE(T2, T) , IS_CONVERTIBLE(T3, T)>
    bool nearly_equal(const Quaternion<T> &x, const std::complex<T2> &y, T3 eps)

    template<typename T , typename T2 , typename T3 , IS_CONVERTIBLE(T2, T) , IS_CONVERTIBLE(T3, T)>
    bool nearly_equal(const std::complex<T2> &y, const Quaternion<T> &x, T3 eps)

    template<typename T>
    bool operator==(const Quaternion<T> &x, const Quaternion<T> &y)

    template<typename T>
    bool operator!=(const Quaternion<T> &x, const Quaternion<T> &y)

    template<typename T , typename T2 , typename T3 , IS_CONVERTIBLE(T2, T) , IS_CONVERTIBLE(T3, T)>
    bool nearly_equal(const Quaternion<T> &x, const Quaternion<T2> &y, T3 eps)

#### Arithmetic
    template<typename T , typename T1>
    Quaternion<T> operator+(const Quaternion<T> &x, T1 y)

    template<typename T , typename T1>
    Quaternion<T> operator+(T1 y, const Quaternion<T> &x)

    template<typename T>
    Quaternion<T> operator+(const Quaternion<T> &x, std::complex<T> &y)

    template<typename T>
    Quaternion<T> operator+(std::complex<T> &y, const Quaternion<T> &x)

    template<typename T>
    Quaternion<T> operator+(const Quaternion<T> &x, const Quaternion<T> &y)

    template<typename T , typename T1>
    Quaternion<T> operator-(const Quaternion<T> &x, T1 y)

    template<typename T , typename T1>
    Quaternion<T> operator-(T1 y, const Quaternion<T> &x)

    template<typename T>
    Quaternion<T> operator-(const Quaternion<T> &x, std::complex<T> &y)

    template<typename T>
    Quaternion<T> operator-(std::complex<T> &y, const Quaternion<T> &x)

    template<typename T>
    Quaternion<T> operator-(const Quaternion<T> &x, const Quaternion<T> &y)

    template<typename T , typename T1>
    Quaternion<T> operator*(const Quaternion<T> &x, T1 y)

    template<typename T , typename T1>
    Quaternion<T> operator*(T1 y, const Quaternion<T> &x)

    template<typename T>
    Quaternion<T> operator*(const Quaternion<T> &x, std::complex<T> &y)

    template<typename T>
    Quaternion<T> operator*(std::complex<T> &y, const Quaternion<T> &x)

    template<typename T>
    Quaternion<T> operator*(const Quaternion<T> &x, const Quaternion<T> &y)

    template<typename T>
    Quaternion<T> inverse(const Quaternion<T> &x)

    template<typename T , typename T1>
    Quaternion<T> operator/(const Quaternion<T> &x, T1 y)

    template<typename T , typename T1>
    Quaternion<T> operator/(T1 y, const Quaternion<T> &x)

    template<typename T>
    Quaternion<T> operator/(const Quaternion<T> &x, std::complex<T> &y)

    template<typename T>
    Quaternion<T> operator/(std::complex<T> &y, const Quaternion<T> &x)

    template<typename T>
    Quaternion<T> operator/(const Quaternion<T> &x, const Quaternion<T> &y)

#### Dot, cross product, commutator
    template<typename T>
    T dot(const Quaternion<T> &x, const Quaternion<T> &y)

    template<typename T>
    Quaternion<T> cross(const Quaternion<T> &x, const Quaternion<T> &y)

    template<typename T>
    Quaternion<T> commutator(const Quaternion<T> &x, const Quaternion<T> &y)

#### Transcendentals
    template<typename T>
    Quaternion<T> exp(const Quaternion<T> &x)

    template<typename T>
    Quaternion<T> log(const Quaternion<T> &x)

    template<typename T>
    Quaternion<T> pow2(const Quaternion<T> &x)

    template<typename T>
    Quaternion<T> pow3(const Quaternion<T> &x)

    template<typename T>
    Quaternion<T> pow4(const Quaternion<T> &x)

    template<typename T>
    Quaternion<T> pow(const Quaternion<T> &x, int expt)

    template<typename T>
    Quaternion<T> pow(const Quaternion<T> &x, T a)

    template<typename T>
    Quaternion<T> pow(const Quaternion<T> &x, const Quaternion<T> &a)

    template<typename T>
    Quaternion<T> cos(const Quaternion<T> &x)

    template<typename T>
    Quaternion<T> sin(const Quaternion<T> &x)

    template<typename T>
    Quaternion<T> tan(const Quaternion<T> &x)

    template<typename T>
    Quaternion<T> cosh(const Quaternion<T> &x)

    template<typename T>
    Quaternion<T> sinh(const Quaternion<T> &x)

    template<typename T>
    Quaternion<T> tanh(const Quaternion<T> &x)

#### Algorithms
    template<typename T , typename K>
    Quaternion<T> axby(K k1, const Quaternion<T> &x, K k2, const Quaternion<T> &y)