#include "Quaternion.h"

#include <set>
#include <random>
#include <chrono>
#include <array>
#include <iomanip>

#include <boost/math/quaternion.hpp>
#include <boost/rational.hpp>

using namespace std;
using namespace boost;
using namespace boost::math;

//----------------------------------------------------------------------------------------------------------------------
// Test utilities
//----------------------------------------------------------------------------------------------------------------------
/**
 * Prints out arrays to streams
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
 * Prints out a vector of elements of type T
 */
template <typename T>
inline ostream& operator<<(ostream& out, const vector<T>& x) {
  for (T v : x)
    out << v << " ";
  return out;
}

/**
 * Simple test function for the unit tests
 */
void EXPECT(bool v) {
  if (!v) {
    cerr << "Test failed" << endl;
    exit(EXIT_FAILURE);
  }
}

/**
 * Compare a boost quaternion to Quaternion, relatively, with epsilon.
 */
template <typename T, typename T1>
inline bool equals(const Quaternion<T>& us, const quaternion<T>& them, T1 eps) {
  T teps = static_cast<T>(eps);
  return std::abs((us.a() - them.R_component_1())) < teps
         && std::abs((us.b() - them.R_component_2())) < teps
         && std::abs((us.c() - them.R_component_3())) < teps
         && std::abs((us.d() - them.R_component_4())) < teps;
}

/**
 * A random number generator.
 */
std::default_random_engine generator;
std::uniform_real_distribution<float> distribution(0.0,1.0);
auto rng = bind(distribution, generator);

/**
 * A function to generate random quaternions.
 */
template <typename T, typename G>
inline Quaternion<T> random_quaternion(G& g) {
  return Quaternion<T>(g(), g(), g(), g());
}

/**
 * This method useful in unit tests, to compare against boost, which has
 * a (supposedly) well tested quaternion library.
 * NOTE: hard to compare large values, requires floating point comparison.
 */
bool operator==(const Quaternion<float>& x, const quaternion<float>& boost_y) {
  return x.a() == boost_y.R_component_1()
         && x.b() == boost_y.R_component_2()
         && x.c() == boost_y.R_component_3()
         && x.d() == boost_y.R_component_4();
}

bool operator==(const quaternion<float>& boost_y, const Quaternion<float>& x) {
  return x == boost_y;
}

bool operator==(const Quaternion<double>& x, const quaternion<double>& boost_y) {
  return x.a() == boost_y.R_component_1()
         && x.b() == boost_y.R_component_2()
         && x.c() == boost_y.R_component_3()
         && x.d() == boost_y.R_component_4();
}

bool operator==(const quaternion<double>& boost_y, const Quaternion<double>& x) {
  return x == boost_y;
}

bool operator==(const Quaternion<long double>& x, const quaternion<long double>& boost_y) {
  return x.a() == boost_y.R_component_1()
         && x.b() == boost_y.R_component_2()
         && x.c() == boost_y.R_component_3()
         && x.d() == boost_y.R_component_4();
}

bool operator==(const quaternion<long double>& boost_y, const Quaternion<long double>& x) {
  return x == boost_y;
}

//----------------------------------------------------------------------------------------------------------------------
// Unit tests
//----------------------------------------------------------------------------------------------------------------------
void test_constructors() {
  cout << "Testing constructors" << endl;
  {
    Qf x;
    EXPECT(x.a() == 0 && x.b() == 0 && x.c() == 0 && x.d() == 0);
    EXPECT(x.real() == 0);
    EXPECT(x.unreal() == 0);
  }

  {
    Qf x((float)0, (float)0, (float)0, (float)0);
    EXPECT(x.a() == 0 && x.b() == 0 && x.c() == 0 && x.d() == 0);
  }

  {
    Qf x((float)1);
    EXPECT(x.a() == 1 && x.b() == 0 && x.c() == 0 && x.d() == 0);
  }

  {
    Qf x((float)1, (float)2);
    EXPECT(x.a() == 1 && x.b() == 2 && x.c() == 0 && x.d() == 0);
  }

  {
    Qf x((float)1, (float)2, (float)3);
    EXPECT(x.a() == 1 && x.b() == 2 && x.c() == 3 && x.d() == 0);
  }

  {
    Qf x((float)1, (float)2, (float)3, (float)4);
    EXPECT(x.a() == 1 && x.b() == 2 && x.c() == 3 && x.d() == 4);

  }

  {
    Qf x((int)1,(int)2);
    EXPECT(x.a() == 1 && x.b() == 2 && x.c() == 0 && x.d() == 0);
  }

  {
    Qf x((int)1,(int)2,(int)3,(int)4);
    EXPECT(x.a() == 1 && x.b() == 2 && x.c() == 3 && x.d() == 4);
  }

  {
    Qf x(complex<float>(1,2));
    EXPECT(x.a() == 1 && x.b() == 2 && x.c() == 0 && x.d() == 0);
  }

  {
    Qf x(complex<float>(1,2), complex<float>(3,4));
    EXPECT(x.a() == 1 && x.b() == 2 && x.c() == 3 && x.d() == 4);
    EXPECT(x.c1() == complex<float>(1,2));
    EXPECT(x.c2() == complex<float>(3,4));
  }

  {
    Qd x(complex<float>(1,2), complex<float>(3,4));
    EXPECT(x.a() == 1 && x.b() == 2 && x.c() == 3 && x.d() == 4);
  }

  {
    Qd y(complex<float>(1,2), complex<float>(3,4));
    Qd x(y);
    EXPECT(x.a() == 1 && x.b() == 2 && x.c() == 3 && x.d() == 4);
  }

  {
    Qf y((int)1,(int)2,(int)3,(int)4);
    Qd x(y);
    EXPECT(x.a() == 1 && x.b() == 2 && x.c() == 3 && x.d() == 4);
  }

  {
    Qd y(complex<float>(1,2), complex<float>(3,4));
    Qd x = y;
    EXPECT(x.a() == 1 && x.b() == 2 && x.c() == 3 && x.d() == 4);
  }

  {
    Qf y((int)1,(int)2,(int)3,(int)4);
    Qd x = y;
    EXPECT(x.a() == 1 && x.b() == 2 && x.c() == 3 && x.d() == 4);
    EXPECT(x.c1() == complex<double>(1,2));
    EXPECT(x.c2() == complex<double>(3,4));
  }

  {
    int d[4] = {1,2,3,4};
    Qf x(d);
    EXPECT(x.a() == 1 && x.b() == 2 && x.c() == 3 && x.d() == 4);
  }

  {
    vector<int> d{1,2,3,4};
    Qf x(d.begin());
    EXPECT(x.a() == 1 && x.b() == 2 && x.c() == 3 && x.d() == 4);
  }

  {
    set<int> d{1,2,3,4};
    Qf x(d.begin());
    EXPECT(x.a() == 1 && x.b() == 2 && x.c() == 3 && x.d() == 4);
  }

  {
    Qf x{1,0,2,0};
    EXPECT(x.a() == 1 && x.b() == 0 && x.c() == 2 && x.d() == 0);
  }
}

void test_trigonometric_constructors() {
  cout << "Testing trigonometric constructors" << endl;
  {
    Qf x; x.spherical(0,0,0,0);
    EXPECT(x.a() == 0 && x.b() == 0 && x.c() == 0 && x.d() == 0);
  }

  {
    Qd x = Qf::spherical((float)10, (float)3.1415, float(3.1415/2), float(3.1415/4));
    cout << x << endl;
    //EXPECT(x.a() == 0 && x.b() == 0 && x.c() == 0 && x.d() == 0);
  }

  {
    Qd x = Qf::semipolar((float)10, (float)3.1415, float(3.1415/2), float(3.1415/4));
    cout << x << endl;
    //EXPECT(x.a() == 0 && x.b() == 0 && x.c() == 0 && x.d() == 0);
  }

  {
    Qd x = Qf::multipolar((float)10, (float)3.1415, float(3.1415/2), float(3.1415/4));
    cout << x << endl;
    //EXPECT(x.a() == 0 && x.b() == 0 && x.c() == 0 && x.d() == 0);
  }

  {
    Qd x = Qf::cylindrospherical((float)10, (float)3.1415, float(3.1415/2), float(3.1415/4));
    cout << x << endl;
    //EXPECT(x.a() == 0 && x.b() == 0 && x.c() == 0 && x.d() == 0);
  }

  {
    Qd x = Qf::cylindrical((float)10, (float)3.1415, float(3.1415/2), float(3.1415/4));
    cout << x << endl;
    //EXPECT(x.a() == 0 && x.b() == 0 && x.c() == 0 && x.d() == 0);
  }
}

void test_IJK() {
  cout << "Testing IJK" << endl;
  EXPECT(Qf_0 == 0);
  EXPECT(Qf_1 == 1);
  EXPECT(Qf_i * Qf_i == -1);
  EXPECT(Qf_j * Qf_j == -1);
  EXPECT(Qf_k * Qf_k == -1);
  EXPECT(Qf_i * Qf_j * Qf_k == -1);

  EXPECT(Qd_0 == 0);
  EXPECT(Qd_1 == 1);
  EXPECT(Qd_i * Qd_i == -1);
  EXPECT(Qd_j * Qd_j == -1);
  EXPECT(Qd_k * Qd_k == -1);
  EXPECT(Qd_i * Qd_j * Qd_k == -1);

  EXPECT(Qld_0 == 0);
  EXPECT(Qld_1 == 1);
  EXPECT(Qld_i * Qld_i == -1);
  EXPECT(Qld_j * Qld_j == -1);
  EXPECT(Qld_k * Qld_k == -1);
  EXPECT(Qld_i * Qld_j * Qld_k == -1);
}

void test_accessors() {
  cout << "Testing accessors" << endl;
  {
    Qf x(1, 2, 3, 4);
    EXPECT(x.real() == 1);
    EXPECT(x.unreal() != 0); // TODO: check what is happening here exactly
    EXPECT(x.unreal() == Qf(0, 2, 3, 4));
    EXPECT(x.conjugate() == Qf(1, -2, -3, -4));
    EXPECT(x.conjugate().conjugate() == x);
    EXPECT(x.conjugate() + x == Qf(2, 0, 0, 0));
    EXPECT(x - x.conjugate() == Qf(0, 4, 6, 8));
    EXPECT(!x.is_real());
    EXPECT(!x.is_complex());
    EXPECT(!x.is_unreal());
    EXPECT(!x.is_unit());
    EXPECT(-x == Qf(-1,-2,-3,-4));
  }

  {
    Qf x(3.14);
    float a = x;
    EXPECT(a == 3.14f);
  }

  {
    Qf x(3.14, 2.71);
    complex<float> c = x;
    EXPECT(c.real() == 3.14f && c.imag() == 2.71f);
  }

  {
    Qf x(1,2,3,4);
    array<float,4> a = x;
    array<float,4> r{{1,2,3,4}};
    EXPECT(a == r);
  }
}

void test_to_matrix() {
  cout << "Testing to_matrix" << endl;
  typedef std::complex<double> CD;
  Qd x(1,2,3,4);
  Qd::matrix_representation r;
  r[0] = {{CD(1,2),CD(3,4)}};
  r[1] = {{CD(-3,4),CD(1,-2)}};
  EXPECT(x.to_matrix_representation() == r);
}

void test_to_polar_representation() {
  cout << "Testing polar representation" << endl;
  Qd x(1);
  cout << Qd(1,0,0,0).to_polar_representation() << endl;
  cout << Qd(0,1,0,0).to_polar_representation() << endl;
  cout << Qd(0,0,1,0).to_polar_representation() << endl;
  cout << Qd(0,0,0,1).to_polar_representation() << endl;
}

void test_norms() {
  cout << "Testing norms" << endl;
  {
    Qf x{1,2,3,4};
    EXPECT(norm2(x) == 1+4+9+16);
    EXPECT(std::abs(norm(x) - std::sqrt(1+4+9+16)) < 1e-6);
    EXPECT(normalize(x).is_unit(1e-6));
    Qf::scalar_zero_threshold = 1e-6f;
    EXPECT(normalize(x) == x/std::sqrt(1+4+9+16));
  }
}

void test_equality() {
  cout << "Testing equality" << endl;
  EXPECT(Qf(1,2,3,4) == Qf(1, 2, 3, 4));
  EXPECT(Qf(1,2,3,4) != Qf(4,3,2,1));
  //TODO: refine for precision
}

void test_unaries() {
  cout << "Testing unaries" << endl;
  EXPECT(-Qf(1,2,3,4) == Qf(-1,-2,-3,-4));
  EXPECT(- -Qf(1,2,3,4) == Qf(1,2,3,4));
  EXPECT(- - -Qf(1,2,3,4) == Qf(-1,-2,-3,-4));
  EXPECT(+Qf(1,2,3,4) == Qf(1,2,3,4));
  EXPECT(+ +Qf(1,2,3,4) == Qf(1,2,3,4));
  EXPECT(+ + +Qf(1,2,3,4) == Qf(1,2,3,4));
}

void test_unary_w_scalar() {
  cout << "Testing unary operations with scalar" << endl;
  {
    Qd x(1,2,3,4);
    x += (long double) 3;
    EXPECT(x == Qd(4,2,3,4));
  }

  {
    Qd x(1,2,3,4);
    x -= 3;
    EXPECT(x == Qd(-2,2,3,4));
  }

  {
    Qd x(1,2,3,4);
    x *= -3.5;
    EXPECT(x == Qd(-3.5,-3.5*2,-3.5*3,-3.5*4));
  }

  {
    Qd x(1,2,3,4);
    x /= 2;
    EXPECT(x == Qd(1.0/2,2.0/2,3.0/2,4.0/2));
  }

  {
    Qd x(1,2);
    x /= 3;
    EXPECT(x == complex<double>(1.0/3, 2.0/3));
  }
}

// TODO: test unaries with complex

void test_pow2() {
  cout << "Testing pow2" << endl;
  EXPECT(pow2(Qf_0) == 0);
  EXPECT(pow2(Qf_1) == 1);
  EXPECT(pow2(Qf_1) == Qf_1 * Qf_1);
  EXPECT(pow2(Qf_i) == -1);
  EXPECT(pow2(Qf_j) == -1);
  EXPECT(pow2(Qf_k) == -1);
  EXPECT(pow2(Qf(1,2,3,4)) == Qf(1,2,3,4) * Qf(1,2,3,4));
  EXPECT(pow2(Qf(1,2)) == complex<float>(1,2) * complex<float>(1,2));
}

void test_pow3() {
  cout << "Testing pow3" << endl;
  EXPECT(pow3(Qd_0) == 0);
  EXPECT(pow3(Qd_1) == 1);
  EXPECT(pow3(Qd_i) == -Qd_i);
  EXPECT(pow3(Qd_j) == -Qd_j);
  EXPECT(pow3(Qd_k) == -Qd_k);

  EXPECT(pow3(Qld_0) == 0);
  EXPECT(pow3(Qld_1) == 1);
  EXPECT(pow3(Qld_i) == -Qld_i);
  EXPECT(pow3(Qld_j) == -Qld_j);
  EXPECT(pow3(Qld_k) == -Qld_k);
  EXPECT(pow3(Qf(1,2,3,4)) == Qf(1,2,3,4) * Qf(1,2,3,4) * Qf(1,2,3,4));
}

void test_pow() {
  cout << "Testing pow" << endl;

  for (size_t i = 0; i < 1000; ++i) {
    int n = (int) random() % 20;
    Qld x(rand()%5,rand()%5,rand()%5,rand()%5);
    Qld y = Qld_1;
    for (int j = 0; j < n; ++j)
      y *= x;
    EXPECT(norm2(y - pow(x,n)) < 1e-10);
  }

  for (size_t i = 0; i < 1000; ++i) {
    double n = double(random() % 20)/(1 + double(random() % 10));
    Qd x(rand()%5,rand()%5,rand()%5,rand()%5);
    EXPECT(std::abs(norm(pow(x,n)) - norm(exp(n * log(x)))) < 1e-6);
  }
}

void test_q_pow() {
  cout << "Testing q pow" << endl;

  for (size_t i = 0; i < 10; ++i) {
    Qld x(rand()%5,rand()%5,rand()%5,rand()%5);
    Qld y(rand()%5,rand()%5,rand()%5,rand()%5);
    EXPECT(std::abs(norm(pow(x,y)) - norm(exp(y * log(x)))) < 1e-6);
  }
}

void test_log() {
  cout << "Testing log" << endl;

  { // Make sure it works for reals
    Qld x(2,0,0,0);
    EXPECT(std::abs(log(x).real() - std::log(2)) < 1e-10);
  }

  size_t N = 10;
  Qld::scalar_zero_threshold = 1e-6;
  for (size_t i = 0; i < N; ++i) {
    Qld x(rand() % 5, rand() % 5, rand() % 5, rand() % 5);
    Qld::scalar_zero_threshold = 1e-12;
    EXPECT(x == exp(log(x))); // but not the other way around!
  }
}

void test_exp() {
  cout << "Testing exp" << endl;

  { // Make sure it works for reals
    Qld x(2,0,0,0);
    EXPECT(std::abs(exp(x).real() - std::exp(2)) < 1e-10);
  }

  {
    size_t N = 1000;
    for (size_t i = 0; i < N; ++i) {
      Qld x(rand() % 5, rand() % 5, rand() % 5, rand() % 5);
      quaternion<long double> bx(x.a(), x.b(), x.c(), x.d());
      Qld y = exp(x);
      quaternion<long double> by = exp(bx);
      EXPECT(equals(y, by, 1e-6));
    }
  }
}

void test_dot() {
  cout << "Testing dot product" << endl;
  EXPECT(dot(Qf(1, 2, 3, 4), Qf(2, 2, 2, 2)) == 20);
  EXPECT(dot(Qf(-1, 2, -3, 4), Qf(-1, 2, -3, 4)) == norm2(Qf(1,2,3,4)));
  // TODO: verify this works in general
  Qf a(1,2,3,4), b(5,6,7,8);
  Qf d1 = dot(a,b);
  Qf d2 = .5*(conjugate(b) * a + conjugate(a) * b);
  EXPECT(d1 == d2);
}

void test_cross() {
  cout << "Testing cross product" << endl;
  // TODO: verify
  EXPECT(Qf(0,-2,4,-2) == cross(Qf(1, 2, 3, 4), Qf(2, 2, 2, 2)));
  // TODO: verify this works in general
  Qf a(0,2,3,4), b(0,6,7,8);
  Qf p1 = cross(a,b);
  Qf p2 = .5*(a * b - conjugate(b) * conjugate(a));
  cout << p1 << endl;
  cout << p2 << endl;
  EXPECT(p1 == p2);
}

void test_commutator() {
  cout << "Testing commutator" << endl;
  Qf c = commutator(Qf(1,2,3,4), Qf(2,2,2,2));
  EXPECT(c == Qf(0,-4,8,-4));
  EXPECT(c == Qf(1,2,3,4) * Qf(2,2,2,2) - Qf(2,2,2,2) * Qf(1,2,3,4));
  // TODO: verify this works in general
  Qf a(1,2,3,4), b(5,6,7,8);
  Qf c1 = commutator(a,b);
  Qf c2 = 2 * cross(a,b);
  EXPECT(c1 == c2);
}

void test_io_eps() {
  cout << "Testing io/eps" << endl;
  Qd x(1e-1, 1e-2, 1e-3, 1e-4);
  cout << set_eps(1e-6) << x << endl;
  cout << set_eps(1e-4) << x << endl;
  cout << set_eps(1e-3) << x << endl;
  cout << set_eps(1e-2) << x << endl;
  cout << set_eps(1e-1) << x << endl;
  cout << set_eps(1) << x << endl;
}

void test_io_style() {
  cout << "Testing io/style" << endl;
  Qd x(1,2,3,4);
  cout << set_style_nice<double>() << x << endl;
  cout << set_style_compact<double>() << x << endl;
}

void test_stl() {
  cout << "Testing STL" << endl;
  vector<Qf> qs{{1,2,3,4},{5,6,7,8},{1,3,5,7},{2,4,6,8}};
  auto v = accumulate(qs.begin(), qs.end(), Qf_1, multiplies<Qf>());
  EXPECT(v == Qf_1 * Qf(1,2,3,4) * Qf(5,6,7,8) * Qf(1,3,5,7) * Qf(2,4,6,8));
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
    EXPECT(z.a() == Rl(-554,315));
    EXPECT(z.b() == Rl(481,540));
    EXPECT(z.c() == Rl(641,630));
    EXPECT(z.d() == Rl(253,252));
  }
}

/**
 * Between standard C++, boost and vectorclass which uses intrinsics, no difference.
 * The compiler is probably optimizing well enough, or the intrinsics are not used properly.
 * Whoever is in first position (boot or Quaternion) in this micro-benchmark, is twice as
 * fast as the second.... But only on the Mac/clang or gcc. That doesn't happen on Linux/gcc.
 */
void test_multiplication_speed() {
  cout << "Testing multiplication speed" << endl;
  size_t N = 100000;

  Quaternion<float> q1 = random_quaternion<float>(rng), q2 = random_quaternion<float>(rng);

  {
    float certificate = 0.0;
    auto start = std::chrono::system_clock::now();
    for (size_t i = 0; i < N; ++i) {
      Quaternion<float> r = q1 * (q2 + (float)i);
      certificate += r.a() + r.b() + r.c() + r.d();
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::nanoseconds diff = end - start;
    cout << "Quaternion: " << (diff.count() / N) << "ns" << endl;
    cout << "Certificate=" << certificate << endl;
  }

  { // With Boost
    quaternion<float> a(q1.a(),q1.b(),q1.c(),q1.d()), b(q2.a(),q2.b(),q2.c(),q2.d());
    float certificate = 0.0;
    auto start = std::chrono::system_clock::now();
    for (size_t i = 0; i < N; ++i) {
      quaternion<float> r = a * (b + (float)i);
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
  Quaternion<float> q1 = random_quaternion<float>(rng);


  { // With Boost
    quaternion<float> a(q1.a(),q1.b(),q1.c(),q1.d());
    float certificate = 0.0;
    auto start = std::chrono::system_clock::now();
    for (size_t i = 0; i < N; ++i) {
      quaternion<float> r = pow(a, 15);
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
      Quaternion<float> r = pow(q1, 15);
      certificate += r.a() + r.b() + r.c() + r.d();
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::nanoseconds diff = end - start;
    cout << "Quaternion: " << (diff.count() / N) << "ns" << endl;
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
    cout << "Quaternion: " << (diff.count() / N) << "ns" << endl;
    cout << "Certificate=" << certificate << endl;
  }
}

void test_axby_speed() {
  cout << "Testing axby speed" << endl;
  size_t N = 100000;

  Quaternion<float> q1 = random_quaternion<float>(rng), q2 = random_quaternion<float>(rng);

  { // With Boost
    quaternion<float> a(q1.a(),q1.b(),q1.c(),q1.d()), b(q2.a(),q2.b(),q2.c(),q2.d());
    float certificate = 0.0;
    auto start = std::chrono::system_clock::now();
    for (size_t i = 0; i < N; ++i) {
      quaternion<float> r = ((float) i) * a + ((float) i+1) * b;
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
      Quaternion<float> r = axby(i, q1, i+1, q2);
      certificate += r.a() + r.b() + r.c() + r.d();
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::nanoseconds diff = end - start;
    cout << "Quaternion: " << (diff.count() / N) << "ns" << endl;
    cout << "Certificate=" << certificate << endl;
  }
}

int main() {
  test_constructors();
  test_trigonometric_constructors();
  test_IJK();
  test_accessors();
  test_to_matrix();
  test_to_polar_representation();
  test_norms();
  test_equality();
  test_unaries();
  test_unary_w_scalar();
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
  test_io_eps();
  test_io_style();

//  test_exp_speed();
//  test_multiplication_speed();
//  test_pow_speed();
//  test_axby_speed();

  return 0;
}
