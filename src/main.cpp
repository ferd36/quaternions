#include "Quaternion.h"

#include <random>
#include <chrono>
#include <boost/math/quaternion.hpp>

using namespace std;
using namespace boost::math;


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

void test_IJK() {
  cout << "Testing IJK" << endl;
  assert(Qf_0 == 0);
  assert(Qf_1 == 1);
  assert(Qf_i * Qf_i == -1);
  assert(Qf_j * Qf_j == -1);
  assert(Qf_k * Qf_k == -1);
  assert(Qf_i * Qf_j * Qf_k == -1);
}

void test_pow2() {
  cout << "Testing pow2" << endl;
  assert(pow2(Qf_0) == 0);
  assert(pow2(Qf_1) == 1);
  assert(pow2(Qf_1) == Qf_1 * Qf_1);
  assert(pow2(Qf_i) == -1);
  assert(pow2(Qf_j) == -1);
  assert(pow2(Qf_k) == -1);
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
}

void test_pow() {
  cout << "Testing pow" << endl;
  for (size_t i = 0; i < 100; ++i) {
    int n = (int) random() % 20;
    Qld x(rand()%5,rand()%5,rand()%5,rand()%5);
    Qld y = Qld_1;
    for (int j = 0; j < n; ++j)
      y *= x;
    assert(pow(x, n) == y);
    // Test against boost
    quaternion<long double> bx(x.a(),x.b(),x.c(),x.d());
    assert(pow(bx, n) == y);
  }
}

void test_to_matrix() {
  Qld x(rand()%5,rand()%5,rand()%5,rand()%5);
  auto X = x.to_matrix();
  //cout << X << endl;
}

void test_exp() {
  cout << "Testing exp" << endl;
  Qld x(rand()%5,rand()%5,rand()%5,rand()%5);
  quaternion<long double> bx(x.a(),x.b(),x.c(),x.d());
  cout << exp(x) << endl;
  cout << exp(bx) << endl;
  assert(exp(x) == exp(bx));
}

void test_addition() {
  cout << Qf(1, 2, 3, 4) + Qf(4, 3, 2, 1) << endl;
}

void test_multiplication() {

  cout << Qf(1, 2, 3, 4) * Qf(4, 3, 2, 1) << endl;
  cout << Qf(4, 3, 2, 1) * Qf(1, 2, 3, 4) << endl;
  cout << Qf(4, 3, 2, 1) * 1.5 << endl;
  cout << Qf(4, 3, 2, 1) * 0.0 << endl;
  cout << Qf_i * Qf_i << endl;
}

void test_equality() {
  cout << (Qf(4, 3, 2, 1) == Qf(1, 2, 3, 4)) << endl;
  cout << (Qf(4, 3, 2, 1) == Qf(4, 3, 2, 1)) << endl;
}

void test_unary_minus() {
  cout << (Qf(4, 3, 2, 1) + -Qf(4, 3, 2, 1)) << endl;
}

void test_conjugate() {
  cout << conjugate(Qf(1, 2, 3, 4)) << endl;
}

void test_norm() {
  cout << norm(Qf(1, 2, 3, 4)) << " " << sqrt(30) << endl;
}

void test_is_unit() {
  cout << is_unit(Qf(0, 1, 0, 0)) << endl;
}

void test_inverse() {
  cout << inverse(Qf(0, 1, 0, 0)) << endl;
  cout << (inverse(Qf(0, 1, 0, 0)) == Qf_i) << endl;
}

void test_division() {
  cout << (Qf(1, 2, 3, 4) / Qf(2, 2, 2, 2)) << endl;
}

void test_dot() {
  cout << dot(Qf(1, 2, 3, 4), Qf(2, 2, 2, 2)) << endl;
  cout << .5*(Qf(1,2,3,4) * conjugate(Qf(2,2,2,2)) + Qf(2,2,2,2) * conjugate(Qf(1,2,3,4))) << endl;
}

void test_cross() {
  cout << cross(Qf(1, 2, 3, 4), Qf(2, 2, 2, 2)) << endl;
}

void test_commutator() {
  cout << commutator(Qf(1,2,3,4), Qf(2,2,2,2)) << endl;
}

void test_normalize() {
  assert(norm(normalize(Qf(1,2,3,4))) == 1);
}

/**
 * Between standard C++, boost and vectorclass which uses intrinsics, no difference.
 * The compiler is probably optimizing well enough, or the intrinsics are not used properly.
 */
void test_multiplication_speed() {
  cout << "Testing multiplication speed" << endl;
  size_t N = 10000;

  Quaternion<float> q1 = random_quaternion<float>(rng), q2 = random_quaternion<float>(rng);

  { // With Boost
    quaternion<float> a(q1.a(),q1.b(),q1.c(),q1.d()), b(q2.a(),q2.b(),q2.c(),q2.d());
    float certificate = 0.0;
    auto start = std::chrono::system_clock::now();
    for (size_t i = 0; i < N; ++i) {
      quaternion<float> r = a * b;
      certificate += r.R_component_1() + r.R_component_2() + r.R_component_3() + r.R_component_4();
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<float> diff = end - start;
    cout << "Boost: " << (1e9 * diff.count() / N) << "ns" << endl;
    cout << "Certificate=" << certificate << endl;
  }

//  { // With vectorclass, which uses intrinsics - didn't turn out to be faster
//    Quaternion4f a(q1.a(),q1.b(),q1.c(),q1.d()), b(q2.a(),q2.b(),q2.c(),q2.d());
//    float certificate = 0.0;
//    auto start = std::chrono::system_clock::now();
//    for (size_t i = 0; i < N; ++i) {
//      Quaternion4f r = a * b;
//      certificate += r.extract(0) + r.extract(1) + r.extract(2) + r.extract(3);
//    }
//    auto end = std::chrono::system_clock::now();
//    std::chrono::duration<float> diff = end - start;
//    cout << "vectorclass: " << (1e9 * diff.count() / N) << "ns" << endl;
//    cout << "Certificate=" << certificate << endl;
//  }

  {
    float certificate = 0.0;
    auto start = std::chrono::system_clock::now();
    for (size_t i = 0; i < N; ++i) {
      Quaternion<float> r = q1 * q2;
      certificate += r.a() + r.b() + r.c() + r.d();
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<float> diff = end - start;
    cout << "Quaternion: " << (1e9 * diff.count() / N) << "ns" << endl;
    cout << "Certificate=" << certificate << endl;
  }
  // TODO: match all 3 results
}

void test_pow_speed() {
  cout << "Testing pow speed" << endl;
  size_t N = 10000;

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
    std::chrono::duration<float> diff = end - start;
    cout << "Boost: " << (1e9 * diff.count() / N) << "ns" << endl;
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
    std::chrono::duration<float> diff = end - start;
    cout << "Quaternion: " << (1e9 * diff.count() / N) << "ns" << endl;
    cout << "Certificate=" << certificate << endl;
  }

  // TODO: match results
}

int main() {
  test_IJK();
  test_pow2();
  test_pow3();
  test_pow();
  test_exp();
  test_multiplication_speed();
  test_pow_speed();
  test_to_matrix();
  return 0;
}