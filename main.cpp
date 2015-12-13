#include <iostream>
#include <cmath>
#include <array>

using namespace std;

template <typename T =double> // assert operations for numeric
class Quaternion {
public:
  Quaternion(T a =0, T b =0, T c =0, T d =0)
      : _a(a), _b(b), _c(c), _d(d)
  {}

  Quaternion<T> operator-() const {
    return Quaternion<T>(-a(), -b(), -c(), -d());
  }

  T a() const { return _a; }
  T b() const { return _b; }
  T c() const { return _c; }
  T d() const { return _d; }

private:
  T _a,_b,_c,_d;
};

template <typename T=double, int i = 0, int j = 1, int k = 2, int l = 3>
class CQuaternion {
public:
  CQuaternion(T a =0, T b =0, T c =0, T d =0)
  {
    if (i >= 0) _val[i] = a;
    if (j >= 0) _val[j] = b;
    if (k >= 0) _val[k] = c;
    if (l >= 0) _val[l] = d;
  }

  T a() const { return i == -1 ? 0 : _val[i]; }
  T b() const { return j == -1 ? 0 : _val[j]; }
  T c() const { return k == -1 ? 0 : _val[k]; }
  T d() const { return l == -1 ? 0 : _val[l]; }

private:
  array<T, (i >= 0) + (j >= 0) + (k >=0) + (l >= 0)> _val;
};

typedef Quaternion<> Q;

template <typename T=double>
const Quaternion<T> Q_zero = Quaternion<T>(0,0,0,0);

template <typename T=double>
const Quaternion<T> Q_1 = Quaternion<T>(1,0,0,0);

template <typename T=double>
const Quaternion<T> Q_i = Quaternion<T>(0,1,0,0);

template <typename T=double>
const Quaternion<T> Q_j = Quaternion<T>(0,0,1,0);

template <typename T=double>
const Quaternion<T> Q_k = Quaternion<T>(1,0,0,0);

template <typename T>
inline ostream& operator<<(ostream& out, const Quaternion<T>& q) {
  if (q == Q_zero<T>)
    return out << 0;
  if (q == Q_1<T>)
    return out << 1;
  if (q == -Q_1<T>)
    return out << -1;
  if (q == Q_i<T>)
    return out << "i";
  if (q == -Q_i<T>)
    return out << "-i";
  if (q == Q_j<T>)
    return out << "j";
  if (q == -Q_j<T>)
    return out << "-j";
  if (q == Q_k<T>)
    return out << "k";
  if (q == -Q_k<T>)
    return out << "-k";
  auto s = [](T x) { return x < 0 ? "" : "+"; };
  if (std::abs(q.a()) > numeric_limits<T>::epsilon())
    out << q.a();
  if (std::abs(q.b()) > numeric_limits<T>::epsilon())
    out << s(q.b()) << q.b() << "i";
  if (std::abs(q.c()) > numeric_limits<T>::epsilon())
    out << s(q.c()) << q.c() << "j";
  if (std::abs(q.d()) > numeric_limits<T>::epsilon())
    out << s(q.d()) << q.d() << "k";
  return out;
}

template <typename T, typename T1>
inline Quaternion<T> operator*(T1 k, const Quaternion<T>& x) {
  if (std::abs(k) < numeric_limits<T>::epsilon())
    return Q_zero<T>;
  if (std::abs(k - 1) < numeric_limits<T>::epsilon())
    return x;
  return Quaternion<T>(k*x.a(), k*x.b(), k*x.c(), k*x.d());
}

template <typename T, typename T1>
inline Quaternion<T> operator*(const Quaternion<T>& x, T1 k) {
  if (std::abs(k) < numeric_limits<T>::epsilon())
    return Q_zero<T>;
  if (std::abs(k - 1) < numeric_limits<T>::epsilon())
    return x;
  return k * x;
}

template <typename T, typename T1>
inline Quaternion<T> operator/(const Quaternion<T>& x, T1 k) {
  if (std::abs(k - 1) < numeric_limits<T>::epsilon())
    return x;
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
inline bool is_unit_q(const Quaternion<T>& x) {
  return std::abs(norm(x) - 1) < numeric_limits<T>::epsilon();
}

template <typename T>
inline bool operator==(const Quaternion<T>& x, const Quaternion<T>& y) {
  return norm2(x - y) < numeric_limits<T>::epsilon();
}

template <typename T>
inline Quaternion<T> operator+(const Quaternion<T>& x, const Quaternion<T>& y) {
  return Quaternion<T>(x.a()+y.a(),x.b()+y.b(),x.c()+y.c(),x.d()+y.d());
}

template <typename T>
inline Quaternion<T> operator-(const Quaternion<T>& x, const Quaternion<T>& y) {
  return Quaternion<T>(x.a()-y.a(),x.b()-y.b(),x.c()-y.c(),x.d()-y.d());
}

template <typename T>
inline Quaternion<T> operator*(const Quaternion<T>& x, const Quaternion<T>& y) {
  return Quaternion<T>(x.a()*y.a() - x.b()*y.b() - x.c()*y.c() - x.d()*y.d(),
                       x.a()*y.b() + x.b()*y.a() + x.c()*y.d() - x.d()*y.c(),
                       x.a()*y.c() - x.b()*y.d() + x.c()*y.a() + x.d()*y.b(),
                       x.a()*y.d() + x.b()*y.c() - x.c()*y.b() + x.d()*y.a());
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


template <typename T, int i1, int j1, int k1, int l1,
    int i2, int j2, int k2, int l2,
    int i3, int j3, int k3, int l3>
inline CQuaternion<T,i3,j3,k3,l3>
operator+(const CQuaternion<T,i1,j1,k2,l1>& x, const CQuaternion<T,i2,j2,k2,l2>& y) {
  T a = (i1 == -1 || i2 == -1) ? 0 : x.a() + y.a();
  T b = (j1 == -1 || j2 == -1) ? 0 : x.b() * y.b();
  T c = (k1 == -1 || k2 == -1) ? 0 : x.c() * y.c();
  T d = (l1 == -1 || l2 == -1) ? 0 : x.d() * y.d();
  if (std::abs(a) < numeric_limits<T>::epsilon
      && std::abs(b) < numeric_limits<T>::epsilon
      && std::abs(c) < numeric_limits<T>::epsilon
      && std::abs(d) < numeric_limits<T>::epsilon)
    return CQuaternion<T,-1,-1,-1>();
}

int main() {
  cout << numeric_limits<double>::epsilon() << endl;
  cout << Q(1,2,3,4) + Q(4,3,2,1) << endl;
  cout << Q(1,2,3,4) * Q(4,3,2,1) << endl;
  cout << Q(4,3,2,1) * Q(1,2,3,4) << endl;
  cout << Q(4,3,2,1) * 1.5 << endl;
  cout << Q(4,3,2,1) * 0.0 << endl;
  cout << (Q(4,3,2,1) == Q(1,2,3,4)) << endl;
  cout << (Q(4,3,2,1) == Q(4,3,2,1)) << endl;
  cout << (Q(4,3,2,1) + -Q(4,3,2,1)) << endl;
  cout << Q_i<double> * Q_i<double> << endl;
  cout << conjugate(Q(1,2,3,4)) << endl;
  cout << norm(Q(1,2,3,4)) << " " << sqrt(30) << endl;
  cout << Q(0,1,0,0) << endl;
  cout << is_unit_q(Q(0,1,0,0)) << endl;
  cout << inverse(Q(0,1,0,0)) << endl;
  cout << (inverse(Q(0,1,0,0)) == Q_i<>) << endl;
  cout << (Q(1,2,3,4) / Q(2,2,2,2)) << endl;
  cout << (Q(1,2,3,4) * Q(2,2,2,2)) << endl;
  cout << .5*(Q(1,2,3,4) * conjugate(Q(2,2,2,2)) + Q(2,2,2,2) * conjugate(Q(1,2,3,4))) << endl;
  cout << dot(Q(1,2,3,4), Q(2,2,2,2)) << endl;
  cout << cross(Q(1,2,3,4), Q(2,2,2,2)) << endl;
  cout << commutator(Q(1,2,3,4), Q(2,2,2,2)) << endl;
  cout << sizeof(CQuaternion<double,0,1>(2,3)) << endl;
  cout << sizeof(CQuaternion<double,0,1,-1,-1>(2)) << endl;
  return 0;
}