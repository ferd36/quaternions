//
// Created by Frank Astier on 2015/12/13.
//

#ifndef QUATERNIONS_CQUATERNION_H
#define QUATERNIONS_CQUATERNION_H

/**
 * This is a "compressed" quaternion class, where we try to allocate memory just
 * for the non-zero components of a quaternion.
 * However, this might not be a good idea, because:
 * 1. there might not be many zero components, unless in rare cases,
 * 2. the cost of maintaining the minimum storage throughout calculations might
 *    be too much: it requires a lot of if-statements that derail the flow of
 *    calculations ; it requires a lot of boilerplate code to do right with
 *    templates ; and zero/non-zero components can creep up at anytime in e.g.
 *    the product of 2 quaternions. In short, it's a mess that might not be
 *    worth optimizing for.
 */
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

// Example mess - incomplete.
//template <typename T, int i1, int j1, int k1, int l1,
//    int i2, int j2, int k2, int l2,
//    int i3, int j3, int k3, int l3>
//inline CQuaternion<T,i3,j3,k3,l3>
//operator+(const CQuaternion<T,i1,j1,k2,l1>& x, const CQuaternion<T,i2,j2,k2,l2>& y) {
//  T a = (i1 == -1 || i2 == -1) ? 0 : x.a() + y.a();
//  T b = (j1 == -1 || j2 == -1) ? 0 : x.b() * y.b();
//  T c = (k1 == -1 || k2 == -1) ? 0 : x.c() * y.c();
//  T d = (l1 == -1 || l2 == -1) ? 0 : x.d() * y.d();
//  if (std::abs(a) < numeric_limits<T>::epsilon
//      && std::abs(b) < numeric_limits<T>::epsilon
//      && std::abs(c) < numeric_limits<T>::epsilon
//      && std::abs(d) < numeric_limits<T>::epsilon)
//    return CQuaternion<T,-1,-1,-1>();
//}


#endif //QUATERNIONS_CQUATERNION_H
