//
// Created by Frank Astier on 2015/12/13.
//

#ifndef QUATERNIONS_INTRINSICS_H
#define QUATERNIONS_INTRINSICS_H

#include <xmmintrin.h>
#include <pmmintrin.h>

/**
 * This code is from StackOverflow. It it probably wrong, the description
 * of the product of two quaternions seems to be incorrect.
 * In any case, this happens to not be faster than a naive implementation.
 */
inline __m128 _mm_cross4_ps(const __m128 xyzw, const __m128 abcd)
{
  /* The product of two quaternions is:                                 */
  /* (X,Y,Z,W) = (xd+yc-zb+wa, -xc+yd+za+wb, xb-ya+zd+wc, -xa-yb-zc+wd) */

  __m128 wzyx = _mm_shuffle_ps(xyzw, xyzw, _MM_SHUFFLE(0,1,2,3));
  __m128 baba = _mm_shuffle_ps(abcd, abcd, _MM_SHUFFLE(0,1,0,1));
  __m128 dcdc = _mm_shuffle_ps(abcd, abcd, _MM_SHUFFLE(2,3,2,3));

  /* variable names below are for parts of componens of result (X,Y,Z,W) */
  /* nX stands for -X and similarly for the other components             */

  /* znxwy  = (xb - ya, zb - wa, wd - zc, yd - xc) */
  __m128 ZnXWY = _mm_hsub_ps(_mm_mul_ps(xyzw, baba), _mm_mul_ps(wzyx, dcdc));

  /* xzynw  = (xd + yc, zd + wc, wb + za, yb + xa) */
  __m128 XZYnW = _mm_hadd_ps(_mm_mul_ps(xyzw, dcdc), _mm_mul_ps(wzyx, baba));

  /* _mm_shuffle_ps(XZYnW, ZnXWY, _MM_SHUFFLE(3,2,1,0)) */
  /*      = (xd + yc, zd + wc, wd - zc, yd - xc)        */
  /* _mm_shuffle_ps(ZnXWY, XZYnW, _MM_SHUFFLE(2,3,0,1)) */
  /*      = (zb - wa, xb - ya, yb + xa, wb + za)        */

  /* _mm_addsub_ps adds elements 1 and 3 and subtracts elements 0 and 2, so we get: */
  /* _mm_addsub_ps(*, *) = (xd+yc-zb+wa, xb-ya+zd+wc, wd-zc+yb+xa, yd-xc+wb+za)     */

  __m128 XZWY = _mm_addsub_ps(_mm_shuffle_ps(XZYnW, ZnXWY, _MM_SHUFFLE(3,2,1,0)),
                              _mm_shuffle_ps(ZnXWY, XZYnW, _MM_SHUFFLE(2,3,0,1)));

  /* now we only need to shuffle the components in place and return the result      */
  return XZWY; //_mm_shuffle_ps(XZWY, XZWY, _MM_SHUFFLE(2,1,3,0));

  /* operations: 6 shuffles, 4 multiplications, 3 compound additions/subtractions   */
}

// Not faster than regular C++ for now
//template <>
//inline Quaternion<float> operator*(const Quaternion<float>& x, const Quaternion<float>& y) {
//  __m128 r = _mm_cross4_ps(*reinterpret_cast<const __m128*>(&x), *reinterpret_cast<const __m128*>(&y));
//  float* p = reinterpret_cast<float*>(&r);
//  return Quaternion<float>(*(p+2), *(p+1), *(p+3), *p);
//}

#endif //QUATERNIONS_INTRINSICS_H
