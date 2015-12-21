//
// Created by Frank Astier on 2015/12/18.
//

#ifndef QUATERNIONS_UTILS_H
#define QUATERNIONS_UTILS_H

#include <cmath>
#include <iterator>

#define IS_ITERATOR(X) typename std::enable_if<!std::is_same<typename std::iterator_traits<X>::value_type, void>::value>::type* =nullptr
#define IS_NOT_ITERATOR(X) typename std::enable_if<std::is_same<typename std::iterator_traits<X>::value_type, void>::value>::type* =nullptr

#define IS_NOT_COMPLEX(X) typename std::enable_if<std::is_same<

/**
 * Utility function to work with numbers approximately equal to zero.
 * If eps == 0, does a "hard" comparison to 0.
 * Otherwise, uses a ball of radius eps around 0. If the scalar is inside
 * that ball, it is equivalent to 0.
 * TODO: refine for large floating point numbers
 */
template<typename T>
inline bool is_scalar_zero(T x, T eps = 0) {
  if (eps > 0)
    return std::fabs(x) < eps;
  else
    return x == 0;
}

#endif //QUATERNIONS_UTILS_H
