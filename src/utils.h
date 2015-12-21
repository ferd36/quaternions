/**
 * The MIT License (MIT)
 *
 * Copyright (c) 2015 Frank Astier
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/**
 * Utilities for working with Quaternion.
 */
#ifndef QUATERNIONS_UTILS_H
#define QUATERNIONS_UTILS_H

#include <cmath>
#include <iterator>

/**
 * A couple of macros to use in template declarations in order to select
 * the correct function.
 */
#define IS_ITERATOR(X) typename std::enable_if<!std::is_same<typename std::iterator_traits<X>::value_type, void>::value>::type* =nullptr
#define IS_NOT_ITERATOR(X) typename std::enable_if<std::is_same<typename std::iterator_traits<X>::value_type, void>::value>::type* =nullptr

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
