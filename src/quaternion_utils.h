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
 * Utilities for working with quaternion.
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
 */
template<typename T, typename T1>
inline bool is_scalar_zero(T x, T1 eps = 0) {
  return std::abs(x) <= eps;
}

/**
 * Compares 2 floating point numbers "relatively": if the numbers are
 * very large, differences are still "small" if they are "small"
 * relative to the magnitudes of the quantities.
 */
template<typename T, typename T2>
inline bool is_near_equal_relative(T x, T y, T2 eps = 0) {
  if (x == 0)
    return is_scalar_zero(y, eps);
  else if (y == 0)
    return is_scalar_zero(x, eps);
  else
    return is_scalar_zero((x-y)/std::min(x,y), eps);
}

/**
 * TODO: try this "fast" cos
 *
 * template<typename T>
inline T cos(T x) noexcept
{
    constexpr T tp = 1./(2.*M_PI);
    x *= tp;
    x -= T(.25) + std::floor(x + T(.25));
    x *= T(16.) * (std::abs(x) - T(.5));
    #if EXTRA_PRECISION
    x += T(.225) * x * (std::abs(x) - T(1.));
    #endif
    return x;
}

 http://forum.devmaster.net/t/fast-and-accurate-sine-cosine/9648
 http://stackoverflow.com/questions/345085/how-do-trigonometric-functions-work/394512#394512
 https://en.wikipedia.org/wiki/CORDIC
 https://en.wikipedia.org/wiki/Approximation_theory
 https://en.wikipedia.org/wiki/Methods_of_computing_square_roots#Approximations_that_depend_on_the_floating_point_representation
 // cos(x) = sin(x + pi/2)
addps xmm0, PI_2
movaps xmm1, xmm0
cmpnltps xmm1, PI
andps xmm1, PIx2
subps xmm0, xmm1

// Parabola
movaps xmm1, xmm0
andps xmm1, abs
mulps xmm1, xmm0
mulps xmm0, B
mulps xmm1, C
addps xmm0, xmm1

// Extra precision
movaps xmm1, xmm0
andps xmm1, abs
mulps xmm1, xmm0
subps xmm1, xmm0
mulps xmm1, P
addps xmm0, xmm1
 */

#endif //QUATERNIONS_UTILS_H
