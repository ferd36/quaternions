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

#ifndef QUATERNIONS_QUATERNION_IO_H
#define QUATERNIONS_QUATERNION_IO_H

#include <iostream>

#include "quaternion.h"
#include "quaternion_utils.h"

namespace quaternion {

/**
 * Facilities to display quaternions for humans as well as facilities to read/write them with C++ streams.
 */
struct QuaternionIO {
  /**
   * Print format control flags.
   */
  static long double scalar_zero_threshold;
  // if 0, does "hard" equality tests for zero
  static int print_style;

  /**
   * Print a quaternion to a stream in various formats.
   * TODO: introduce eps and make faster with constants?
   */
  template<typename T>
  static std::ostream& print(std::ostream& out, const Quaternion <T>& q) {
    if (print_style == 0) {
      if (q == 0)
        return out << 0;
      if (q == Quaternion<T>(1))
        return out << 1;
      if (q == Quaternion<T>(-1))
        return out << -1;
      if (q == Quaternion<T>(0, 1))
        return out << "i";
      if (q == Quaternion<T>(0, -1))
        return out << "-i";
      if (q == Quaternion<T>(0, 0, 1))
        return out << "j";
      if (q == Quaternion<T>(0, 0, -1))
        return out << "-j";
      if (q == Quaternion<T>(0, 0, 0, 1))
        return out << "k";
      if (q == Quaternion<T>(0, 0, 0, -1))
        return out << "-k";
      auto s = [](T x) { return x < 0 ? "" : "+"; }; // print out the sign correctly
      if (!is_scalar_zero(q.a(), QuaternionIO::scalar_zero_threshold))
        out << q.a();
      if (!is_scalar_zero(q.b(), QuaternionIO::scalar_zero_threshold))
        out << s(q.b()) << q.b() << "i";
      if (!is_scalar_zero(q.c(), QuaternionIO::scalar_zero_threshold))
        out << s(q.c()) << q.c() << "j";
      if (!is_scalar_zero(q.d(), QuaternionIO::scalar_zero_threshold))
        out << s(q.d()) << q.d() << "k";
    } else if (print_style == 1) {
      out << "{" << q.a() << "," << q.b() << "," << q.c() << "," << q.d() << "}";
    }
    return out;
  }

  // TODO: read from stream
};

/**
 * IO manipulators to control the format when printing quaternions out to a stream.
 */
long double QuaternionIO::scalar_zero_threshold = 0;
int QuaternionIO::print_style;

struct SetScalarZeroThreshold {
  long double eps = 0;
};
inline SetScalarZeroThreshold set_display_eps(long double eps) {
  SetScalarZeroThreshold sszt;
  sszt.eps = eps;
  return sszt;
}

inline std::ostream& operator<<(std::ostream& out, SetScalarZeroThreshold sszt) {
  QuaternionIO::scalar_zero_threshold = sszt.eps;
  return out;
}

enum DisplayStyle {
  q_nice, q_compact
};
struct SetDisplayStyle {
  DisplayStyle style;
};
inline SetDisplayStyle set_display_style(DisplayStyle ds) {
  SetDisplayStyle sds;
  sds.style = ds;
  return sds;
}

inline std::ostream& operator<<(std::ostream& out, SetDisplayStyle sds) {
  QuaternionIO::print_style = sds.style;
  return out;
}

/**
 * This streaming operator made me wonder if I should sneak "smart" code
 * in the quaternion arithmetic, in order to optimize it for space, but that
 * turned out not worthwhile (see CQuaternion).
 * TODO: control format for file or human readable. Also write operator>>
 */
template<typename T>
inline std::ostream& operator<<(std::ostream& out, const Quaternion <T>& q) {
  return QuaternionIO::print(out, q);
}

} // end namespace quaternion
#endif //QUATERNIONS_QUATERNION_IO_H
