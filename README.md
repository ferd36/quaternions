# Quaternions
A library to work with quaternions.

## Notes
- Boost provides quaternions at: http://www.boost.org/doc/libs/1_59_0/libs/math/doc/html/quaternions.html
  This implementation is as fast as or faster than Boost. In particular pow is much faster than Boost. The speedups
  in this implementation were obtained by factorizing the a few low powers, and re-using those factorizations for
  higher powers.
- I have tried to use intrinsics (SSE), but didn't find it to be faster than "naive" code.

## TODO
- Check if the multiplication in ASM uses SSE instructions correctly.

## Others
- *vectorclass* from Agner Fog, can be found at: http://www.agner.org/optimize/, and provides quaternion classes.