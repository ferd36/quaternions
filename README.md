# Quaternions  [![Build status](https://travis-ci.org/FrankAstier/quaternions.svg?branch=master)](https://travis-ci.org/FrankAstier/quaternions) [![Coverage Status](https://coveralls.io/repos/FrankAstier/quaternions/badge.svg?branch=master&service=github)](https://coveralls.io/github/FrankAstier/quaternions?branch=master)

A library to work with quaternions, as a single header file.

## Requirements
- Requires boost (http://www.boost.org) for the unit tests, but not for the Quaternion class itself.

## Notes
- Boost provides quaternions at: http://www.boost.org/doc/libs/1_59_0/libs/math/doc/html/quaternions.html
  This implementation is as fast as or faster than Boost. In particular pow is much faster than Boost. The speedups
  in this implementation were obtained by factorizing the a few low powers, and re-using those factorizations for
  higher powers.
- I have tried to use intrinsics (SSE), but didn't find it to be faster than "naive" code.

## Others
- *vectorclass* from Agner Fog, can be found at: http://www.agner.org/optimize/, and provides quaternion classes.