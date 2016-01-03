# Quaternions  [![Build status](https://travis-ci.org/FrankAstier/quaternions.svg?branch=master)](https://travis-ci.org/FrankAstier/quaternions) [![Coverage Status](https://coveralls.io/repos/FrankAstier/quaternions/badge.svg?branch=master&service=github&bust=1)](https://coveralls.io/github/FrankAstier/quaternions?branch=master)

A library to work with quaternions, as a single header file.

## Design objectives
- This library was designed to be simple, fast and convenient.
- For simplicity, there is a single (parametric) class, quaternion::Quaternion, that lives in a single header file.
- Computing the power of a quaternion in particular has been optimized to be significantly faster than boost.
- The methods provided have been designed to make it as natural as possible to use quaternion::Quaternion.

## Notes
- Boost provides quaternions at: http://www.boost.org/doc/libs/1_59_0/libs/math/doc/html/quaternions.html
  This implementation is as fast as or faster than Boost. In particular pow is much faster than Boost. The speedups
  in this implementation were obtained by factorizing the a few low powers, and re-using those factorizations for
  higher powers.
- I have tried to use intrinsics (SSE), but didn't find it to be faster than "naive" code.
- Expression templates: unless the expression templates do some serious work, gcc can optimize the code to get very
  good performance, essentially equal to expression templates. So expression templates seem to be an older technique
  that's no longer required by modern compilers (although clang seems to lag compared to gcc in terms of optimization).

## Requirements
- Boost is required for the unit tests, but not for the Quaternion class itself.

## Others/References
- http://www.boost.org/doc/libs/1_59_0/libs/math/doc/html/quaternions.html
- *vectorclass* from Agner Fog, can be found at: http://www.agner.org/optimize/, and provides quaternion classes.
- http://www.geometrictools.com/index.html