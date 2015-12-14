//
// Created by Frank Astier on 2015/12/13.
//

#ifndef QUATERNIONS_GENERATORS_H_H
#define QUATERNIONS_GENERATORS_H_H

#include "Quaternion.h"

#include <random>

using namespace std;

template <typename T, typename G>
inline Quaternion<T> random_quaternion(G& g) {
  return Quaternion<T>(g(),g(),g(),g());
}

#endif //QUATERNIONS_GENERATORS_H_H
