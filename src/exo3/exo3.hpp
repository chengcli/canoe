#pragma once

#include <snap/stride_iterator.hpp>

// covariant velocity to contravariant velocity
template <typename T>
std::array<Real, 3> vec_raise(StrideIterator<T*> w, StrideIterator<T*> m) {
  std::array<Real, 3> v{w[IVX], w[IVY], w[IVZ]};
  return v;
}

template <typename T>
void vec_raise_inplace(StrideIterator<T*> w, StrideIterator<T*> m) {}

template <typename T>
std::array<Real, 3> vec_lower(StrideIterator<T*> w, StrideIterator<T*> m) {
  std::array<Real, 3> v{w[IVX], w[IVY], w[IVZ]};
  return v;
}

template <typename T>
void vec_lower_inplace(StrideIterator<T*> w, StrideIterator<T*> m) {}
