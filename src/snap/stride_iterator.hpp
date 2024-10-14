#pragma once

// C/C++
#include <cassert>
#include <iterator>

template <typename T>
class StrideIterator {
 public:
  // public typedefs
  typedef typename std::iterator_traits<T>::value_type value_type;
  typedef typename std::iterator_traits<T>::reference reference;
  typedef typename std::iterator_traits<T>::difference_type difference_type;
  typedef typename std::iterator_traits<T>::pointer pointer;
  typedef std::random_access_iterator_tag iterator_category;

  // constructors
  StrideIterator() : data(NULL), step(0){};
  StrideIterator(const StrideIterator& x) : data(x.data), step(x.step) {}
  StrideIterator(T x, difference_type n = 1) : data(x), step(n) {}
  StrideIterator(std::vector<T> x) : data(x.data()), step(1) {}

  difference_type stride() const { return step; }

  // operators
  StrideIterator& operator++() {
    data += step;
    return *this;
  }

  StrideIterator operator++(int) {
    StrideIterator tmp = *this;
    data += step;
    return tmp;
  }

  StrideIterator& operator+=(difference_type x) {
    data += x * step;
    return *this;
  }

  StrideIterator& operator--() {
    data -= step;
    return *this;
  }

  StrideIterator operator--(int) {
    StrideIterator tmp = *this;
    data -= step;
    return tmp;
  }

  StrideIterator& operator-=(difference_type x) {
    data -= x * step;
    return *this;
  }

  reference operator[](difference_type n) const { return data[n * step]; }

  reference operator*() const { return *data; }

  // friend operators
  friend bool operator==(const StrideIterator& x, const StrideIterator& y) {
    assert(x.step == y.step);
    return x.data == y.data;
  }

  friend bool operator!=(const StrideIterator& x, const StrideIterator& y) {
    assert(x.step == y.step);
    return x.data != y.data;
  }

  friend bool operator<(const StrideIterator& x, const StrideIterator& y) {
    assert(x.step == y.step);
    return x.data < y.data;
  }

  friend difference_type operator-(const StrideIterator& x,
                                   const StrideIterator& y) {
    assert(x.step == y.step);
    return (x.data - y.data) / x.step;
  }

  friend StrideIterator operator+(const StrideIterator& x, difference_type y) {
    return x += y * x.step;
  }

  friend StrideIterator operator+(difference_type x, const StrideIterator& y) {
    return y += x * y.step;
  }

 private:
  T data;
  difference_type step;
};
