#pragma once

#include <cstddef>
#include <utility>
#include <vector>

#if defined(GRAPHOS_HAVE_CHAI)
#include "chai/ManagedArray.hpp"
#endif

namespace graphos::exec {

// Persistent storage for frozen complexes — the CHAI seam. With CHAI
// enabled, data lives in a chai::ManagedArray whose value-copies migrate
// between memory spaces when captured into RAJA kernels; the fallback is
// host memory behind the same interface. Distinct from Buffer, which is
// kernel-scratch: Array outlives kernels and is what device execution
// policies will read.
//
// Move-only RAII owner. host access through data()/operator[]; view()
// yields the handle to capture BY VALUE in kernels (ManagedArray with CHAI,
// raw pointer on the host fallback).
template <typename T>
class Array {
 public:
#if defined(GRAPHOS_HAVE_CHAI)
  using view_type = chai::ManagedArray<T>;
#else
  using view_type = T*;
#endif

  Array() = default;

  explicit Array(std::size_t n) { allocate(n); }

  Array(const T* src, std::size_t n) {
    allocate(n);
    T* d = data();
    for (std::size_t i = 0; i < n; ++i) d[i] = src[i];
  }

  explicit Array(const std::vector<T>& v) : Array(v.data(), v.size()) {}

  Array(const Array&) = delete;
  Array& operator=(const Array&) = delete;

  Array(Array&& other) noexcept { swap_with(other); }

  Array& operator=(Array&& other) noexcept {
    if (this != &other) {
      release();
      swap_with(other);
    }
    return *this;
  }

  ~Array() { release(); }

  T& operator[](std::size_t i) { return data()[i]; }
  const T& operator[](std::size_t i) const { return data()[i]; }

#if defined(GRAPHOS_HAVE_CHAI)

  T* data() { return size_ == 0 ? nullptr : arr_.data(chai::CPU); }
  const T* data() const {
    return size_ == 0 ? nullptr : const_cast<chai::ManagedArray<T>&>(arr_).data(chai::CPU);
  }
  std::size_t size() const { return size_; }
  view_type view() { return arr_; }

 private:
  void allocate(std::size_t n) {
    size_ = n;
    if (n != 0) arr_.allocate(n);
  }
  void release() {
    if (size_ != 0) arr_.free();
    size_ = 0;
  }
  void swap_with(Array& other) noexcept {
    std::swap(arr_, other.arr_);
    std::swap(size_, other.size_);
  }

  chai::ManagedArray<T> arr_;
  std::size_t size_ = 0;

#else

  T* data() { return data_.data(); }
  const T* data() const { return data_.data(); }
  std::size_t size() const { return data_.size(); }
  view_type view() { return data_.data(); }

 private:
  void allocate(std::size_t n) { data_.resize(n); }
  void release() { data_.clear(); }
  void swap_with(Array& other) noexcept { std::swap(data_, other.data_); }

  std::vector<T> data_;

#endif
};

}  // namespace graphos::exec
