#pragma once

#include <cstddef>
#include <utility>

#if defined(GRAPHOS_HAVE_UMPIRE)
#include "umpire/Allocator.hpp"
#include "umpire/ResourceManager.hpp"
#endif

namespace graphos::exec {

// Kernel scratch. Under Umpire the storage comes from the shared pools; the
// fallback is heap memory behind the same interface, so operation code has
// one path.
//
// T must be trivially constructible and destructible: pool allocations are
// raw bytes.
template <typename T>
class Buffer {
 public:
  Buffer() = default;
  explicit Buffer(std::size_t n) { allocate(n); }
  Buffer(const Buffer&) = delete;
  Buffer& operator=(const Buffer&) = delete;
  Buffer(Buffer&& other) noexcept { swap(other); }
  Buffer& operator=(Buffer&& other) noexcept {
    release();
    swap(other);
    return *this;
  }
  ~Buffer() { release(); }

  T* data() { return data_; }
  const T* data() const { return data_; }
  std::size_t size() const { return size_; }
  T& operator[](std::size_t i) { return data_[i]; }
  const T& operator[](std::size_t i) const { return data_[i]; }

 private:
  void allocate(std::size_t n) {
    size_ = n;
    if (n == 0) return;
#if defined(GRAPHOS_HAVE_UMPIRE)
    auto& rm = umpire::ResourceManager::getInstance();
    data_ = static_cast<T*>(rm.getAllocator("HOST").allocate(n * sizeof(T)));
#else
    data_ = new T[n];
#endif
  }

  void release() {
    if (data_ == nullptr) return;
#if defined(GRAPHOS_HAVE_UMPIRE)
    umpire::ResourceManager::getInstance().getAllocator("HOST").deallocate(data_);
#else
    delete[] data_;
#endif
    data_ = nullptr;
    size_ = 0;
  }

  void swap(Buffer& other) noexcept {
    std::swap(data_, other.data_);
    std::swap(size_, other.size_);
  }

  T* data_ = nullptr;
  std::size_t size_ = 0;
};

}  // namespace graphos::exec
