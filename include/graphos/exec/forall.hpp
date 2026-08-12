#pragma once

#include <utility>

#include "graphos/core/types.hpp"

#if defined(GRAPHOS_HAVE_RAJA)
#include "RAJA/RAJA.hpp"
#endif

namespace graphos::exec {

// The execution seam: every bulk kernel phase in graphos ops goes through
// these primitives. With RAJA enabled they dispatch to RAJA policies
// (OpenMP when RAJA was built with it, sequential otherwise); without RAJA
// they are plain loops with identical semantics.
//
// Kernel bodies must therefore be data-parallel over their index range:
// writes to distinct indices only, no cross-iteration ordering assumptions.
// Device policies (CUDA/HIP/SYCL) arrive once persistent storage moves to
// CHAI-managed arrays; the op code is written to survive that unchanged.

#if defined(GRAPHOS_HAVE_RAJA)

#if defined(RAJA_ENABLE_OPENMP)
using for_policy = RAJA::omp_parallel_for_exec;
#else
using for_policy = RAJA::seq_exec;
#endif

template <typename F>
void forall(Index n, F&& body) {
  RAJA::forall<for_policy>(RAJA::TypedRangeSegment<Index>(0, n), std::forward<F>(body));
}

inline Index exclusive_scan(const Index* in, Index* out, Index n) {
  if (n == 0) return 0;
  RAJA::exclusive_scan<for_policy>(RAJA::make_span(in, n), RAJA::make_span(out, n));
  return out[n - 1] + in[n - 1];
}

inline void inclusive_scan_inplace(Index* a, Index n) {
  if (n == 0) return;
  RAJA::inclusive_scan_inplace<for_policy>(RAJA::make_span(a, n));
}

#else

template <typename F>
void forall(Index n, F&& body) {
  for (Index i = 0; i < n; ++i) body(i);
}

inline Index exclusive_scan(const Index* in, Index* out, Index n) {
  Index running = 0;
  for (Index i = 0; i < n; ++i) {
    out[i] = running;
    running += in[i];
  }
  return running;
}

inline void inclusive_scan_inplace(Index* a, Index n) {
  for (Index i = 1; i < n; ++i) a[i] += a[i - 1];
}

#endif

}  // namespace graphos::exec
