#pragma once

#include <musa_runtime.h>
#include <math.h>
#include <stdio.h>

#include "code_head.cl"

static_assert(sizeof(0.0) == sizeof(float),
              "Use -Xclang -cl-single-precision-constant for the Vina device code");

namespace {

__device__ inline int get_global_id(int axis) {
    if (axis == 0) return blockIdx.x * blockDim.x + threadIdx.x;
    if (axis == 1) return blockIdx.y * blockDim.y + threadIdx.y;
    return blockIdx.z * blockDim.z + threadIdx.z;
}

__device__ inline int get_global_size(int axis) {
    if (axis == 0) return gridDim.x * blockDim.x;
    if (axis == 1) return gridDim.y * blockDim.y;
    return gridDim.z * blockDim.z;
}

__device__ inline int get_global_linear_id() {
    return get_global_id(0) + get_global_size(0) *
        (get_global_id(1) + get_global_size(1) * get_global_id(2));
}

template <int exponent>
__device__ inline float opencl_pown(float value) {
    static_assert(exponent == 2, "Vina only uses pown with exponent 2");
    return value * value;
}

} // namespace

#define __global
#define __constant const
// The upstream files put their device helpers before the entry point.
#define __kernel _Pragma("clang attribute pop") __global__

#ifdef M_PI
#undef M_PI
#endif
#define M_PI 3.14159265358979323846f
#define pown(value, exponent) opencl_pown<exponent>(value)
#define sqrt sqrtf
#define cos cosf
#define sin sinf
#define exp expf
#define fabs fabsf
#define ceil ceilf
