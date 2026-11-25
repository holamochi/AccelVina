#include "vina_musa_launch.h"
#include "opencl_device_compat.h"

namespace {

#pragma clang attribute push(__attribute__((device)), apply_to = function)
#include "kernel1.cl"

} // namespace

extern "C" int vina_musa_launch_kernel1(const void* pre, const void* pa,
                                         const void* gb, const void* ar,
                                         void* grids, const void* mis,
                                         const void* needed, int atu, int nat,
                                         void* stream) {
    kernel1<<<dim3(32, 32, 32), dim3(4, 4, 4), 0, static_cast<musaStream_t>(stream)>>>(
        static_cast<const pre_cl*>(pre), static_cast<const pa_cl*>(pa),
        static_cast<const gb_cl*>(gb), static_cast<const ar_cl*>(ar),
        static_cast<grids_cl*>(grids), static_cast<const mis_cl*>(mis),
        static_cast<const float*>(needed), atu, nat);
    return musaGetLastError();
}
