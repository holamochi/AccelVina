#include "vina_musa_launch.h"
#include "opencl_device_compat.h"

namespace {

#pragma clang attribute push(__attribute__((device)), apply_to = function)
#include "mutate_conf.cl"
#include "matrix.cl"
#include "quasi_newton.cl"
#include "kernel2.cl"

} // namespace

extern "C" int vina_musa_launch_kernel2(const void* ric, void* model,
                                         const void* pre, const void* grids,
                                         const void* random_maps, void* coords,
                                         void* results, const void* mis,
                                         int torsion_size, int search_depth,
                                         int max_bfgs_steps, int rilc_bfgs,
                                         void* stream) {
    kernel2<<<dim3(32, 16), dim3(16, 2), 0, static_cast<musaStream_t>(stream)>>>(
        static_cast<const output_type_cl*>(ric), static_cast<m_cl*>(model),
        static_cast<const pre_cl*>(pre), static_cast<const grids_cl*>(grids),
        static_cast<const ::random_maps*>(random_maps),
        static_cast<ligand_atom_coords_cl*>(coords), static_cast<output_type_cl*>(results),
        static_cast<const mis_cl*>(mis), torsion_size, search_depth, max_bfgs_steps,
        rilc_bfgs);
    return musaGetLastError();
}
