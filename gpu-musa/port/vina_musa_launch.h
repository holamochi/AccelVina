#pragma once

#ifdef __cplusplus
extern "C" {
#endif

int vina_musa_launch_kernel1(const void* pre, const void* pa, const void* gb,
                             const void* ar, void* grids, const void* mis,
                             const void* needed, int atu, int nat, void* stream);

int vina_musa_launch_kernel2(const void* ric, void* model, const void* pre,
                             const void* grids, const void* random_maps,
                             void* coords, void* results, const void* mis,
                             int torsion_size, int search_depth, int max_bfgs_steps,
                             int rilc_bfgs, void* stream);

#ifdef __cplusplus
}
#endif
