/* SPDX-License-Identifier: Apache-2.0 */
#ifndef VINA_C_KERNEL_H
#define VINA_C_KERNEL_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

enum { VINA_CK_PAIR_WIDTH = 4 };

typedef struct {
    size_t a;
    size_t begin;
    size_t end;
} vina_ck_pair_row;

typedef struct {
    size_t b;
    int64_t table[VINA_CK_PAIR_WIDTH]; /* offsets in doubles */
    unsigned mask;
} vina_ck_pair_chunk;

/* forces contains tightly packed XYZ triples. x, y and z are padded by width-1. */
double vina_ck_pairs_deriv(const vina_ck_pair_row *rows, size_t row_count,
                          const vina_ck_pair_chunk *chunks, const double *smooth,
                          const double *x, const double *y, const double *z,
                          void *forces, double factor, double cutoff_sqr, double cap);

void vina_ck_ad4_deriv(const double *const maps[3], const size_t dims[3],
                      const double origin[3], const double factor[3],
                      const double factor_inv[3], const double last[3],
                      const double location[3], double slope, double cap,
                      double values[3], double *const gradients[3]);

#ifdef __cplusplus
}
#endif

#endif
